import bz2
import json
import os

import pandas as pd
import pytest

import pyexocross.api as api
from pyexocross.base.large_file import read_trans_chunks
from pyexocross.config import Config
from pyexocross.database.data import (
    filterlinelist,
    loaddata,
    memorysource,
    parquetsource,
    readlinecache,
    writelinecache,
)
from pyexocross.database.load_exomol import get_part_transfiles, get_statesfile
from pyexocross.process.stick_xsec_filepath import cross_section_filepath, crosssectiondetails


def makedata(tmp_path):
    folder = tmp_path / 'CO' / '12C-16O' / 'TEST'
    folder.mkdir(parents=True)
    base = folder / '12C-16O__TEST'
    with bz2.open(str(base) + '.states.bz2', 'wt') as stream:
        stream.write('1 0.0 1 0\n')
    with open(str(base) + '.states', 'w', encoding='utf-8') as stream:
        stream.write('1 0.0 1 0\n')
    with bz2.open(str(base) + '.trans.bz2', 'wt') as stream:
        stream.write('2 1 1.0e-3\n')
    with open(str(base) + '.trans', 'w', encoding='utf-8') as stream:
        stream.write('2 1 2.0e-3\n3 1 3.0e-3\n')
    return folder, str(base)


def makerangedata(tmp_path):
    folder = tmp_path / 'NH3' / '14N-1H3' / 'CoYuTe'
    folder.mkdir(parents=True)
    states = folder / '14N-1H3__CoYuTe.states'
    transitions = folder / '14N-1H3__CoYuTe__03300-03400.trans'
    with open(states, 'w', encoding='utf-8') as stream:
        stream.write('1 0.0 1 0\n2 3335.0 1 1\n3 3395.0 1 1\n')
    with open(transitions, 'w', encoding='utf-8') as stream:
        stream.write('2 1 2.0e-3\n3 1 3.0e-3\n')
    return transitions


def test_uncompressed_files_are_preferred(tmp_path):
    _, base = makedata(tmp_path)
    info = ['CO', '12C-16O', 'TEST']
    read_path = str(tmp_path) + os.sep

    assert get_statesfile(read_path, info) == base + '.states'
    assert get_part_transfiles(
        read_path,
        info,
        0,
        100,
        prepare=False,
    ) == [base + '.trans']


def test_memory_and_parquet_sources_select_columns(tmp_path):
    _, base = makedata(tmp_path)
    filepath = base + '.trans'

    memory = memorysource(filepath, 1)
    memory_chunks = list(read_trans_chunks(memory, [0, 2], ['uid', 'A'], 1))
    assert pd.concat(memory_chunks, ignore_index=True).to_dict('list') == {
        'uid': [2, 3],
        'A': [2.0e-3, 3.0e-3],
    }

    parquet = parquetsource(filepath, str(tmp_path / 'cache'), 1)
    parquet_chunks = list(read_trans_chunks(parquet, [0, 2], ['uid', 'A'], 1))
    assert pd.concat(parquet_chunks, ignore_index=True).to_dict('list') == {
        'uid': [2, 3],
        'A': [2.0e-3, 3.0e-3],
    }


def test_device_and_run_mode_are_compatible_aliases(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    common = {
        'database': 'ExoMol',
        'molecule': 'CO',
        'isotopologue': '12C-16O',
        'dataset': 'TEST',
    }

    config = Config(device='GPU', **common)
    assert config.device == 'GPU'
    assert config.run_mode == 'GPU'
    assert config.cache == 'auto'
    assert config.max_memory == 512

    with pytest.raises(ValueError, match='same CPU/GPU mode'):
        Config(device='CPU', run_mode='GPU', **common)


def test_api_abundance_defaults_to_one_and_accepts_override(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    common = {
        'database': 'ExoMol',
        'molecule': 'CO',
        'isotopologue': '12C-16O',
        'dataset': 'TEST',
    }

    assert Config(**common).abundance == 1.0
    assert Config(abundance=0.75, **common).abundance == 0.75


def test_transition_apis_forward_loaded_data(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    monkeypatch.setattr(api, '_ensure_logging', lambda *args, **kwargs: None)
    calls = []
    monkeypatch.setattr(api, 'get_results', lambda config, data=None: calls.append((config, data)))
    config = Config(
        database='ExoMol',
        molecule='CO',
        isotopologue='12C-16O',
        dataset='TEST',
    )
    loaded = api.LoadedData(
        config=config,
        states=pd.DataFrame(),
        preparedstates=pd.DataFrame(),
        transitions=(),
        cache='memory',
        storage='memory',
        cachedir='',
        alltrans=True,
    )

    functions = (
        api.conversion,
        api.lifetimes,
        api.cooling_functions,
        api.oscillator_strengths,
    )
    for function in functions:
        function(data=loaded)

    assert [data for _, data in calls] == [loaded] * len(functions)
    assert calls[0][0].conversion == 1
    assert calls[1][0].lifetimes == 1
    assert calls[2][0].cooling_functions == 1
    assert calls[3][0].oscillator_strengths == 1


def test_whole_list_calculations_expand_loaded_data_once(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    config = Config(
        database='ExoMol',
        molecule='CO',
        isotopologue='12C-16O',
        dataset='TEST',
        cache='parquet',
    )
    ranged = api.LoadedData(
        config=config,
        states=pd.DataFrame({'id': [1]}),
        preparedstates=pd.DataFrame({'id': [1]}),
        transitions=('range',),
        cache='parquet',
        storage='parquet',
        cachedir='cache',
        alltrans=False,
    )
    expanded = api.LoadedData(
        config=config,
        states=pd.DataFrame({'id': [1]}),
        preparedstates=pd.DataFrame({'id': [1]}),
        transitions=('all',),
        cache='parquet',
        storage='parquet',
        cachedir='cache',
        alltrans=True,
    )
    loads = []

    def expand(*args, **kwargs):
        loads.append(kwargs)
        return expanded

    monkeypatch.setattr(api, 'loaddata', expand)

    assert api.requirealltransitions(ranged) is ranged
    assert api.requirealltransitions(ranged) is ranged
    assert ranged.alltrans is True
    assert ranged.transitions == ('all',)
    assert len(loads) == 1
    assert loads[0]['alltrans'] is True


def test_load_prints_timing(monkeypatch, capsys):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    marker = object()
    monkeypatch.setattr(api, '_ensure_logging', lambda *args, **kwargs: None)
    monkeypatch.setattr(api, 'loaddata', lambda *args, **kwargs: marker)

    result = api.load(
        database='ExoMol',
        molecule='CO',
        isotopologue='12C-16O',
        dataset='TEST',
    )

    output = capsys.readouterr().out
    assert result is marker
    assert 'Loading reusable line-list data' in output
    assert 'Finished loading reusable line-list data' in output
    assert 'Running time on CPU' in output
    assert 'Running time on system' in output


def test_default_parquet_cache_is_stored_with_input_data(tmp_path):
    makedata(tmp_path)

    class CacheConfig:
        database = 'ExoMol'
        read_path = str(tmp_path) + os.sep
        data_info = ['CO', '12C-16O', 'TEST']
        min_wn = 0
        max_wn = 100
        chunk_size = 1
        check_uncertainty = False
        states_col = []
        states_fmt = []

    loaded = loaddata(
        CacheConfig(),
        cache='parquet',
        alltrans=True,
        preparestates=False,
    )

    expected = tmp_path / 'CO' / '12C-16O' / 'TEST' / '.pyexocross_cache'
    assert loaded.cachedir == str(expected)
    assert sorted(path.suffix for path in expected.iterdir()) == ['.json', '.parquet', '.parquet']
    statespath = expected / '12C-16O__TEST.states.parquet'
    assert statespath.exists()
    statesmtime = statespath.stat().st_mtime_ns

    loaddata(
        CacheConfig(),
        cache='parquet',
        alltrans=True,
        preparestates=False,
    )
    assert statespath.stat().st_mtime_ns == statesmtime


def test_cross_section_filename_formats_binsize_and_cutoff(tmp_path):
    filepath = cross_section_filepath(
        str(tmp_path) + os.sep,
        ['NO', '14N-16O', 'XABC'],
        300,
        1,
        0,
        [],
        [],
        '1000',
        '1100',
        'cm-1__',
        'WN',
        'None',
        'None',
        'ExoMol',
        'Ab',
        0.01,
        25.0,
        'Voigt',
        '__LTE',
        '',
        'L',
        True,
    )

    assert '__BinSize0.0100cm-1__Cutoff25.0__Voigt__' in filepath
    assert crosssectiondetails(0.01, 'cm-1__', 25.0, 'Voigt') in filepath


def test_range_cache_uses_cutoff_without_adding_adjacent_files(tmp_path):
    source = makerangedata(tmp_path)

    class RangeConfig:
        database = 'ExoMol'
        read_path = str(tmp_path) + os.sep
        data_info = ['NH3', '14N-1H3', 'CoYuTe']
        min_wn = 3380
        max_wn = 3390
        cutoff = 30
        chunk_size = 1
        check_uncertainty = False
        states_col = []
        states_fmt = []

    loaded = loaddata(
        RangeConfig(),
        cache='parquet',
        preparestates=False,
    )

    rangesdir = source.parent / '.pyexocross_cache' / 'ranges'
    expected = rangesdir / '14N-1H3__CoYuTe__03350-03400.trans.parquet'
    assert loaded.transitions[0].path == str(expected)
    assert loaded.transitions[0].minwn == 3350
    assert loaded.transitions[0].maxwn == 3400

    chunks = list(read_trans_chunks(loaded.transitions[0], [0, 2], ['uid', 'A'], 10))
    assert pd.concat(chunks, ignore_index=True).to_dict('list') == {
        'uid': [3],
        'A': [3.0e-3],
    }

    manifestpath = source.parent / '.pyexocross_cache' / 'manifest.json'
    with open(manifestpath, 'r', encoding='utf-8') as stream:
        manifest = json.load(stream)
    assert manifest['ranges'][0]['range'] == [3350, 3400]
    assert manifest['ranges'][0]['source']['path'] == str(source)

    cachemtime = expected.stat().st_mtime_ns
    RangeConfig.min_wn = 3382
    RangeConfig.max_wn = 3388
    RangeConfig.cutoff = 5
    reused = loaddata(
        RangeConfig(),
        cache='parquet',
        preparestates=False,
    )
    assert reused.transitions[0].path == str(expected)
    assert expected.stat().st_mtime_ns == cachemtime
    with open(manifestpath, 'r', encoding='utf-8') as stream:
        manifest = json.load(stream)
    assert len(manifest['ranges']) == 1

    narrower = loaddata(
        RangeConfig(),
        cache='parquet',
        preparestates=False,
        refresh=True,
    )
    narrowpath = rangesdir / '14N-1H3__CoYuTe__03377-03393.trans.parquet'
    assert narrower.transitions[0].path == str(narrowpath)

    RangeConfig.min_wn = 3383
    RangeConfig.max_wn = 3387
    RangeConfig.cutoff = 4
    smallest = loaddata(
        RangeConfig(),
        cache='parquet',
        preparestates=False,
    )
    assert smallest.transitions[0].path == str(narrowpath)


def test_normalized_line_cache_range_and_database_filters(tmp_path):
    frame = pd.DataFrame({
        'v': [10.0, 20.0, 30.0],
        'S': [1e-25, 1e-20, 1e-15],
        'Unc': [1, 5, 9],
        'A': [1.0, 2.0, 3.0],
    })
    cachepath = tmp_path / 'HITRAN__NO__14N-16O.linelist.parquet'
    writelinecache(frame, str(cachepath))

    ranged = readlinecache(str(cachepath), 15, 30)
    filtered = filterlinelist(
        ranged,
        'HITRAN',
        15,
        30,
        uncfilter=1e-5,
        threshold=1e-18,
    )
    assert filtered['v'].tolist() == [30.0]
