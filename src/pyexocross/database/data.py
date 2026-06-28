"""Reusable database loading and line-list cache support."""

import bz2
import copy
import json
import math
import os
import re
import time
from dataclasses import dataclass

import pandas as pd

from pyexocross.base.large_file import TransSource
from pyexocross.database.load_exomol import (
    get_part_transfiles,
    get_statesfile,
    get_transfiles,
    read_all_states,
    read_part_states,
)


PREPROCESS_KEYS = {
    'database',
    'molecule',
    'isotopologue',
    'dataset',
    'atom',
    'read_path',
    'unc_filter',
    'nlte_method',
    'nlte_path',
    'qnslabel_list',
    'qnsformat_list',
    'qns_filter',
    'qns_label',
    'qns_value',
    'vib_label',
    'rot_label',
    'states_col',
    'states_fmt',
    'check_uncertainty',
    'check_lifetime',
    'check_gfactor',
    'cache',
    'cache_dir',
    'max_memory',
    'refresh',
    'refresh_cache',
}


def uncvalue(value):
    """Normalize a disabled or numeric uncertainty filter."""
    if value in (None, 'None'):
        return None
    return float(value)


def qnconstraints(config):
    """Return effective state-level QN constraints, excluding wildcard labels."""
    constraints = {}
    labels = getattr(config, 'qns_label', [])
    valueslist = getattr(config, 'qns_value', [])
    for label, values in zip(labels, valueslist):
        allowed = []
        wildcard = False
        for value in values:
            parts = [part.strip() for part in str(value).split(',')]
            if '' in parts:
                wildcard = True
                break
            allowed.extend(parts)
        if not wildcard and allowed:
            constraints[label] = set(allowed)
    return constraints


def qnrefines(loadedconfig, requestedconfig):
    """Return whether requested QN constraints only narrow loaded constraints."""
    loaded = qnconstraints(loadedconfig)
    requested = qnconstraints(requestedconfig)
    for label, loadedvalues in loaded.items():
        requestedvalues = requested.get(label)
        if requestedvalues is None or not requestedvalues.issubset(loadedvalues):
            return False
    return True


def fileidentity(filepath):
    """Return source metadata used to validate persistent caches."""
    stat = os.stat(filepath)
    return {
        'path': os.path.abspath(filepath),
        'size': stat.st_size,
        'mtime_ns': stat.st_mtime_ns,
    }


def loadmanifest(cachedir):
    """Read the cache manifest or create an empty manifest structure."""
    manifestpath = os.path.join(cachedir, 'manifest.json')
    if not os.path.exists(manifestpath):
        return {'version': 1, 'ranges': []}
    with open(manifestpath, 'r', encoding='utf-8') as stream:
        manifest = json.load(stream)
    if manifest.get('version') != 1:
        return {'version': 1, 'ranges': []}
    manifest.setdefault('ranges', [])
    return manifest


def savemanifest(cachedir, manifest):
    """Write the cache manifest atomically."""
    os.makedirs(cachedir, exist_ok=True)
    manifestpath = os.path.join(cachedir, 'manifest.json')
    temppath = manifestpath + '.tmp'
    with open(temppath, 'w', encoding='utf-8') as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
    os.replace(temppath, manifestpath)


def sourcebounds(filepath):
    """Return wavenumber bounds encoded in a partitioned transition filename."""
    name = os.path.basename(filepath)
    match = re.search(r'__(\d+)-(\d+)\.trans(?:\.bz2)?$', name)
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def rangebounds(filepath, minwn, maxwn, cutoff):
    """Return integer cache bounds clipped to the initially selected file."""
    expansion = 0.0 if cutoff in (None, 'None') else float(cutoff)
    lower = max(0, math.floor(float(minwn) - expansion))
    upper = math.ceil(float(maxwn) + expansion)
    physical = sourcebounds(filepath)
    if physical is not None:
        lower = max(lower, physical[0])
        upper = min(upper, physical[1])
    return lower, upper


def rangefilename(data_info, lower, upper, filepath=None):
    """Return an ExoMol-style range cache filename."""
    prefix = '__'.join(data_info[-2:])
    if filepath is not None:
        stem = os.path.basename(filepath)
        if stem.endswith('.bz2'):
            stem = stem[:-4]
        if stem.endswith('.trans'):
            stem = stem[:-6]
        prefix = re.sub(r'__\d+-\d+$', '', stem) or prefix
    return f'{prefix}__{lower:05d}-{upper:05d}.trans.parquet'


def fullrangebounds(filepath, states):
    """Return complete source bounds for an all-transition range cache."""
    physical = sourcebounds(filepath)
    if physical is not None:
        return physical
    energies = pd.to_numeric(states['E'], errors='coerce').dropna()
    if energies.empty:
        raise ValueError('Cannot determine transition cache range from empty states energies.')
    return 0, math.ceil(float(energies.max() - energies.min()))


def textcolumns(filepath):
    """Count columns in the first non-empty transition row."""
    opener = bz2.open if filepath.endswith('.bz2') else open
    with opener(filepath, 'rt', encoding='utf-8') as stream:
        for line in stream:
            values = line.split()
            if values:
                return len(values)
    return 0


def transitionchunks(filepath, chunksize):
    """Read all numeric transition columns in bounded chunks."""
    count = textcolumns(filepath)
    if count == 0:
        return iter(())
    dtypes = {0: 'int64', 1: 'int64'}
    dtypes.update({index: 'float64' for index in range(2, count)})
    kwargs = {
        'sep': r'\s+',
        'header': None,
        'chunksize': chunksize,
        'iterator': True,
        'low_memory': False,
        'dtype': dtypes,
    }
    if filepath.endswith('.bz2'):
        kwargs['compression'] = 'bz2'
    return pd.read_csv(filepath, **kwargs)


def memorysource(filepath, chunksize):
    """Load a transition file into one reusable in-memory frame."""
    chunks = list(transitionchunks(filepath, chunksize))
    frame = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    return TransSource(
        path=filepath,
        origin=filepath,
        frame=frame,
        size=int(frame.memory_usage(index=True, deep=True).sum()),
    )


def rangecandidate(manifest, cachedir, filepath, statesfile, lower, upper):
    """Return the smallest valid range cache covering the requested interval."""
    source = fileidentity(filepath)
    state_source = fileidentity(statesfile)
    candidates = []
    for entry in manifest['ranges']:
        cachepath = os.path.join(cachedir, entry['path'])
        cachedrange = entry['range']
        if (
            entry.get('source') == source
            and entry.get('states') == state_source
            and os.path.exists(cachepath)
            and cachedrange[0] <= lower
            and cachedrange[1] >= upper
        ):
            candidates.append((cachedrange[1] - cachedrange[0], cachepath, cachedrange))
    if not candidates:
        return None
    _, cachepath, cachedrange = min(candidates, key=lambda item: item[0])
    return cachepath, cachedrange[0], cachedrange[1]


def writerangecache(filepath, states, cachepath, chunksize, lower, upper):
    """Create a transition Parquet file filtered by a wavenumber interval."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    energies = states[['id', 'E']].drop_duplicates('id').set_index('id')['E']
    temppath = cachepath + '.tmp'
    writer = None
    count = textcolumns(filepath)
    for chunk in transitionchunks(filepath, chunksize):
        upper_energy = chunk[0].map(energies)
        lower_energy = chunk[1].map(energies)
        wavenumber = upper_energy - lower_energy
        selected = chunk[wavenumber.between(lower, upper)].copy()
        if selected.empty:
            continue
        selected.columns = [f'c{index}' for index in range(selected.shape[1])]
        table = pa.Table.from_pandas(selected, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(temppath, table.schema, compression='zstd')
        writer.write_table(table)
    if writer is None:
        empty = {
            f'c{index}': pd.Series(dtype='int64' if index < 2 else 'float64')
            for index in range(count)
        }
        table = pa.Table.from_pandas(pd.DataFrame(empty), preserve_index=False)
        writer = pq.ParquetWriter(temppath, table.schema, compression='zstd')
    writer.close()
    os.replace(temppath, cachepath)


def rangecachesource(
    filepath,
    states,
    statesfile,
    data_info,
    cachedir,
    chunksize,
    lower,
    upper,
    manifest,
    refresh=False,
):
    """Reuse or create a range-filtered transition cache."""
    candidate = None if refresh else rangecandidate(
        manifest,
        cachedir,
        filepath,
        statesfile,
        lower,
        upper,
    )
    if candidate is not None:
        cachepath, cachedlower, cachedupper = candidate
        print('Using range Parquet cache:', cachepath)
        return TransSource(
            path=cachepath,
            origin=filepath,
            size=os.path.getsize(cachepath),
            minwn=cachedlower,
            maxwn=cachedupper,
        )

    rangesdir = os.path.join(cachedir, 'ranges')
    os.makedirs(rangesdir, exist_ok=True)
    cachepath = os.path.join(
        rangesdir,
        rangefilename(data_info, lower, upper, filepath=filepath),
    )
    print('Creating range Parquet cache:', cachepath)
    writerangecache(filepath, states, cachepath, chunksize, lower, upper)

    relativepath = os.path.relpath(cachepath, cachedir)
    manifest['ranges'] = [
        entry for entry in manifest['ranges']
        if entry.get('path') != relativepath
    ]
    manifest['ranges'].append({
        'path': relativepath,
        'range': [lower, upper],
        'source': fileidentity(filepath),
        'states': fileidentity(statesfile),
        'created': int(time.time()),
    })
    return TransSource(
        path=cachepath,
        origin=filepath,
        size=os.path.getsize(cachepath),
        minwn=lower,
        maxwn=upper,
    )


def printcacheinfo(cachedir, manifest):
    """Print range cache count, ranges, size, and a cleanup reminder."""
    entries = []
    for entry in manifest['ranges']:
        cachepath = os.path.join(cachedir, entry['path'])
        if os.path.exists(cachepath):
            entries.append((entry, cachepath))
    print('{:45s} : {}'.format('Number of range cache files', len(entries)))
    if not entries:
        return
    ranges = [f"{entry['range'][0]}-{entry['range'][1]} cm-1" for entry, _ in entries]
    totalsize = sum(os.path.getsize(cachepath) for _, cachepath in entries)
    print('{:45s} : {}'.format('Range cache intervals', ', '.join(ranges)))
    print('{:45s} : {:.2f} MB'.format('Total range cache size', totalsize / 1024**2))
    print('Cache reminder: remove overlapping or unused range caches when they are no longer needed.')


def hitranpath(config):
    """Resolve a HITRAN/HITEMP parameter file without relying on core globals."""
    expected = f'{config.data_info[0]}__{config.data_info[1]}.par'
    cleanpath = config.read_path.rstrip('/')
    if os.path.isdir(cleanpath):
        filepath = os.path.join(
            cleanpath,
            config.data_info[0],
            config.data_info[1],
            expected,
        )
    else:
        filepath = cleanpath
    if not os.path.exists(filepath):
        raise ValueError(f'Missing HITRAN line list file: {filepath}')
    if os.path.basename(filepath) != expected:
        raise ValueError(f'Invalid HITRAN line list filename: expected {expected}')
    return filepath


def linesource(config):
    """Return the primary line-list source for non-ExoMol databases."""
    if config.database in ('HITRAN', 'HITEMP'):
        return hitranpath(config)
    if config.database == 'ExoMolHR':
        from pyexocross.database.load_exomolhr import resolve_exomolhr_filepaths

        csvpath, _, _ = resolve_exomolhr_filepaths(
            config.read_path,
            config.data_info[0],
            config.data_info[1],
        )
        return csvpath
    raise ValueError(f'Unsupported line-list cache database: {config.database}')


def parselinelist(config):
    """Parse a complete normalized HITRAN/HITEMP or ExoMolHR line list."""
    if config.database in ('HITRAN', 'HITEMP'):
        from pyexocross.database.load_hitran import read_hitran_parfile, read_parfile

        raw = read_parfile(config.read_path)
        return read_hitran_parfile(
            config.read_path,
            raw,
            -float('inf'),
            float('inf'),
            'None',
            'None',
        ).reset_index(drop=True)
    from pyexocross.database.load_exomolhr import read_exomolhr_df

    return read_exomolhr_df(
        config.read_path,
        config.data_info,
        -float('inf'),
        float('inf'),
        'None',
    ).reset_index(drop=True)


def linecacheentry(manifest, database, source):
    """Return a valid canonical line-list cache entry."""
    identity = fileidentity(source)
    for entry in manifest.setdefault('linelists', []):
        path = entry.get('path', '')
        if (
            entry.get('database') == database
            and entry.get('source') == identity
            and entry.get('range') is not None
            and re.search(r'__\d+-\d+\.linelist\.parquet$', path)
        ):
            return entry
    return None


def linecachefilename(database, source, frame):
    """Return a range-named normalized line-list cache filename."""
    wavenumbers = pd.to_numeric(frame['v'], errors='coerce').dropna()
    if wavenumbers.empty:
        raise ValueError('Cannot determine line-list cache range from an empty v column.')
    lower = max(0, math.floor(float(wavenumbers.min())))
    upper = math.ceil(float(wavenumbers.max()))
    stem = os.path.splitext(os.path.basename(source))[0]
    return f'{database}__{stem}__{lower:05d}-{upper:05d}.linelist.parquet', lower, upper


def writelinecache(frame, cachepath):
    """Write a normalized line list sorted by wavenumber with bounded row groups."""
    temppath = cachepath + '.tmp'
    ordered = frame.sort_values('v').reset_index(drop=True)
    ordered.to_parquet(
        temppath,
        engine='pyarrow',
        compression='zstd',
        index=False,
        row_group_size=100_000,
    )
    os.replace(temppath, cachepath)


def readlinecache(cachepath, minwn, maxwn):
    """Read only Parquet row groups whose statistics overlap a range."""
    import pyarrow.parquet as pq

    parquet = pq.ParquetFile(cachepath)
    names = parquet.schema_arrow.names
    vindex = names.index('v')
    frames = []
    for index in range(parquet.num_row_groups):
        column = parquet.metadata.row_group(index).column(vindex)
        stats = column.statistics
        if stats is not None and (stats.max < minwn or stats.min > maxwn):
            continue
        frame = parquet.read_row_group(index).to_pandas()
        frames.append(frame[frame['v'].between(minwn, maxwn)])
    if frames:
        return pd.concat(frames, ignore_index=True)
    return parquet.schema_arrow.empty_table().to_pandas()


def filterlinelist(frame, database, minwn, maxwn, uncfilter='None', threshold='None'):
    """Apply database-specific range, uncertainty, and intensity filters."""
    filtered = frame[frame['v'].between(minwn, maxwn)].copy()
    if database in ('HITRAN', 'HITEMP'):
        if uncfilter != 'None':
            code = int(('%e' % uncfilter)[-1])
            filtered = filtered[filtered['Unc'] >= code]
        if threshold != 'None':
            filtered = filtered[filtered['S'] >= threshold]
    elif uncfilter != 'None':
        if 'unc' not in filtered.columns:
            raise ValueError('No uncertainties in ExoMolHR line list.')
        filtered = filtered[filtered['unc'] <= uncfilter]
    return filtered.reset_index(drop=True)


def loadlinedata(config, cache, cachedir, maxmemory, refresh):
    """Load reusable HITRAN/HITEMP or ExoMolHR line-list data."""
    source = linesource(config)
    if cachedir is None:
        cachedir = os.path.join(os.path.dirname(source), '.pyexocross_cache')
    estimate = os.path.getsize(source) * 2
    memorylimit = maxmemory * 1024**2
    resolved = cache
    if cache == 'auto':
        resolved = 'none' if maxmemory <= 0 else ('memory' if estimate <= memorylimit else 'parquet')

    lineframe = None
    linepath = None
    if resolved == 'parquet':
        os.makedirs(cachedir, exist_ok=True)
        manifest = loadmanifest(cachedir)
        entry = None if refresh else linecacheentry(manifest, config.database, source)
        if entry is not None and os.path.exists(os.path.join(cachedir, entry['path'])):
            linepath = os.path.join(cachedir, entry['path'])
            print('Using line-list Parquet cache:', linepath)
        else:
            lineframe = parselinelist(config)
            filename, lower, upper = linecachefilename(
                config.database,
                source,
                lineframe,
            )
            linepath = os.path.join(
                cachedir,
                filename,
            )
            print('Creating line-list Parquet cache:', linepath)
            writelinecache(lineframe, linepath)
            relative = os.path.relpath(linepath, cachedir)
            manifest.setdefault('linelists', [])
            manifest['linelists'] = [
                item for item in manifest['linelists']
                if not (
                    item.get('database') == config.database
                    and item.get('path') == relative
                )
            ]
            manifest['linelists'].append({
                'database': config.database,
                'path': relative,
                'range': [lower, upper],
                'source': fileidentity(source),
                'created': int(time.time()),
            })
            savemanifest(cachedir, manifest)
            lineframe = None
    else:
        lineframe = parselinelist(config)

    print('Loaded data cache mode:', cache, f'({resolved})' if cache == 'auto' else '')
    return LoadedData(
        config=config,
        states=None,
        preparedstates=None,
        transitions=(),
        cache=cache,
        storage=resolved,
        cachedir=cachedir,
        lineframe=lineframe,
        linepath=linepath,
        linedatabase=config.database,
    )


def parquetsources(config, selected, states, cachedir, alltrans, refresh=False):
    """Create canonical or range-filtered Parquet transition sources."""
    statesfile = get_statesfile(config.read_path, config.data_info)
    manifest = loadmanifest(cachedir)
    sources = []
    for path in selected:
        if alltrans:
            lower, upper = fullrangebounds(path, states)
        else:
            lower, upper = rangebounds(
                path,
                config.min_wn,
                config.max_wn,
                config.cutoff,
            )
        sources.append(rangecachesource(
            path,
            states,
            statesfile,
            config.data_info,
            cachedir,
            config.chunk_size,
            lower,
            upper,
            manifest,
            refresh=refresh,
        ))
    savemanifest(cachedir, manifest)
    return tuple(sources)


def statesdata(config, cache, cachedir, refresh=False):
    """Read states from the preferred source, optionally reusing Parquet."""
    statesfile = get_statesfile(config.read_path, config.data_info)
    if cache != 'parquet':
        return read_all_states(
            config.read_path,
            config.data_info,
            config.check_uncertainty,
            config.states_col,
            config.states_fmt,
        )
    import pyarrow  # noqa: F401

    os.makedirs(cachedir, exist_ok=True)
    stem = os.path.basename(statesfile)
    if stem.endswith('.bz2'):
        stem = stem[:-4]
    parquetpath = os.path.join(cachedir, f'{stem}.parquet')
    manifest = loadmanifest(cachedir)
    stateentry = manifest.get('states', {})
    validcache = (
        stateentry.get('source') == fileidentity(statesfile)
        and stateentry.get('path') == os.path.basename(parquetpath)
        and os.path.exists(parquetpath)
    )
    if not refresh and validcache:
        import pyarrow.parquet as pq

        print('Reading states from Parquet cache:', parquetpath)
        return pq.ParquetFile(parquetpath).read().to_pandas()
    states = read_all_states(
        config.read_path,
        config.data_info,
        config.check_uncertainty,
        config.states_col,
        config.states_fmt,
    )
    temppath = parquetpath + '.tmp'
    print('Creating states Parquet cache:', parquetpath)
    states.to_parquet(temppath, engine='pyarrow', compression='zstd', index=False)
    os.replace(temppath, parquetpath)
    manifest['states'] = {
        'path': os.path.basename(parquetpath),
        'source': fileidentity(statesfile),
        'created': time.time(),
    }
    savemanifest(cachedir, manifest)
    return states


@dataclass
class LoadedData:
    """Prepared states and reusable transition sources for repeated calculations."""

    config: object
    states: pd.DataFrame
    preparedstates: pd.DataFrame
    transitions: tuple
    cache: str
    storage: str
    cachedir: str
    alltrans: bool = False
    lineframe: object = None
    linepath: object = None
    linedatabase: object = None

    def configfor(self, **kwargs):
        """Create a calculation config while protecting preprocessing identity."""
        changed = PREPROCESS_KEYS.intersection(kwargs)
        if self.linedatabase in ('HITRAN', 'HITEMP', 'ExoMolHR'):
            changed -= {
                'nlte_method',
                'nlte_path',
                'qnsformat_list',
                'qns_filter',
                'qns_label',
                'qnslabel_list',
                'qns_value',
                'vib_label',
                'rot_label',
            }
        config = copy.deepcopy(self.config)
        config._load_from_kwargs(**kwargs)
        if 'unc_filter' in changed:
            loadedunc = uncvalue(self.config.unc_filter)
            requestedunc = uncvalue(kwargs['unc_filter'])
            if (
                requestedunc is None and loadedunc is None
            ) or (
                requestedunc is not None
                and (loadedunc is None or requestedunc <= loadedunc)
            ):
                changed.remove('unc_filter')
            else:
                raise ValueError(
                    f'unc_filter={kwargs["unc_filter"]} includes states excluded by '
                    f'px.load(unc_filter={self.config.unc_filter}). Call px.load again '
                    f'with unc_filter={kwargs["unc_filter"]}.'
                )
        if 'qns_filter' in changed:
            if qnrefines(self.config, config):
                changed.remove('qns_filter')
            else:
                raise ValueError(
                    f'qns_filter={kwargs["qns_filter"]} includes states excluded by '
                    f'px.load(qns_filter={self.config.qns_filter}). Call px.load '
                    f'again with qns_filter={kwargs["qns_filter"]}.'
                )
        if changed:
            names = ', '.join(sorted(changed))
            raise ValueError(
                f'{names} affect loaded data. Call px.load again instead of overriding them.'
            )
        config.conversion = 0
        config.partition_functions = 0
        config.specific_heats = 0
        config.lifetimes = 0
        config.cooling_functions = 0
        config.oscillator_strengths = 0
        config.stick_spectra = 0
        config.cross_sections = 0
        return config

    def prepared(self, config):
        """Return prepared states with any valid stricter uncertainty filter applied."""
        if self.preparedstates is None:
            raise ValueError('Loaded data does not contain prepared states.')
        states = self.preparedstates.copy()
        loadedunc = uncvalue(self.config.unc_filter)
        requestedunc = uncvalue(config.unc_filter)
        if requestedunc is not None and (
            loadedunc is None or requestedunc < loadedunc
        ):
            if 'unc' not in states.columns:
                raise ValueError(
                    'No uncertainties are available in the loaded states. '
                    'Call px.load again without unc_filter.'
                )
            states = states[states['unc'].astype(float) <= requestedunc].copy()
        for label, allowed in qnconstraints(config).items():
            if label not in states.columns:
                raise ValueError(
                    f'Quantum number {label} is not available in loaded states. '
                    'Call px.load again with the required QN metadata.'
                )
            values = states[label].astype(str).str.strip()
            states = states[values.isin(allowed)].copy()
        return states

    def sources(self, min_wn, max_wn):
        """Select already-loaded transition sources for a calculation range."""
        selected = get_part_transfiles(
            self.config.read_path,
            self.config.data_info,
            min_wn,
            max_wn,
            prepare=False,
        )
        source_map = {os.path.abspath(source.origin): source for source in self.transitions}
        missing = [path for path in selected if os.path.abspath(path) not in source_map]
        if missing:
            raise ValueError(
                'The requested range needs transition files not included by px.load. '
                'Load the wider range before calculating.'
            )
        sources = [source_map[os.path.abspath(path)] for path in selected]
        uncovered = []
        for source in sources:
            if source.minwn is None or source.maxwn is None:
                continue
            lower, upper = rangebounds(source.origin, min_wn, max_wn, None)
            if source.minwn > lower or source.maxwn < upper:
                uncovered.append((source.path, lower, upper))
        if uncovered:
            raise ValueError(
                'The requested range is not covered by the loaded range cache. '
                'Call px.load again with the wider calculation range.'
            )
        return sources

    def lines(self, min_wn, max_wn, uncfilter='None', threshold='None'):
        """Return a filtered normalized line list for a calculation."""
        if self.linedatabase is None:
            raise ValueError('Loaded data does not contain a normalized line list.')
        if self.linepath is not None:
            frame = readlinecache(self.linepath, min_wn, max_wn)
        else:
            frame = self.lineframe
        return filterlinelist(
            frame,
            self.linedatabase,
            min_wn,
            max_wn,
            uncfilter,
            threshold,
        )


def loaddata(
    config,
    cache='auto',
    cachedir=None,
    maxmemory=512,
    refresh=False,
    alltrans=False,
    preparestates=True,
):
    """Create reusable data for a supported line-list database."""
    if config.database in ('HITRAN', 'HITEMP', 'ExoMolHR'):
        return loadlinedata(config, cache, cachedir, maxmemory, refresh)
    if config.database not in ('ExoMol', 'ExoAtom'):
        raise ValueError('Unsupported database for px.load.')
    if cache not in ('auto', 'parquet', 'none'):
        raise ValueError("cache must be 'auto', 'parquet', or 'none'.")

    if alltrans:
        selected = get_transfiles(
            config.read_path,
            config.data_info,
            prepare=False,
        )
    else:
        selected = get_part_transfiles(
            config.read_path,
            config.data_info,
            config.min_wn,
            config.max_wn,
            prepare=False,
        )
    if not selected:
        raise ValueError('No transition files found for the requested range.')

    if cachedir is None:
        cachedir = os.path.join(
            config.read_path,
            *config.data_info,
            '.pyexocross_cache',
        )

    estimate = sum(
        os.path.getsize(path) * (10 if path.endswith('.bz2') else 2)
        for path in selected
    )
    memorylimit = maxmemory * 1024**2
    resolved = cache
    if cache == 'auto':
        resolved = 'none' if maxmemory <= 0 else ('memory' if estimate <= memorylimit else 'parquet')

    statescache = 'parquet' if resolved == 'parquet' else 'none'
    states = statesdata(config, statescache, cachedir, refresh=refresh)
    prepared = None
    if preparestates:
        prepared = read_part_states(
            states,
            config.unc_filter,
            config.nlte_method,
            config.nlte_path,
            config.check_uncertainty,
            config.check_lifetime,
            config.check_gfactor,
            config.qnslabel_list,
            config.states_col,
            config.states_fmt,
            config.qns_label,
            config.qns_filter,
            config.qns_value,
            config.vib_label,
            config.rot_label,
        )

    if resolved == 'memory':
        memory_sources = []
        actual_memory = 0
        for path in selected:
            source = memorysource(path, config.chunk_size)
            memory_sources.append(source)
            actual_memory += source.size
            if actual_memory > memorylimit:
                break
        if actual_memory <= memorylimit and len(memory_sources) == len(selected):
            sources = tuple(memory_sources)
        else:
            print('Memory estimate exceeded max_memory; switching to Parquet cache.')
            resolved = 'parquet'
            memory_sources = None
            sources = parquetsources(
                config,
                selected,
                states,
                cachedir,
                alltrans,
                refresh=refresh,
            )
    elif resolved == 'parquet':
        sources = parquetsources(
            config,
            selected,
            states,
            cachedir,
            alltrans,
            refresh=refresh,
        )
    else:
        sources = tuple(
            TransSource(
                path=path,
                origin=path,
                size=os.path.getsize(path) * (10 if path.endswith('.bz2') else 2),
            )
            for path in selected
        )

    printcacheinfo(cachedir, loadmanifest(cachedir))
    print('Loaded data cache mode:', cache, f'({resolved})' if cache == 'auto' else '')
    return LoadedData(
        config,
        states,
        prepared,
        sources,
        cache,
        resolved,
        cachedir,
        alltrans,
    )
