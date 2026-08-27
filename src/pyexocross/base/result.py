"""In-memory results and timing for interactive PyExoCross calculations."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class Condition:
    T: float | None = None
    P: float | None = None
    Tvib: float | None = None
    Trot: float | None = None


@dataclass
class Parameters:
    values: dict

    def __getitem__(self, key):
        return self.values[key]

    def to_dict(self):
        return dict(self.values)

    def keys(self):
        return self.values.keys()

    def items(self):
        return self.values.items()

    def get(self, key, default=None):
        return self.values.get(key, default)

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return 'Calculation parameters\n' + '\n'.join(
            f'  {key:<18}: {value}' for key, value in sorted(self.values.items())
            if not key.startswith('_')
        )


@dataclass
class Result:
    kind: str
    coords: dict = field(default_factory=dict)
    units: dict = field(default_factory=dict)
    data: dict = field(default_factory=dict)
    params: Parameters | None = None
    timing: dict = field(default_factory=dict)
    _product: str | None = field(default=None, repr=False)

    @property
    def conditions(self):
        if self._product is None or self._product == self.kind:
            return list(self.data)
        return {
            product: list(values)
            for product, values in self.data.items()
        }

    def select(self, T=None, P=None, Tvib=None, Trot=None, product=None, squeeze=True):
        source = self.data
        if self._product == 'combined':
            if product is None:
                raise ValueError(
                    f'product is required. Available products: {list(self.data)}'
                )
            source = self.data[product]
        wanted = {'T': T, 'P': P, 'Tvib': Tvib, 'Trot': Trot}
        found = {condition: values for condition, values in source.items()
                 if all(value is None or getattr(condition, name) == value
                        for name, value in wanted.items())}
        if not found:
            raise KeyError(f'No result matches {wanted}. Available: {self.conditions}')
        if squeeze and len(found) == 1:
            return next(iter(found.values()))
        if squeeze:
            raise ValueError(f'More than one result matches {wanted}; add conditions or use squeeze=False.')
        return found

    def summary(self):
        print(self)

    def help(self):
        """Print the public fields and selection syntax for this result."""
        print(
            'Result fields:\n'
            '  result.data        calculated arrays grouped by condition\n'
            '  result.coords      coordinate grids\n'
            '  result.units       units for coordinates and data\n'
            '  result.conditions  available T, P, Tvib and Trot combinations\n'
            '  result.params      resolved parameters for this calculation\n'
            '  result.timing      load, prepare, calculate, save and total times\n'
            '  result.select(...) select an array by conditions\n'
            '  result.summary()   print a result summary'
        )
        if self._product == 'combined':
            print(f'\nAvailable products: {list(self.data)}')
            print('Example: result.select(T=300, P=1.0, product="cross_section")')
        else:
            print(f'\nAvailable conditions: {self.conditions}')
            condition = self.conditions[0] if self.conditions else Condition()
            selectors = ', '.join(
                f'{name}={getattr(condition, name)!r}'
                for name in ('T', 'P', 'Tvib', 'Trot')
                if getattr(condition, name) is not None
            )
            print(f'Example: result.select({selectors})')

    def __repr__(self):
        axis = ', '.join(
            f'{name} ({len(values)} condition grids)'
            if isinstance(values, dict)
            else f'{name} ({len(values)} points)'
            for name, values in self.coords.items()
        )
        timing = '\n'.join(
            f'  {name:<18}: {values["system"]:.6f} s'
            for name, values in self.timing.items()
        )
        text = (
            f'{self.kind}Result\n'
            f'  coordinates: {axis or "none"}\n'
            f'  conditions: {self.conditions}'
        )
        if self.params is not None:
            text += f'\n{self.params}'
        if timing:
            text += f'\nTiming\n{timing}'
        return text


COMMON_PARAMETERS = {
    'database', 'molecule', 'isotopologue', 'dataset', 'atom',
    'read_path', 'save_path', 'logs_path', 'output', 'cache', 'cache_dir',
    'max_memory', 'refresh_cache',
}
COMPUTE_PARAMETERS = {
    'ncputrans', 'ncpufiles', 'chunk_size', 'run_mode', 'device',
    'gpu_backend', 'gpu_batch_lines', 'gpu_batch_grid',
}
TEMPERATURE_PARAMETERS = {'ntemp', 'tmax'}
SPECTRUM_PARAMETERS = {
    'T_list', 'wn_wl', 'wn_wl_unit', 'min_wnl', 'max_wnl', 'threshold',
    'unc_filter', 'nlte_method', 'tvib_list', 'trot_list', 'vib_label',
    'rot_label', 'nlte_path', 'abs_emi', 'abundance', 'qns_label',
    'qns_value', 'qns_filter',
}
CROSS_SECTION_PARAMETERS = {
    'P_list', 'N_point', 'bin_size', 'profile', 'cutoff', 'predissoc_yn',
    'broadeners', 'ratios', 'gamma_broadening', 'n_broadening',
    'doppler_hwhm_yn', 'lorentzian_hwhm_yn',
    'alpha_hwhm', 'gamma_hwhm',
}
PARAMETERS_BY_FUNCTION = {
    'partition_functions': TEMPERATURE_PARAMETERS,
    'specific_heats': TEMPERATURE_PARAMETERS,
    'cooling_functions': TEMPERATURE_PARAMETERS | COMPUTE_PARAMETERS,
    'lifetimes': {'compress_yn'} | COMPUTE_PARAMETERS,
    'oscillator_strengths': {
        'gf_or_f', 'wn_wl', 'wn_wl_unit', 'min_wnl', 'max_wnl',
        'unc_filter', 'qns_label', 'qns_value', 'qns_filter',
        'plot_oscillator_strength_yn', 'plot_oscillator_strength_method',
        'plot_oscillator_strength_wn_wl', 'plot_oscillator_strength_unit',
        'limit_yaxis_os',
    } | COMPUTE_PARAMETERS,
    'stick_spectra': SPECTRUM_PARAMETERS | COMPUTE_PARAMETERS | {
        'plot_stick_spectra_yn', 'plot_stick_spectra_method',
        'plot_stick_spectra_wn_wl', 'plot_stick_spectra_unit',
        'limit_yaxis_stick_spectra',
    },
    'cross_sections': SPECTRUM_PARAMETERS | CROSS_SECTION_PARAMETERS
    | COMPUTE_PARAMETERS | {
        'plot_cross_section_yn', 'plot_cross_section_method',
        'plot_cross_section_wn_wl', 'plot_cross_section_unit',
        'limit_yaxis_xsec',
    },
}


def relevant_parameters(config):
    """Return resolved parameters used by the enabled calculation."""
    enabled = [
        name for name in PARAMETERS_BY_FUNCTION
        if getattr(config, name, 0) == 1
    ]
    keys = set(COMMON_PARAMETERS)
    for name in enabled:
        keys.update(PARAMETERS_BY_FUNCTION[name])
    values = vars(config)
    return {
        name: values[name] for name in keys
        if name in values and not (
            name in {'molecule', 'isotopologue', 'dataset', 'atom'}
            and values[name] is None
        )
    }


def record(kind, values, coords, units=None, condition=None):
    """Add an in-memory product when the current run requested one."""
    from pyexocross import core
    if core.active_result is None:
        return
    result = core.active_result
    condition = condition or Condition()
    if result._product == 'combined':
        previous_conditions = list(result.data.get(kind, {}))
    else:
        previous_conditions = list(result.data)
    if result._product is None:
        result.kind = kind
        result._product = kind
    elif result._product != kind and result._product != 'combined':
        previous = result._product
        result.data = {previous: result.data}
        result.kind = 'Combined'
        result._product = 'combined'
    if result._product == 'combined':
        result.data.setdefault(kind, {})[condition] = values.copy()
    else:
        result.data[condition] = values.copy()
    for name, coordinate in coords.items():
        coordinate = np.asarray(coordinate)
        current = result.coords.get(name)
        if current is None:
            result.coords[name] = coordinate.copy()
        elif isinstance(current, dict):
            current[condition] = coordinate.copy()
        elif not np.array_equal(current, coordinate):
            first = previous_conditions[0] if previous_conditions else Condition()
            result.coords[name] = {
                first: current,
                condition: coordinate.copy(),
            }
    if units:
        result.units.update(units)


def saving_enabled():
    from pyexocross import core
    return core.output != 'memory'


def add_timing(name, cpu, system):
    """Accumulate a measured phase in the active result."""
    from pyexocross import core
    if not core.timing_active:
        return
    values = core.timing_summary.setdefault(name, {'cpu': 0.0, 'system': 0.0})
    values['cpu'] += cpu
    values['system'] += system
    if core.active_result is not None:
        core.active_result.timing = core.timing_summary


def timing_value(name):
    from pyexocross import core
    values = core.timing_summary.get(name, {})
    return values.get('cpu', 0.0), values.get('system', 0.0)


def end_calculation(timer, save_start=None):
    """Record elapsed calculation time without nested file-saving time."""
    timer.end()
    if save_start is None:
        save_start = (timer.start_save_CPU, timer.start_save_sys)
    save_end = timing_value('save')
    add_timing(
        'calculate',
        max(0.0, timer.interval_CPU - (save_end[0] - save_start[0])),
        max(0.0, timer.interval_sys - (save_end[1] - save_start[1])),
    )
