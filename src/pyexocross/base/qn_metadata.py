"""Helpers for deriving quantum-number metadata from states definitions."""


def short_states_column_label(label):
    """Return the unqualified state-column label used inside PyExoCross."""
    label = str(label).split(':')[-1]
    normalized = label.lower().replace('-', '').replace('_', '')
    if normalized in ('unc', 'uncertainty'):
        return 'unc'
    if normalized in ('gfactor', 'gfactorlande', 'landeg', 'landegfactor'):
        return 'gfac'
    if normalized in ('tau', 'lifetime'):
        return 'tau'
    return label


def normalized_states_columns(raw_columns, ncols=None):
    """Normalize states metadata column names and avoid duplicate labels."""
    main_columns = ['id', 'E', 'g', 'J']
    columns = list(raw_columns if ncols is None else raw_columns[:ncols])
    names = []
    seen = {}
    for idx, raw_column in enumerate(columns):
        if idx < len(main_columns):
            name = main_columns[idx]
        else:
            name = short_states_column_label(raw_column)
        if name in seen:
            seen[name] += 1
            name = f'{name}.{seen[name]}'
        else:
            seen[name] = 0
        names.append(name)
    return names


def qn_labels_formats_from_states(states_col, states_fmt):
    """Derive QN labels/formats from states metadata."""
    labels = []
    formats = []
    for raw_label, fmt, label in zip(
        list(states_col)[4:],
        list(states_fmt)[4:],
        normalized_states_columns(states_col)[4:],
    ):
        if label in ('unc', 'tau', 'gfac'):
            continue
        if '.' in label:
            continue
        if str(raw_label).startswith('Auxiliary:'):
            continue
        labels.append(label)
        formats.append(fmt)
    return labels, formats


def qn_format_from_states(label, states_col, states_fmt):
    """Return the format for a normalized state QN label."""
    normalized = normalized_states_columns(states_col)
    if label in normalized:
        return list(states_fmt)[normalized.index(label)]
    return None


def qn_formats_for_labels(labels, qnslabel_list, qnsformat_list, states_col=None, states_fmt=None):
    """Resolve formats for requested QN labels."""
    formats = []
    missing = []
    for label in labels:
        fmt = None
        if label in qnslabel_list:
            fmt = qnsformat_list[qnslabel_list.index(label)]
        elif states_col is not None and states_fmt is not None:
            fmt = qn_format_from_states(label, states_col, states_fmt)
        if fmt is None:
            missing.append(label)
        else:
            formats.append(fmt)
    if missing:
        raise ValueError(
            'Requested quantum-number labels are missing from states metadata: '
            + ', '.join(missing)
        )
    return formats
