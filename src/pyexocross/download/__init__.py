"""Download helpers for external molecular databases."""

from .download_exoatom import download_exoatom, get_exoatom_urls
from .download_exomol import download_exomol, get_exomol_urls
from .download_exomolhr import download_exomolhr, get_exomolhr_urls
from .download_hitran import download_hitran, get_hitran_targets


def download(
    database,
    file_path=None,
    species_info=None,
    molecules=None,
    **kwargs,
):
    """
    Download database files using a database-specific downloader.

    Parameters
    ----------
    database : str
        One of ``ExoMol``, ``ExoAtom``, ``ExoMolHR``, or ``HITRAN``.
    file_path : str
        Folder where downloaded files are saved.
    species_info : dict, optional
        Molecule/isotopologue or atom configuration.
    **kwargs
        Additional downloader-specific options.
    """
    database_name = database.lower()
    if database_name == 'exomol':
        if species_info is None:
            raise ValueError('species_info must be provided for ExoMol downloads.')
        return download_exomol(
            species_info=species_info,
            file_path=file_path,
            **kwargs,
        )
    if database_name == 'exoatom':
        if species_info is None:
            raise ValueError('species_info must be provided for ExoAtom downloads.')
        return download_exoatom(
            species_info=species_info,
            file_path=file_path,
            **kwargs,
        )
    if database_name == 'exomolhr':
        if species_info is None:
            raise ValueError('species_info must be provided for ExoMolHR downloads.')
        return download_exomolhr(
            species_info=species_info,
            file_path=file_path,
            **kwargs,
        )
    if database_name == 'hitran':
        if species_info is None:
            raise ValueError('species_info must be provided for HITRAN downloads.')
        return download_hitran(
            species_info=species_info,
            file_path=file_path,
            **kwargs,
        )

    raise ValueError(
        f"Unsupported download database: {database}. "
        "Supported databases are ExoMol, ExoAtom, ExoMolHR, and HITRAN."
    )


__all__ = [
    'download',
    'download_exomol',
    'get_exomol_urls',
    'download_exoatom',
    'get_exoatom_urls',
    'download_exomolhr',
    'get_exomolhr_urls',
    'download_hitran',
    'get_hitran_targets',
]

