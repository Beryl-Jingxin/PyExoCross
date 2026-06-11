"""Shared helpers for database download modules."""

import os
import requests
import subprocess
import zipfile
from io import BytesIO


def resolve_save_path(save_path=None, file_path=None):
    """Resolve ``save_path`` with ``file_path`` as a public API alias."""
    if save_path is None:
        save_path = file_path
    if save_path is None:
        raise ValueError('save_path or file_path must be provided.')
    return save_path


def write_urls(urls, url_path):
    """Write a URL list to disk."""
    url_dir = os.path.dirname(url_path)
    if url_dir:
        os.makedirs(url_dir, exist_ok=True)
    with open(url_path, 'w', encoding='utf-8') as file_obj:
        for url in urls:
            file_obj.write(f'{url}\n')


def default_url_path(save_path, database_name):
    """Return the default URL-list path for a database."""
    return os.path.join(save_path, 'url', f'{database_name.lower()}__urls.txt')


def download_url_file(url_path, save_path):
    """Download files listed in ``url_path`` using wget."""
    command = ['wget', '-r', '-nH', '--cut-dirs=1', '-P', save_path, '-i', url_path]
    subprocess.run(command, check=True)


def _looks_like_html(content):
    prefix = content[:512].lstrip().lower()
    return prefix.startswith(b'<!doctype html') or prefix.startswith(b'<html')


def _validate_download_content(url, target_path, content, expected=None):
    if len(content) == 0:
        raise ValueError(f'Downloaded empty file from {url}')
    if _looks_like_html(content):
        raise ValueError(
            f'Downloaded HTML instead of data from {url}. '
            'The URL may be a landing/login/results page rather than a direct data file.'
        )
    if expected == 'bz2' and not content.startswith(b'BZh'):
        raise ValueError(f'Downloaded file is not a bzip2 file: {url}')
    if expected == 'text' and b'\x00' in content[:1024]:
        raise ValueError(f'Downloaded file is binary, expected text: {url}')
    if expected == 'csv':
        first_line = content.splitlines()[0] if content.splitlines() else b''
        if b',' not in first_line:
            raise ValueError(f'Downloaded file does not look like CSV: {url}')


def download_url_to_path_with_wget(url, target_path, timeout=60, expected=None):
    """Download one URL to an exact local filepath using wget."""
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    command = ['wget', f'--timeout={timeout}', f'--tries=1', '-O', target_path, url]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        if os.path.exists(target_path):
            os.remove(target_path)
        raise RuntimeError(f'wget failed for {url}') from exc
    with open(target_path, 'rb') as file_obj:
        sample = file_obj.read(2048)
    try:
        _validate_download_content(url, target_path, sample, expected=expected)
    except ValueError:
        if os.path.exists(target_path):
            os.remove(target_path)
        raise
    print(f'Downloaded {url}')
    print(f'Saved to {target_path} ({os.path.getsize(target_path)} bytes)')


def download_url_to_path(url, target_path, timeout=60, expected=None):
    """Download one URL to an exact local filepath and reject HTML/empty files."""
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    content = response.content
    _validate_download_content(url, target_path, content, expected=expected)
    with open(target_path, 'wb') as file_obj:
        file_obj.write(content)
    print(f'Downloaded {url}')
    print(f'Saved to {target_path} ({len(content)} bytes)')


def download_zip_csv_to_path(url, target_path, timeout=60):
    """Download a ZIP archive, extract the first CSV member, and save it."""
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    content = response.content
    _validate_download_content(url, target_path, content)
    if not content.startswith(b'PK'):
        raise ValueError(f'Downloaded file is not a ZIP archive: {url}')
    with zipfile.ZipFile(BytesIO(content)) as archive:
        csv_members = [name for name in archive.namelist() if name.lower().endswith('.csv')]
        if len(csv_members) == 0:
            raise ValueError(f'No CSV file found inside ExoMolHR ZIP archive: {url}')
        csv_content = archive.read(csv_members[0])
    _validate_download_content(url, target_path, csv_content, expected='csv')
    with open(target_path, 'wb') as file_obj:
        file_obj.write(csv_content)
    print(f'Downloaded {url}')
    print(f'Extracted {csv_members[0]}')
    print(f'Saved to {target_path} ({len(csv_content)} bytes)')
