"""
Test PyExoCross API function: downloading databases

"""
import sys
import os

# Ensure src layout is importable in local runs
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

import pyexocross as px


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------
def test_download_exomol():
    """Test download ExoMol files."""
    print('\n' + '='*70)
    print('TEST: px.download()')
    print('='*70)
    px.download(    
        file_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/',    # Write that file_path or save_path are the same in download
        database='ExoMol',
        molecule_isotopologues={
            'MgH': {
                '24Mg-1H': {'wn_range': None},
                '25Mg-1H': {'wn_range': None},
            },
            'H2O': {
                '1H2-16O': {'wn_range': [41000, 41200]},
            },
        },
        download=True,
    )
    print('PASSED: download()')
    
    
def test_download_exomolhr():
    """Test download ExoMolHR files."""
    print('\n' + '='*70)
    print('TEST: px.download()')
    print('='*70)
    px.download(    
        file_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMolHR/',    # Write that file_path or save_path are the same in download
        database='ExoMolHR',
        molecule_isotopologues={
            'MgH': {
                '24Mg-1H': None,
                '25Mg-1H': {'T': 1000, 'wn_range': [0, 500], 'threshold': 1e-30},
            },
            'AlH': {
                '27Al-1H': {'T': 500, 'wn_range': [0, 500], 'threshold': 1e-30}
            },
        },
        download=True,
    )
    print('PASSED: download()')
    
    
def test_download_exoatom():
    """Test download ExoAtom files."""
    print('\n' + '='*70)
    print('TEST: px.download()')
    print('='*70)
    px.download(    
        file_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoAtom/',    # Write that file_path or save_path are the same in download
        database='ExoAtom',
        atom_datasets={
            'He': {'dataset': 'NIST'},
            'He_p': {'3He_p': {'dataset': 'NIST'}},
            'Ar': {'dataset': ['NIST', 'Kurucz']},
        },
        download=True,
    )
    print('PASSED: download()')
    
    
def test_download_hitran():
    """Test download HITRAN files."""
    print('\n' + '='*70)
    print('TEST: px.download()')
    print('='*70)
    px.download(    
        file_path='/Users/beryl/Academic/UCL/PhD/Data/database/HITRAN/',    # Write that file_path or save_path are the same in download
        database='HITRAN',
        molecule_isotopologues={
            'NO': {
                '14N-16O': {'wn_range': [0,100]},
                '15N-16O': {'wn_range': [100,150]},
            },
            'H2O': {
                '1H2-16O': {'wn_range': [100, 110]},
            },
        },
        download=True,
    )
    print('PASSED: download()')    
 

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('PyExoCross API Test — Download Databases')
    print(f'pyexocross version: {px.__version__}')

    tests = [
        test_download_exomol,
        test_download_exoatom,
        test_download_exomolhr,
        test_download_hitran,
    ]

    passed = 0
    failed = 0
    for test_fn in tests:
        try:
            test_fn()
            passed += 1
        except Exception as exc:
            failed += 1
            print(f'FAILED: {test_fn.__name__}  —  {exc}')
            import traceback
            traceback.print_exc()

    print('\n' + '='*70)
    print(f'ExoMol API test results: {passed} passed, {failed} failed out of {len(tests)}')
    print('='*70)
    sys.exit(1 if failed else 0)
