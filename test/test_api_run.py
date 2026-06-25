import sys
import os

# Ensure src layout is importable in local runs
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

import pyexocross as px

def main():
    px.run('/Users/beryl/Academic/UCL/PhD/Code/PyExoCross/input/MgH_ExoMol_recommended.inp')
    
    px.run('/Users/beryl/Academic/UCL/PhD/Code/PyExoCross/input/NO_ExoMolHR.inp')

    px.run('/Users/beryl/Academic/UCL/PhD/Code/PyExoCross/input/Ar_QDB_nlte.inp')

    px.run('/Users/beryl/Academic/UCL/PhD/Code/PyExoCross/input/NO_HITRAN.inp')


if __name__ == '__main__':
    main()
