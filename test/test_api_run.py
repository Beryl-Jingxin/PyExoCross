import sys
import os

# Ensure the project root is on sys.path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import pyexocross as px

px.run('/home/jingxin/LHD/Program/PyExoCross/.input/MgH_ExoMol.inp')

px.run('/home/jingxin/LHD/Program/PyExoCross/.input/Ar_QDB_nlte.inp')

px.run('/home/jingxin/LHD/Program/PyExoCross/.input/NO_HITRAN.inp')
