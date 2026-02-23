import os
import sys
import multiprocessing as mp
import concurrent.futures as cf

project_root = os.path.dirname(os.path.abspath(__file__))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

try:
    mp.set_start_method('fork')
except (RuntimeError, ValueError):
    pass

# macOS stability: avoid nested process-pool deadlocks in long workflows.
if sys.platform == 'darwin':
    cf.ProcessPoolExecutor = cf.ThreadPoolExecutor

from pyexocross.config import Config
from pyexocross.core import get_results
from pyexocross.base.log import parse_logging_info, setup_logging
from pyexocross.base.input import parse_args

if __name__ == '__main__':
    # inp_path = "/home/jingxin/LHD/Program/PyExoCross/.input/MgH_ExoMol.inp"
    inp_path = parse_args()

    setup_logging(parse_logging_info(inp_path))
    cfg = Config(inp_filepath=inp_path)
    get_results(cfg)
