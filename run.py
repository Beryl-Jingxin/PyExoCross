from pyexocross.config import Config
from pyexocross.core import get_results
from src.base.log import parse_logging_info, setup_logging
from src.base.input import parse_args

# inp_path = "/home/jingxin/LHD/Program/PyExoCross/.input/MgH_ExoMol.inp"
inp_path = parse_args()

setup_logging(parse_logging_info(inp_path))
cfg = Config(inp_filepath=inp_path)
get_results(cfg)
