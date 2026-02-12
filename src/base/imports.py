## Import all what we need.
# encoding: utf-8
import os
import re
import bz2
import sys
import glob
import time
import atexit
import requests
import argparse
import datetime
import subprocess
import numpy as np
import pandas as pd
import numexpr as ne
import dask.array as da
import astropy.units as au
import dask.dataframe as dd
import astropy.constants as ac
import matplotlib.pyplot as plt
from tqdm import tqdm
from io import StringIO
from itertools import chain
from tabulate import tabulate
from pandarallel import pandarallel
from matplotlib.collections import LineCollection
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED, as_completed
from functools import partial
from multiprocessing import Pool, freeze_support
from scipy.special import voigt_profile, wofz, erf, roots_hermite

import warnings
warnings.simplefilter("ignore", np.ComplexWarning)
pd.options.mode.chained_assignment = None

import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

import multiprocessing as mp
freeze_support()
