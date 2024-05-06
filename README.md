# PyExoCross

See the PyExoCross ***manual*** at [https://pyexocross.readthedocs.io](https://pyexocross.readthedocs.io).

See the PyExoCross &nbsp; ***paper***&nbsp; at [https://doi.org/10.1093/rasti/rzae016](https://doi.org/10.1093/rasti/rzae016).


## Download

Download PyExoCross program by command:

```bash
git clone https://github.com/Beryl-Jingxin/PyExoCross.git
```

## Install Python packages

```bash
pip install -r requirements.txt
```

Python packages version

| Python packages     | Version  |
| :------------------ | -------- |
| argparse            | 1.1      |
| astropy             | 6.0.0    |
| dask                | 2024.1.0 |
| indexed_bzip2       | 1.5.0    |
| ipython             | 8.12.3   |
| matplotlib          | 3.8.2    |
| numexpr             | 2.8.8    |
| numpy               | 1.22.3   |
| pandarallel         | 1.6.5    |
| pandas              | 2.0.3    |
| python_version      | >"3.8"   |
| requests            | 2.31.0   |
| scipy               | 1.11.4   |
| tqdm                | 4.66.1   |
| urllib3             | 1.26.13  |
| version-information | 1.0.4    |
| Watermark           | 2.4.3    |

## Run PyExoCross

In the terminal, use the following commands to run PyExoCross:

```bash
python3 pyexocross.py -p input_filepath
```

If the input filepath is `/home/username/PyExoCross/input/H2O_exomol.inp`

```bash
python3 pyexocross.py -p ./input/H2O_exomol.inp
# OR 
python3 pyexocross.py -p /home/username/PyExoCross/input/H2O_exomol.inp
```

If you want to run program in conda environment which is named as 'exomol', please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python pyexocross.py -p ./input/H2O_exomol.inp
```

If you need to run program in background, please use command:

```bash
nohup python3 -u pyexocross.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
# OR 
nohup /home/username/anaconda3/envs/exomol/bin/python -u pyexocross.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
```

## Notes for input file

1. All information can be written in the input file. Just change the information you will use.You don't need to change any other unnecessary information.Please do not change the first column strings.
2. If you met problems, jupyter notebook `.ipynb` code is stored for checking and testing.
