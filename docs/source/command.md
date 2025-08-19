# Getting PyExoCross

## Download PyExoCross

Download PyExoCross from [GitHub](https://github.com/Beryl-Jingxin/PyExoCross.git "GitHub").

```bash
git clone https://github.com/Beryl-Jingxin/PyExoCross.git
```

## Install Python packages

```bash
pip install -r requirements.txt
```

## Run PyExoCross

Prepare an input file *filename.inp* (see examples in the 'input' folder on GitHub) and run the program with command:

```bash
python pyexocross.py -p input_filepath
# OR
python3 pyexocross.py -p input_filepath
```

*Example*

If the input filepath is `/home/username/PyExoCross/input/H2O_exomol.inp`.

```bash
python pyexocross.py -p ./input/H2O_exomol.inp
# OR 
python pyexocross.py -p /home/username/PyExoCross/input/H2O_exomol.inp
```

If you want to run program in conda environment which is named as 'exomol', please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python pyexocross.py -p ./input/H2O_exomol.inp
```

If you need to run program in background, please use command:

```bash
nohup python -u pyexocross.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
# OR 
nohup /home/username/anaconda3/envs/exomol/bin/python -u pyexocross.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
```

## Notes for input file

1. All information can be written in the input file. Just change the information you will use.You don't need to change any other unnecessary information.Please do not change the first column strings.
2. If you met problems, jupyter notebook `.ipynb` code is stored for checking and testing.
