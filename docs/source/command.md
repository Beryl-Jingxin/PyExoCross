# Getting PyExoCross

## Download PyExoCross

Download PyExoCross from [GitHub](https://github.com/Beryl-Jingxin/PyExoCross.git "GitHub").

```
git clone https://github.com/Beryl-Jingxin/PyExoCross.git
```

## Run PyExoCross

Prepare an input file *filename.inp* (see examples in the 'input' folder on GitHub) and run the program with command:

```bash
python pyexocross.py -p input_filepath
```

*Example*

If the input filepath is `/home/username/PyExoCross/input/MgH_exomol.inp`.

```bash
python pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp
```

If you want to run program in conda environment, please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp
```

## Notes for input file

1. All information can be written in the input file. Just change the information you will use, please do not change the first column strings or delete the other unnecessary information.
2. If you met problems, jupyter notebook `.ipynb` code is stored for checking and testing.
