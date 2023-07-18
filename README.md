# PyExoCross

See the PyExoCross manual at [https://pyexocross.readthedocs.io](https://pyexocross.readthedocs.io)

## Download

Download PyExoCross program by command:

```bash
git clone https://github.com/Beryl-Jingxin/PyExoCross.git
```

## Run PyExoCross

In the terminal, use the following commands to run PyExoCross:

```bash
python pyexocross.py -p input_filepath
```

If the input filepath is `/home/username/PyExoCross/input/MgH_exomol.inp`

```bash
python pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp
```

If you want to run program in conda environment, please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp
```

If you need to run program in background, please use command:

```bash
nohup python -u pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp > MgH_exomol.log 2>&1 &
```

```bash
nohup /home/username/anaconda3/envs/exomol/bin/python -u pyexocross.py -p /home/username/PyExoCross/input/MgH_exomol.inp > MgH_exomol.log 2>&1 &
```

All information could be written in the input file. Just change the information you will use, please do not change the other unnecessary information.
