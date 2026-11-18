# Getting PyExoCross Program

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
python run.py -p input_filepath
```

*Example*

If the input filepath is `/home/username/PyExoCross/input/H2O_exomol.inp`.

```bash
python run.py -p ./input/H2O_exomol.inp
# OR 
python run.py -p /home/username/PyExoCross/input/H2O_exomol.inp
```

If you want to run program in conda environment which is named as 'exomol', please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python run.py -p ./input/H2O_exomol.inp
```

If you need to run program in background, please use command:

```bash
# Discard redirected terminal output. LogFilePath still controls log-file output.
nohup python -u run.py -p ./input/H2O_exomol.inp > /dev/null 2>&1
# Save redirected terminal output separately
nohup python -u run.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
# OR 
nohup /home/username/anaconda3/envs/exomol/bin/python -u run.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
```

To disable PyExoCross log-file output, set this in the input file:

```text
LogFilePath                             None
```

This is independent of terminal output. Use `Verbose False` to hide normal
terminal output and progress bars.

## Notes for input file

All information can be written in the input file. Just change the information you will use.You don't need to change any other unnecessary information.Please do not change the first column strings.
