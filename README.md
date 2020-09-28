# Prototype_Analysis
An analysis program for JUNO prototype data.
## Compile Need
CERN ROOT
Cable map in current directory
## Compile Command
```bash
make
```

And clean compile

```bash
make clean
```

## Command
input arguments:

1. "./getGain"
2. Input data type (prtJUNO_externalTrigger or prtJUNO_physics).
3. Run Number.
4. The number of files.
5. Output txt file name (for test).
6. Output root file name.
7. Input calibration gain file.

### An example
./getGain prtJUNO_physics 3211 5 PENum.txt Calib_3211.root gain_3211.txt

