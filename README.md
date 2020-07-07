# ModIsing

Simulation of phase transition in modified Ising model with anisotropic force field.

## Prerequisites

- GFortran 95
- Make
- Bash

## Instruction

The parameters of the anisotropic force field can be modified in driver/const.f. To create executables, navigate to folder driver/ and run with
```
make all
```
The simulation program can be run with
```bash
./isingModel.exe $lSze $temp $stpNum $calStp
```
where lSze and temp are the size and temperature of the system, stpNum and calStp are the total number of steps and the number of steps per saving of the simulation. The result then needed to be processed by 2 different ways
```bash
./isingAnal.exe $lSze $temp $stpNum
```
which analyzes the physical properties that depend only on one "snapshot", while
```bash
./dataAnal.exe $lSze $temp $stpNum
```
evalutes the time-series variables, meaning the analysis is done over multiple "snapshots". Because the experiment must be done over different values of size and temperature, bash scripts are provided in folder driver/loop/ to faciliate job submission. Finally, run statPlot.py to see the graphed results.
