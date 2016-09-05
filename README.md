#Binary Lennard-Jones parallel Molecular Dynamics

This is yet another parallel molecular dynamics simulation code to compute mean square displacement, velocity autocorrelation function, Static structure factor, phonon density of state, pair distribution function and coordination function during simulation. It also comes with several utility tools including a shell script to setup computing environment at the HPC cluster and an iPython Notebook to plot the aformentioned properties into 2D graph.

The original parallel L-J MD code is taken from Prof. Aiichiro Nakano's class material. See his class web sites for more information. 

CSCI 596: SCIENTIFIC COMPUTING & VISUALIZATION

http://cacs.usc.edu/education/cs596.html

CS 653: High-Performance Computing and Visualization. 

http://cacs.usc.edu/education/cs653-code.html

This code also supports binary LJ interatomic potential. The used parameters are taken from a glass forming material study by Stillinger et al. 
Signatures of distinct dynamical regimes in the energy landscape of a glass-forming liquid.
http://www.nature.com/nature/journal/v393/n6685/abs/393554a0.html

## Basic steps to perform MD simulation

Move to working directory:
> cd binary-LJ-pmd/

clean up working directory:
> make clean

Setup enviroment:
> source ./util/setup.sh

Compile:
> make

Run:
> mpirun -np 4 ./pmd (4 MPIrank case)

## Molecualr Dynamics Parameters
All necessary parameters to perform MD simulation are stored in pmd.in file. 'pmd.in' file is read from MPI processes at the begging of each run. Typical 'pmd.in' file would look like below.

```
1                   // mdmode: 0=initial run, 1=normal run, 4=T control
2 2 1               // vproc[0],vproc[1],vproc[2]
8 8 18              // InitUcell[0],InitUcell[1],InitUcell[2]
1.0                 // Density
0.5                 // InitTemp
0.005               // DeltaT
2000                // StepLimit
10                  // StepAvg
```

## Solid vs Liquid phase
Here are some system properties at different temperatures. 

The system is thermalized at T = 0.8. All properties indicate that the system is solid. 

<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp0.8-Solid.png" width="400">

The system is thermalized at T = 1.2 where the plots show liquid like behavior. 

<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp1.6-Melt.png" width="400">

