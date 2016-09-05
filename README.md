#Binary Lennard-Jones parallel Molecular Dynamics Code

The original code is taken from Prof. Aiichiro Nakano's class material; CS653 High-Performance Computing and Visualization. 
http://cacs.usc.edu/education/cs653-code.html

In addition to the original parallel L-J potential MD code, mean square displacement, velocity autocorrelation function, Static structure factor, phonon density of state, pair distribution function and coordination function have been added, along with several utility tools to setup computing environment at the HPC cluster and a iPython Notebook script to plot the aformentioned properties into 2D graph.

This code also supports binary LJ interatomic potential whose parameters are used in a glass forming material study by Stillinger et al. 
Signatures of distinct dynamical regimes in the energy landscape of a glass-forming liquid.
http://www.nature.com/nature/journal/v393/n6685/abs/393554a0.html


## Basics
clean working directory 
> make clean

Setup enviroment
> source ./util/setup.sh

Compile
> make

Run
> mpirun -np 4 ./pmd (4 MPIrank case)

## Molecualr Dynamics Parameters. 
All necessary parameters to run this code is stored in pmd.in file, which is read from MPI processes at the begging of each run. pmd.in would look like, 

```
1
2 2 1
8 8 18
1.0
0.5
0.005
2000
10
--------------------------------------------
mdmode: 0=initial run, 1=normal run, 4=T control, 5=MSD 
vproc[0],vproc[1],vproc[2]
InitUcell[0],InitUcell[1],InitUcell[2]
Density
InitTemp
DeltaT
StepLimit
StepAvg
```

## Solid vs Liquid phase
Here are some system properties at different temperatures. 
The system is thermalized at T = 0.8. All properties indicate that the system is solid. 
<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp0.8-Solid.png" width="400">

The system is thermalized at T = 1.2 where the plots show liquid like behavior. 
<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp1.6-Melt.png" width="400">

