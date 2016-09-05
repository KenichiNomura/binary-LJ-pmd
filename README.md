#Binary Lennard-Jones parallel Molecular Dynamics Code

The original code is taken from Prof. Aiichiro Nakano's class material; CS653 High-Performance Computing and Visualization. 
http://cacs.usc.edu/education/cs653-code.html

In addition to the original parallel L-J potential MD code, mean square displacement, velocity autocorrelation function, Static structure factor, phonon density of state, pair distribution function and coordination function have been added, along with several utility tools to setup computing environment at the HPC cluster and a iPython Notebook script to plot the aformentioned properties into 2D graph.

## To clean working directory, setup enviroment, compile the code, and run the executable. 
> make clean

> source ./util/setup.sh

> make

> mpirun -np 4 ./pmd (4 MPIrank case)

## Here are some system properties at different temperatures. 
This system is thermalized at T = 0.8. All properties indicate that the system is solid. 
<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp0.8-Solid.png" width="400">

This system is thermalized at T = 1.2 where the plots show liquid like behavior. 
<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp1.6-Melt.png" width="400">

