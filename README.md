#Binary Lennard-Jones parallel Molecular Dynamics Code

The original code is taken from Prof. Aiichiro Nakano's class material; CS653 High-Performance Computing and Visualization. 
http://cacs.usc.edu/education/cs653-code.html

In addition to the original parallel L-J potential MD code, mean square displacement, velocity autocorrelation function, Static structure factor, phonon density of state, pair distribution function and coordination function have been added, along with several utility tools to setup computing environment at the HPC cluster and a iPython Notebook script to plot the aformentioned properties into 2D graph.

## To clean working directory, setup enviroment, compile the code, and run the executable. 
> make clean

> source ./util/setup.sh

> make

> mpirun -np 4 ./pmd (4 MPIrank case)


