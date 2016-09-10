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

## Steps to perform MD simulation
### Download Source Code
If you are not faimilier with github, go this page first **https://github.com/KenichiNomura/binary-LJ-pmd**. 
You can download the source code to your compure by clicking **clone or download** button, then select **Download Zip** tab. You will have a file called **binary-LJ-pmd-master.zip** Aftre unzipped it, change directory to **binary-LJ-pmd**. 
```
cd binary-LJ-pmd/
```
You will be working in this directory most of time. 

### Setup Your Environment
If you are using HPCC cluster, it is exteremely important to keep in mind that you **CAN NOT** run any program on the file servers. The file server is a shared resource among all students, staff and researchers and serving disk access to thousands of computing nodes. If your processe occupies CPU, memory, network traffic etc on the file server, it would impact to everyone using the HPCC cluster. Instead of using the shared machine you are allowed to grab some compute nodes from the cluster and work (compile/debug/run) there by using  **PBS interactive mode**.
```
qsub -I -d .
```

Next step is set up your environment so that you can use compiler and Message Passing Interface (MPI) library to generate exectuable file. You need to **source** a shell script under **util** directory. 
```
source ./util/setup.sh
```
Please note that this **setup.sh** script works for HPCC environment only. 

Now you are ready to compile the code. To do so simply type **make** like below. 
```
make
```
You should see a file called **pmd** (* mark after the filename indicates this is an executable file)
```
$ ls
Makefile   data/      pmd*       pmd.h      utils/
README.md  docs/      pmd.c      pmd.in
```

Now you should be able to perform simulation. A standard way to execute a command is like typing 
```
./a.out
```
but this is a *parallel* MD code so you have to call a MPI command, **mpirun**, to let the program and OS know how many MPI processes you will create. To do so, you give **-np** option and a interger before **pmd**. Example below creates 2 MPI processes whose ranks are 0 and 1. If the system you are on has enouch resources, e.g. 2 CPU cores, your program will run at fastest speed. 

```
mpirun -np 2 ./pmd
```

## Input Parameters
Okay, hopefully you are able to run the code successfully now. Next step is to understand input parameters so that you have full control of your simulations. All input parameters in **pmd.in** file, which typically looks like,

```
1                   // mdmode: 0=initial run, 1=normal run, 4=T control
2 1 1               // vproc[0],vproc[1],vproc[2]
2 4 4              // InitUcell[0],InitUcell[1],InitUcell[2]
1.0                 // Density
0.5                 // InitTemp
0.005               // DeltaT
2000                // StepLimit
10                  // StepAvg
```

## Analysis

### 1. Mean Square Displacement (MSD)
output file: msd.d

### 2. Velocity Autocorrelation Function (VAF)
output file: vac.d

### 3. Pair Distribution Function
output file: gr.d

### 4. Coordination Number
output file: gr.d

### 5. Static Strcture Factor
output file: sq.d

### 6. Phonon Density of State
output file: ft.d

## Solid vs Liquid phase
Here are some results of systems at different temperatures. 

1) Thermalized at T = 0.8. The flat MSD and oscillating VAF with rather long correlation indiate that the system is solid. 

<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp0.8-Solid.png" width="400">

2) Thermalized at T = 1.2. MSD increases linearly wrt time, VAF dies out quickly and PDF stays around 1 after the first peak. These are typical behaviors of liqud.  

<img src="https://github.com/KenichiNomura/binary-LJ-pmd/blob/master/docs/Temp1.6-Melt.png" width="400">

