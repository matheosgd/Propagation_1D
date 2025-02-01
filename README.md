# Propagation_1D :  
*Developed by*:  
Mathéo Segaud $^{\dagger}$  
*Under the supervision of:*  
Research Professor : David Lauvergnat $^{\dagger,1,2}$  
Ph.D. Student : Rabiou Issa $^{1}$  

${^{\dagger}}$ Université Paris-Saclay, 3 rue Joliot Curie, Bâtiment Breguet, 91190 Gif-sur-Yvette  
${^1}$ Institut de Chimie Physique, 310 Rue Michel Magat, 91400 Orsay  
${^2}$ CNRS  

## Contents
- [Quick introduction and Objectives](#Quick-introduction-and-Objectives) 
- [Setting up environment](#Setting-up-environment)
- [Propagation_1D Installation](#Propagation_1D-Installation)
- [General view and Structure](#General-view-and-Structure)
    - [```APP```](#APP)
    - [```DATA```](#DATA)
    - [```Ext_Lib```](#Ext_Lib)
    - [```GP```](#GP)
    - [```OBJ```](#OBJ)
    - [```OUT```](#OUT)
    - [```SRC```](#SRC)
        - [```P1D_basis_m.f90```](#P1D_basis_m.f90)
        - [```P1D_lanczos_m.f90```](#P1D_lanczos_m.f90)
        - [```P1D_operators_m.f90```](#P1D_operators_m.f90)
        - [```P1D_perturbations_m.f90```](#P1D_perturbations_m.f90)
        - [```P1D_propagation_m.f90```](#P1D_propagation_m.f90)
        - [```P1D_wavepacket_m.f90```](#P1D_wavepacket_m.f90)
    - [```TESTS```](#TESTS)
- [Compiling the library](#Compiling-the-library)
- [Running the tests](#Running-the-tests)
- [Running the application](#Running-the-application)
- [Authors](#authors)
- [Reference](#reference)


# Quick Introduction and Objectives
The code allows to determine eigenvalues and eigenvectors of a systeme then to propagate a wavepacket on a potential in a 'BoxAB' simulation grid : it simulate the propagation of a particule in box of bounds A and B.  

Four types of potential are possible to choose in the input files ("data_*"). The parameters of the wavepacket and the choosed potential can be change in the same input file. Then one can execute the code via the script "compilation.sh".  

The eigenvalues and the firts eigenvectors are printed in the file "results" as well as the initial energy, position, and impulsion of the wavepacket. The values are calculated in atomic units but energies are also printed in cm-1 and eV.  

Four propagators are available and can be choosed in the same input file for each propagation :  
  The 'Spectral' propagator which use the evolution operator.  
  The 'Euler' propagator which use the Euler method for solve differential equations (1st ordre Taylor developpement).  
  The 'Taylor_4' propagator which use a 4th order Taylor developpement.  
  And the Short Iterative Lanczos ('SIL') propagator which use the Krylov vector space and the Lanczos diagonalisation method.  

The evolution over time of the wavepacket properties (energy, position, impulsion, autocorrelation function) can be monitored in the "prop" file, then the position over time can be traced via the script gnuplot "trace_positions.gp".  

The final vector of the wavepacket on the simulation grid is printed in the "psif" file and can be traced via gnuplot whith the "trace_psif.gp" script.  

More detailed information about the model will be found in the futur Manual "Propagation_1D_Manual.pdf".  


# Setting up environment
This library is designed to work within a GNU/Linux operating system and the compilation and execution files (makefile and run.sh) were written to use the gfortran compiler for Fortran90. It uses the QDUtilLib Fortran library [[1]](#reference) designed by David Lauvergnat (cf. Authors).  


# Propagation_1D installation
One can either install the library by downloading its files (if they are not interested into having a Git control of their version) or by Cloning it (otherwise) from its GitHub repository.  


## Downloading the files
    1. Choose directory on your computer where the Propagation_1D directory is to be created.  
    2. Open the GitHub repository of Propagation_1D and click on "<> Code" on top right of the page > Download zip  
    3. Place it into the chosed directory and extract the files  


## Cloning the files
    1. Choose directory on your computer where the Propagation_1D directory is to be created.  
    2. Open the GitHub repository of Propagation_1D and click on "<> Code" at the top right of the project page  
    3. In the menue that appeared, select HTTPS to make the URL appear and copy it.  
    4. Open the terminal of the computer and ```bash cd``` to the chosed directory  
    5. execute the command : ```bash git clone <the copied link>```   


# General view and Structure
The Propagation_1D library is structured as follows :  

```.
Propagation_1D
├── APP
│   ├── App_Propagation_1D.f90
│   └── Main.f90
│
├── DATA
│   ├── data_app.nml
│   ├── data_db_pts.nml
│   ├── data_escalier.nml
│   ├── data_harmo.nml
│   ├── data_morse.nml
│   └── data_tests.nml
│
├── Ext_Lib
│   ├── QDUtilLib
│   ├── QDUtilLib_loc
│   ├── QDUtilLib-main
│   ├── .gitignore
│   ├── cleanlib
│   └── get_Lib.sh
│
├── GP
│   ├── trace_gif.gp
│   ├── trace_harmo.gp
│   ├── trace_morse.gp
│   ├── trace_positions.gp
│   ├── trace_psif.gp
│   └── trace.gp
│
├── OBJ
│   └── obj
│       ├── *.o
│       └── *.mod
│    
├── OUT
│   └── *.log
│
├── RESULTS
│   └── output.gif
│
├── SRC
│   ├── P1D_basis_m.f90
│   ├── P1D_lanczos_m.f90
│   ├── P1D_operators_m.f90
│   ├── P1D_perturbations_m.f90
│   ├── P1D_propagation_m.f90
│   └── P1D_wavepacket_m.f90
│
├── TESTS
│
├── .gitignore
├── README.md
└── run.sh
```  

More details will be found about each of these files and directories in the futur manual. Here will be only a general presentation.  


## APP
This directory contains the source ```.f90``` file of the application program. It is written to provide an illustration of what the library can realise and how to set-up the simulations, as well as an additional test.  


## DATA
This directory contains the data files used to feed the library with the boxAB and propagation parameters. They are written as Fortran namelists. Different one are provided to be used as convenient. Please do not modify the two dedicated respectively to the execution of the tests and the application (data_tests.nml and data_app.nml). The others can be modified and used by the user for specific uses.  


## Ext_Lib
This directory contains the external libraries used by Propagation_1D. Up to now, the QDUtilLib [[1]](#reference) library designed by David Lauvergnat.  


## GP
This directory contains several gnuplot scripts that can be used to plot the data from the output files, or as inspiration to build the user's own scripts.  


## OBJ
This directory contains all the object ```.o``` and ```.mod``` files generated as the library is compiled.  


## OUT
This directory contains all the output ```.log``` files generated as the tests or the application are executed.  


## SRC
This directory contains all the source ```.f90``` files of the modules that composes the library.  


### P1D_basis_m.f90
This module contains the procedures to build and describe the boxAB and the basis set.  


### P1D_lanczos_m.f90
This module contains the procedures needed for the propagation using the Short Iterative Lanczos method.  


### P1D_operators_m.f90
This module contains the procedures to build the quantum mechanics operators relative to the system, and to compute their action on its vector.  


### P1D_perturbations_m.f90
This module will contain, in future versions of this library, the procedures to compute more realistic system by perturbating the initial one with quantum perturbation theory.  


### P1D_propagation_m.f90
This module contains the procedures to propagate the system along time.  


### P1D_wavepacket_m.f90
This module contains the procedures to build and describe the wavepacket representing the system.  


## TESTS
This directory will contain all the source ```.f90``` files of the programs designed to test automatically that the library is working as expected. There is none for now.  


# Compiling the library
The library can be built using the run.sh shell script. It supports the following commands to manage the library :  

    - "getlib" : install the external library if not already downloaded.  
    - "lib" : Compiles the external library's modules and build the .a static library file if needed, then does the same for the Propagation_1D library.  
    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and link them with the static library file into executables.  
    - "clean" : Deletes all the object, executable and output files from Propagation_1D (not those of the external libraries).  
    - "cleanall" : Performs the same cleaning as "clean" in addition also to delete the .mod and static library files, and run the cleanlib shell script from the external QDUtilLib library that takes the same actions upon its files.  

They can be executed by commad-line :  
```bash
make <name of the command>
./run.sh <name of the commands> (up to four arguments)
```


# Running the tests
The tests can be compiled and executed using the run.sh shell script. The shell script allows to choose to compile and execute only one of the tests and to choose the data file to read if one is not interested in the default one. It uses these commands :  

    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and applications and link them with the static library file into executables.  
    - "ut" : Performs the same actions as "all" in addition also to execute the tests, using the indicated data file and direct the output into the corresponding .log output file. Then it grabs the final sentence of these files for each test and display them on the screen. It is the default command if no arguments are passed when they are called.  


# Running the application
The application can be compiled and executed using either the makefile or the run.sh shell script. The shell script allows to choose data file to read. It uses these commands on both the makefile and the script :  

    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and applications and link them with the static library file into executables.  
    - "app" : Performs the same actions as "all" in addition also to compile the application source file, link it with the library into an executable and execute it using the indicated data file and direct the several output files into the corresponding P1D_<name of the app>_<name of the data file>_<name of the output file kind>.log output file.  

The different output files are :  
    - "results" the standard output of the code.  
    - "matrix" that displays the elements of some matrices for checking purposes.  
    - "potential" that contains the values of the potential on each grid point.  
    - "propagation" that contains the wavepacket values for the main observables along the propagation. For each time point are available : the averaged position, momentum and energy, the norm and the autocorrelation function value.  
    - "psif" that contains, for each grid point (lines), the squared value of several wavepackets at several propagation time (columns). It is used to plot the animation of the propagation.  
    - "psit" that contains, for each time point, the wavepacket vector's coefficients of its decomposition on the basis.  

## Authors
*Developed by*:  
Mathéo Segaud $^{\dagger}$  
*Under the supervision of:*  
Research Professor : David Lauvergnat $^{\dagger,1,2}$  
Ph.D. Student : Rabiou Issa $^{1}$  

${^{\dagger}}$ Université Paris-Saclay, 3 rue Joliot Curie, Bâtiment Breguet, 91190 Gif-sur-Yvette  
${^1}$ Institut de Chimie Physique, 310 Rue Michel Magat, 91400 Orsay  
${^2}$ CNRS  

## Reference
[1] [https://github.com/lauvergn/QDUtilLib)  
cf. manual (to come)  
