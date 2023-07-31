# Propagation_1D
Quantum dynamic code which aims to resolve 1D quatum chemistry problems, developped within the scope of a bachelor's interneship in the Institut de Chimie Physique of Saclay.

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

The code needs the library QDUtil to work. It can be downloaded thanks to the script 'get_QDUtilLib.sh' or in the Lauvergn's GitHub: https://github.com/lauvergn/QDUtilLib.