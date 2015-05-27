PoissonSolver
===============================

This program solves the discrete poisson equation in the unit square in two dimensions.

INSTALLATION
-------------------------------

1.) Download the source code
2.) Run the Makefile


RUN THE PROGRAM
-------------------------------

1.) Execute ./poissonSolver program
2.) add three more command line arguments:
..* m: the parameter which specifies h = 1/m (you observer N = m - 1 grid points in x and y direction respectively)
..* x: the chosen algorithm (CG-Method, PCG-Method with incomplete LU decomposition, Jacobi-Method, Multigrid-Method)
..* k: the chose function -> default: laplace(u) = -4 in Omega and u(x) = x^2 + y^2 in the bound of Omega

FURTHER INFORMATION
-------------------------------

1.) Please, visit http://homepages.uni-regensburg.de/~bam23651/PoissonSolver and try the website!