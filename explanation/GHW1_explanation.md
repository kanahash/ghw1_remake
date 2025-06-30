=============================================================================
                                     GHW1 
=============================================================================

2D delta-f gyrofluid (modified) Hasegawa-Wakatani Code with periodic b.c.
(by Alexander Kendl, Innsbruck; ghw1 version V1.25 of 10.01.2025)

( changes to V1.24: - speedup of FFT by using r2c/c2r instead of complex FFT
                    - put initial conditions into functions for clarity )

Physics: 
Simulation of 2D drift wave turbulence and vortex dynamics in magnetised plasma.
In cold ion (taui=0) limit: comparable to results by Numata in arXiv:0708.4317.

Examples for usage (V1.0) : arXiv:1708.06213
Normalisation and numerical methods: in similar full-f version see arXiv:1502.05494.
More background on gyrofluid theory and models see references therein.
Parameters are set in 'ghw1.inp' input file. Explanations given within the code below.

Dependencies:
Needs some libraries (e.g. fftw3, omp) available in every Linux distribution.

For command line execution under Linux first apply these settings:
  ulimit -s unlimited
  export OMP_STACKSIZE=8G
Compilation with attached makefile:
  make ghw1
Execution with run time control:
  time ./ghw1

Parallelisation with OpenMP is not extremely efficient, but good speed-up up to np=8.
Run time depends on problem size: a few minutes for 64^2, many hours for 1024^2.

In this code version most subroutines are inlined instead of calling functions:
allthough not elegant in coding, this considerably speeds up the code execution.

Time step and hyperviscosity have to be refined with the grid resolution
to prevent crashing through numerical instability. (Set "dt" and "diff" by trial and error.)
Time step (and typical time scales) also sensitive to choice of "chat" and "delta".

Running in a RAM disc folder (here: at location '/ramdisc/') saves some output time:
'mkdir /ramdisc', 'sudo mount -t tmpfs none /ramdisc -o size=2048m'.

Visualisation of 1D graphs with any plotting programme (e.g. xmgrace),
visualisation of 2D contour plots with attached gnuplot script file "gpl-ctr.plt".
Watch the simulation live by uncommenting "reread". In gnuplot: "load 'gpl-ctr.plt'".

Makefile example:
exec = ghw1
objs = ghw1.o 
LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3 
CFLAGS = -I.
LIBS = /usr/lib64 
all: $(exec)
$(exec): $(objs)
	g++  $(exec).C -o $@  $(LFLAGS) -L$(LIBS)
	rm *.o
ghw.o: ghw1.C 
	g++  -c $< -o $@

==============================================================================