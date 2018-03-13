# CFD-Godunov-scheme-1D-with-MPI

The one-dimensional Godunov scheme of the first order for the equations of fluid dynamics
has been implemented via using OpenMP + MPI technology for multi-threading.

The two versions of the scheme are enabled: with nonlinear (1959) and linear (2018) solution of the Riemann problem on the boundaries.
The user can switch them in the code by itself (see definitions.h)

For better performance, please compile the code with Intel Compiler:

mpiicpc -ipo -O2 -qopenmp -o a.out main.cpp functions.cpp arrays.cpp

Run the appilcation as

mpirun -n NUMBER_OF_PROCESSORS ./a.out

Important: if you run the application at 1 cluster node, please disable OpenMP threading for better performance.

If not, set the number (n) of OpenMP threads with command:

set OMP_NUM_THREADS=n (for Windows)

export OMP_NUM_THREADS=n (for Linux)

The measured performance on Intel(R) Xeon(R) CPU E5-2630 v4 @ 2.20GHz:

24300 points of the grid. 1 node. OpenMP is disabled.

1 proc - 56.7 sec

2 proc - 22.9 sec

4 proc - 12.1 sec

8 proc - 7.2 sec


Please, ask your questions at dmitriy_klyuchinskiy@mail.ru
