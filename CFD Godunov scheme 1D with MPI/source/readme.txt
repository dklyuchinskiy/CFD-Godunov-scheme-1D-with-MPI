The one-dimensional Godunov scheme of the first order for the equations of fluid dynamics
has been implemented via using OpenMP + MPI technology for multi-threading.

The two versions of the scheme are enabled: with nonlinear (1959) and linear (2018) solution of the Riemann problem on the boundaries.
The user can switch them in the code by itself (see definitions.h)

For better performance, please compile the code with Intel Compiler:
mpiicpc -ipo -O2 -qopenmp -o a.out main.cpp functions.cpp arrays.cpp

Run the appilcation as
mpirun -n NUMBER_OF_PROCESSORS ./a.out

Before running, set the number (n) of OpenMP threads with command:
set OMP_NUM_THREADS=n (for Windows)
export OMP_NUM_THREADS=n (for Linux)


Please, ask your questions at dmitriy_klyuchinskiy@mail.ru
