to run: 
1. make
2. mpiexec --oversubscribe  -np \<number of processes\> ./new_mpi \<number of cities\>

make sure the number of processes is a power of 2. if the number of cities is not at least twice the number of processes, the program will make it so.