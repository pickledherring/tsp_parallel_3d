todo: mpi
mpi: mpi.cpp
	mpicxx mpi.cpp -o mpi -std=c++17
clean:
	rm mpi