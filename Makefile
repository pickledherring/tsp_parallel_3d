todo: new_mpi
new_mpi: new_mpi.cpp
	mpicxx new_mpi.cpp -o new_mpi -std=c++17
clean:
	rm new_mpi