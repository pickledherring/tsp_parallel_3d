todo: serial pthread mpi
serial: serial.cpp
	g++ -o serial serial.cpp -std=c++17
pthread: pthread.cpp
	g++ -pthread -o pthread pthread.cpp -std=c++17
mpi: mpi.cpp
	mpicxx mpi.cpp -o mpi -std=c++17
clean:
	rm serial pthread mpi