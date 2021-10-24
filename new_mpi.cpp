#include <vector>
#include <map>
#include <limits>
#include <string>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <random>

struct Vertex {
    float weight;
    int name;
    float coords[2];
    Vertex* parent;
    std::vector<Vertex*> children;
};

// struct Edge {
//     float weight;
//     Vertex* start, * end;
// };

void gen_verts(int rank, int n_cities, int n_processes, std::vector<Vertex> &verts)
// creates cities for each process, equally divides n_cities among them
{
    int block_length = 500 / sqrt(n_processes);
    for (int i = 0; i < (n_cities / n_processes); i++) {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist_block(0, block_length);
        std::uniform_int_distribution<std::mt19937::result_type> dist_infect(0, 100);
        float x_rand = dist_block(rng);
        float y_rand = dist_block(rng);
        float infection = (dist_infect(rng))/100.0;

        Vertex new_vert = {infection, rank * n_cities / n_processes + i,
            {block_length * (rank % (int)sqrt(n_processes)) + x_rand,
            block_length * (rank / (int)sqrt(n_processes)) + y_rand
            }
        };
        verts.push_back(new_vert);
    }
}

void MST(std::vector<Vertex> &verts) 
// make a tree from the vertexes, optimal but might cross edges, then make cycle
{
    // find maximums to normalize weights
    float inf_max, dist_max = 0;
    for (int i = 0; i < verts.size(); i++) {
        for (int j = 0; j < verts.size(); j++) {
            if (verts[i].weight * verts[j].weight > inf_max) {
                inf_max = verts[i].weight * verts[j].weight;
            }
            float x = verts[i].coords[0] - verts[j].coords[0];
            float y = verts[i].coords[1] - verts[j].coords[1];
            float distance = sqrt(x * x + y * y);
            if (distance > dist_max) {
                dist_max = distance;
            }
        }
    }
    // build the tree
    std::vector<bool> selected(verts.size(), false);
    selected[0] = true;
    int n_edges = 0;
    while (n_edges < verts.size() - 1) {
        float weight_min = 3; // maximum value for the normalized weight is 2
        int start, end;
        for (int i = 0; i < selected.size(); i++) {
            if (selected[i]) {
                for (int j = 0; j < selected.size(); j++) {
                    if ((i != j) && !selected[j]) {
                        float x = verts[i].coords[0] - verts[j].coords[0];
                        float y = verts[i].coords[1] - verts[j].coords[1];
                        float distance = sqrt(x * x + y * y);

                        float infection = verts[i].weight * verts[j].weight;

                        float norm = distance / dist_max + infection / inf_max;

                        if (norm < weight_min) {
                            std::cout<<i<<" and "<<j<<" have new normal: "<<norm<<std::endl;
                            weight_min = norm;
                            start = i;
                            end = j;
                        }
                    }
                }
            }
        }
                                    
        std::cout<<start<<" and "<<end<<" have new edge c -> p"<<std::endl;
        verts[start].children.push_back(&verts[end]);
        verts[end].parent = &verts[start];
        selected[end] = true;
        n_edges++;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int n_processes, rank, n_cities;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (sqrt((float)n_processes) > (int)sqrt(n_processes) && rank == 0) {
        std::cout<<"non-square (1, 4, 9, 16, ...) # of processes, use "<<
                (int)sqrt(n_processes) * (int)sqrt(n_processes)<<" processes instead."<<std::endl;
        MPI_Finalize();
        exit(0);
    }

    n_cities = std::stoi(argv[1], nullptr);
    // small discrepancy from equal division - increases n_cities if necessary
    if (n_cities % n_processes != 0) {
        std::cout<<"Number of cities increased to divide equally among processes\n"<<std::endl;
        n_cities += n_cities % n_processes;
        std::cout<<"New number of cities is "<<n_cities<<std::endl;
    }

    std::vector<Vertex> verts;
    gen_verts(rank, n_cities, n_processes, verts);
    for (int i = 0; i < n_processes; i++) {
        if (rank == i) {
            std::cout<<"rank is "<<i<<", verts: "<<std::endl;
            for (int j = 0; j < (n_cities / n_processes); j++) {
                std::cout<<"\t"<<verts[j].name<<"\t"<<verts[j].weight<<std::endl;
                std::cout<<"\t"<<verts[j].coords[0]<<", "<<verts[j].coords[1]<<std::endl;
            }
        }
    }

    MST(verts);
    MPI_Finalize();
    exit(0);
}