#include <vector>
#include <map>
#include <limits>
#include <string>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <random>
#include <float.h>

struct Vertex {
    float weight;
    int name;
    float coords[2];
    Vertex* parent;
    std::vector<Vertex*> children;
};

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

float edge_weight(Vertex &left, Vertex &right, float inf_max, float dist_max) {
    float x = left.coords[0] - right.coords[0];
    float y = left.coords[1] - right.coords[1];
    float distance = sqrt(x * x + y * y);

    float infection = left.weight * right.weight;

    return distance / dist_max + infection / inf_max;
}

void MST(std::vector<Vertex> &verts, float inf_max, float dist_max) {
    // make a tree from the vertexes, optimal but might cross edges, then make cycle
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
                        float norm = edge_weight(verts[i], verts[j], inf_max, dist_max);

                        if (norm < weight_min) {
                            // std::cout<<i<<" and "<<j<<" have new normal: "<<norm<<std::endl;
                            weight_min = norm;
                            start = i;
                            end = j;
                        }
                    }
                }
            }
        }
                                    
        // std::cout<<start<<" and "<<end<<" have new edge c -> p"<<std::endl;
        verts[start].children.push_back(&verts[end]);
        verts[end].parent = &verts[start];
        selected[end] = true;
        n_edges++;
    }
}

void traverse(std::vector<Vertex> &cycle, Vertex vert) {
    // std::cout<<vert.name<<std::endl;
    cycle.push_back(vert);
    if (vert.children.size() > 0) {
        for (Vertex* child : vert.children) {traverse(cycle, *child);}
    }
}

std::vector<Vertex> build_cycle(std::vector<Vertex> &verts) {
    // helper for recursive traverse
    std::vector<Vertex> cycle;
    traverse(cycle, verts[0]);
    return cycle;
}

void row_join(int &rank, int &n_cities, int &n_processes, std::vector<Vertex> &cycle,
            float inf_max, float dist_max, int start, int end) {
    int width = (int)sqrt(n_processes);
    int row = rank / width;
    int col = rank % width;
    int range = end - start;
    int mid = start + range / 2;
    // stop criteria
    if (range > 1) {
        if (col > mid) {row_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, mid + 1, end);}
        if (col <= mid) {row_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, start, mid);}
    }

    // receive if start, compute lowest cost for new edges between cycles, join
    if (col == start) {
        int msg_size = log2(range + 1) * (n_cities / n_processes) * 4;
        float rec_cycle[msg_size];
        int source = row * width + mid + 1;
        MPI_Status Stat;
        MPI_Recv(&rec_cycle, msg_size, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &Stat);

        // create vector of vertexes out of received data
        std::vector<Vertex> vec_rec_cycle;
        for (int i = 0; i < (msg_size / 4); i++) {
            Vertex new_vert = {
                rec_cycle[i * 4],
                (int)rec_cycle[i * 4 + 1],
                {rec_cycle[i * 4 + 2], rec_cycle[i * 4 + 3]}
            };
            vec_rec_cycle.push_back(new_vert);
        }

        // minimum cost edges to switch
        float min = FLT_MAX;
        int l_index = 0;
        int r_index = 0;
        for (int i = 0; i < cycle.size(); i++) {
            float dist_left = edge_weight(cycle[i], cycle[i + 1], inf_max, dist_max);

            for (int j = 0; j < vec_rec_cycle.size(); j++) {
                float dist_right = edge_weight(vec_rec_cycle[j], vec_rec_cycle[j + 1], inf_max, dist_max);
                float dist_upper = edge_weight(cycle[i], vec_rec_cycle[j + 1], inf_max, dist_max);
                float dist_lower = edge_weight(cycle[i + 1], vec_rec_cycle[j], inf_max, dist_max);
                float cost = dist_upper + dist_lower - dist_left - dist_right;

                if (cost < min) {
                    min = cost;
                    l_index = i;
                    r_index = j;
                }
            }
        }
        std::cout<<"rank "<<rank<<" has vertexes: "<<std::endl;
        std::cout<<"from cycle: "<<std::endl;
        for (int i = 0; i < cycle.size(); i++) {
            std::cout<<"\t"<<cycle[i].name;
        }
        std::cout<<std::endl;
        std::cout<<"from vec_rec_cycle: "<<std::endl;
        for (int i = 0; i < vec_rec_cycle.size(); i++) {
            std::cout<<"\t"<<vec_rec_cycle[i].name;
        }
        std::cout<<std::endl;
        // the joining
        if (r_index > 0) {
            cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin(), vec_rec_cycle.begin() + r_index + 1);
        }
        cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin() + r_index + 1, vec_rec_cycle.end());
        
        std::cout<<"rank "<<rank<<" has vertexes: "<<std::endl;
        for (int i = 0; i < cycle.size(); i++) {
            std::cout<<"\t"<<cycle[i].name;
        }
        std::cout<<std::endl;
    }

    // send if mid + 1
    if (col ==  mid + 1) {
        // convert cycle to floats for message passing
        float float_cycle[4 * cycle.size()];
        for (int i = 0; i < cycle.size(); i++) {
            float_cycle[i * 4] = cycle[i].weight;
            float_cycle[i * 4 + 1] = (float)cycle[i].name;
            float_cycle[i * 4 + 2] = cycle[i].coords[0];
            float_cycle[i * 4 + 3] = cycle[i].coords[1];
        }
        int msg_size = log2(range + 1) * (n_cities / n_processes) * 4;
        int dest = row * width + start;
        MPI_Send(&float_cycle, msg_size, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    }
}

void col_join(int &rank, int &n_cities, int &n_processes, std::vector<Vertex> &cycle,
            float inf_max, float dist_max, int start, int end) {
    int width = (int)sqrt(n_processes);
    int row = rank / width;
    int range = end - start;
    int mid = start + range / 2;
    // stop criteria
    if (range > 1) {
        if (row > mid) {col_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, mid + 1, end);}
        if (row <= mid) {col_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, start, mid);}
    }

    // receive if start, compute lowest cost for new edges between cycles, join
    if (row == start) {
        int msg_size = log2(range + 1) * (n_cities / n_processes) * width * 4;
        float rec_cycle[msg_size];
        int source = width * (mid + 1);
        MPI_Status Stat;
        MPI_Recv(&rec_cycle, msg_size, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &Stat);

        // create vector of vertexes out of received data
        std::vector<Vertex> vec_rec_cycle;
        for (int i = 0; i < (msg_size / 4); i++) {
            Vertex new_vert = {
                rec_cycle[i * 4],
                (int)rec_cycle[i * 4 + 1],
                {rec_cycle[i * 4 + 2], rec_cycle[i * 4 + 3]}
            };
            vec_rec_cycle.push_back(new_vert);
        }

        // minimum cost edges to switch
        float min = FLT_MAX;
        int l_index = 0;
        int r_index = 0;
        for (int i = 0; i < cycle.size(); i++) {
            float dist_left = edge_weight(cycle[i], cycle[i + 1], inf_max, dist_max);

            for (int j = 0; j < vec_rec_cycle.size(); j++) {
                float dist_right = edge_weight(vec_rec_cycle[j], vec_rec_cycle[j + 1], inf_max, dist_max);
                float dist_upper = edge_weight(cycle[i], vec_rec_cycle[j + 1], inf_max, dist_max);
                float dist_lower = edge_weight(cycle[i + 1], vec_rec_cycle[j], inf_max, dist_max);
                float cost = dist_upper + dist_lower - dist_left - dist_right;

                if (cost < min) {
                    min = cost;
                    l_index = i;
                    r_index = j;
                }
            }
        }
        std::cout<<"rank "<<rank<<" has vertexes: "<<std::endl;
        std::cout<<"from cycle: "<<std::endl;
        for (int i = 0; i < cycle.size(); i++) {
            std::cout<<"\t"<<cycle[i].name;
        }
        std::cout<<std::endl;
        std::cout<<"from vec_rec_cycle: "<<std::endl;
        for (int i = 0; i < vec_rec_cycle.size(); i++) {
            std::cout<<"\t"<<vec_rec_cycle[i].name;
        }
        std::cout<<std::endl;
        // the joining
        if (r_index > 0) {
            cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin(), vec_rec_cycle.begin() + r_index + 1);
        }
        cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin() + r_index + 1, vec_rec_cycle.end());
        
        std::cout<<"rank "<<rank<<" has vertexes: "<<std::endl;
        for (int i = 0; i < cycle.size(); i++) {
            std::cout<<"\t"<<cycle[i].name;
        }
        std::cout<<std::endl;
    }

    // send if mid + 1
    if (row ==  mid + 1) {
        // convert cycle to floats for message passing
        float float_cycle[4 * cycle.size()];
        for (int i = 0; i < cycle.size(); i++) {
            float_cycle[i * 4] = cycle[i].weight;
            float_cycle[i * 4 + 1] = (float)cycle[i].name;
            float_cycle[i * 4 + 2] = cycle[i].coords[0];
            float_cycle[i * 4 + 3] = cycle[i].coords[1];
        }
        int msg_size = log2(range + 1) * (n_cities / n_processes) * width * 4;
        int dest = width * start;
        MPI_Send(&float_cycle, msg_size, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int n_processes, rank, n_cities;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (sqrt((float)n_processes) > (int)sqrt(n_processes) && rank == 0) {
        std::cout<<"Non-square (1, 4, 9, 16, ...) # of processes, use "<<
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

    // create our cities, compute the max infection and distance edge weights
    std::vector<Vertex> verts;
    gen_verts(rank, n_cities, n_processes, verts);
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

    float all_inf_max;
    MPI_Reduce(&inf_max, &all_inf_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    float all_dist_max;
    MPI_Reduce(&dist_max, &all_dist_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    // start timing


    // build the MST for each process, make cycle out of it
    MST(verts, inf_max, dist_max);
    std::vector<Vertex> cycle = build_cycle(verts);

    // stitch by row, then column
    row_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, 0, (int)sqrt(n_processes) - 1);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank % (int)sqrt(n_processes) == 0) {
        col_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, 0, (int)sqrt(n_processes) - 1);
    }

    // end timing


    // for (int i = 0; i < n_processes; i++) {
    //     if (rank == i) {
    //         std::cout<<"rank is "<<i<<", cycle: "<<std::endl;
    //         for (int j = 0; j < (n_cities / n_processes); j++) {
    //             std::cout<<"\t"<<cycle[j].name<<"\t"<<cycle[j].weight<<std::endl;
    //             std::cout<<"\t"<<cycle[j].coords[0]<<", "<<cycle[j].coords[1]<<std::endl;
    //         }
    //     }
    // }

    MPI_Finalize();
    exit(0);
}