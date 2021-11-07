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
    // build a minimum spanning tree
    std::vector<bool> selected(verts.size(), false);  // which vertexes are in the tree
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
                            weight_min = norm;
                            start = i;
                            end = j;
                        }
                    }
                }
            }
        }
                                    
        verts[start].children.push_back(&verts[end]);
        verts[end].parent = &verts[start];
        selected[end] = true;
        n_edges++;
    }
}

void traverse_build(std::vector<Vertex> &cycle, Vertex vert) {
    cycle.push_back(vert);
    if (vert.children.size() > 0) {
        for (Vertex* child : vert.children) {traverse_build(cycle, *child);}
    }
}

std::vector<Vertex> build_cycle(std::vector<Vertex> &tree) {
    // helper for recursive traverse, builds one cycle out of a tree
    std::vector<Vertex> cycle;
    traverse_build(cycle, tree[0]);
    return cycle;
}

int orientation(float p[], float q[], float r[]) {
    /*
    based on https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/?ref=lbp
    p, q, and r are points
    is (q[1] - p[1]) / (q[0] - p[0]) > (r[1] - q[1]) / (r[0] - q[0]) ?
    val is the cross product of that inequality
    */
    float val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    if (val > 0) {
        return 1; // Clockwise orientation
    } else if (val < 0) {
        return 2;// Counterclockwise orientation
    } else {
        return 0; // Collinear orientation
    }
}

// collinear points p, q, r, checks if point r lies on line segment 'pq'
bool on_segment(float p[], float q[], float r[]) {
    if (
         (r[0] <= std::max(p[0], q[0])) &
         (r[0] >= std::min(p[0], q[0])) &
         (r[1] <= std::max(p[1], q[1])) &
         (r[1] >= std::min(p[1], q[1]))
        ) {
        return true;
    }
    return false;
}

// returns true if the line segments 'pq' and 'rs' intersect
bool do_intersect(float p[], float q[], float r[], float s[]) {
    int pqr = orientation(p, q, r);
    int pqs = orientation(p, q, s);
    int rsp = orientation(r, s, p);
    int rsq = orientation(r, s, q);

    // general case
    if ((pqr != pqs) & (rsp != rsq)) {
        return true;
    }

    // if points are collinear and the third lies on the segment formed
    // by the first two
    if (
        ((pqr == 0) & on_segment(p, q, r)) |
        ((pqs == 0) & on_segment(p, q, s)) |
        ((rsp == 0) & on_segment(r, p, s)) |
        ((rsq == 0) & on_segment(r, s, q))
        ) {
        return true;
    }
    return false;
}

void inversion_handle(std::vector<Vertex> &cycle) {
    // uncrosses any crossed edges by switching the first two vertexes from each
    if (cycle.size() > 3) {
        // compare each edge in our cycle to the subsequent edges
        for (int j = 0; j < cycle.size() - 3; j++) {
            for (int k = 2; k < cycle.size() - 1; k++) {
                if (do_intersect(cycle[j].coords, cycle[j + 1].coords, cycle[k].coords, cycle[k + 1].coords)) {
                    Vertex temp = {
                        cycle[j].weight,
                        cycle[j].name,
                        {cycle[j].coords[0], cycle[j].coords[1]},
                        cycle[j].parent,
                        cycle[j].children,
                    };
                    cycle[j] = {
                        cycle[k].weight,
                        cycle[k].name,
                        {cycle[k].coords[0], cycle[k].coords[1]},
                        cycle[k].parent,
                        cycle[k].children,
                    };
                    cycle[k] = {
                        temp.weight,
                        temp.name,
                        {temp.coords[0], temp.coords[1]},
                        temp.parent,
                        temp.children,
                    };
                }
            }
        }
    }
}

void row_join(int &rank, int &n_cities, int &n_processes, std::vector<Vertex> &cycle,
            float inf_max, float dist_max, int start, int end) {
    // stitches sectors in a row into one big cycle on the leftmost process of the row
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

        // the joining
        if (r_index > 0) {
            cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin(), vec_rec_cycle.begin() + r_index + 1);
        }
        cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin() + r_index + 1, vec_rec_cycle.end());
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
    // stitches column together into one big cycle on the uppermost process
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
        
        // the joining
        if (r_index > 0) {
            cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin(), vec_rec_cycle.begin() + r_index + 1);
        }
        cycle.insert(cycle.begin() + l_index + 1, vec_rec_cycle.begin() + r_index + 1, vec_rec_cycle.end());
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
    
    // currently only supports a square # of processes
    if (sqrt((float)n_processes) > (int)sqrt(n_processes) && rank == 0) {
        std::cout<<"Non-square (1, 4, 9, 16, ...) # of processes, use "<<
                (int)sqrt(n_processes) * (int)sqrt(n_processes)<<" processes instead."<<std::endl;
        MPI_Finalize();
        exit(0);
    }

    n_cities = std::stoi(argv[1], nullptr);
    // small discrepancy from equal division - increases n_cities if necessary
    
    if (n_cities % n_processes != 0) {
        if (rank == 0) {std::cout<<"Number of cities increased to divide equally among processes\n"<<std::endl;}
        n_cities += n_processes - n_cities % n_processes;
        if (rank == 0) {std::cout<<"New number of cities is "<<n_cities<<std::endl;}
    }
    // if (n_cities <  2 * n_processes) {
    //     if (rank == 0) {std::cout<<"Number of cities increased to two per process. Run with fewer processes"<<
    //                                                     " if you want to route through fewer cities\n"<<std::endl;}
    //     n_cities =  2 * n_processes;
    //     if (rank == 0) {std::cout<<"New number of cities is "<<n_cities<<std::endl;}
    // }

    // create our cities, compute the max infection probability and max distance
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
    struct timespec start, end;
    if (rank == 0) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    }

    // build the MST for each process, make cycle out of it
    MST(verts, inf_max, dist_max);
    std::vector<Vertex> cycle = build_cycle(verts);
    inversion_handle(cycle);

    // stitch by row, then column
    row_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, 0, (int)sqrt(n_processes) - 1);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank % (int)sqrt(n_processes) == 0) {
        col_join(rank, n_cities, n_processes, cycle, inf_max, dist_max, 0, (int)sqrt(n_processes) - 1);
    }

    if (rank == 0) {
        inversion_handle(cycle);

        // end timing
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;
        std::cout<<"Finished in "<<diff<<" ms!"<<std::endl;

        // uncomment to show final path
        // std::cout<<"start:"<<std::endl;
        // for (Vertex v: cycle) {
        //     std::cout<<v.name<<std::endl;
        // }
        // std::cout<<"finish"<<std::endl;
    }

    MPI_Finalize();
    return 0;
}