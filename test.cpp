#include <vector>
#include <map>
#include <limits>
#include <string>
#include <iostream>
#include <math.h>

struct Vertex {
    float weight;
    int name;
    float coords[2];
    Vertex* parent;
    std::vector<Vertex*> children;
};

struct Edge {
    float weight;
    Vertex* start, * end;
};

void traverse(Vertex* vert) { // traverse and print from passed vertex through tree
    std::cout<<vert->name<<std::endl;
    if (vert->children.size() > 0) {
        for (Vertex* child : vert->children) {traverse(child);}
    }
}

void hamiltonian(Vertex* vert) { // helper for traverse
    traverse(vert);
    std::cout<<vert->name<<std::endl;
}

int main(int argc, char *argv[]) {
    std::vector<Vertex> vertexes;
    for (int i = 1; i < argc - 1; i+=3) {
        Vertex v = {std::stof(argv[i], nullptr),
                            i/3 - 1,
                            {std::stof(argv[i + 1], nullptr), std::stof(argv[ i + 2], nullptr)}
                            };
        vertexes.push_back(v);
    }

    int n_threads = strtol(argv[argc - 1], NULL, 10);
    std::vector<Vertex> vert_cuts[n_threads];
    int rows = (int)sqrt(n_threads);
    float size = 500 / rows; // hard coded coordinate ranges, would be better with a range finding function
    // start threaded code
    for (int i = 0; i < n_threads; i++) {
        for (Vertex v: vertexes) {
            if ((v.coords[0] > i * (size) - 1)
            & (v.coords[0] < (i + 1) * size)
            & (v.coords[1] > (i % rows) * size - 1)
            & (v.coords[1] <  ((i + 1) % rows) * size - 1)) {
                vert_cuts[i].push_back(v);
            }
        }
    }
    