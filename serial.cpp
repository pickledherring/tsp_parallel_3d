#include <vector>
#include <map>
#include <limits>
#include <string>
#include <iostream>

struct Vertex {
    float weight;
    int name;
    Vertex* parent;
    std::vector<Vertex*> children;
};

struct Edge {
    float weight;
    Vertex* start, * end;
};

void traverse(Vertex* vert) { // traverse and print from passed vertex through tree
    // std::cout<<vert->name<<std::endl; // uncomment to print city names
    if (vert->children.size() > 0) {
        for (Vertex* child : vert->children) {traverse(child);}
    }
}

void hamiltonian(Vertex* vert) { // helper for traverse
    // std::cout<<"Start traversal"<<std::endl;
    traverse(vert);
}

int main(int argc, char *argv[]) {
    std::vector<Vertex> vertexes;
    for (int i = 3; i < argc; i+=3) {
        Vertex v = {std::stof(argv[i], nullptr), i/3 - 1};
        vertexes.push_back(v);
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    std::vector<Edge> edges;
    std::vector<bool> selected(vertexes.size());
    std::fill(selected.begin(), selected.end(), false);
    selected[0] = true;
    //building our MST
    while (edges.size() < vertexes.size() - 1) {
        float minimum = 1.1; // each probability product's max is 1
        Edge to_add;

        for (int i = 0; i < selected.size(); i++) {
            if (selected[i]) {
                for (int j = 0; j < selected.size(); j++) {
                    if ((i != j) && !selected[j]) {
                        float weight = vertexes[i].weight * vertexes[j].weight;

                        if (weight < minimum) {
                            minimum = weight;
                            to_add.weight = weight;
                            to_add.start = &vertexes[i];
                            to_add.end = &vertexes[j];}
                    }
                }
            }
        }

        vertexes[to_add.start->name].children.push_back(to_add.end);
        vertexes[to_add.end->name].parent = to_add.start;
        edges.push_back(to_add);
        selected[to_add.end->name] = true;
    }

    hamiltonian(&vertexes[0]);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

    std::cout<<"Finished in "<<diff<<" ms!"<<std::endl;
}