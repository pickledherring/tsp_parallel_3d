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
    std::vector<Vertex*> vertexes;
    for (int i = 3; i < argc; i+=3) {
        Vertex v = {std::stof(argv[i], nullptr), i/3 - 1};
    }
    