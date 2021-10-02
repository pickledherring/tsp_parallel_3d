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
    std::vector<Vertex> vertexes;
    for (int i = 3; i < argc; i+=3) {
        Vertex v = {std::stof(argv[i], nullptr), i/3 - 1};
        vertexes.push_back(v);
    }

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
                            std::cout<<"minimum: "<<minimum<<std::endl;
                            std::cout<<"i: "<<vertexes[i].name<<std::endl;
                            std::cout<<"j: "<<vertexes[j].name<<std::endl;
                            to_add.weight = weight;
                            to_add.start = &vertexes[i];
                            to_add.end = &vertexes[j];}
                    }
                }
            }
        }
        std::cout<<"weight: "<<to_add.weight<<std::endl;
        std::cout<<"start: "<<to_add.start->name<<std::endl;
        std::cout<<"end: "<<to_add.end->name<<std::endl;
        vertexes[to_add.start->name].children.push_back(to_add.end);
        vertexes[to_add.end->name].parent = to_add.start;
        edges.push_back(to_add);
        selected[to_add.end->name] = true;
    }
    // for (int i = 0; i < edges.size(); i++) {
    //     std::cout<<"edges["<<i<<"]: "<<std::endl;
    //     std::cout<<"weight: "<<edges[i]->weight<<std::endl;
    //     std::cout<<"start: "<<edges[i]->start->name<<std::endl;
    //     std::cout<<"end: "<<edges[i]->end->name<<std::endl;
    // }
    for (int i = 0; i < vertexes.size(); i++) {
        std::cout<<"vertexes["<<i<<"]: "<<std::endl;
        std::cout<<"name: "<<vertexes[i].name<<std::endl;
        std::cout<<"weight: "<<vertexes[i].weight<<std::endl;
    }
    
    // hamiltonian(vertexes[0]);
}