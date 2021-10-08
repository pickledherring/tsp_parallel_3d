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

void* MST(void* vertexes) {
    std::vector<Vertex*> verts = *reinterpret_cast<std::vector<Vertex*>*>(vertexes);
    std::vector<Edge> edges;
    std::vector<bool> selected(verts.size(), false);
    // std::cout<<"verts size: "<<verts.size()<<"\tselected size: "<<selected.size()<<std::endl;
    selected[0] = true;
    // std::cout<<"made it to e\tselected 0, 1: "<<selected[0]<<"\t"<<selected[1]<<std::endl;
    //building our MST

    // for (int i = 0; i <verts.size(); i++) {
    //     std::cout<<"in MST: name: "<<verts[i]->name<<std::endl;
    // }
    // std::cout<<"made it to e"<<std::endl;
    while (edges.size() < verts.size() - 1) {
        float minimum = 1.1; // each probability product's max is 1
        Edge to_add;

        for (int i = 0; i < selected.size(); i++) {
            if (selected[i]) {
                for (int j = 0; j < selected.size(); j++) {
                    if ((i != j) && !selected[j]) {
                        float weight = verts[i]->weight * verts[j]->weight;

                        if (weight < minimum) {
                            minimum = weight;
                            // std::cout<<"minimum: "<<minimum<<std::endl;
                            // std::cout<<"i: "<<verts[i]->name<<std::endl;
                            // std::cout<<"j: "<<verts[j]->name<<std::endl;
                            to_add.weight = weight;
                            to_add.start = verts[i];
                            to_add.end = verts[j];}
                    }
                }
            }
        }
        // std::cout<<"weight: "<<to_add.weight<<std::endl;
        // std::cout<<"start: "<<to_add.start->name<<std::endl;
        // std::cout<<"end: "<<to_add.end->name<<std::endl;
        for (int i = 0; i < verts.size(); i++) {
            if (verts[i]->name == to_add.start->name) {
                verts[i]->children.push_back(to_add.end);
            }
            if (verts[i]->name == to_add.end->name) {
                verts[i]->parent = to_add.start;
                selected[i] = true;
            }
        }
        edges.push_back(to_add);
    }
    // // std::cout<<"size: "<<verts.size()<<std::endl;
    // for (int i = 0; i <verts.size(); i++) {
    //     std::cout<<"name: "<<verts[i]->name<<std::endl;
    //     if (verts[i]->parent) {std::cout<<"\tparent: "<<verts[i]->parent->name<<std::endl;}
    //     for (int j = 0; j < verts[i]->children.size(); j++) {
    //         std::cout<<"\tchildren["<<j<<"]: "<<verts[i]->children[j]->name<<std::endl;
    //     }
    // }
    // std::cout<<"made it to f"<<std::endl;
    pthread_exit(0);
}

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
        Vertex v = {std::stof(argv[i + 2], nullptr),
                            i / 3,
                            {std::stof(argv[i], nullptr), std::stof(argv[ i + 1], nullptr)}
                            };
        vertexes.push_back(v);
    }
    int n_threads = strtol(argv[argc - 1], NULL, 10);
    std::vector<Vertex*> vert_cuts[n_threads];
    int rows = (int)sqrt(n_threads);
    if ((float)rows < sqrt(n_threads)) {
        std::cout<<"non-square (1, 4, 9, 16, ...) # of threads, using "<<
                                rows * rows<<" threads instead."<<std::endl;
        n_threads = rows * rows;
    }
    float size = 500 / rows; // hard coded coordinate ranges, would be better with a range finding function
    // start threaded code
    for (int i = 0; i < n_threads; i++) {
        for (int j = 0; j < vertexes.size(); j++) {
            // std::cout<<"["<<i<<"]["<<vs<<"] vertex: "<<vertexes[j].name<<
            //                 "\tcoords: "<<vertexes[j].coords[0]<<", "<<vertexes[j].coords[1]<<std::endl;
            if ((vertexes[j].coords[0] > i / rows * size)
                && (vertexes[j].coords[0] < (i / rows + 1) * size)
                && (vertexes[j].coords[1] > i % rows * size)
                && (vertexes[j].coords[1] <  ((i % rows + 1) * size))) {
                    vert_cuts[i].push_back(&vertexes[j]);
                    // std::cout<<"vert_cuts["<<i<<"] got vertex "<<vertexes[j].name<<std::endl;
            }
        }
    }
    pthread_t* threads;
    threads = (pthread_t*)malloc(n_threads * sizeof(pthread_t));
    // std::cout<<"vert_cuts[0][1].name: "<<vert_cuts[0][1]->name<<std::endl;

    for (int i = 0; i < n_threads; i++) {
        pthread_create(&threads[i], NULL, MST, (void*) &vert_cuts[i]);
    }
    for (int i = 0; i < n_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    // std::cout<<"made it to g"<<std::endl;
    free(threads);

    for (int i = 0; i < vertexes.size(); i++) {
        std::cout<<"vertexes["<<i<<"]: "<<std::endl;
        std::cout<<"\tname: "<<vertexes[i].name<<std::endl;
        if (vertexes[i].parent) {std::cout<<"\tparent: "<<vertexes[i].parent->name<<std::endl;}
        for (int j = 0; j <vertexes[i].children.size(); j++) {
            std::cout<<"\tchild["<<j<<"]: "<<vertexes[i].children[j]->name<<std::endl;
        }
    }
    // std::cout<<"made it to h"<<std::endl;
    hamiltonian(&vertexes[0]);
}