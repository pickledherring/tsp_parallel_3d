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
    selected[0] = true;
    //building our MST
    float min_in_sector = 1.1; // minimum weight for this thread
    int min_index;
    for (int i = 0; i <verts.size(); i++) {
        if (verts[i]->weight < min_in_sector) {
            min_in_sector = verts[i]->weight;
            min_index = i;
        }
    }

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
                            to_add.weight = weight;
                            to_add.start = verts[i];
                            to_add.end = verts[j];}
                    }
                }
            }
        }
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
    
    return (void*) min_index;
}

void traverse(Vertex* vert) { // traverse and print from passed vertex through tree
    std::cout<<vert->name<<std::endl;
    if (vert->children.size() > 0) {
        for (Vertex* child : vert->children) {traverse(child);}
    }
}

void hamiltonian(std::vector<Vertex*>* vert_cuts, std::vector<int> mins, int n_sectors) {
    for (int i = 1; i < n_sectors; i++) { // linking sectors by their minimum-weighted vertexes
        for (int j = 0; j < vert_cuts[i].size(); j++) {
            if (vert_cuts[i][j]->parent == nullptr) {
                if (i % (int)sqrt(n_sectors) > 0 ) { // trying to link proximal sectors
                    vert_cuts[i][j]->parent = vert_cuts[i - 1][mins[i - 1]];
                    vert_cuts[i - 1][mins[i - 1]]->children.push_back(vert_cuts[i][j]);
                } else {
                    vert_cuts[i][j]->parent = vert_cuts[i - (int)sqrt(n_sectors)][mins[i - (int)sqrt(n_sectors)]];
                    vert_cuts[i - (int)sqrt(n_sectors)][mins[i - (int)sqrt(n_sectors)]]->children.push_back(vert_cuts[i][j]);
                }
                break;
            }
        }
    }

    std::cout<<"Start traversal"<<std::endl;
    // vert_cuts[0][0] has no parent, must start there
    traverse(vert_cuts[0][0]);
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
    // hard coded coordinate ranges, would be better with a range finding function
    float size = 500 / rows;

    // start threaded code
    for (int i = 0; i < n_threads; i++) {
        for (int j = 0; j < vertexes.size(); j++) {
            if ((vertexes[j].coords[0] > i / rows * size)
                && (vertexes[j].coords[0] < (i / rows + 1) * size)
                && (vertexes[j].coords[1] > i % rows * size)
                && (vertexes[j].coords[1] <  ((i % rows + 1) * size))) {
                    vert_cuts[i].push_back(&vertexes[j]);
            }
        }
    }
    std::vector<int> min_from_sector;
    void* min;
    pthread_t* threads;
    threads = (pthread_t*)malloc(n_threads * sizeof(pthread_t));

    for (int i = 0; i < n_threads; i++) {
        pthread_create(&threads[i], NULL, MST, (void*) &vert_cuts[i]);
    }
    for (int i = 0; i < n_threads; i++) {
        pthread_join(threads[i], &min);
        min_from_sector.push_back(reinterpret_cast<long>(min));
    }

    free(threads);
    hamiltonian(vert_cuts, min_from_sector, n_threads);
    std::cout<<"Finished!"<<std::endl;
    return 0;
}