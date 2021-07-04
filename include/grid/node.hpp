#ifndef Node_hpp
#define Node_hpp

#include <set>

#include "vector.hpp"

class Mesh;

class Node {
public:
    long id;
    bool used;
    bool is_bound;

    Vector coords;

    std::set<long> cell_ids;
    std::set<long> edge_ids;

    Node(Mesh &mesh, double x, double y, bool used, bool is_bound);
};

#endif
