#ifndef Node_hpp
#define Node_hpp

#include "Data.hpp"
#include "Vector.hpp"

class Node {
public:
    long id;

    Data data;
    Vector normal;
};

#endif