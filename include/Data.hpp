#ifndef Data_hpp
#define Data_hpp

#include "Vector.hpp"

#include <vector>
#include <unordered_map>

class Data {
public:
    Vector coords;

    // list of scalar variables on n layer of time
    std::unordered_map<int, double> s0;
    // list of scalar variables on n+1/2 layer of time
    std::unordered_map<int, double> s1;
    // list of scalar variables on n+1 layer of time
    std::unordered_map<int, double> s2;

    // list of vector variables on n layer of time
    std::unordered_map<int, Vector> v0;
    // list of vector variables on n+1/2 layer of time
    std::unordered_map<int, Vector> v1;
    // list of vector variables on n+1 layer of time
    std::unordered_map<int, Vector> v2;

    // list of stationary vector variables
    std::unordered_map<int, Vector> vector;
};

#endif
