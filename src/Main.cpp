#include "Data.hpp"
#include "Mesh.hpp"

#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello world, CABARET inda house" << std::endl;

    Data data;
    data.u0 = 1;
    data.u1 = 2;
    data.u2 = 3;

    std::cout << data.u0 << std::endl;
    std::cout << data.u1 << std::endl;
    std::cout << data.u2 << std::endl;
}
