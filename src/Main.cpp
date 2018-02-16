#include "Data.hpp"
#include "Mesh.hpp"

#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello world, CABARET inda house" << std::endl;

	Mesh *grid = new Mesh();
	
	grid->InitMesh(grid);
}	
