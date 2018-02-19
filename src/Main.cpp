#include "Data.hpp"
#include "Mesh.hpp"
#include "Cell.hpp"
#include "OutputUtils.hpp"

#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello world, CABARET inda house" << std::endl;

	Mesh *mesh = new Mesh();

	mesh->InitMesh(mesh);
	// OutputUtils::OutputParaview(mesh);

	return 1;
}	
