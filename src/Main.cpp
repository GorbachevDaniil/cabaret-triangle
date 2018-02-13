#include <sys/stat.h>
#include "unistd.h"
#include <iostream>
#include <stdio.h>
#include "Mesh.cpp"

int main(int argc, char **argv) {

	printf("Hello world, CABARET inda house\n");

	Mesh *grid = new Mesh();
	
	InitMesh(grid);
}
