#include <stdio.h>
#include <iostream>
#include <fstream>

#define L_x 1.0
#define L_y 1.0
#define N_x 5.0
#define N_y 5.0

using namespace std;

int main()
{
	int node_index = 1;
	double h_x = 2 * L_x / N_x;
	double h_y = 2 * L_y / N_y;
	ofstream polyfile("Mesh.poly");
	polyfile << (N_x-1)*(N_y-1)+4 << " 2 0 0" << endl;

	for (int i = 0; i < N_y ; i++) {
		polyfile << node_index << " " <<-L_x << " " << -L_y + i * h_y << endl;
		node_index++;
	}
	for (int i = 0; i < N_x; i++) {
		polyfile << node_index << " " << -L_x + i * h_x << " " << L_y << endl;
		node_index++;
	}
	for (int i = N_y; i > 0; i--) {
		polyfile << node_index << " " << L_x << " " << -L_y + i * h_x << endl;
		node_index++;
	}
	for (int i = N_x; i > 0; i--) {
		polyfile << node_index << " " << -L_x + i * h_x << " " << -L_y << endl;
		node_index++;
	}
	polyfile << (N_x - 1)*(N_y - 1) + 4 << " 0" << endl;
	for (int i = 1; i < node_index-1; i++) {
		polyfile << i << " " << i << " " << i+1 << endl;
	}
	polyfile << node_index - 1 << " " << node_index - 1 << " 1" << endl;
	polyfile <<"0" << endl;

    return 0;
}

