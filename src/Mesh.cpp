#include "Mesh.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

//Make vector of values from string words
vector<string> split(const string& s, char delimiter)
{
	vector<string> values;
	char *pch;
	char *cstr = new char[s.length() + 1];
	strcpy(cstr, s.c_str());
	pch = strtok(cstr," ");
	
	while (pch != NULL){
		values.push_back(pch);
		pch = strtok(NULL, " ");
	}

	delete [] cstr;
	return values;
}


int LoadNodes(Mesh *grid){
	string line;
	ifstream nodefile("triangle/Mesh.1.node");

	if(!nodefile){
		cout << "Error opening node file!" << endl;
		return 1;
	}

	vector<string> values;
	int i=0, nNodes;

	while (getline(nodefile, line) && (i != (nNodes+1))){
		if (i==nNodes+1) continue; // miss last string 

		values = split(line, ' ');
		if (i==0){
			i++;
			nNodes = atoi(values[0].c_str()); //init number of Nodes
			continue;
		} 

		Vector *vec = new Vector();
		vec->x = atof(values[1].c_str());
		vec->y = atof(values[2].c_str());
		Data *data = new Data();
		data->coords = *vec;
		Node *node = new Node();
		node->data = *data;
		grid->nodes.push_back(*node);

		i++;
	}

	return 0;
}

int LoadEdges(Mesh *grid){
	string line;
	ifstream edgefile("triangle/Mesh.1.edge");

	if(!edgefile){
		cout << "Error opening edge file!" << endl;
		return 1;
	}

	vector<string> values;
	int i=0, nEdges;

	while (getline(edgefile, line) && (i != (nEdges+1))){
		if (i==nEdges+1) continue; // miss last string 

		values = split(line, ' ');
		if (i==0){
			i++;
			nEdges = atoi(values[0].c_str()); //init number of Nodes
			continue;
		} 

		Edge *edge = new Edge();
		edge->ID = atol(values[0].c_str());
		edge->nodeIDs.push_back(atol(values[1].c_str())-1);
		edge->nodeIDs.push_back(atol(values[2].c_str())-1);
		edge->boundEdge = atoi(values[3].c_str());
		grid->edges.push_back(*edge);

		i++;
	}

	return 0;
}

int LoadCells(Mesh *grid){
	string line;
	ifstream cellfile("triangle/Mesh.1.ele");

	if(!cellfile){
		cout << "Error opening edge file!" << endl;
		return 1;
	}

	vector<string> values;
	int i=0, nCells;

	while (getline(cellfile, line) && (i != (nCells+1))){
		if (i==nCells+1) continue; // miss last string 

		values = split(line, ' ');
		if (i==0){
			i++;
			nCells = atoi(values[0].c_str()); //init number of Nodes
			continue;
		} 

		Cell *cell = new Cell();
		cell->ID = atol(values[0].c_str());
		cell->nodeIDs.push_back(atol(values[1].c_str())-1);
		cell->nodeIDs.push_back(atol(values[2].c_str())-1);
		cell->nodeIDs.push_back(atol(values[3].c_str())-1);
		grid->cells.push_back(*cell);

		i++;
	}

	return 0;
}


int Mesh::InitMesh(Mesh *GRID){
	cout << "Loading mesh" << endl;

	LoadNodes(GRID);
	LoadEdges(GRID);
	LoadCells(GRID);

	for(int i=0; i<GRID->nodes.size(); i++){
		cout << i << " " << GRID->nodes[i].data.coords.x << " " << GRID->nodes[i].data.coords.y << endl; 
	}

	for(int i=0; i<GRID->edges.size(); i++){
		cout << i << " " << GRID->edges[i].nodeIDs[0] << " " << GRID->edges[i].nodeIDs[1] << endl; 
	}

	for(int i=0; i<GRID->cells.size(); i++){
		cout << i << " " << GRID->cells[i].nodeIDs[0] << " " << GRID->cells[i].nodeIDs[1] << " " << GRID->cells[i].nodeIDs[2] << endl; 
	}
	
	return 0;
}
