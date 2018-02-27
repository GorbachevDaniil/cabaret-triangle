#include "Mesh.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//Make vector of values from string words
vector<string> split(const string& s, const char *delimiter) {
    vector<string> values;
    char *pch;
    char *cstr = new char[s.length() + 1];
    strcpy(cstr, s.c_str());
    pch = std::strtok(cstr, delimiter);

    while (pch != NULL) {
        values.push_back(pch);
        pch = std::strtok(NULL, delimiter);
    }

    delete[] cstr;
    return values;
}

int LoadNodes(Mesh *grid) {
    string line;
    ifstream nodefile("triangle/Mesh.1.node");

    if (!nodefile) {
        cout << "Error opening node file!" << endl;
        return 1;
    }

    vector<string> values;
    int i = 0, nNodes;

    while (getline(nodefile, line) && (i != (nNodes + 1))) {
        if (i == nNodes + 1)
            continue; // miss last string

        values = split(line, " ");
        if (i == 0) {
            i++;
            nNodes = atoi(values[0].c_str()); //init number of Nodes
            continue;
        }

        Node *node = new Node(*grid, atof(values[1].c_str()), atof(values[2].c_str()));
        grid->nodes.push_back(*node);

        i++;
    }

    return 0;
}

int LoadEdges(Mesh *grid) {
    string line;
    ifstream edgefile("triangle/Mesh.1.edge");

    if (!edgefile) {
        cout << "Error opening edge file!" << endl;
        return 1;
    }

    vector<string> values;
    int i = 0, nEdges;

    while (getline(edgefile, line) && (i != (nEdges + 1))) {
        if (i == nEdges + 1)
            continue; // miss last string

        values = split(line, " ");
        if (i == 0) {
            i++;
            nEdges = atoi(values[0].c_str()); //init number of Nodes
            continue;
        }

        Edge *edge = new Edge(*grid, atol(values[0].c_str()) - 1, atol(values[1].c_str()) - 1,
                atol(values[2].c_str()) - 1, atoi(values[3].c_str()));
        grid->edges.push_back(*edge);

        i++;
    }

    return 0;
}

int LoadCells(Mesh *grid) {
    string line;
    ifstream cellfile("triangle/Mesh.1.ele");

    if (!cellfile) {
        cout << "Error opening ele file!" << endl;
        return 1;
    }

    vector<string> values;
    int i = 0, nCells;

    while (getline(cellfile, line) && (i != (nCells + 1))) {
        if (i == nCells + 1)
            continue; // miss last string

        values = split(line, " ");
        if (i == 0) {
            i++;
            nCells = atoi(values[0].c_str()); //init number of Nodes
            continue;
        }

        Cell *cell = new Cell(*grid, atol(values[0].c_str()) - 1, atol(values[1].c_str()) - 1,
                atol(values[2].c_str()) - 1, atol(values[3].c_str()) - 1);
        grid->cells.push_back(*cell);

        i++;
    }

    return 0;
}

int Mesh::InitMesh(Mesh *mesh) {
    cout << "Loading mesh" << endl;

    LoadNodes(mesh);
    LoadEdges(mesh);
    LoadCells(mesh);

    return 0;
}
