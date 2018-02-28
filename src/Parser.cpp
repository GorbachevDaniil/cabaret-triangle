#include "Mesh.hpp"
#include "Parser.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>

//Make vector of values from string words
std::vector<std::string> split(const std::string& s, const char *delimiter) {
    std::vector<std::string> values;
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

int Parser::LoadNodes(Mesh *mesh, std::string nodeFile) {
    std::string line;
    std::ifstream nodefile(nodeFile);

    if (!nodefile) {
        std::cout << "Error opening node file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
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

        Node *node = new Node(*mesh, atof(values[1].c_str()), atof(values[2].c_str()));
        mesh->nodes.push_back(*node);

        i++;
    }

    return 0;
}

int Parser::LoadEdges(Mesh *mesh, std::string edgeFile) {
    std::string line;
    std::ifstream edgefile(edgeFile);

    if (!edgefile) {
        std::cout << "Error opening edge file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
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

        Edge *edge = new Edge(*mesh, atol(values[0].c_str()) - 1, atol(values[1].c_str()) - 1,
                atol(values[2].c_str()) - 1, atoi(values[3].c_str()));
        mesh->edges.push_back(*edge);

        i++;
    }

    return 0;
}

int Parser::LoadCells(Mesh *mesh, std::string cellFile) {
    std::string line;
    std::ifstream cellfile(cellFile);

    if (!cellfile) {
        std::cout << "Error opening ele file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
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

        Cell *cell = new Cell(*mesh, atol(values[0].c_str()) - 1, atol(values[1].c_str()) - 1,
                atol(values[2].c_str()) - 1, atol(values[3].c_str()) - 1);
        mesh->cells.push_back(*cell);

        i++;
    }

    return 0;
}