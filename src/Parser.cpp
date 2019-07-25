#include "Parser.hpp"

#include "Mesh.hpp"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Make vector of values from string words
std::vector<std::string> split(const std::string &s, const char *delimiter) {
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

int Parser::LoadNodes(Mesh *mesh, std::string fileName) {
    std::string line;
    std::ifstream nodeFile(fileName);

    if (!nodeFile) {
        std::cout << "Error opening node file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
    int nNodes;
    getline(nodeFile, line);
    values = split(line, " ");
    nNodes = atoi(values[0].c_str());

    int i = 1;
    while (getline(nodeFile, line) && (i != (nNodes + 1))) {
        values = split(line, " ");
        double x = atof(values[1].c_str());
        double y = atof(values[2].c_str());
        bool boundNode = atoi(values[3].c_str()) == 1 ? true : false;
        Node *node = new Node(*mesh, x, y, mesh->apexNodesUsed, boundNode, false, true);
        mesh->nodes.push_back(*node);

        i++;
    }

    return 0;
}

int Parser::LoadEdges(Mesh *mesh, std::string fileName) {
    std::string line;
    std::ifstream edgeFile(fileName);

    if (!edgeFile) {
        std::cout << "Error opening edge file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
    int nEdges;
    getline(edgeFile, line);
    values = split(line, " ");
    nEdges = atoi(values[0].c_str());

    int i = 1;
    while (getline(edgeFile, line) && (i != (nEdges + 1))) {
        values = split(line, " ");
        unsigned long ID = atol(values[0].c_str()) - 1;
        unsigned long nodeID1 = atol(values[1].c_str()) - 1;
        unsigned long nodeID2 = atol(values[2].c_str()) - 1;
        int bound = atoi(values[3].c_str());
        Edge *edge = new Edge(*mesh, ID, nodeID1, nodeID2, bound, mesh->edgeInnerNodesNumber);
        mesh->edges.push_back(*edge);

        i++;
    }

    return 0;
}

int Parser::LoadCells(Mesh *mesh, std::string fileName) {
    std::string line;
    std::ifstream cellFile(fileName);

    if (!cellFile) {
        std::cout << "Error opening ele file!" << std::endl;
        return 1;
    }

    std::vector<std::string> values;
    int nCells;
    getline(cellFile, line);
    values = split(line, " ");
    nCells = atoi(values[0].c_str());

    int i = 1;
    while (getline(cellFile, line) && (i != (nCells + 1))) {
        values = split(line, " ");
        unsigned long ID = atol(values[0].c_str()) - 1;
        unsigned long nodeID1 = atol(values[1].c_str()) - 1;
        unsigned long nodeID2 = atol(values[2].c_str()) - 1;
        unsigned long nodeID3 = atol(values[3].c_str()) - 1;
        Cell *cell = new Cell(*mesh, ID, nodeID1, nodeID2, nodeID3);
        mesh->cells.push_back(*cell);

        i++;
    }

    return 0;
}
