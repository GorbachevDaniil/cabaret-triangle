#include "Data.hpp"
#include "Mesh.hpp"

double Data::getS0(int pos) {
    return mesh->s0[nodeID][pos];
}

double Data::getS1(int pos) {
    return mesh->s1[nodeID][pos];
}

double Data::getS2(int pos) {
    return mesh->s2[nodeID][pos];
}

Vector Data::getV0(int pos) {
    return mesh->v0[nodeID][pos];
}

Vector Data::getV1(int pos) {
    return mesh->v1[nodeID][pos];
}

Vector Data::getV2(int pos) {
    return mesh->v2[nodeID][pos];
}