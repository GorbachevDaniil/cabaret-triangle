#ifndef MeshUtils_hpp
#define MeshUtils_hpp

#include "Mesh.hpp"

class MeshUtils {
public:
    static void calculateNodeNormals(Mesh &mesh);
    static void calculateVectorsFromCenterToEdges(Mesh &mesh);
};

#endif
