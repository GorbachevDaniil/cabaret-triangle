#ifndef MeshUtils_hpp
#define MeshUtils_hpp

#include "mesh.hpp"

class MeshUtils {
public:
    static void calculateEdgesNormals(Mesh &mesh);
    static void calculateTransferVectors(Mesh &mesh);
};

#endif
