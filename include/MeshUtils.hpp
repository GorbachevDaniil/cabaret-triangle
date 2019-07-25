#ifndef MeshUtils_hpp
#define MeshUtils_hpp

#include "Mesh.hpp"

class MeshUtils {
   public:
    static void calculateEdgesNormals(Mesh &mesh);
    static void calculateTransferVectors(Mesh &mesh);
};

#endif
