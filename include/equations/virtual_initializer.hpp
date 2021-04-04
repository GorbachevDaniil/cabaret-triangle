#ifndef AbstractInitializer_hpp
#define AbstractInitializer_hpp

#include "mesh.hpp"

class Initializer {
   public:
    virtual void initialize(Mesh &mesh) = 0;
};

#endif
