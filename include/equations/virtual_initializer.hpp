#ifndef AbstractInitializer_hpp
#define AbstractInitializer_hpp

#include "grid/mesh.hpp"

class Initializer {
public:
    virtual void initialize(Mesh &mesh) = 0;
};

#endif
