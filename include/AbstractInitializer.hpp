#ifndef AbstractInitializer_hpp
#define AbstractInitializer_hpp

#include "Mesh.hpp"

class AbstractInitializer {
   public:
    virtual void initialize(Mesh &mesh) = 0;
};

#endif
