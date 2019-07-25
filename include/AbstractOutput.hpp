#ifndef AbstractOutput_hpp
#define AbstractOutput_hpp

#include "Mesh.hpp"

class AbstractOutput {
   protected:
    int writePeriod;

   public:
    virtual void writeParaview(Mesh *mesh, int step) = 0;
};

#endif
