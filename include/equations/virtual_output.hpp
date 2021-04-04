#ifndef AbstractOutput_hpp
#define AbstractOutput_hpp

#include "mesh.hpp"

class Output {
   protected:
    int writePeriod;
    bool writeConservative;
    bool writeFlux;

   public:
    virtual void writeParaview(Mesh *mesh, double time, int step) = 0;
};

#endif
