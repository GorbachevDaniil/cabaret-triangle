#ifndef ShallowWaterOutput_hpp
#define ShallowWaterOutput_hpp

#include "AbstractOutput.hpp"

class ShallowWaterOutput : public AbstractOutput {
   public:
    ShallowWaterOutput(int writePeriod, bool writeConservative, bool writeFlux) {
        this->writePeriod = writePeriod;
        this->writeConservative = writeConservative;
        this->writeFlux = writeFlux;
    };

    void writeParaview(Mesh *mesh, double time, int step);
};

#endif
