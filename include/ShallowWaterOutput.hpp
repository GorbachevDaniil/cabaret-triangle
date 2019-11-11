#ifndef ShallowWaterOutput_hpp
#define ShallowWaterOutput_hpp

#include "AbstractOutput.hpp"

class ShallowWaterOutput : public AbstractOutput {
   public:
    ShallowWaterOutput(int writePeriod) { this->writePeriod = writePeriod; };

    void writeParaview(Mesh *mesh, double time, int step);
};

#endif
