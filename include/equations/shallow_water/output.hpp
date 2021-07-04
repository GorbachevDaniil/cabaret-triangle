#ifndef ShallowWaterOutput_hpp
#define ShallowWaterOutput_hpp

#include "../virtual_output.hpp"

class ShallowWaterOutput : public Output {
public:
    ShallowWaterOutput(bool write_conservative,
                       bool write_flux) :
        Output(write_conservative,
               write_flux) {}

    void write_paraview(Mesh& mesh, double time, int step);
};

#endif
