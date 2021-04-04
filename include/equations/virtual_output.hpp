#ifndef AbstractOutput_hpp
#define AbstractOutput_hpp

#include "grid/mesh.hpp"

class Output {
public:
    Output(int write_period,
           bool write_conservative,
           bool write_flux) :
        write_period_(write_period),
        write_conservative_(write_conservative),
        write_flux_(write_flux) {}

    virtual void write_paraview(Mesh& mesh, double time, int step) = 0;

protected:
    const int write_period_;
    const bool write_conservative_;
    const bool write_flux_;
};

#endif
