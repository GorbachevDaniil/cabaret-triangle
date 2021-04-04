#ifndef TransferOutput_hpp
#define TransferOutput_hpp

#include "../virtual_output.hpp"

class TransferOutput : public AbstractOutput {
public:
    TransferOutput(int writePeriod) { this->writePeriod = writePeriod; };

    void writeParaview(Mesh *mesh, double time, int step);
};

#endif
