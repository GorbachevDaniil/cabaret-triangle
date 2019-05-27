#ifndef TransferInitializer_hpp
#define TransferInitializer_hpp

#include "AbstractInitializer.hpp"

class TransferInitializer : public AbstractInitializer {
public:
    void initialize(Mesh &mesh);
};

#endif
