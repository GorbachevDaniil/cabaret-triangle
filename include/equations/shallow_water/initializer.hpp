#ifndef ShallowWaterInitializer_hpp
#define ShallowWaterInitializer_hpp

#include "../virtual_initializer.hpp"

class ShallowWaterInitializer : public Initializer {
   public:
    void initialize(Mesh &mesh);
};

#endif
