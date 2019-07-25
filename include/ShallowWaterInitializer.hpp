#ifndef ShallowWaterInitializer_hpp
#define ShallowWaterInitializer_hpp

#include "AbstractInitializer.hpp"

class ShallowWaterInitializer : public AbstractInitializer {
   public:
    void initialize(Mesh &mesh);
};

#endif
