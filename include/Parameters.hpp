#ifndef Parameters_hpp
#define Parameters_hpp

class Parameters {
public:
    static constexpr double CFL = 0.3;
    static constexpr int    WRITE_STEP_PERIOD = 50;
    static constexpr int    HARMONIZATION_STEP = 1000000000;
};

#endif
