#ifndef Parameters_hpp
#define Parameters_hpp

class Parameters {
public:
    static constexpr double CFL = 0.3;
    static constexpr int    STEPS = 400;
    static constexpr int    WRITE_STEP_PERIOD = 5;
    static constexpr int    EDGE_INNER_NODES_NUMBER = 2;
};

#endif
