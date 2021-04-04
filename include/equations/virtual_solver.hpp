#ifndef AbstractSolver_hpp
#define AbstractSolver_hpp

#include "mesh.hpp"

class Solver {
public:
    virtual double calc_tau() = 0;
    virtual void process_phase_1(double tau) = 0;
    virtual void process_phase_2(double tau) = 0;
    virtual void process_phase_3(double tau) = 0;
    virtual void prepare_next_step() = 0;

private:
    virtual void processPhase2BoundEdge(Edge *edge, double tau) = 0;
    virtual void processPhase2InnerEdge(Edge *edge, double tau) = 0;
};

#endif
