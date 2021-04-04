#ifndef AbstractSolver_hpp
#define AbstractSolver_hpp

#include "grid/mesh.hpp"

class Solver {
public:
    Solver(Mesh& mesh,
           double cfl) :
        mesh_(mesh),
        cfl_(cfl)
    {}

    virtual void calc_tau() = 0;
    virtual void process_phase_1() = 0;
    virtual void process_phase_2() = 0;
    virtual void process_phase_3() = 0;
    virtual void prepare_next_step() = 0;

    double get_tau() const {
        return tau_;
    }

protected:
    Mesh mesh_;
    double cfl_;
    double tau_;

private:
    virtual void process_phase_2_bound(Edge *edge) = 0;
    virtual void process_phase_2_inner(Edge *edge) = 0;
};

#endif
