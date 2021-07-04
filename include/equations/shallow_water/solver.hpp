#ifndef ShallowWaterSolver_hpp
#define ShallowWaterSolver_hpp

#include <functional>

#include "../virtual_solver.hpp"

class ShallowWaterSolver : public Solver {
public:
    ShallowWaterSolver(Mesh& mesh,
                       double cfl,
                       double g) :
        Solver(mesh, cfl),
        g_(g)
    {}

    void calc_tau() final;
    void process_phase_1() final;
    void process_phase_2() final;
    void process_phase_3() final;
    void prepare_next_step() final;

private:
    double g_;

    void process_phase_2_bound(const Edge& edge) final;
    void process_phase_2_inner(const Edge& edge) final;

    [[nodiscard]] double calc_div_1(double h, double u_x, double u_y, Vector div_coef) const;
    [[nodiscard]] double calc_div_2(double h, double u_x, double u_y, Vector div_coef) const;
    [[nodiscard]] double calc_div_3(double h, double u_x, double u_y, Vector div_coef) const;

    [[nodiscard]] double calc_inv_r(double G, double h, const Vector& u, const Vector& n) const;
    [[nodiscard]] double calc_inv_q(double G, double h, const Vector& u, const Vector& n) const;
    [[nodiscard]] double calc_inv_s(const Vector& u, const Vector& n) const;

    [[nodiscard]] double calc_lambda_r(double h, const Vector& u, const Vector& n) const;
    [[nodiscard]] double calc_lambda_q(double h, const Vector& u, const Vector& n) const;
    [[nodiscard]] double calc_lambda_s(const Vector& u, const Vector& n) const;

    std::array<std::array<double, 3>, 2> calc_invs(const double& G, const Vector& n,
                                                   Cell& cell, long node_id);
};

#endif
