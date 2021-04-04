#include <chrono>
#include <iostream>

//#include <libconfig.h++>

#include "mesh.hpp"
#include "equations/shallow_water/initializer.hpp"
#include "equations/shallow_water/output.hpp"
#include "equations/shallow_water/solver.hpp"

int main() {
    std::cout.precision(8);

//    libconfig::Config config;
//    config.readFile("config/config.cfg");

//    int steps = config.lookup("steps");
//    double timeToStop = config.lookup("time_to_stop");
//
//    int writePeriod = config.lookup("write_period");
//    bool writeConservative = config.lookup("write_conservative");
//    bool writeFlux = config.lookup("write_flux");
//
//    int edgeInnerNodesNumber = config.lookup("edge_inner_nodes_number");
//    bool apexNodesUsed = config.lookup("apex_nodes_used");

    int steps = 10;
    double timeToStop = 0.0;

    int writePeriod = 1;
    bool writeConservative = true;
    bool writeFlux = true;

    int edgeInnerNodesNumber = 2;
    bool apexNodesUsed = false;

    Mesh *mesh = new Mesh(edgeInnerNodesNumber, apexNodesUsed);
    mesh->init_mesh(
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/Mesh.node",
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/Mesh.edge",
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/Mesh.ele"
    );
    std::cout << "number of cells = " << mesh->cells.size() << std::endl;

    mesh->calculate_edges_normals();
    mesh->calculate_transfer_vectors();

    int task = 2;
    switch (task) {
        // case 1: {
        //     TransferInitializer initializer = TransferInitializer();
        //     initializer.initialize(*mesh);

        //     TransferOutput output = TransferOutput(writePeriod);
        //     output.writeParaview(mesh, 0);

        //     double cfl = config.lookup("solver.cfl");
        //     TransferSolver solver = TransferSolver(cfl, mesh);
        //     double time = 0;
        //     double tau = solver.calcTau();

        //     for (int i = 0; i < steps; i++) {
        //         if (timeToStop > 0 && time >= timeToStop) {
        //             output.writeParaview(mesh, i + 1);
        //             break;
        //         }

        //         std::cout << "step = " << i << " tau = " << tau << " time = " << time << std::endl;
        //         time += tau;

        //         solver.processPhase1(tau);
        //         solver.processPhase2(tau);
        //         solver.processPhase3(tau);
        //         solver.prepareNextStep();

        //         output.writeParaview(mesh, i + 1);
        //     }
        //     break;
        // }

        case 2: {
            ShallowWaterInitializer initializer = ShallowWaterInitializer();
            initializer.initialize(*mesh);

            ShallowWaterOutput output = ShallowWaterOutput(writePeriod,
                                                           writeConservative,
                                                           writeFlux);
            output.write_paraview(mesh, 0, 0);

//            double cfl = config.lookup("solver.cfl");
//            double g = config.lookup("solver.shallow_water.g");
            double cfl = 0.3;
            double g = 1.0;
            ShallowWaterSolver solver = ShallowWaterSolver(cfl, g, mesh);

            std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

            std::chrono::duration<double> durCalcTau = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase1 = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase2 = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase3 = std::chrono::duration<double>(0);

            double time = 0;
            int real_steps = 0;
            for (int i = 0; i < steps; i++) {
                std::chrono::system_clock::time_point startCalcTau = std::chrono::system_clock::now();
                double tau = solver.calc_tau();
                std::cout << "step = " << i + 1 << " tau = " << tau << " time = " << time << std::endl;
                time += tau;
                std::chrono::system_clock::time_point stopCalcTau = std::chrono::system_clock::now();
                durCalcTau +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopCalcTau - startCalcTau);

                std::chrono::system_clock::time_point startPhase1 = std::chrono::system_clock::now();
                solver.process_phase_1(tau);
                std::chrono::system_clock::time_point stopPhase1 = std::chrono::system_clock::now();
                durPhase1 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase1 - startPhase1);

                std::chrono::system_clock::time_point startPhase2 = std::chrono::system_clock::now();
                solver.process_phase_2(tau);
                std::chrono::system_clock::time_point stopPhase2 = std::chrono::system_clock::now();
                durPhase2 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase2 - startPhase2);

                std::chrono::system_clock::time_point startPhase3 = std::chrono::system_clock::now();
                solver.process_phase_3(tau);
                std::chrono::system_clock::time_point stopPhase3 = std::chrono::system_clock::now();
                durPhase3 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase3 - startPhase3);

                solver.prepare_next_step();

                real_steps += 1;

                if (timeToStop > 0 && time >= timeToStop) {
                    output.write_paraview(mesh, time, i + 1);
                    break;
                }
                output.write_paraview(mesh, time, i + 1);
            }
            std::chrono::system_clock::time_point stop = std::chrono::system_clock::now();
            std::chrono::duration<double> dur =
                std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);

            std::cout << "it took " << dur.count() / real_steps << " sec to process one step" << std::endl;
            std::cout << "it took " << durCalcTau.count() / real_steps << " sec to calculate tau" << std::endl;
            std::cout << "it took " << durPhase1.count() / real_steps << " sec to process phase 1" << std::endl;
            std::cout << "it took " << durPhase2.count() / real_steps << " sec to process phase 2" << std::endl;
            std::cout << "it took " << durPhase3.count() / real_steps << " sec to process phase 3" << std::endl;
            break;
        }
    }

    return 0;
}
