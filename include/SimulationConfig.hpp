#pragma once

#include "SimTypes.hpp"
#include "Constants.hpp"

#include <string>
#include <vector>

struct SimulationConfig
{
    sim::CellMode mode = sim::CellMode::HALF;
    sim::Electrode half_electrode = sim::Electrode::ANODE;

    const char* config_file = "../inputs/run_config.txt";

    const char* mesh_file  = "../inputs/colored_labels_labels.tif";   ///< Default mesh file path.
    const char* dsF_file_A = "../inputs/dummy.gf";  ///< Distance-function file for anode.
    const char* dsF_file_C = "../inputs/dummy.gf";  ///< Distance-function file for cathode.
    const char* mesh_type  = "v";

    int order = Constants::order;
    int num_timesteps = 1000;

    bool combine_particle_groups = false;

    std::vector<double> init_anode_particles;
    std::vector<double> init_cathode_particles;

    std::vector<sim::MaterialType> anode_materials;
    std::vector<sim::MaterialType> cathode_materials;

    double init_CnA = 0.95;
    double init_CnC = 0.30;
    double init_CnE = 0.001;
    double init_BvA = -0.01;
    double init_BvC = 3.30;
    double init_BvE = -0.10;
};

SimulationConfig ParseSimulationArgs(int argc, char *argv[]);
void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[]);
void PrintAvailableSimulationOptions();