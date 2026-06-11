#include "mfem.hpp"
#include "../include/BESFEM_All.hpp"
#include "../include/SimulationConfig.hpp"


Adjust::Adjust(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), cfg(cfg)

{}


// Full Cell
void Adjust::AdjustConstantCurrent(double current_A, double current_C, ElectrodePotential &anode_potential, ElectrodePotential &cathode_potential,
    mfem::ParGridFunction &phA_gf, mfem::ParGridFunction &phC_gf, double &VCell)
{
    double Vsr;
    double dCrnt = std::abs(current_A - domain_parameters.gTrgI);

    // --- Current deviation scaling ---
    if (dCrnt < std::abs(domain_parameters.gTrgI) * 0.05)
        Vsr = 0.025 * cfg.Vsr0;
    else if (dCrnt < std::abs(domain_parameters.gTrgI) * 0.10)
        Vsr = 0.25 * cfg.Vsr0;
    else
        Vsr = 1.0 * cfg.Vsr0;

    // --- Anode voltage correction ---
    double sgnA = std::copysign(1.0, domain_parameters.gTrgI - std::abs(current_A));
    double dV_A = cfg.dt * Vsr * sgnA * 0.10;
    anode_potential.AddBoundaryVoltage(dV_A);
    phA_gf += dV_A;

    // --- Cathode voltage correction ---
    double sgnC = std::copysign(1.0, domain_parameters.gTrgI - current_C);
    double dV_C = cfg.dt * Vsr * sgnC * 2.0;
    cathode_potential.AddBoundaryVoltage(-dV_C);
    phC_gf -= dV_C;

    // --- Compute overall cell voltage ---
    VCell = cathode_potential.GetBoundaryVoltage() - anode_potential.GetBoundaryVoltage();
}



