#ifndef ELECTRODE_POTENTIAL_HPP
#define ELECTRODE_POTENTIAL_HPP

#include "Utils.hpp"
#include "Potentials_Base.hpp"
#include "SimTypes.hpp"

class Initialize_Geometry;
class Domain_Parameters;
class BoundaryConditions;

class ElectrodePotential : public PotentialBase
{
public:
    ElectrodePotential(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::Electrode electrode, sim::MaterialType material);

    void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx);
    void AssembleSystem(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, mfem::ParGridFunction &potential);
    void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

    void AssembleSystem(mfem::ParGridFunction &Cn,mfem::ParGridFunction &psx,mfem::ParGridFunction &potential) override;

    double GetBoundaryVoltage() const override { return Bv; }
    void AddBoundaryVoltage(double dV){Bv += dV;}

private:

    Initialize_Geometry  &geometry;            ///< Geometry and mesh infrastructure
    Domain_Parameters    &domain_parameters;   ///< Material/phase-field parameters
    BoundaryConditions   &boundary_conditions; ///< Electrode boundary configuration

    FEMOperators fem; ///< FEM operator assembly utilities
    Utils utils;      ///< Utility helpers (clamping, reductions, etc.)

    sim::Electrode electrode;
    sim::MaterialType material;

    mfem::ParGridFunction kap; // Conductivity field
    mfem::GridFunctionCoefficient cKp; // Coefficient for conductivity
    mfem::ParGridFunction RpP; // Reaction term scaled by Frd
    mfem::GridFunctionCoefficient cRp; // Coefficient for reaction term

    std::unique_ptr<mfem::ParLinearForm> Bp2; // Linear form for the reaction term
    std::unique_ptr<mfem::ParBilinearForm> Kp2; // Bilinear form for the conductivity

    mfem::ParGridFunction pP0; // Previous potential field

    mfem::ParLinearForm B1t; // Temporary linear form for assembling
    mfem::HypreParVector X1v; // Temporary solution vector
    mfem::HypreParVector B1v; // Temporary RHS vector
    mfem::HypreParVector Fpb; // Reaction contribution vector
    mfem::HypreParVector Xs0; // Previous solution vector (true DO
    mfem::ParLinearForm Fpt; // Linear form for the reaction term in the potential update
    mfem::HypreParMatrix KmP; // System matrix for the potential update
    std::unique_ptr<mfem::HypreBoomerAMG> Mpp; // AMG preconditioner for the potential update

    mfem::Array<int> ess_tdof_list; // List of essential DOFs for boundary conditions

    double Bv; // Boundary value for potential
    double gtPs; // Global scaling factor for conductivity

    mfem::CGSolver cgPP_solver; // Conjugate gradient solver for potential update

    void ParticleConductivityMulti(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups);
};

#endif // ELECTRODE_POTENTIAL_HPP