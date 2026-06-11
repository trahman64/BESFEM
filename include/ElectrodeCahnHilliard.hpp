#ifndef ELECTRODE_CAHN_HILLIARD_HPP
#define ELECTRODE_CAHN_HILLIARD_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"

class Initialize_Geometry;
class Domain_Parameters;

class ElectrodeCahnHilliard : public ConcentrationBase
{
public:
    ElectrodeCahnHilliard(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat, const SimulationConfig &cfg);

    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                             double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms); 

    const SimulationConfig& cfg;

private:
    mfem::ParGridFunction Mub; 
    mfem::ParGridFunction Mob; 
    mfem::ParGridFunction RxA; 

    mfem::ParMesh *pmesh;

    mfem::HypreParVector Lp1;
    mfem::HypreParVector Lp2;
    mfem::HypreParVector MuV;

    mfem::HypreParVector PsVc; 
    mfem::HypreParVector CpV0; 
    mfem::HypreParVector RHCp; 
    mfem::HypreParVector CpVn; 

    std::unique_ptr<mfem::ParBilinearForm> M_init; 
    mfem::HypreParMatrix MmatCH;                 

    std::unique_ptr<mfem::ParBilinearForm> Grad_MForm; 
    std::unique_ptr<mfem::ParBilinearForm> Grad_EForm; 

    mfem::HypreParMatrix Grad_MM;
    mfem::HypreParMatrix Grad_EM;

    std::unique_ptr<mfem::ParLinearForm> B_init;   

    mfem::HypreParVector Fcb; 
    mfem::ParLinearForm Fct; 

    mfem::GridFunctionCoefficient cAp; 
    mfem::GridFunctionCoefficient cDp; 

    std::unique_ptr<mfem::HypreParMatrix> TmatCH;

    double gtPsA = 0.0;
    double gtPsi = 0.0;

    mfem::CGSolver       MCH_solver;            
    mfem::HypreSmoother  MCH_prec;  
    
    bool combine_particle_groups = false;

};

#endif // ELECTRODE_CAHN_HILLIARD_HPP