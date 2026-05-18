#ifndef ELECTRODE_DIFFUSION_HPP
#define ELECTRODE_DIFFUSION_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"

class Initialize_Geometry;
class Domain_Parameters;

class ElectrodeDiffusion : public ConcentrationBase
{
public:
    ElectrodeDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat);

    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                             double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms);

private:

    mfem::ParGridFunction Dp;
    mfem::ParGridFunction Rxn;

    mfem::HypreParVector PsVc; 
    mfem::HypreParVector CpV0; 
    mfem::HypreParVector RHCp; 
    mfem::HypreParVector CpVn; 

    std::unique_ptr<mfem::ParBilinearForm> Mt; 
    mfem::HypreParMatrix Mmatp;                

    mfem::CGSolver       Mp_solver;            
    mfem::HypreSmoother  Mp_prec;              

    std::unique_ptr<mfem::ParBilinearForm> Kc2; 
    mfem::HypreParMatrix Kmatp;                 

    std::unique_ptr<mfem::ParLinearForm> Bc2;   

    mfem::HypreParVector Fcb; 
    mfem::ParLinearForm Fct; 

    mfem::GridFunctionCoefficient cAp; 
    mfem::GridFunctionCoefficient cDp; 

    std::unique_ptr<mfem::HypreParMatrix> Tmatp;

    double gtPsC = 0.0;
    double gtPsi = 0.0;
};

#endif // ELECTRODE_DIFFUSION_HPP