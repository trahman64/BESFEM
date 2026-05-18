#ifndef ELECTROLYTE_DIFFUSION_HPP
#define ELECTROLYTE_DIFFUSION_HPP

#include "Utils.hpp"
#include "Concentrations_Base.hpp"

class Initialize_Geometry;
class Domain_Parameters;

class ElectrolyteDiffusion : public ConcentrationBase
{ 
public:
    ElectrolyteDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode, sim::MaterialType mat);

    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                             double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms);

    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) override;

    BoundaryConditions  &boundary_conditions;
    sim::CellMode mode_;

    mfem::HypreParVector CeVn; 
    mfem::ParLinearForm Fet;                 
    mfem::Array<int>    nbc_bdr;             
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; 
    mfem::ConstantCoefficient nbcCoef;       
    mfem::GridFunctionCoefficient matCoef_R; 
    mfem::Array<int> boundary_dofs;          


private:

    std::unique_ptr<mfem::ParBilinearForm> Me_init; 
    std::unique_ptr<mfem::ParLinearForm>   Be_init; 

    mfem::HypreParMatrix Mmate; 
    mfem::HypreParMatrix Kmate; 
    std::unique_ptr<mfem::ParBilinearForm> Ke2; 

    mfem::CGSolver      Me_solver; 
    mfem::HypreSmoother Me_prec;   

    mfem::GridFunctionCoefficient cAe; 
    mfem::GridFunctionCoefficient cDe; 

    mfem::ParGridFunction De;  
    mfem::ParGridFunction Rxe; 
    mfem::ParGridFunction PeR; 

    double eCrnt = 0.0; 
    double infx  = 0.0; 
    double L_w   = 0.0; 

    mfem::HypreParVector Feb; 
    mfem::HypreParVector X1v; 

    std::unique_ptr<mfem::HypreParMatrix> TmatR; 
    std::unique_ptr<mfem::HypreParMatrix> TmatL; 

    mfem::HypreParVector CeV0; 
    mfem::HypreParVector RHCe; 
    mfem::HypreParVector PsVc; 

};

#endif // ELECTROLYTE_DIFFUSION_HPP