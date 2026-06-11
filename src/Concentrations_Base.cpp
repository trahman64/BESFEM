#include "../include/Concentrations_Base.hpp"
#include "../include/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/SimulationConfig.hpp"
#include "mfem.hpp"

ConcentrationBase::ConcentrationBase(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat, const SimulationConfig &cfg)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), cfg(cfg), EVol(para.EVol), gtPsi(para.gtPsi), gtPse(para.gtPse),
    nE(geo.nE), nC(geo.nC), nV(geo.nV), VtxVal(geo.nC), EAvg(geo.nE), gmesh(geo.globalMesh.get()), material(mat), fem(geo.parfespace), utils(geo, para, cfg)

{
    // Allocate once for reuse
    VtxVal.SetSize(nC); // Set size for nodal values
    EAvg.SetSize(nE);   // Set size for average element contributions
    TmpF = std::make_unique<mfem::ParGridFunction>(fespace.get());

}
