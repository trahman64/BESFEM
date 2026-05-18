// MaterialProperties.hpp
#pragma once

#include "SimTypes.hpp"

namespace MaterialProperties
{
    double CathodeOCV(sim::MaterialType material, double c);
    double ChemicalPotential(sim::MaterialType material, double c);
    double CathodeExchangeCurrentDensity(sim::MaterialType material, double c);
    double AnodeChemicalPotential(sim::MaterialType material, double c); 
    double CathodeKfw(sim::MaterialType material, double c);
    double CathodeKbw(sim::MaterialType material, double c);
    bool UsesDirectReactionTables(sim::MaterialType material); 
    double LFP_ChpValue(double c); 
    double Diffusivity(sim::MaterialType material, double c);
    double Mobility(sim::MaterialType material, double c);
}