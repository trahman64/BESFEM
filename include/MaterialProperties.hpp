// MaterialProperties.hpp
#pragma once

#include "SimTypes.hpp"

namespace MaterialProperties
{
    double OCV(sim::MaterialType material, double c);
    double ChemicalPotential(sim::MaterialType material, double c);
    double ExchangeCurrentDensity(sim::MaterialType material, double c);
    double LFP_ChpValue(double c); 
    double Diffusivity(sim::MaterialType material, double c);
    double Mobility(sim::MaterialType material, double c);
    double Conductivity(sim::MaterialType material, double c);
    double SiteDensity(sim::MaterialType material);

}