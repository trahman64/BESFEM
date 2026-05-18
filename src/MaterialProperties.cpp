// MaterialProperties.cpp
#include "../include/MaterialProperties.hpp"
#include "../include/Constants.hpp"
#include "mfem.hpp"
#include <cmath>
#include <fstream>

namespace MaterialProperties
{
    static double GetTableValues(double x, const mfem::Vector &ticks, const mfem::Vector &data)
    {
        const int n = ticks.Size();

        if (x <= ticks(0)){return data(0);}
        if (x >= ticks(n - 1)){return data(n - 1);}

        int idx = 0;
        while (idx < n - 2 && ticks(idx + 1) < x) { idx++;}

        const double dx = ticks(idx + 1) - ticks(idx);
        return data(idx) + (x - ticks(idx)) / dx * (data(idx + 1) - data(idx));
    }
    
    static double NMC_OCV(double c)
    {
        return (1.095 * c * c) - (8.234e-7 * std::exp(14.31 * c)) + (4.692 * std::exp(-0.5389 * c));
    }

    static double NMC_mu(double c)
    {
        double val = -Constants::Frd * NMC_OCV(c);
        return val;
    }

    static double NMC_i0(double c)
    {
        double val = -0.2 * (c - 0.37) - 1.559 - 0.9376 * std::tanh(8.961 * c - 3.195);
        return std::pow(10.0, val) * 1.0e-3;
    }

    static double NMC_diff(double c)
    {
        double val = (0.0277 - 0.084 * c + 0.1003 * c * c) * 1.0e-8;
        return val;
    }

    static double LFP_i0(double c)
    {
        double val = 1.2e-5 * std::pow(c, 0.35) * std::pow(1.0 - c, 2.5);
        return val;
    }

    static double LFP_mu(double c)
    {
        static mfem::Vector Ticks(201);
        static mfem::Vector chmPot(201);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/LFP_Chm_Pot_Ticks.txt");
            std::ifstream mydFfile("../inputs/LFP_Chm_Pot.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open LFP chemical potential input files.");
            }

            for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

            loaded = true;
        }

        double val = (-1 * GetTableValues(c, Ticks, chmPot) +3.4) * -Constants::Frd;
        return val;
    }

    static double LFP_OCV(double c)
    {
        static mfem::Vector Ticks(201);
        static mfem::Vector chmPot(201);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/LFP_Chm_Pot_Ticks.txt");
            std::ifstream mydFfile("../inputs/LFP_Chm_Pot.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open LFP chemical potential input files.");
            }

            for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

            loaded = true;
        }
        
        return (-1*GetTableValues(c, Ticks, chmPot)) + 3.4;    
    }

    double LFP_ChpValue(double c)
    {
        static mfem::Vector Ticks(201);
        static mfem::Vector chmPot(201);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/LFP_Chm_Pot_Ticks.txt");
            std::ifstream mydFfile("../inputs/LFP_Chm_Pot.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open LFP chemical potential input files.");
            }

            for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, chmPot) ;
    }

    double CathodeOCV(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_OCV(c);

            case sim::MaterialType::LFP:
                return LFP_OCV(c);

            default:
                mfem::mfem_error("Unknown cathode material in CathodeOCV.");
                return 0.0;
        }
    }

    double CathodeChemicalPotential(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_mu(c);

            case sim::MaterialType::LFP:
                return LFP_mu(c);

            default:
                mfem::mfem_error("Unknown cathode material in CathodeChemicalPotential.");
                return 0.0;
        }
    }

    double CathodeExchangeCurrentDensity(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_i0(c);

            case sim::MaterialType::LFP:
                return LFP_i0(c);

            default:
                mfem::mfem_error("Unknown cathode material in CathodeExchangeCurrentDensity.");
                return 0.0;
        }
    }

    static double Electrolyte_diff(double c)
    {
        return Constants::D0 * std::exp(-7.02 - 830 * c + 50000 * c * c);
    }

    double Diffusivity(sim::MaterialType material, double c)
    {
        switch (material)
        {
            // case sim::MaterialType::Graphite:
            //     return Graphite_diff(c);

            case sim::MaterialType::NMC:
                return NMC_diff(c);

            case sim::MaterialType::LFP:
                return NMC_diff(c);

            case sim::MaterialType::Electrolyte:
                return Electrolyte_diff(c);

            default:
                mfem::mfem_error("Unknown material in Diffusivity.");
                return 0.0;
        }
    }

    static double Graphite_mu(double c)
    {
        static mfem::Vector Ticks(101);
        static mfem::Vector chmPot(101);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/C_Li_X_101.txt");
            std::ifstream mydFfile("../inputs/C_Li_M6_101.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open graphite chemical potential input files.");
            }

            for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, chmPot);
    }

    double AnodeChemicalPotential(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::Graphite:
                return Graphite_mu(c);
            default:
                mfem::mfem_error("Unknown anode material in AnodeChemicalPotential.");
                return 0.0;
        }
    }
}