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

    bool UsesDirectReactionTables(sim::MaterialType material)
    {
        // return material == sim::MaterialType::LFP;
        return false;
    }
    
    static double NMC_OCV(double c)
    {
        return (1.095 * c * c) - (8.234e-7 * std::exp(14.31 * c)) + (4.692 * std::exp(-0.5389 * c));
    }

    static double NMC_mu(double c)
    {
        return -Constants::Frd * NMC_OCV(c);
    }

    static double NMC_i0(double c)
    {
        double val = -0.2 * (c - 0.37) - 1.559 - 0.9376 * std::tanh(8.961 * c - 3.195);
        return std::pow(10.0, val) * 1.0e-3;
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

        // return -Constants::Frd * (GetTableValues(c, Ticks, chmPot) + 3.4);
        return (-1 * GetTableValues(c, Ticks, chmPot) +3.4) * -Constants::Frd;
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


    static double LFP_Kfw(double c)
    {
        static mfem::Vector WmuAxis(201);
        static mfem::Vector KfLog(201);
        static bool loaded = false;
        static bool printed = false;

        if (!loaded)
        {
            std::ifstream file("../inputs/LFP_RCAKfL_201_A.txt");

            if (!file)
            {
                mfem::mfem_error("Could not open LFP_RCAKfL_201_A.txt");
            }

            for (int i = 0; i < 201; i++)
            {
                file >> WmuAxis(i);
                file >> KfLog(i);
            }

            loaded = true;
        }

        // ===== Fortran equivalent =====
        const double dfdX = LFP_ChpValue(c);

        const double LpC = dfdX * Constants::Cst1;

        // equivalent to:
        // Wmu = (LpC - RxTbl(1,1)) / (RxTbl(1,RxI)-RxTbl(1,1))

        const double RxMin = WmuAxis(0);
        const double RxMax = WmuAxis(WmuAxis.Size() - 1);

        // std::cout << WmuAxis(0) << " " << WmuAxis(200) << std::endl;

        double Wmu = (LpC - RxMin) / (RxMax - RxMin);

        // optional clamp
        // Wmu = std::max(0.0, std::min(1.0, Wmu));

        // interpolate ln(Kfw)
        const double logK =
            GetTableValues(Wmu, WmuAxis, KfLog);

        if (!printed && mfem::Mpi::WorldRank() == 0)
        {
            std::cout
                << "LFP DEBUG: "
                << "c = " << c
                << " dfdX = " << dfdX
                << " LpC = " << LpC
                << " Wmu = " << Wmu
                << " logK = " << logK
                << std::endl;

            printed = true;
        }

        // std::cout << "Kfw Wmu Axis: "<< WmuAxis(0) << " " << WmuAxis(200) << std::endl;


        return std::exp(logK);
    }

    static double LFP_Kbw(double c)
    {
        static mfem::Vector WmuAxis(201);
        static mfem::Vector KbLog(201);
        static bool loaded = false;
        static bool printed = false;

        if (!loaded)
        {
            std::ifstream file("../inputs/LFP_RCAKbL_201_A.txt");

            if (!file)
            {
                mfem::mfem_error("Could not open LFP_RCAKbL_201_A.txt");
            }

            for (int i = 0; i < 201; i++)
            {
                file >> WmuAxis(i);
                file >> KbLog(i);
            }

            loaded = true;
        }

        const double dfdX = LFP_ChpValue(c);
        const double LpC = dfdX * Constants::Cst1;

        const double RxMin = WmuAxis(0);
        const double RxMax = WmuAxis(WmuAxis.Size() - 1);

        double Wmu = (LpC - RxMin) / (RxMax - RxMin);

        const double logK =
            GetTableValues(Wmu, WmuAxis, KbLog);

        if (!printed && mfem::Mpi::WorldRank() == 0)
        {
            std::cout
                << "LFP DEBUG: "
                << "c = " << c
                << " dfdX = " << dfdX
                << " LpC = " << LpC
                << " Wmu = " << Wmu
                << " logK = " << logK
                << std::endl;

            printed = true;
        }

        // std::cout << "Kbw Wmu Axis: "<< WmuAxis(0) << " " << WmuAxis(200) << std::endl;

        return std::exp(logK);
    }

    double CathodeKfw(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::LFP:
                return LFP_Kfw(c);

            default:
                mfem::mfem_error("CathodeKfw only implemented for table-based materials.");
                return 0.0;
        }
    }

    double CathodeKbw(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::LFP:
                return LFP_Kbw(c);

            default:
                mfem::mfem_error("CathodeKbw only implemented for table-based materials.");
                return 0.0;
        }
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

            // case sim::MaterialType::LFP:
            //     // placeholder until LFP kinetics are added
            //     return NMC_i0(c);

            default:
                mfem::mfem_error("Unknown cathode material in CathodeExchangeCurrentDensity.");
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