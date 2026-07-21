#include "../include/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <sstream>

static inline void GlobalMinMax(const mfem::ParGridFunction& gf,
                                double& gmin, double& gmax,
                                MPI_Comm comm = MPI_COMM_WORLD)
{
    double lmin =  std::numeric_limits<double>::infinity();
    double lmax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < gf.Size(); ++i) {
        const double v = gf(i);
        if (v < lmin) lmin = v;
        if (v > lmax) lmax = v;
    }
    MPI_Allreduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, comm);
}

double gTrgI = 0.0;

Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo, const SimulationConfig &cfg)
    : geometry(geo), cfg(cfg), nV(geo.nV), nE(geo.nE), nC(geo.nC), 
    pmesh(geo.parallelMesh.get()), fespace(geo.parfespace),
    particle_labels(geo.particle_labels)
{}

// Destructor
Domain_Parameters::~Domain_Parameters() {}

void Domain_Parameters::SetupDomainParameters(){

    InitializeGridFunctions();
    InterpolateDomainParameters();
    CalculatePhasePotentialsAndTargetCurrent();

    psi->SaveAsOne("psi");
    pse->SaveAsOne("pse");
    AvP->SaveAsOne("AvP");
    AvE->SaveAsOne("AvE");
    AvB->SaveAsOne("AvB");
    pmesh->SaveAsOne("pmesh");

    PrintInfo();
}

void Domain_Parameters::InitializeGridFunctions() {

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }

    psi = std::make_unique<mfem::ParGridFunction>(fespace.get());
    pse = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvE = std::make_unique<mfem::ParGridFunction>(fespace.get());
    denom = std::make_unique<mfem::ParGridFunction>(fespace.get());

    ps.clear();
    ps.resize(particle_labels.size());

    AvPs.clear();
    AvPs.resize(particle_labels.size());

    AvEs.clear();
    AvEs.resize(particle_labels.size());

    WeightEs.clear();
    WeightEs.resize(particle_labels.size());

    for (int k = 0; k < (int)particle_labels.size(); ++k)
    {
        ps[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvPs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvEs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        WeightEs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

    AvP_Pairs.clear();
    AvP_Pairs.resize(particle_labels.size());
    psi_Pairs.clear();
    psi_Pairs.resize(particle_labels.size());
    WeightPairs.clear();
    WeightPairs.resize(particle_labels.size());

    for (int j = 0; j < (int)particle_labels.size(); ++j)
    {
        AvP_Pairs[j].resize(particle_labels.size());
        psi_Pairs[j].resize(particle_labels.size());
        WeightPairs[j].resize(particle_labels.size());

        for (int k = 0; k < (int)particle_labels.size(); ++k)
        {
            if (k != j)
            {
                AvP_Pairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
                psi_Pairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
                WeightPairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            }
        }
    }

    tPs.clear();
    tPs.resize(particle_labels.size());

    gtPs.clear();
    gtPs.resize(particle_labels.size());

    gTrgPs.clear();
    gTrgPs.resize(particle_labels.size());

    for (int k = 0; k < (int)particle_labels.size(); ++k)
    { 
        tPs[k] = 0.0; 
        gtPs[k] = 0.0;
        gTrgPs[k] = 0.0;
    }

}

void Domain_Parameters::InterpolateDomainParameters() {

    nV = pmesh->GetNV();

    MFEM_VERIFY(pmesh, "Parallel mesh is not initialized.");
    MFEM_VERIFY(geometry.MaskFilterPse, "Electrolyte mask is not initialized.");
    MFEM_VERIFY(static_cast<int>(geometry.MaskFilters.size()) == static_cast<int>(ps.size()),
                "Particle mask count does not match particle field count.");

    *pse = *geometry.MaskFilterPse;
    *psi = 0.0;

    for (int k = 0; k < (int)ps.size(); ++k)
    {
        MFEM_VERIFY(geometry.MaskFilters[k], "A particle mask is not initialized.");
        *ps[k] = *geometry.MaskFilters[k];
        *psi += *ps[k];
    }

    for (int i = 0; i < psi->Size(); i++)
    {
        if ((*psi)(i) < 0.0) { (*psi)(i) = 0.0; }
        if ((*psi)(i) > 1.0) { (*psi)(i) = 1.0; }

        if ((*pse)(i) < 0.0) { (*pse)(i) = 0.0; }
        if ((*pse)(i) > 1.0) { (*pse)(i) = 1.0; }

        (*psi)(i) += 1.0e-6;
        (*pse)(i) += 1.0e-6;
    }

    for (int k = 0; k < (int)ps.size(); ++k)
    {
        for (int i = 0; i < ps[k]->Size(); ++i)
        {
            if ((*ps[k])(i) < 0.0) { (*ps[k])(i) = 0.0; }
            if ((*ps[k])(i) > 1.0) { (*ps[k])(i) = 1.0; }
            (*ps[k])(i) += 1.0e-6;
        }
    }
    // -------------------------------------------------
    // AMR SECTION
    // -------------------------------------------------
    if (cfg.amr_levels > 0)
    {
        const double outer_half_width = 0.45; // half width of band

        for (int lev = 0; lev < cfg.amr_levels; ++lev)
        {
            mfem::Array<int> refinement_list;

            const double band_fraction = static_cast<double>(cfg.amr_levels - lev) / static_cast<double>(cfg.amr_levels);
            const double half_width = outer_half_width * band_fraction;

            // widest possible band is 0.05 to 0.95
            const double psi_lower = 0.5 - half_width;
            const double psi_upper = 0.5 + half_width;

            for (int ei = 0; ei < pmesh->GetNE(); ++ei)
            {
                mfem::Array<double> psi_vals;
                psi->GetNodalValues(ei, psi_vals);

                double psi_avg = 0.0;

                for (int j = 0; j < psi_vals.Size(); ++j)
                {
                    psi_avg += psi_vals[j];
                }

                psi_avg /= static_cast<double>(psi_vals.Size());

                if (psi_avg > psi_lower && psi_avg < psi_upper)
                {
                    refinement_list.Append(ei);
                }
            }

            const int local_marked = refinement_list.Size();
            const int local_elements = pmesh->GetNE();

            int global_marked = 0;
            int global_elements = 0;

            MPI_Allreduce(&local_marked, &global_marked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout << "[AMR] band " << lev + 1
                    << ": psi range = (" << psi_lower
                    << ", " << psi_upper << ")"
                    << ", marked " << global_marked
                    << " / " << global_elements
                    << " elements globally"
                    << std::endl;
            }

            if (global_marked == 0)
            {
                if (mfem::Mpi::WorldRank() == 0)
                {
                    std::cout
                        << "[AMR] No elements found in band "
                        << lev + 1 << ". Stopping refinement."
                        << std::endl;
                }

                break;
            }

            pmesh->GeneralRefinement(refinement_list, 1);
            fespace->Update();

            // Find the element-size range after refining this band.
            double local_hmin = std::numeric_limits<double>::max();
            double local_hmax = 0.0;

            for (int ei = 0; ei < pmesh->GetNE(); ++ei)
            {
                const double h = pmesh->GetElementSize(ei);

                local_hmin = std::min(local_hmin, h);
                local_hmax = std::max(local_hmax, h);
            }

            double global_hmin = 0.0;
            double global_hmax = 0.0;

            MPI_Allreduce(&local_hmin, &global_hmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&local_hmax, &global_hmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout
                    << "[AMR] band " << lev + 1
                    << " element size after refinement: "
                    << global_hmin << std::endl;
            }

            auto UpdateGF = [](std::unique_ptr<mfem::ParGridFunction> &gf)
            {
                if (gf)
                {
                    gf->Update();
                }
            };

            UpdateGF(psi);
            UpdateGF(pse);
            UpdateGF(AvP);
            UpdateGF(AvB);
            UpdateGF(AvE);
            UpdateGF(denom);

            for (int k = 0; k < static_cast<int>(ps.size()); ++k)
            {
                UpdateGF(ps[k]);
                UpdateGF(AvPs[k]);
                UpdateGF(AvEs[k]);
                UpdateGF(WeightEs[k]);
            }

            for (int j = 0; j < static_cast<int>(ps.size()); ++j)
            {
                for (int k = 0; k < static_cast<int>(ps.size()); ++k)
                {
                    if (j != k)
                    {
                        UpdateGF(AvP_Pairs[j][k]);
                        UpdateGF(psi_Pairs[j][k]);
                        UpdateGF(WeightPairs[j][k]);
                    }
                }
            }

            nV = pmesh->GetNV();
            nE = pmesh->GetNE();
            nC = pmesh->GetElement(0)->GetNVertices();

            geometry.nV = nV;
            geometry.nE = nE;
            geometry.nC = nC;

            int local_elements_after = pmesh->GetNE();
            int global_elements_after = 0;

            MPI_Allreduce(&local_elements_after, &global_elements_after, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (mfem::Mpi::WorldRank() == 0)
            {
                std::cout
                    << "[AMR] band " << lev + 1
                    << " complete: "
                    << global_elements_after
                    << " total elements"
                    << std::endl;
            }
        }
    }

    // -------------------------------------------------
    // MULTI PARTICLE AvP
    // -------------------------------------------------

    auto ComputeGradMagnitude = [&](const mfem::ParGridFunction &phase_in,
                                    mfem::ParGridFunction &AvP_out)
    {
        const int dim = pmesh->Dimension();
        mfem::ParGridFunction dphase(fespace.get());

        AvP_out = 0.0;
        for (int d = 0; d < dim; ++d)
        {
            dphase = 0.0;
            mfem::ParGridFunction phase_tmp(phase_in);
            phase_tmp.GetDerivative(1, d, dphase);

            for (int vi = 0; vi < nV; ++vi)
            {
                const double v = dphase(vi);
                AvP_out(vi) += v * v;
            }
        }

        for (int vi = 0; vi < nV; ++vi)
        {
            AvP_out(vi) = std::sqrt(AvP_out(vi));
        }
    };



    auto BuildPairInterface = [&](mfem::ParGridFunction &out,
                        const mfem::ParGridFunction &psa,
                        const mfem::ParGridFunction &psb,
                        const mfem::ParGridFunction &AvPa,
                        const mfem::ParGridFunction &AvPb)
    {
        out = psa;      // psa
        out *= AvPb;    // psa * AvPb

        mfem::ParGridFunction tmp(fespace.get());
        tmp = psb;      // psb
        tmp *= AvPa;    // psb * AvPa
        out += tmp;     // psa*AvPb + psb*AvPa

        mfem::ParGridFunction overlap(fespace.get());
        overlap = psa;
        overlap *= psb; // psa * psb

        out *= overlap;
        out *= 4.0;

        for (int vi = 0; vi < out.Size(); ++vi)
        {
            if (out(vi) > 9000.0)
            {
                out(vi) = 1.4e4;
            }
        }
    };

    auto BuildElectrolyteInterface = [&](mfem::ParGridFunction &out,
                                const mfem::ParGridFunction &pse_in,
                                const mfem::ParGridFunction &AvPk)
    {
        out = pse_in;
        out *= AvPk;
    };

    
    auto ComputeWeight = [&](mfem::ParGridFunction &weight_out,
                    const mfem::ParGridFunction &num_in,
                    const mfem::ParGridFunction *mask_in = nullptr)
    {
        weight_out = 0.0;

        const double beta = 0.8;
        const double eps  = 1e-30;

        for (int vi = 0; vi < nV; ++vi)
        {
            const double num = num_in(vi);
            const double den = (*denom)(vi);

            double ratio = 0.0;
            if (den > eps)
            {
                ratio = num / den;
                if (ratio < 0.0) ratio = 0.0;
            }

            weight_out(vi) = std::pow(ratio, beta);
        }

        if (mask_in) { weight_out *= *mask_in; }
    };

    auto BuildPsiPairs = [&](mfem::ParGridFunction &out,
                        const mfem::ParGridFunction &psa,
                        const mfem::ParGridFunction &psb)
    {
        out = psa;
        out += psb;

        for (int vi = 0; vi < nV; ++vi)
        {
            if (out(vi) > 1.0) { out(vi) = 1.0; }
        }
    };


    // -------------------------------------------------
    //  Solving
    // -------------------------------------------------

    // Building AvPs
    ComputeGradMagnitude(*psi, *AvP);      // total solid
    ComputeGradMagnitude(*pse, *AvE);      // electrolyte

    for (int k = 0; k < (int)ps.size(); ++k)
    {
        ComputeGradMagnitude(*ps[k], *AvPs[k]);
    }

    // Building Pairs
    for (int j = 0; j < (int)ps.size(); ++j)
    {
        for (int k = j + 1; k < (int)ps.size(); ++k)
        {
            BuildPairInterface(*AvP_Pairs[j][k], *ps[j], *ps[k], *AvPs[j], *AvPs[k]);
            BuildPsiPairs(*psi_Pairs[j][k], *ps[j], *ps[k]);
        }
    }

    // Building AvEs 
    for (int k = 0; k < (int)ps.size(); ++k)
    {
        BuildElectrolyteInterface(*AvEs[k], *pse, *AvPs[k]);
    }

    // Building denominator for weights
    *denom = 0.0;
    for (int j = 0; j < (int)ps.size(); ++j)
    {
        for (int k = j + 1; k < (int)ps.size(); ++k)
        {
            *denom += *AvP_Pairs[j][k];
        }
    }
    for (int k = 0; k < (int)ps.size(); ++k)
    {
        *denom += *AvEs[k];
    }

    // Building weights
    for (int k = 0; k < (int)ps.size(); ++k)
    {
        ComputeWeight(*WeightEs[k], *AvEs[k]);
    }

    for (int j = 0; j < (int)ps.size(); ++j)
    {
        for (int k = j + 1; k < (int)ps.size(); ++k)
        {
            ComputeWeight(*WeightPairs[j][k], *AvP_Pairs[j][k], psi_Pairs[j][k].get());
        }
    }

    // ---- GLOBAL checks for psi -------------------------------------------
    double psi_min = 0.0, psi_max = 0.0;
    GlobalMinMax(*psi, psi_min, psi_max);

    // Basic bounds check
    if (mfem::Mpi::WorldRank() == 0) {std::cout << "[Psi Check] min = " << psi_min 
            << ", max = " << psi_max << " (expected min = 1e-06, max = 1)" << std::endl;}


    if (psi_min < 0.0 || psi_max > 1.0 + 1e-6) {
        std::cerr << "[Psi Check] ERROR: psi values out of [0,1]!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (psi_min > 0.1) {
        std::cerr << "[Psi Check] ERROR: psi_min not near 0." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (psi_max < 0.9) {
        std::cerr << "[Psi Check] ERROR: psi_max not near 1." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void Domain_Parameters::CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total) {
    local_total = 0.0;

    // Average value for each element
    mfem::Vector element_avg_values(nE);

    for (int ei = 0; ei < pmesh->GetNE(); ei++) {

        mfem::Array<double> vertex_values;
        grid_function.GetNodalValues(ei, vertex_values);

        double avg_value = 0.0;
        for (int vt = 0; vt < vertex_values.Size(); vt++) {
            avg_value += vertex_values[vt];
        }
        avg_value /= vertex_values.Size();

        element_avg_values(ei) = avg_value;

        // Accumulate the weighted value for the current element
        local_total += avg_value * element_volumes(ei);
    }

    // Perform global MPI reduction to sum values across all processes
    MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void Domain_Parameters::CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total) {

    EVol.SetSize(nE);
	for (int ei = 0; ei < nE; ei++){
        EVol(ei) = pmesh->GetElementVolume(ei);
	}

    // Call the general `CalculateTotals` method for any field
    CalculateTotals(grid_function, EVol, total, global_total);

}

void Domain_Parameters::CalculatePhasePotentialsAndTargetCurrent() {

    // Calculate totals for Psi and Pse fields
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    CalculateTotalPhaseField(*pse, tPse, gtPse);
    
    const std::vector<sim::MaterialType>& active_materials =
    (cfg.half_electrode == sim::Electrode::CATHODE)
        ? cfg.cathode_materials
        : cfg.anode_materials;

    
    gTrgI = 0.0;

    for (int k = 0; k < ps.size(); ++k)
    {
        CalculateTotalPhaseField(*ps[k], tPs[k], gtPs[k]);
        CalculateTargetCurrent(tPs[k], gTrgPs[k], active_materials[k]);

        gTrgI += gTrgPs[k];
    }             
}

void Domain_Parameters::CalculateTargetCurrent(double total_psi, double &global_total, sim::MaterialType material) {

    double rho = MaterialProperties::SiteDensity(material);

    // Compute target current based on total Psi, rho, Cr, and constants
    trgI = total_psi * rho * (0.95 - 0.3) / (3600.0 / cfg.Cr); // bounds of cathode 

    // Perform global MPI reduction to get the total target current
    MPI_Allreduce(&trgI, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Domain_Parameters::PrintInfo() {

    if (mfem::Mpi::WorldRank() == 0) // only print on rank 0
    {
        cout << "Total Psi: " << gtPsi << endl;
        cout << "Total Pse: " << gtPse << endl;
        cout << "Target Current: " << gTrgI << endl;
    }
}

