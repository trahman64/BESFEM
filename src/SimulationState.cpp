#include "../include/SimulationState.hpp"
#include "../include/MaterialProperties.hpp"

static double GetInitialValue(
    const std::vector<double>& values,
    int k,
    double fallback)
{
    if (k < static_cast<int>(values.size()))
    {
        return values[k];
    }
    return fallback;
}

static void InitializePairWorkspaces(SimulationState& state, Initialize_Geometry& geometry, int np)
{
    state.mu_pair_a.clear();
    state.mu_pair_b.clear();
    state.sum_pairs.clear();

    state.mu_pair_a.resize(np);
    state.mu_pair_b.resize(np);
    state.sum_pairs.resize(np);

    for (int j = 0; j < np; ++j)
    {
        state.mu_pair_a[j].resize(np);
        state.mu_pair_b[j].resize(np);
        state.sum_pairs[j].resize(np);

        for (int k = 0; k < np; ++k)
        {
            if (j < k)
            {
                state.mu_pair_a[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                state.mu_pair_b[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
                state.sum_pairs[j][k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }
        }
    }

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "[DEBUG] Initialized pair workspaces for np = "
                  << np << std::endl;
    }
}


// static inline double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data)
// {
//     const double dx = ticks(1) - ticks(0);

//     if (cn < 1.0e-6) cn = 1.0e-6;
//     if (cn > 0.999999) cn = 0.999999;

//     int idx = std::floor((cn - ticks(0)) / dx);
//     if (idx < 0) idx = 0;

//     const int max_idx = ticks.Size() - 2;
//     if (idx > max_idx) idx = max_idx;

//     return data(idx) + (cn - ticks(idx)) / dx * (data(idx + 1) - data(idx));
// }

// static inline double NMC_mu(double c)
// {
//     return -Constants::Frd * ((1.095 * c * c) - (8.234e-7 * std::exp(14.31 * c)) + (4.692 * std::exp(-0.5389 * c)));
// }

// static inline double LFP_mu(double c)
// {
//     static mfem::Vector Ticks(201);
//     static mfem::Vector chmPot(201);
//     static bool loaded = false;

//     if (!loaded)
//     {
//         std::ifstream myXfile("../inputs/LFP_Chm_Pot_Ticks.txt");
//         std::ifstream mydFfile("../inputs/LFP_Chm_Pot.txt");

//         if (!myXfile || !mydFfile)
//         {
//             mfem::mfem_error("Could not open LFP chemical potential input files.");
//         }

//         for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
//         for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

//         loaded = true;
//     }

//     return -Constants::Frd * (GetTableValues(c, Ticks, chmPot) + 3.4);
// }

// static inline double Graphite_mu(double c)
// {
//     static mfem::Vector Ticks(101);
//     static mfem::Vector chmPot(101);
//     static bool loaded = false;

//     if (!loaded)
//     {
//         std::ifstream myXfile("../inputs/C_Li_X_101.txt");
//         std::ifstream mydFfile("../inputs/C_Li_M6_101.txt");

//         if (!myXfile || !mydFfile)
//         {
//             mfem::mfem_error("Could not open graphite chemical potential input files.");
//         }

//         for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
//         for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);

//         loaded = true;
//     }

//     return GetTableValues(c, Ticks, chmPot);
// }

// static inline double CathodeChemicalPotential(sim::MaterialType material, double c)
// {
//     switch (material)
//     {
//         case sim::MaterialType::NMC:
//             return NMC_mu(c);

//         case sim::MaterialType::LFP:
//             return LFP_mu(c);
        
//         default:
//             mfem::mfem_error("Unknown cathode material type.");
//             return 0.0;
//     }
// }

void UpdateCathodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters)
{
    const int np = static_cast<int>(state.cathode_particles.size());

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *state.cathode_particles[j].Cn_gf;
            auto& Ck = *state.cathode_particles[k].Cn_gf;

            const auto mat_j = state.cathode_particles[j].material;
            const auto mat_k = state.cathode_particles[k].material;

            auto& mu_j = *state.mu_pair_a[j][k];
            auto& mu_k = *state.mu_pair_b[j][k];
            auto& AvP_pair = *domain_parameters.AvP_Pairs[j][k];

            mu_j = 0.0;
            mu_k = 0.0;

            // bool printed_pair = false;
            // int interface_count = 0;

            for (int vi = 0; vi < geometry.nV; ++vi)
            {
                if (AvP_pair(vi) > 1000.0)
                {
                    // interface_count++;

                    mu_j(vi) = MaterialProperties::CathodeChemicalPotential(mat_j, Cj(vi));
                    mu_k(vi) = MaterialProperties::CathodeChemicalPotential(mat_k, Ck(vi));

                    // if (!printed_pair && mfem::Mpi::WorldRank() == 0)
                    // {
                    //     auto MaterialName = [](sim::MaterialType mat)
                    //     {
                    //         switch (mat)
                    //         {
                    //             case sim::MaterialType::NMC: return "NMC";
                    //             case sim::MaterialType::LFP: return "LFP";
                    //             default: return "Unknown";
                    //         }
                    //     };

                        // std::cout << "\n[DEBUG MU TEST]" << std::endl;
                        // std::cout << "Pair (" << j << "," << k << ")" << std::endl;

                        // std::cout << "Particle " << j
                        //         << " material = " << MaterialName(mat_j)
                        //         << ", C = " << Cj(vi)
                        //         << ", mu = " << mu_j(vi)
                        //         << std::endl;

                        // std::cout << "Particle " << k
                        //         << " material = " << MaterialName(mat_k)
                        //         << ", C = " << Ck(vi)
                        //         << ", mu = " << mu_k(vi)
                        //         << std::endl;

                        // printed_pair = true;
                    // }
                }
            }

            // if (mfem::Mpi::WorldRank() == 0)
            // {
            //     std::cout << "[DEBUG MU] Pair (" << j << "," << k
            //             << ") interface_count = " << interface_count
            //             << std::endl;
            // }
        }
    }
}

void UpdateAnodePairChemicalPotentials(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters)
{
    const int np = static_cast<int>(state.anode_particles.size());

    for (int j = 0; j < np; ++j)
    {
        for (int k = j + 1; k < np; ++k)
        {
            auto& Cj = *state.anode_particles[j].Cn_gf;
            auto& Ck = *state.anode_particles[k].Cn_gf;

            const auto mat_j = state.anode_particles[j].material;
            const auto mat_k = state.anode_particles[k].material;

            auto& mu_j = *state.mu_pair_a[j][k];
            auto& mu_k = *state.mu_pair_b[j][k];
            auto& AvP_pair = *domain_parameters.AvP_Pairs[j][k];

            mu_j = 0.0;
            mu_k = 0.0;

            for (int vi = 0; vi < geometry.nV; ++vi)
            {
                if (AvP_pair(vi) > 1000.0)
                {
                    mu_j(vi) = MaterialProperties::AnodeChemicalPotential(mat_j, Cj(vi));
                    mu_k(vi) = MaterialProperties::AnodeChemicalPotential(mat_k, Ck(vi));
                }
            }
        }
    }
}


static void InitializeAnodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const SimulationConfig& cfg, BoundaryConditions& bc)
{
    const int np = static_cast<int>(domain_parameters.ps.size());
    state.anode_particles.clear();
    state.anode_particles.resize(np);

    const std::vector<double>& init_values = cfg.init_anode_particles;

    if (np == 0)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] No different anode particles defined in the configuration." << std::endl;
        }
        return;
    }

    for (int k = 0; k < np; ++k)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Creating Anode Particle " << k
                    << " (label = " << domain_parameters.particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.anode_particles[k];
        p.label = domain_parameters.particle_labels[k];
        p.material = cfg.anode_materials[k];

        p.concentration = std::make_unique<CnA>(geometry, domain_parameters);
        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.potential     = std::make_unique<PotA>(geometry, domain_parameters, bc);
        p.ph_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        
        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(init_values, k, Constants::init_CnA);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *domain_parameters.ps[k], domain_parameters.gtPs[k]);
        p.potential->SetupField(*p.ph_gf, Constants::init_BvA, *domain_parameters.ps[k]);
    }
}

static void InitializeCathodeParticles(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    const SimulationConfig& cfg, BoundaryConditions& bc)
{
    const int np = static_cast<int>(domain_parameters.ps.size());
    state.cathode_particles.clear();
    state.cathode_particles.resize(np);

    const std::vector<double>& init_values = cfg.init_cathode_particles;

    if (np == 0)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] No different cathode particles defined in the configuration." << std::endl;
        }
        return;
    }

    for (int k = 0; k < np; ++k)
    {
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Creating Cathode Particle " << k
                    << " (label = " << domain_parameters.particle_labels[k] << ")"
                    << std::endl;
        }

        auto& p = state.cathode_particles[k];
        p.label = domain_parameters.particle_labels[k];
        p.material = cfg.cathode_materials[k];

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG] Cathode Particle " << k
                    << " assigned material = ";

            switch (p.material)
            {
                case sim::MaterialType::NMC:
                    std::cout << "NMC";
                    break;

                case sim::MaterialType::LFP:
                    std::cout << "LFP";
                    break;

                case sim::MaterialType::Graphite:
                    std::cout << "Graphite";
                    break;
            }

            std::cout << std::endl;
        }

        p.concentration = std::make_unique<CnC>(geometry, domain_parameters);
        p.Cn_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Cn_gf_psi     = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction      = std::make_unique<Reaction>(geometry, domain_parameters);
        p.Rxn_gf        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        p.Rx_src        = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.potential     = std::make_unique<PotC>(geometry, domain_parameters, bc);
        p.ph_gf         = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        p.reaction->Initialize(*p.Rxn_gf, Constants::init_Rxn);

        const double init_cn = GetInitialValue(init_values, k, Constants::init_CnC);

        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[DEBUG]   Initial concentration = " << init_cn << std::endl;
        }

        p.concentration->SetupField(*p.Cn_gf, init_cn, *domain_parameters.ps[k], domain_parameters.gtPs[k]);
        p.potential->SetupField(*p.ph_gf, Constants::init_BvC, *domain_parameters.ps[k]);
    }
}

void InitializeFields(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
    BoundaryConditions& bc, const SimulationConfig& cfg)
{
    state.CnP_together = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.CnE_gf_psi   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

    state.electrolyte_concentration = std::make_unique<CnE>(geometry, domain_parameters, bc, cfg.mode);
    state.CnE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_concentration->SetupField(*state.CnE_gf, Constants::init_CnE, *domain_parameters.pse, domain_parameters.gtPse);

    state.electrolyte_potential = std::make_unique<PotE>(geometry, domain_parameters, bc, cfg.mode);
    state.phE_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.electrolyte_potential->SetupField(*state.phE_gf, Constants::init_BvE, *domain_parameters.pse);

    state.reaction = std::make_unique<Reaction>(geometry, domain_parameters);
    state.Rxn_gf   = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
    state.reaction->Initialize(*state.Rxn_gf, Constants::init_Rxn);

    if (cfg.mode == sim::CellMode::HALF)
    {
        if (cfg.half_electrode == sim::Electrode::ANODE)
        {
            const int np = static_cast<int>(domain_parameters.ps.size());

            state.CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            state.anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
            state.CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.anode_concentration->SetupField(*state.CnA_gf, Constants::init_CnA, *domain_parameters.psi, domain_parameters.gtPsi);

            state.anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
            state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.anode_potential->SetupField(*state.phA_gf, Constants::init_BvA, *domain_parameters.psi);

            InitializeAnodeParticles(state, geometry, domain_parameters, cfg, bc);
            InitializePairWorkspaces(state, geometry, static_cast<int>(state.anode_particles.size()));

            state.anode_out.clear();
            state.anode_out.resize(np);

            for (int k = 0; k < np; ++k)
            {
                state.anode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }

        }
        else
        {
            const int np = static_cast<int>(domain_parameters.ps.size());

            state.CnC_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

            state.cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
            state.CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.cathode_concentration->SetupField(*state.CnC_gf, Constants::init_CnC, *domain_parameters.psi, domain_parameters.gtPsi);

            state.cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
            state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            state.cathode_potential->SetupField(*state.phC_gf, Constants::init_BvC, *domain_parameters.psi);

            InitializeCathodeParticles(state, geometry, domain_parameters, cfg, bc);
            InitializePairWorkspaces(state, geometry,static_cast<int>(state.cathode_particles.size()));

            state.cathode_out.clear();
            state.cathode_out.resize(np);

            for (int k = 0; k < np; ++k)
            {
                state.cathode_out[k] = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
            }
        }
    }
    else
    {
        state.CnA_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.CnC_gf_psi = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());

        state.anode_concentration = std::make_unique<CnA>(geometry, domain_parameters);
        state.CnA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.anode_concentration->SetupField(*state.CnA_gf, Constants::init_CnA, *domain_parameters.psA, domain_parameters.gtPsA);

        state.anode_potential = std::make_unique<PotA>(geometry, domain_parameters, bc);
        state.phA_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.anode_potential->SetupField(*state.phA_gf, Constants::init_BvA, *domain_parameters.psA);

        state.cathode_concentration = std::make_unique<CnC>(geometry, domain_parameters);
        state.CnC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.cathode_concentration->SetupField(*state.CnC_gf, Constants::init_CnC, *domain_parameters.psC, domain_parameters.gtPsC);

        state.cathode_potential = std::make_unique<PotC>(geometry, domain_parameters, bc);
        state.phC_gf = std::make_unique<mfem::ParGridFunction>(geometry.parfespace.get());
        state.cathode_potential->SetupField(*state.phC_gf, Constants::init_BvC, *domain_parameters.psC);
    }

    if (mfem::Mpi::WorldRank() == 0 && (state.anode_particles.size() > 0 || state.cathode_particles.size() > 0))
    {
        std::cout << "[DEBUG] Finished InitializeFields()" << std::endl;

        std::cout << "    Anode particles:   "
                << state.anode_particles.size() << std::endl;

        std::cout << "    Cathode particles: "
                << state.cathode_particles.size() << std::endl;
    }
}

void Pairs(SimulationState& state, Initialize_Geometry& geometry, Domain_Parameters& domain_parameters, int j, std::vector<ConcentrationBase::PairCoupling>& pair_terms, int np, int t)
{                        
    for (int k = 0; k < np; ++k)
    {
        if (j == k) { continue; }

        const int a = std::min(j, k);
        const int b = std::max(j, k);

        ConcentrationBase::PairCoupling pair;
        pair.sum_part = state.sum_pairs[a][b].get();
        pair.weight   = domain_parameters.WeightPairs[a][b].get();
        pair.grad_psi = domain_parameters.AvP_Pairs[a][b].get();

        if (j < k)
        {
            pair.mu_self = state.mu_pair_a[a][b].get();
            pair.mu_nbr  = state.mu_pair_b[a][b].get();
        }
        else
        {
            pair.mu_self = state.mu_pair_b[a][b].get();
            pair.mu_nbr  = state.mu_pair_a[a][b].get();
        }

        pair_terms.push_back(pair);

        if (mfem::Mpi::WorldRank() == 0 && t == 1)
        {
            std::cout << "[DEBUG] Pair (j,k) = (" << j << "," << k << ")"; 
            std::cout << std::endl;
        }

    }
}
