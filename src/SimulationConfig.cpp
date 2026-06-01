#include "../include/SimulationConfig.hpp"
#include "mfem.hpp"

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

static std::string Trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";

    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> SplitString(const std::string& text, char delimiter)
{
    std::vector<std::string> tokens;
    std::stringstream ss(text);
    std::string item;

    while (std::getline(ss, item, delimiter))
    {
        item = Trim(item);
        if (!item.empty())
        {
            tokens.push_back(item);
        }
    }

    return tokens;
}

static std::unordered_map<std::string, std::string>
ReadConfigFile(const std::string& filename)
{
    std::unordered_map<std::string, std::string> data;

    std::ifstream file(filename);
    if (!file)
    {
        mfem::mfem_error(("Could not open config file: " + filename).c_str());
    }

    std::string line;

    while (std::getline(file, line))
    {
        line = Trim(line);

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
        {
            line = Trim(line.substr(0, comment_pos));
        }

        const auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = Trim(line.substr(0, eq_pos));
        std::string val = Trim(line.substr(eq_pos + 1));

        data[key] = val;
    }

    return data;
}

static bool HasKey(
    const std::unordered_map<std::string, std::string>& data,
    const std::string& key)
{
    return data.find(key) != data.end();
}

static std::string GetValue(
    const std::unordered_map<std::string, std::string>& data,
    const std::string& key)
{
    auto it = data.find(key);

    if (it == data.end())
    {
        mfem::mfem_error(("Missing required config key: " + key).c_str());
    }

    return it->second;
}

static bool ParseBool(const std::string& value)
{
    if (value == "true" || value == "1" || value == "yes") return true;
    if (value == "false" || value == "0" || value == "no") return false;

    mfem::mfem_error(("Invalid boolean value: " + value +
                      ". Use true/false.").c_str());

    return false;
}

static sim::MaterialType ParseMaterial(const std::string& name)
{
    if (name == "Graphite" || name == "graphite")
        return sim::MaterialType::Graphite;

    if (name == "NMC" || name == "nmc")
        return sim::MaterialType::NMC;

    if (name == "LFP" || name == "lfp")
        return sim::MaterialType::LFP;

    mfem::mfem_error(("Invalid material: " + name +
                      ". Use Graphite, NMC, or LFP.").c_str());

    return sim::MaterialType::Electrolyte;
}

static std::vector<sim::MaterialType> ParseMaterialList(const std::string& text)
{
    std::vector<sim::MaterialType> materials;

    for (const auto& token : SplitString(text, ','))
    {
        materials.push_back(ParseMaterial(token));
    }

    return materials;
}

static std::vector<double> ParseDoubleList(const std::string& text)
{
    std::vector<double> values;

    for (const auto& token : SplitString(text, ','))
    {
        values.push_back(std::stod(token));
    }

    return values;
}

static void ApplyConfigFile(SimulationConfig& cfg)
{
    auto data = ReadConfigFile(cfg.config_file);

    if (HasKey(data, "mode"))
    {
        std::string mode = GetValue(data, "mode");

        if (mode == "half")
            cfg.mode = sim::CellMode::HALF;
        else if (mode == "full")
            cfg.mode = sim::CellMode::FULL;
        else
            mfem::mfem_error("Invalid config mode. Use: half | full.");
    }

    if (HasKey(data, "electrode"))
    {
        std::string electrode = GetValue(data, "electrode");

        if (electrode == "anode")
            cfg.half_electrode = sim::Electrode::ANODE;
        else if (electrode == "cathode")
            cfg.half_electrode = sim::Electrode::CATHODE;
        else
            mfem::mfem_error("Invalid config electrode. Use: anode | cathode.");
    }

    if (HasKey(data, "mesh_file"))
    {
        cfg.mesh_file =
            strdup(GetValue(data, "mesh_file").c_str());
    }

    if (HasKey(data, "mesh_type"))
    {
        cfg.mesh_type =
            strdup(GetValue(data, "mesh_type").c_str());
    }

    if (HasKey(data, "anode_distance"))
    {
        cfg.dsF_file_A =
            strdup(GetValue(data, "anode_distance").c_str());
    }

    if (HasKey(data, "cathode_distance"))
    {
        cfg.dsF_file_C =
            strdup(GetValue(data, "cathode_distance").c_str());
    }

    if (HasKey(data, "num_steps"))
        cfg.num_timesteps = std::stoi(GetValue(data, "num_steps"));

    if (HasKey(data, "combine_particles"))
        cfg.combine_particle_groups = ParseBool(GetValue(data, "combine_particles"));

    if (HasKey(data, "cathode_materials"))
        cfg.cathode_materials = ParseMaterialList(GetValue(data, "cathode_materials"));

    if (HasKey(data, "anode_materials"))
        cfg.anode_materials = ParseMaterialList(GetValue(data, "anode_materials"));

    if (HasKey(data, "init_cathode_particles"))
        cfg.init_cathode_particles = ParseDoubleList(GetValue(data, "init_cathode_particles"));

    if (HasKey(data, "init_anode_particles"))
        cfg.init_anode_particles = ParseDoubleList(GetValue(data, "init_anode_particles"));

    if (HasKey(data, "init_CnE"))
        cfg.init_CnE = std::stod(GetValue(data, "init_CnE"));

    if (HasKey(data, "init_BvA"))
        cfg.init_BvA = std::stod(GetValue(data, "init_BvA"));

    if (HasKey(data, "init_BvC"))
        cfg.init_BvC = std::stod(GetValue(data, "init_BvC"));

    if (HasKey(data, "init_BvE"))
        cfg.init_BvE = std::stod(GetValue(data, "init_BvE"));
}

SimulationConfig ParseSimulationArgs(int argc, char *argv[])
{
    SimulationConfig cfg;

    const char* config_file = "../inputs/run_config.txt";
    bool list_options = false;

    mfem::OptionsParser args(argc, argv);

    args.AddOption(&config_file,
                   "-cfg", "--config",
                   "User-editable BESFEM config file.");

    args.AddOption(&list_options,
                   "-opts", "--list-options",
                   "-no-opts", "--no-list-options",
                   "Print available BESFEM option choices.");

    args.ParseCheck();

    cfg.config_file = config_file;

    if (list_options)
    {
        PrintAvailableSimulationOptions();
        std::exit(0);
    }

    ApplyConfigFile(cfg);

    return cfg;
}

void ValidateConfig(const SimulationConfig &cfg, int argc, char *argv[])
{
    if (cfg.mode == sim::CellMode::FULL)
    {
        if (!cfg.dsF_file_A || !cfg.dsF_file_C)
        {
            mfem::mfem_error("FULL mode requires both anode_distance and cathode_distance.");
        }
    }
    else
    {
        const bool cathode = cfg.half_electrode == sim::Electrode::CATHODE;

        if (cathode && !cfg.dsF_file_C)
        {
            mfem::mfem_error("HALF-CATHODE requires cathode_distance.");
        }

        if (!cathode && !cfg.dsF_file_A)
        {
            mfem::mfem_error("HALF-ANODE requires anode_distance.");
        }
    }

    if (!cfg.mesh_file)
        mfem::mfem_error("mesh_file cannot be empty.");

    if (std::strcmp(cfg.mesh_type, "ml") != 0 &&
        std::strcmp(cfg.mesh_type, "v")  != 0)
    {
        mfem::mfem_error("mesh_type must be ml or v.");
    }

    if (cfg.num_timesteps <= 0)
        mfem::mfem_error("num_steps must be positive.");

    if (cfg.cathode_materials.empty())
        mfem::mfem_error("cathode_materials cannot be empty.");

    if (cfg.anode_materials.empty())
        mfem::mfem_error("anode_materials cannot be empty.");

    if (cfg.init_cathode_particles.empty())
        mfem::mfem_error("init_cathode_particles cannot be empty.");

    if (cfg.init_anode_particles.empty())
        mfem::mfem_error("init_anode_particles cannot be empty.");

    if (cfg.cathode_materials.size() != cfg.init_cathode_particles.size())
    {
        mfem::mfem_error(
            "cathode_materials and init_cathode_particles must have the same length.");
    }

    if (cfg.anode_materials.size() != cfg.init_anode_particles.size())
    {
        mfem::mfem_error(
            "anode_materials and init_anode_particles must have the same length.");
    }
}

void PrintAvailableSimulationOptions()
{
    std::cout << "\nAvailable BESFEM command-line choices to edit in ../inputs/run_config.txt:\n\n";

    std::cout << "  Run config:\n";
    std::cout << "    -cfg ../inputs/run_config.txt\n\n";

    std::cout << "  Cell modes:\n";
    std::cout << "    half\n";
    std::cout << "    full\n\n";

    std::cout << "  Half-cell electrodes:\n";
    std::cout << "    anode\n";
    std::cout << "    cathode\n\n";

    std::cout << "  Materials:\n";
    std::cout << "    Cathode: NMC, LFP\n";
    std::cout << "    Anode:   Graphite\n\n";

    std::cout << "  Mixed cathode example:\n";
    std::cout << "    cathode_materials = LFP,LFP,NMC\n\n";

    std::cout << "  Mesh types:\n";
    std::cout << "    ml   MATLAB mesh\n";
    std::cout << "    v    voxel mesh\n\n";

    std::cout << "  Particle grouping:\n";
    std::cout << "    combine_particles = true\n";
    std::cout << "    combine_particles = false\n\n";

    std::cout << "  Initial particle concentrations:\n";
    std::cout << "    init_cathode_particles = 0.2,0.15,0.10\n";
    std::cout << "    init_anode_particles = 0.8,0.85,0.90\n\n";

    std::cout << "  Number of timesteps:\n";
    std::cout << "    num_steps = 1000\n\n";

    std::cout << " Initial boundary values:\n";
    std::cout << "    init_BvA = -0.01\n";
    std::cout << "    init_BvC = 3.30\n";
    std::cout << "    init_BvE = -0.10\n\n";

    std::cout << "  Example command:\n";
    std::cout << "    mpirun -np 4 ./battery_simulation -cfg ../inputs/run_config.txt\n\n";
}