#ifndef INITIALIZE_GEOMETRY_HPP
#define INITIALIZE_GEOMETRY_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>
#include "SimTypes.hpp"
#include "SimulationConfig.hpp"
#include <set>

/**
 * @file Initialize_Geometry.hpp
 * @brief Defines mesh, geometry, voxel, distance-field, and FE-space setup for BESFEM.
 */

using namespace std;

/**
 * @class Initialize_Geometry
 * @brief Central geometry handler for BESFEM simulations.
 *
 * This class is responsible for:
 * - Loading meshes (serial or voxel-based)
 * - Creating global and parallel MFEM meshes
 * - Initializing serial and parallel FE spaces (H1 and DG)
 * - Mapping global fields to parallel fields
 * - Loading distance functions ψ (solid), ψ_e (electrolyte), and voxel data
 * - Setting up boundary-condition markers
 * - Selecting pinned DOFs for potential anchoring
 *
 * It serves as the lowest-level mesh/geometry infrastructure on which all
 * physics components (CnA, CnC, CnE, potentials, reactions) depend.
 */
class Initialize_Geometry {
private:
    bool tiffDataLoaded = false; ///< Whether TIFF voxel data has been loaded.
    const SimulationConfig& cfg;

protected:
    mfem::Vector elementVolumes;     ///< Per-element volumes (global or parallel).
    mfem::Array<int> boundaryMarkers; ///< Boundary attribute markers.
    std::vector<std::vector<std::vector<int>>> data; ///< Raw voxel data container.

public:

    /**
     * @brief Construct the geometry handler.
     *
     * Stores the simulation configuration and prepares geometry-related storage.
     *
     * @param cfg Reference to the simulation configuration.
     */
    Initialize_Geometry(const SimulationConfig& cfg);

    /// Destructor.
    virtual ~Initialize_Geometry();

    bool combine_particle_groups = false; ///< Whether to combine particle groups for performance.


    // -------------------------------------------------------------------------
    // Distance function preprocessing
    // -------------------------------------------------------------------------

    /**
     * @brief Adjust distance function file to account for mesh spacing (dh).
     *
     * Applies smoothing, scaling, or clipping to ensure consistent SBM interface
     * thickness across meshes of different resolutions.
     *
     * @param distanceFile Path to the distance function file.
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     */
    void AdjustDistanceFile(const char* distanceFile, const char* mesh_type);

    // -------------------------------------------------------------------------
    // Mesh initialization (serial + parallel)
    // -------------------------------------------------------------------------

    /**
     * @brief Initialize mesh and distance fields for a half-cell simulation.
     *
     * Creates the global mesh, assigns the distance field ψ or ψ_e from file,
     * builds FE spaces, and initializes parallel mesh/FESpaces.
     *
     * @param meshFile Path to mesh file.
     * @param distanceFile Distance function file.
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFile, const char* mesh_type, MPI_Comm comm, int order);

    /**
     * @brief Initialize mesh and distance fields for a full-cell simulation.
     *
     * Loads separate anode and cathode distance functions.
     *
     * @param meshFile Path to mesh file.
     * @param distanceFileA Anode distance function file (ψ_A).
     * @param distanceFileC Cathode distance function file (ψ_C).
     * @param mesh_type Type of mesh ("ml" = MATLAB, "v" = voxel).
     * @param comm MPI communicator.
     * @param order Polynomial order for FE space.
     */
    void InitializeMesh(const char* meshFile, const char* distanceFileA, const char* distanceFileC, const char* mesh_type, MPI_Comm comm, int order);

    /**
     * @brief Load and construct global serial MFEM mesh.
     *
     * @param meshFile Path to mesh file.
     */
    void InitializeGlobalMesh(const char* meshFile);

    /**
     * @brief Build a distributed (parallel) mesh from the global mesh.
     *
     * @param comm MPI communicator.
     */
    void InitializeParallelMesh(MPI_Comm comm);

    // -------------------------------------------------------------------------
    // Voxel / TIFF mesh support
    // -------------------------------------------------------------------------

    /**
     * @brief Read TIFF voxel data for voxel-mesh construction.
     *
     * @param meshFile TIFF volume file.
     * @return 3D voxel array.
     */
    std::vector<std::vector<std::vector<int>>> ReadTiffFile(const char* meshFile);

    /**
     * @brief Construct global mesh from voxelized TIFF data.
     *
     * @param tiffData 3D voxel grid.
     * @return Newly constructed MFEM mesh.
     */
    std::unique_ptr<mfem::Mesh>
    CreateGlobalMeshFromTiffData(const std::vector<std::vector<std::vector<int>>>& tiffData);

    // -------------------------------------------------------------------------
    // Finite element spaces
    // -------------------------------------------------------------------------

    /**
     * @brief Setup global serial FE space.
     *
     * @param order FE polynomial order.
     */
    void SetupFiniteElementSpace(int order);

    /**
     * @brief Setup parallel FE space.
     *
     * @param order FE polynomial order.
     */
    void SetupParFiniteElementSpace(int order);

    // -------------------------------------------------------------------------
    // Global → parallel mapping and distance functions
    // -------------------------------------------------------------------------

    /**
     * @brief Assign global fields: distance functions, voxel mask, etc.
     *
     * Loads, projects, and assigns global fields onto the global mesh.
     *
     * @param mesh_file Mesh path.
     * @param distanceFile Distance field file.
     * @param gDsF_out Output global distance function ψ.
     */
    void AssignGlobalValues(const char* mesh_file, const char* distanceFile, std::unique_ptr<mfem::GridFunction>& gDsF_out);

    /**
     * @brief Transfer global mesh quantities to parallel fields.
     *
     * @param meshFile Mesh file for consistency checking.
     */
    void MapGlobalToLocal(const char* meshFile);

    // -------------------------------------------------------------------------
    // Boundary conditions and pinning
    // -------------------------------------------------------------------------

    /**
     * @brief Set Neumann/Dirichlet boundary condition markers.
     *
     * Uses boundary attributes + cell mode to mark anode, cathode, electrolyte regions.
     *
     * @param mode HALF or FULL cell.
     * @param electrode ANODE, CATHODE, BOTH.
     */
    void SetupBoundaryConditions(sim::CellMode mode, sim::Electrode electrode);

    /**
     * @brief Print mesh information (elements, vertices, dimensions).
     */
    void PrintMeshInfo();

    /**
     * @brief Select and assign a pinned DOF for potential anchoring.
     *
     * @param fespace FE space in which to anchor a single true DOF.
     */
    void SetupPinnedDOF(mfem::ParFiniteElementSpace& fespace);

    
    
    /**
     * @brief Save TIFF voxel data to a series of PGM files.
     * 
     * @param data voxel data.
     * @param filename Base filename for output PGM files.
     */
    void SaveTiffDataToPGM(const std::vector<std::vector<std::vector<int>>> &data,
                       const std::string &filename);

    /**
     * @brief Use MFEM-based solvers to compute distance function from voxel mask.
     * @param dist Output unsigned distance function.
     * @param filt_gf Output filtered level set function.
     * @param mode 0 = psi (keep 1s to right), 1 = pse (keep 0s to any boundary)
     */
    void ComputePDEFilter(
        mfem::ParGridFunction &dist,            // output unsigned distance
        mfem::ParGridFunction &filt_gf,        // output filtered level set
        int mode                           // 0 = psi (keep 1s to right), 1 = pse (keep 0s to any boundary)
    );

    /**
     * @brief Compute a filtered distance field for a specific voxel label.
     *
     * Builds a label-specific mask from the TIFF data, applies the PDE filter,
     * and computes a corresponding distance field. This is used for individual
     * particle or material-region masks in segmented microstructures.
     *
     * @param dist Output unsigned distance field.
     * @param filt_gf Output filtered level-set field.
     * @param target_label Voxel label to isolate.
     * @param keep_boundary_connected Whether to keep only the boundary-connected region.
     * @param seed_side_or_face Optional boundary side/face used as the seed region.
     */
    void ComputePDEFilterLabel(mfem::ParGridFunction &dist,
                           mfem::ParGridFunction &filt_gf,
                           int target_label,
                           bool keep_boundary_connected,
                           int seed_side_or_face = -1
    );

    /**
     * @brief Return the unique particle/material labels found in the TIFF data.
     *
     * @return Vector of integer labels present in the voxelized geometry.
     */
    std::vector<int> GetParticleLabelsFromTiff() const;

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Access distributed mesh.
     * @return Pointer to parallel mesh.
     */
    mfem::ParMesh *GetParallelMesh() const { return parallelMesh.get(); }

    /**
     * @brief Access parallel FE space.
     * @return Parallel FE space (shared_ptr).
     */
    std::shared_ptr<mfem::ParFiniteElementSpace>
    GetParFiniteElementSpace() const { return parfespace; }

    // -------------------------------------------------------------------------
    // Public geometry/mesh fields
    // -------------------------------------------------------------------------
    int nV = 0; ///< Number of vertices.
    int nE = 0; ///< Number of elements.
    int nC = 0; ///< Corners per element.

    int gei = 0; ///< Global element index.
    int ei  = 0; ///< Local element index.

    // Boundary condition marker arrays
    mfem::Array<int> nbc_w_bdr, nbc_s_bdr, nbc_e_bdr, nbc_n_bdr;
    mfem::Array<int> nbc_bdr, dbc_bdr;
    mfem::Array<int> dbc_w_bdr, dbc_e_bdr; 
    mfem::Array<int> ess_tdof_list_w, ess_tdof_list_e;

    mfem::Array<int> gVTX; ///< Global vertex IDs of current element.
    mfem::Array<int> VTX;  ///< Local vertex IDs of current element.

    std::unique_ptr<mfem::Mesh> globalMesh; ///< Serial/global mesh.
    mfem::Array<HYPRE_BigInt> E_L2G; ///< Local-to-global element mapping.

    double Onm = 0.0; ///< Number of grid function entries.

    // Pinned DOF information
    mfem::Array<int> ess_tdof_potE;
    bool anchor_set = false;
    HYPRE_BigInt global_anchor_potE = -1;
    int anchor_owner_potE = -1;
    bool pin = false;

    mfem::Array<int> ess_tdof_listPinned;
    mfem::Array<int> boundary_dofs;
    mfem::Array<int> ess_tdof_marker;

    int myid = 0; ///< MPI rank.
    int rkpp = -1; ///< Rank that owns the pinned DOF.

    std::shared_ptr<mfem::ParMesh> parallelMesh; ///< Distributed parallel mesh.

    std::shared_ptr<mfem::FiniteElementSpace> feSpace; ///< Serial H1 finite element space.
    std::shared_ptr<mfem::FiniteElementSpace> globalfespace; ///< Global serial finite element space.

    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel H1 finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg; ///< Parallel DG finite element space.
    std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg; ///< Vector-valued parallel DG finite element space.

    std::unique_ptr<mfem::GridFunction> gDsF; ///< Global serial distance field.
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< Parallel distance field.

    std::unique_ptr<mfem::GridFunction> gDsF_A, gDsF_C; ///< Global anode/cathode distance fields.
    std::unique_ptr<mfem::ParGridFunction> dsF_A, dsF_C; ///< Parallel anode/cathode distance fields.

    std::unique_ptr<mfem::GridFunction> gVox; ///< Global voxel-label field.
    std::unique_ptr<mfem::ParGridFunction> Vox; ///< Parallel voxel-label field.

    std::vector<std::vector<std::vector<int>>> tiffData; ///< Raw TIFF voxel labels.

    std::unique_ptr<mfem::H1_FECollection> gfec, pfec; ///< Serial/parallel H1 finite element collections.
    std::unique_ptr<mfem::DG_FECollection> pfec_dg; ///< Parallel DG finite element collection.

    std::unique_ptr<mfem::ParGridFunction> distMask; ///< Unsigned distance-to-mask field.
    std::unique_ptr<mfem::ParGridFunction> distMaskSigned; ///< Optional signed distance-to-mask field.
    std::unique_ptr<mfem::ParGridFunction> MaskFilter; ///< Filtered solid-mask level-set field.
    std::unique_ptr<mfem::ParGridFunction> MaskFilterPse; ///< Filtered electrolyte-mask level-set field.

    std::vector<int> particle_labels; ///< Unique particle/material labels from the TIFF geometry.
    std::vector<std::unique_ptr<mfem::ParGridFunction>> MaskFilters; ///< Per-label filtered mask fields.

};

#endif // INITIALIZE_GEOMETRY_HPP
