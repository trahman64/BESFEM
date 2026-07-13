#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * @file Constants.hpp
 * @brief Declares global constants used throughout BESFEM simulations.
 *
 * Includes mesh/input file names, polynomial order, grid spacing parameters,
 * numerical tolerances, electrochemical constants, and initial conditions.
 */

/**
 * @namespace Constants
 * @brief Global user-adjustable simulation constants.
 *
 * This namespace provides a central location for:
 * - mesh and distance-function file paths,
 * - finite element order,
 * - physical parameters,
 * - model constants,
 * - default initialization values for potentials, concentrations, and reactions.
 *
 * Many of these values are read at runtime or overridden by command-line
 * arguments (see `SimulationConfig`).
 */
namespace Constants {

    // -------------------------------------------------------------------------
    // Discretization parameters
    // -------------------------------------------------------------------------

    extern const int    order; ///< Polynomial order for FE basis.
    extern const double zeta;  ///< Interface thickness parameter for SBM.

    // -------------------------------------------------------------------------
    // Numerical tolerances and thresholds
    // -------------------------------------------------------------------------

    extern const double thres; ///< Threshold for phase-field cutoff.
    extern const double eps;   ///< Small epsilon used to avoid division-by-zero.

    // -------------------------------------------------------------------------
    // Electrochemical constants
    // -------------------------------------------------------------------------

    extern const double t_minus; ///< Cation transference number.
    extern const double D0;      ///< Base diffusivity.
    extern const double Frd;     ///< Faraday constant scaling factor.
    extern const double Cst1;    ///< Constant used in electrolyte potential transport term.
    extern const double alp;     ///< Charge-transfer coefficient (α).

    // -------------------------------------------------------------------------
    // Initial conditions
    // -------------------------------------------------------------------------

    extern const double init_Rxn; ///< Initial reaction rate (global).
    extern const double init_RxA; ///< Initial anode reaction rate.
    extern const double init_RxC; ///< Initial cathode reaction rate.

    extern const double RT;   ///< RT constant at 300K.
    extern const double Perm; ///< Permittivity constant.

} // namespace Constants

#endif // CONSTANTS_HPP
