/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef CHEMGCSOLVER_H
#define CHEMGCSOLVER_H

#include <array>
#include <vector>

#include <opm/simulators/geochemistry/Core/ChemBasVec.h>
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>
#include <opm/simulators/geochemistry/Core/SurfaceChemistry.hpp>

#include <opm/simulators/geochemistry/Common/Enums.hpp>
#include <opm/simulators/geochemistry/Common/ChemGlobal.h>
#include <opm/simulators/geochemistry/Common/Constants.hpp>
#include <opm/simulators/geochemistry/Common/CustomExceptions.hpp>
#include <opm/simulators/geochemistry/Numerical/GaussLegendre.hpp>
#include <opm/simulators/geochemistry/Numerical/SecantMethod.hpp>
//#include <opm/simulators/geochemistry/Numerical/BrentMethod.hpp>
#include <opm/simulators/geochemistry/SplayTree/SplayTree.h>
#include <opm/simulators/geochemistry/SplayTree/TreeNode.h>
#include <opm/simulators/geochemistry/Thermo/hkf.h>
#include <opm/simulators/geochemistry/Numerical/BrentMethod.hpp>
#include <opm/simulators/geochemistry/Numerical/RootSolvers.hpp>

struct GCSolverOptions
{
    int CALC_JACOBI_NUM_{0};

    double MBAL_CONV_CRITERION_{1e-5};
    double MBAL_MAX_STEP_{1.2};  // Max length of log10 step in each iteration for mass balance
    double MBAL_SURF_MAX_STEP_{1.2};  // Max length of log10 step in each iteration for surface potential
    int CHEM_MAXITER_MASS_{90};  // Max number of iterations in the mass balance routine
    double MBAL_SCALE_SURF_{1.0};

    double PH_CONV_CRITERION_{1e-5};
    double DELTA_PH_CONV_CRITERION{1.0e-3};
    double PH_MAX_STEP_{0.5};  // Max length of pH step in each iteration
    double PH_SECANT_METHOD_X1_FAC{1.3};  // Start secant method with x0=current_pH, x1=X1_fac*current_pH
    int CHEM_MAX_PH_ITER_{150};  // Max number of iterations in the charge balance solver
};



class GCSolver
{

  public:

    GCSolver(GeoLogger* logger, BasVec* Vchem, int sizeB, int sizeM, int MaxSizeBasis, int MaxSizeComp, int MaxSizeMin);
    ~GCSolver();

    void setLogger(GeoLogger* logger)
    {
        logger_ = logger;
    }

    int noSolverCalls() const { return solver_state_.noCalls_; }
    bool massBalanceHasConverged() const;

    void updateCurrentTimestep(double dt) { dt_ = dt; };

    bool SolveChem(); // Solve chemistry as read from the input file

    void write_convergence_status(const std::string& name, bool converged);
    // Solver used by splay tree    
    int SolveChem_ST(int key,
                     double dt,
                     const std::array<double, 3>& mass_phase,
                     const double* c_bas,
                     const double* c_ads,
                     double* c_min,
                     double* log_a_min,
                     double* output_values);  // dc, dq, (sigma, psi), pH)

    void writeCurrentBasVecToFile(const std::string& file_name);
    void write_logK_to_file_T(double T_start, double T_final, double dT, double P, const std::string& file_name, const int type);

    // Concentrations at the start of a time step
    RealVec c0_;  // Tot aqueous species concentration
    RealVec c0_ads_;  // TOT adsorbed species
    RealVec c0_min_;  // mineral concentration

    // Concentrations at end of a time step
    RealVec cT_;  // TOT aqueous species concentrations
    RealVec cT_ads_;  // TOT adsorbed species
    RealVec cT_ads_io_;  // adsorbed species - ion exchange
    RealVec cT_ads_sc_;  // TOT adsorbed species - surface complexes
    RealVec cT_DL_;  // TOT conc in Diffusive Layer
    RealVec cT_min_;  // mineral concentration
    RealVec SI_min_;  // log Saturation index minerals
    std::array<double, 3> mass_phase_;  // mass of water, oil, gas in kg

    RealVec alpha_a_gas_phase_;
    double b_gas_phase_{ 0. };
    double alpha_a_gas_phase_tot_{ 0. };
    double Z_factor_gas_{ 1. };
    int sizeB_; // size of basis vector
    int sizeM_; // size of mineral vector
    int MaxSizeBasis_;
    int MaxSizeComp_;
    int MaxSizeMin_;

    // Various options, flags, etc.
    GCSolverOptions options_;

    BasVec* Vchem_{nullptr};
    std::unique_ptr<hkf> HKF_EOS_{nullptr};

    // =================== Helper stuff to keep track of convergence behaviour ===================
    enum class ConvergenceStatus{ NOT_ATTEMPTED_SOLVE, CONVERGED, NOT_CONVERGED };

    struct
    {
        int noCalls_ {0};
        int noIterMassBalance_{0};
        ConvergenceStatus massBalanceStatus_ { ConvergenceStatus::NOT_ATTEMPTED_SOLVE };
        int noIterChargeBalance_{0};
        ConvergenceStatus chargeBalanceStatus_ { ConvergenceStatus::NOT_ATTEMPTED_SOLVE };

        // Go back to this value in case of no convergence
        double cached_best_pH_{ 7.0 };
        double cached_best_f_pH_{ 1.0e8 };

        void reset_mbal()
        {
            noIterMassBalance_ = 0;
            massBalanceStatus_ = ConvergenceStatus::NOT_ATTEMPTED_SOLVE;
        }

        void reset_pH()
        {
            noIterChargeBalance_ = 0;
            chargeBalanceStatus_ = ConvergenceStatus::NOT_ATTEMPTED_SOLVE;
            cached_best_pH_ = 7.0;
            cached_best_f_pH_ = 1.0e8;
        }

    } solver_state_;

    DiffuseLayerProperties diffuse_layer_props_{};

  private:

    GeoLogger* logger_{nullptr};
    static constexpr auto numerical_integrator_ = GaussLegendreIntegrator<double, 10>();

    // Current solver state
    RealVec mpow_;
    RealVec smpow_;
    RealVec Fchem_;
    double** Jacobi_{nullptr};
    double** Jacobi_dum_{nullptr};
    RealVec wksp_linalg_{};

    double dt_{0.0};  // For mineral kinetics

    RealVec x0_;
    RealVec x_init_;
    RealVec delta_calc_;
    std::vector<int> chem_indx_;

    bool solve_chemistry();

    void newton_ah_log(const RealVec& x, RealVec& Fchem, double** jacobi);
    void mass_balance();
    double f_pH_charge_balance(double x);

    void add_F_rock_min(int j, RealVec& Fchem, const RealVec& delta_calc);
    void calc_F_spec(int j, int pos, const RealVec& mpow, const RealVec& smpow, RealVec& Fchem) const;
    void calc_F_sup_min(int j, int pos_j, RealVec& Fchem);
    void calc_F_gas_phase(int j, int pos_j, RealVec& Fchem);
    void calc_F_gas(int j, RealVec& Fchem);
    void calc_F_grahame(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow);
    void calc_F_grahame_simple(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow);
    void calc_F_grahame_symmetrical(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow);

    void calc_jacobi_spec(int j, int pos, int q, int pos_q, double** jacobi, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_rock_min(int j, int q, double** jacobi, double** jac_dum) const;
    void calc_jacobi_sup_min(int j, int pos, int q, int pos_q, double** jacobi) const;
    void calc_jacobi_gas(int j, int q, double** jacobi, double** jac_dum) const;
    void calc_jacobi_gas_phase(int j, int pos, int q, int pos_q, double** jacobi) const;

    void calc_jacobi_grahame_wrt_E0(double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_grahame_wrt_E0_simple(double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_grahame_symmetrical_wrt_E0(double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;

    void calc_jacobi_grahame_wrt_mq(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_grahame_wrt_mq_simple(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_grahame_symmetrical_wrt_mq(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;

    void calc_surface_charge(const RealVec& mpow, const RealVec& smpow);
    double res_grahame_simple();

    // EOS for gases

    void update_fugacity_and_mol_volume();
    void update_fugacity_and_mol_volume_mixtures();

    void calculate_mixing_parameters_EOS();
    double alpha_EOS(double T, double omega, double Tc);
    double A_EOS(double P, double T, double omega, double Pc, double Tc);
    double B_EOS(double P, double T, double Pc, double Tc);
    double small_a_EOS(double Tc, double Pc);
    double small_b_EOS(double Tc, double Pc);


    // Diffuse layer calculations
    RealVec c_dl_excess_;
    void calc_diffuse_layer_conc(double* c_dl) const;
    double calc_kappa() const;
    double update_monovalent_divalent_charge_factor(double E); /** Also returns the updated value. */
    void integrate_diffuse_layer_for_mono_or_divalent_ions();
    void integrate_diffuse_layer_for_mono_or_divalent_ions_using_analytical_expressions(double E);
    void integrate_diffuse_layer_for_mono_or_divalent_ions_using_aksel_hiorth_method(double E);
    void integrate_diffuse_layer_for_mono_or_divalent_ions_using_borkovec_westall_method(double E);
    void integrate_diffuse_layer_for_mono_or_divalent_ions_using_mean_potential_approximation(double E);

    // Functions used to calculate gi-factors (surface excess in diffuse layer)
    double integrate_gi_ah(double a, double b, double charge);
    double integrate_gi_borkovec_westall(double a, double b, double charge);
    double diffuse_layer_integrand(double E, double charge);
    double diffuse_layer_integrand_simple(double E, double charge);

    // Diffuse layer contributions to Jacobian and residual?
    void calc_F_diffuse_layer_cbal_simple(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow);
    void calc_jacobi_diffuse_layer_cbal_simple_wrt_mq(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;
    void calc_jacobi_diffuse_layer_cbal_simple_wrt_E0(double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const;

    void set_Jacobian_to_zero();
    void write_jacobi(const RealVec& xinit, RealVec& Fchem, double** jc, int size, const char* fname);

    // Debug
    void calc_analytical_and_numerical_jacobian(RealVec& x);
    void write_mass_balance_debug_info_to_screen(DebugInfoTypeMBAL type);

};
#endif
