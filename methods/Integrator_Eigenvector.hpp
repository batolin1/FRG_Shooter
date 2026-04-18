#ifndef INTEGRATOR_EIGENVECTOR_HPP
#define INTEGRATOR_EIGENVECTOR_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <boost/numeric/odeint.hpp>
#include <array>

using namespace boost::numeric::odeint;
using integrable_element = std::array<double, 2>;
using trajectory = std::vector<double>;

/**
    This class realizes the integration for the ODE corresponding to the 
    linearized eigenvector equation, given a "guess" for the eigenvalue. 
*/
class Integrator_Eigenvector {

    public:

        /**
            Constructs an instance of this class. For more detail regarding the
            input parameters, see comments on "Eigenperturbation_Method.hpp".
            Crucially, sigma is the index parameter for the particular family  
            of initial conditions that integrates the ODE corresponding to the 
            RG flow of the potential at a fixed point. It thus characterizes
            the particular fixed point being investigated. The wavefunction and
            the potential_prime trajectories correspond to the ODE evolution of
            the potential and derivatives at the fixed point of the RG flow. 
            Those trajectories are required as the ODE for the eigenvector 
            depends on those trajectories. 
            @param dimension              The dimension for the system in question. 
            @param anomalous_dimension    The anomalous dimension. 
            @param s_factor               The s-factor for the model.
            @param symmetry_factor_N      The N factor for an O(N) model.
            @param sigma                  The parameter sigma corresponding 
                                          to the solution of the ODE for the 
                                          potential RG flow at a critical pt.
            @param eigenvalue             The "guess" eigenvalue corresponding
                                          to a potential solution of this ODE.
            @param wavefunction           The trajectory for wavefunctions. 
            @param potential_0prime       The trajectory for the potential
                                          along the solution
                                          to the potential RG flow ODE. 
            @param potential_1prime       The trajectory for the potential's 
                                          first derivative along the solution
                                          to the potential RG flow ODE. 
            @param potential_2prime       The trajectory for the potential's 
                                          second derivative along the solution
                                          to the potential RG flow ODE. 
        */
        Integrator_Eigenvector (
            // Simulation parameters (from input file).
            const double dimension,
            const double anomalous_dimension,
            const double s_factor,
            const double sigma, 
            const double eigenvalue, 
            trajectory wavefunction, 
            trajectory potential_0prime, 
            trajectory potential_1prime, 
            trajectory potential_2prime) {
                this->dimension = dimension;
                this->anomalous_dimension = anomalous_dimension;
                this->s_factor = s_factor;
                this->sigma = sigma;
                this->eigenvalue = eigenvalue;
                this->wavefunction = wavefunction;
                this->potential_0prime = potential_0prime;
                this->potential_1prime = potential_1prime;
                this->potential_2prime = potential_2prime;
                // The wavefunction is initialized at an infinitesimal 
                // perturbation near zero. given input parameters, auxiliary 
                // parameters and the initial condition are calculated.
                asymptotic_wavefunction = WAVEFUNCTION_PERTURBATION;
                update_derived_parameters ();
                set_initial_condition ();
            }

        /**

            This is the implementation of the shooting method for this 
            particular problem. Given a potential solution characterized by the
            particular eigenvalue, this method calculates the asymptotic 
            eigenvector; which should be zero for physical solutions.
            @return         the asymptotic eigenvector for this family.  
        */
        double compute_asymptotic_eigenvector () {

            // Prepares the ODE integrator. 
            typedef runge_kutta_cash_karp54<integrable_element> stepper_type;
            auto stepper = make_controlled(1e-6, 1e-9, stepper_type());
            double integration_time_step = INTEGRATION_TIME_DEFAULT;

            // Sets-up the initial condition. 
            update_derived_parameters ();
            set_initial_condition ();
            asymptotic_wavefunction = WAVEFUNCTION_PERTURBATION;

            // Integrate outwards, that is, from a small wavefunction near zero 
            // towards the threhsold; or up until a divergence occurs.
            while (asymptotic_wavefunction < WAVEFUNCTION_THRESHOLD) {

                // Triggers termination event if eigenvector diverges.
                bool is_termination_event =
                    std::abs (eigenvector[1]) > PRACTICALLY_INFINITY ||
                    std::abs (eigenvector[0]) > PRACTICALLY_INFINITY;
                // If termination even is triggered, this step is completed;
                // return eigenvector at which termination occurred.
                if (is_termination_event) return eigenvector[0];
            
                // Now actually try the ODE step. If step failed, dynamically
                // reduce the timestep. WARNING: in principle this method could
                // freeze the algorithm if non-physical parameters are chosen, 
                // that the stepper cannot resolve the ODE. No error-handling
                // was implemented for this initial version of the algorithm.
                const auto ODE = [this] (
                    const integrable_element &x,
                    integrable_element &dxdt,
                    double t) {
                        ODE_eigenvector(x, dxdt, t);
                };
                const controlled_step_result step_result = stepper.try_step (
                    ODE, eigenvector, asymptotic_wavefunction, integration_time_step);
                if (step_result == success) {
                } else {
                    integration_time_step *= 0.5;
                }
            }
            return eigenvector[0];
        }
        

    private:

        // Fixed constants.
        static constexpr double PRACTICALLY_ZERO = 1e-6;
        static constexpr double PRACTICALLY_INFINITY = 1e+6;
        static constexpr double INTEGRATION_TIME_DEFAULT = 1e-3;
        static constexpr double WAVEFUNCTION_PERTURBATION = 1e-8;
        static constexpr double WAVEFUNCTION_THRESHOLD = 1000;
        static constexpr double PI = 3.14159265358979323846;

        // Simulation parameters (from input file).
        double dimension;
        double anomalous_dimension;
        double s_factor;
        double sigma;
        double eigenvalue;
        trajectory wavefunction;
        trajectory potential_0prime;
        trajectory potential_1prime;
        trajectory potential_2prime;

        // Implied parameters (given input).
        double dimension_factor;
        double anomalous_constant;
        double implied_s_factor;
        double s_constant;

        // Parameters obtained after simulation execution
        double asymptotic_wavefunction;
        integrable_element eigenvector;

        // For more details on meaning of those parameters, read article and
        // additionally "Eigenperturbation_Method.hpp".

        /**
            Given input, this method calculates the remaining parameters that 
            are implied from the input.
        */
        void update_derived_parameters () {
            dimension_factor = 4.0 / dimension / 
                std::pow (2.0, dimension + 1.0) / 
                std::pow (PI, dimension * 0.5) / 
                std::tgamma (0.5 * dimension);
            anomalous_constant =
                1.0 - anomalous_dimension / (2.0 + dimension);
            implied_s_factor = s_factor - anomalous_dimension;
            s_constant = s_factor / 2.0 * anomalous_constant;
        }

        /**
            Given a family labelled by an eigenvalue, calculates the initial
            conditions for the eigenvectors corresponding to this family.
            @return               an array containing the first and second 
                                  derivatives of the potential at the initial
                                  condition and for the family sigma. 
        */
        void set_initial_condition () {
            // The eigenvector at start is a free parameter so we fix to 1.0.
            const double eigenvector = 1.0;
            const double prefactor = 
                (1.0 + potential_0prime[0]) * 
                (1.0 + potential_0prime[0]) /
                 s_constant / dimension_factor;
            const double eigenvector_derivative = 
                -prefactor * (dimension + eigenvalue);
            const double eigenvector_2derivative = 
                0.5 * prefactor * prefactor * 
                (eigenvector + dimension) * 
                (eigenvector * implied_s_factor);
            // Introduces perturbation to the initial condition. 
            const double eigenvector_perturbed = eigenvector + 
                WAVEFUNCTION_PERTURBATION * eigenvector_derivative + 0.5 *
                std::pow (WAVEFUNCTION_PERTURBATION, 2.0) * eigenvector_2derivative;
            const double eigenvector_derivative_perturbed = eigenvector_derivative + 
                WAVEFUNCTION_PERTURBATION * eigenvector_2derivative;

            this->eigenvector = {eigenvector_perturbed, eigenvector_derivative_perturbed};
        }

        /**
            The ordinary differential equation for the eigenvectors of the RG.
            @param eigenvector               Array containing the eigenvector 
                                             and its first derivative. 
            @param eigenvector_derivative    Array for the first and second
                                             derivatives of the eigenvector.
            @param wavefunction              The wavefunction at this step. 
        */
        void ODE_eigenvector 
            (const integrable_element &eigenvector, 
            integrable_element &eigenvector_derivative,
            const double wavefunction) {
                const double potential_contribution = 
                    1.0 + get_potential_1prime (wavefunction) + 
                    2.0 * wavefunction * get_potential_2prime (wavefunction);
                const double prefactor = 
                    std::pow (potential_contribution, 2.0) / 
                    dimension_factor / 
                    s_constant;
                const double eigenvector_2derivative = - 0.5 / wavefunction * (
                    eigenvector[1] + prefactor * (
                        (dimension + eigenvalue) * eigenvector[0] -
                        (dimension - implied_s_factor) * wavefunction * eigenvector[1]
                    )
                );
                eigenvector_derivative[0] = eigenvector [1];
                eigenvector_derivative[1] = eigenvector_2derivative;
        }

        /**
            An auxiliary method to extract (by interpolation) the 
            zeroth derivative of the potential (pre-calculated).
            @param x    The value for the wavefunction at which 
                        the parameter is returned. 
            @return     The zeroth derivative of the potential at 
                        wavefunction valued x.
        */
        double get_potential_0prime (const double x) {
            return get_potential (x, 0);
        }

        /**
            An auxiliary method to extract (by interpolation) the 
            first derivative of the potential (pre-calculated).
            @param x    The value for the wavefunction at which 
                        the parameter is returned. 
            @return     The first derivative of the potential at 
                        wavefunction valued x.
        */
        double get_potential_1prime (const double x) {
            return get_potential (x, 1);
        }

        /**
            An auxiliary method to extract (by interpolation) the 
            second derivative of the potential (pre-calculated).
            @param x    The value for the wavefunction at which 
                        the parameter is returned. 
            @return     The second derivative of the potential at 
                        wavefunction valued x.
        */
        double get_potential_2prime (const double x) {
            return get_potential (x, 2);
        }

        /** 
            An interpolation method to, given a wavefunction, extract the
            potentials and its derivatives. 
            @param x             The wavefunction at which the interpolation
                                 is being evaluated. 
            @param derivative    Whether to extract zeroth, first or 
                                 second derivative. 
            @return              Returns the interpolated value for the 
                                 potential or its derivatives at x.
        */
        double get_potential (const double x, const int derivative) {
            // Firs find index closest to the solution. 
            int i = 0;
            while (i < wavefunction.size () && x > wavefunction[i]) {
                i++;
            }
            // Handles out of boundary issues.
            if (i >= wavefunction.size ()) i = wavefunction.size () - 1;
            if (i==0) i=1;
            // Chooses whether to take the zeroth, first, or second derivative.
            const trajectory* potential = nullptr;
            if (derivative == 0) {
                potential = &potential_0prime;
            } else if (derivative == 1) {
                potential = &potential_1prime;
            } else if (derivative == 2) {
                potential = &potential_2prime;
            }
            // Interpolates between (i-1) and (i). 
            const double overstep = 
                (x - wavefunction[i-1]) / 
                (wavefunction[i] - wavefunction[i-1]);
            const double y = (*potential)[i-1] + 
                overstep * ((*potential)[i] - (*potential)[i-1]);
            return y;
        }
};

#endif