#ifndef INTEGRATOR_POTENTIAL_FLOW_HPP
#define INTEGRATOR_POTENTIAL_FLOW_HPP

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
    An integrator class to perform the integration of the RG flow of the 
    potential at the fixed point, given model parameters and an arbitrary 
    initial condition indexed by sigma. 
*/
class Integrator_Potential_Flow {

    public:

        /**
            Constructs an instance of this class. For more details on the input
            parameters, read comments on "Shooting_Method.hpp".
            @param dimension              The dimension of the model.
            @param anomalous_dimension    The anomalous dimension of the model.
            @param s_factor               The s-factor for the model.
            @param symmetry_factor_N      The N factor for an O(N) model.
            @param sigma                  The parameter sigma which indexes the 
                                          initial condition for this family.
        */
        Integrator_Potential_Flow (
            const double dimension,
            const double anomalous_dimension,
            const double s_factor,
            const double symmetry_factor_N,
            const double sigma) {
                this->dimension = dimension;
                this->anomalous_dimension = anomalous_dimension;
                this->s_factor = s_factor;
                this->symmetry_factor_N = symmetry_factor_N;
                this->sigma = sigma;
                // Initializes wavefunction at a tiny value near zero.
                asymptotic_wavefunction = WAVEFUNCTION_PERTURBATION;
                // The initial condition and other parameters are calculated
                // from the given input.
                update_derived_parameters ();
                set_initial_condition ();
            }

        /**
            This is the implementation of the shooting method. Given some 
            family labelled by sigma, this method calculates the asymptotic
            wavefunction; that is to say, given sigma, the maximum wavefunction
            reachable before divergence (or the threshold itself, if it is 
            reached, that is, if no divergence occurs.).   
            @return         the asymptotic wavefunction for this family.  
        */
        double compute_asymptotic_wavefunction () {

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
                // Here we verify if method triggered an asympthotic event.
                // 1) Compute all the denominators.
                const double symmetry_contribution = 
                    - (symmetry_factor_N - 1.0) * 
                    dimension_factor * s_constant / 
                    (1.0 + potential[0]);
                const double denominator = 
                    symmetry_contribution + dimension * potential[0] - 
                    (dimension - implied_s_factor) * asymptotic_wavefunction * potential[1];
                const double potential_2derivative = 0.5 / asymptotic_wavefunction * (
                    dimension_factor * s_constant / denominator - 1 - potential[1]);
                const double the_real_denominator = 1 + potential[1] + 
                    2.0 * asymptotic_wavefunction * potential_2derivative;
                // 2) Triggers termination if either denominator diverges.
                bool is_termination_event =
                    std::abs (denominator) < PRACTICALLY_ZERO ||
                    std::abs (the_real_denominator) < PRACTICALLY_ZERO ||
                    std::abs (1.0 + potential[1]) < PRACTICALLY_ZERO;
                // 3) If termination even is triggered, this step is completed.
                //    Return wavefunction at which termination occurred.
                if (is_termination_event) return asymptotic_wavefunction;
            
                // Now actually try the ODE step. If step failed, dynamically
                // reduce the timestep. WARNING: in principle this method could
                // freeze the algorithm if non-physical parameters are chosen, 
                // that the stepper cannot resolve the ODE. No error-handling
                // was implemented for this initial version of the algorithm.
                const auto ODE = [this] (
                    const integrable_element &x,
                    integrable_element &dxdt,
                    double t) {
                        ODE_flow_potential(x, dxdt, t);
                };
                const controlled_step_result step_result = stepper.try_step (
                    ODE, potential, asymptotic_wavefunction, integration_time_step);

                // Exclusively for successful results, stores the trajectories
                // for the potential and its derivatives as function of the 
                // wavefunction. 
                if (step_result == success) {
                    wavefunction.push_back (asymptotic_wavefunction);
                    potential_0prime.push_back (potential[0]);
                    potential_1prime.push_back (potential[1]);
                    potential_2prime.push_back (
                        0.5 / asymptotic_wavefunction * (
                        dimension_factor * s_constant / denominator
                        - 1 - potential[1])
                    );
                // Dinamically reduces step if needed.
                } else {
                    integration_time_step *= 0.5;
                }
            }
            // No divergence? - Return the threshold. 
            return asymptotic_wavefunction;
        }

        /**
            Getter method for the wavefunction.
            @return    the wavefunction. 
        */
        trajectory get_wavefunction () {
            return wavefunction;
        }

        /**
            Getter method for the potential.
            @return    the potential. 
        */
        trajectory get_potential_0prime () {
            return potential_0prime;
        }

        /**
            Getter method for the first derivative of potential.
            @return    the first derivative of potential. 
        */
        trajectory get_potential_1prime () {
            return potential_1prime;
        }

        /**
            Getter method for the second derivative of potential.
            @return    the second derivative of potential. 
        */
        trajectory get_potential_2prime () {
            return potential_2prime;
        }

    private:

        // Fixed constants.
        static constexpr double PRACTICALLY_ZERO = 1e-6;
        static constexpr double INTEGRATION_TIME_DEFAULT = 1e-3;
        static constexpr double WAVEFUNCTION_PERTURBATION = 1e-8;
        static constexpr double WAVEFUNCTION_THRESHOLD = 1000;
        static constexpr double PI = 3.14159265358979323846;

        // Simulation parameters (from input file).
        double dimension;
        double anomalous_dimension;
        double s_factor;
        double symmetry_factor_N;
        double sigma;

        // Implied parameters (given input).
        double dimension_factor;
        double anomalous_constant;
        double implied_s_factor;
        double s_constant;

        // Parameters obtained after simulation execution
        double asymptotic_wavefunction;
        integrable_element potential;
        trajectory wavefunction;
        trajectory potential_0prime;
        trajectory potential_1prime;
        trajectory potential_2prime;

        // For more details on meaning of those parameters, read article and
        // additionally "Shooting_Method.hpp".

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
            Given a family labelled by sigma, this method calculates the 
            initial conditions for the shooting method. 
            @return         an array containing the first and second 
                            derivatives of the potential at the initial
                            condition and for the family sigma. 
        */
        void set_initial_condition () {
            const double potential = 
                dimension_factor * s_constant * 
                symmetry_factor_N / dimension / 
                (1.0 + sigma);
            const double potential_derivative = sigma;
            const double potential_2derivative = 
                - implied_s_factor * (1.0 + sigma) * (1.0 + sigma) * 
                sigma / (s_constant * dimension_factor);
            // Introduces perturbation to the initial condition.
            const double potential_perturbed = potential + 
                WAVEFUNCTION_PERTURBATION * potential_derivative + 0.5 *
                std::pow (WAVEFUNCTION_PERTURBATION, 2.0) * potential_2derivative;
            const double potential_derivative_perturbed = potential_derivative + 
                WAVEFUNCTION_PERTURBATION * potential_2derivative;

            this->potential = {potential_perturbed, potential_derivative_perturbed};
        }

        /**
            The ordinary differential equation corresponding to the RG flow
            flow of the potential and its derivatives. In other words, this
            method performs one step along the ODE. 
            @param potential               An array containing the potential 
                                           and its first derivative. 
            @param potential_derivative    An array containing the first and
                                           second derivatives of the potential.
            @param wavefunction            The wavefunction at this step. 
        */
        void ODE_flow_potential 
            (const integrable_element &potential, 
            integrable_element &potential_derivative,
            const double wavefunction) {
                const double symmetry_contribution = - (symmetry_factor_N - 1.0) * 
                    dimension_factor * s_constant / (1.0 + potential[0]);
                const double denominator = 
                    symmetry_contribution + dimension * potential[0] - 
                    (dimension - implied_s_factor) * potential[1] * wavefunction;
                const double potential_2derivative = 0.5 / wavefunction * (
                    dimension_factor * s_constant / denominator - 1.0 - potential[1]);
                potential_derivative[0] = potential [1];
                potential_derivative[1] = potential_2derivative;
        }
};

#endif