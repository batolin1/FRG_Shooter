#ifndef INTEGRATOR_POTENTIAL_HPP
#define INTEGRATOR_POTENTIAL_HPP

#include "Integrator.hpp"

/**
    An integrator class to perform the integration of the RG flow of the 
    potential at the fixed point, given model parameters and an arbitrary 
    initial condition indexed by sigma. 
*/
class Integrator_Potential : public Integrator {

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
        void initialize (
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
        }

        /**
            Getter method for the wavefunction.
            @return    the wavefunction. 
        */
        trajectory& get_wavefunction () {
            return wavefunction;
        }

        /**
            Getter method for the potential.
            @return    the potential. 
        */
        trajectory& get_potential_0prime () {
            return potential_0prime;
        }

        /**
            Getter method for the first derivative of potential.
            @return    the first derivative of potential. 
        */
        trajectory& get_potential_1prime () {
            return potential_1prime;
        }

        /**
            Getter method for the second derivative of potential.
            @return    the second derivative of potential. 
        */
        trajectory& get_potential_2prime () {
            return potential_2prime;
        }


    private:

        /**
            Given a family labelled by sigma, this method calculates the 
            initial conditions for the shooting method. 
            @return         an array containing the first and second 
                            derivatives of the potential at the initial
                            condition and for the family sigma. 
        */
        void set_initial_condition () override {
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

            this->state = {potential_perturbed, potential_derivative_perturbed};
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
        void ODE_step 
            (const integrable_element &state, 
            integrable_element &state_derivative,
            const double wavefunction) override {
                const double symmetry_contribution = - (symmetry_factor_N - 1.0) * 
                    dimension_factor * s_constant / (1.0 + state[0]);
                const double denominator = 
                    symmetry_contribution + dimension * state[0] - 
                    (dimension - implied_s_factor) * state[1] * wavefunction;
                const double state_2derivative = 0.5 / wavefunction * (
                    dimension_factor * s_constant / denominator - 1.0 - state[1]);
                state_derivative[0] = state [1];
                state_derivative[1] = state_2derivative;
        }

        double result() override {
            return asymptotic_wavefunction;
        }

        bool termination_event() override {

            const double symmetry_contribution =
                - (symmetry_factor_N - 1.0) *
                dimension_factor * s_constant / (1.0 + state[0]);
        
            const double denominator =
                symmetry_contribution + dimension * state[0] -
                (dimension - implied_s_factor) * state[1] * asymptotic_wavefunction;
        
            const double state_2derivative =
                0.5 / asymptotic_wavefunction * (
                    dimension_factor * s_constant / denominator - 1.0 - state[1]
                );
        
            const double the_real_denominator =
                1.0 + state[1] +
                2.0 * asymptotic_wavefunction * state_2derivative;
        
            return
                std::abs(denominator) < PRACTICALLY_ZERO ||
                std::abs(the_real_denominator) < PRACTICALLY_ZERO ||
                std::abs(1.0 + state[1]) < PRACTICALLY_ZERO;
        }

        void on_success_step() override {

            const double symmetry_contribution =
                - (symmetry_factor_N - 1.0) *
                dimension_factor * s_constant / (1.0 + state[0]);
        
            const double denominator =
                symmetry_contribution + dimension * state[0] -
                (dimension - implied_s_factor) * state[1] * asymptotic_wavefunction;
        
            const double state_2derivative =
                0.5 / asymptotic_wavefunction * (
                    dimension_factor * s_constant / denominator - 1.0 - state[1]
                );
        
            wavefunction.push_back(asymptotic_wavefunction);
            potential_0prime.push_back(state[0]);
            potential_1prime.push_back(state[1]);
            potential_2prime.push_back(state_2derivative);
        }

};

#endif