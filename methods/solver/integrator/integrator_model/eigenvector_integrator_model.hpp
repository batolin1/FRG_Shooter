#ifndef EIGENVECTOR_INTEGRATOR_HPP
#define EIGENVECTOR_INTEGRATOR_HPP

#include <vector>
#include <cmath>

#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"

/**
    @brief A class representing an integrator model for the integral of the differential equation 
           for eigenvalues and eigenvectors near fixed point.
        
    @see Integrator_Model_Concept, Solver_Concept, integrator
*/
class Eigenvector_Integrator_Model {

    public:

        // The definition of a parameter. 
        using Parameter = Eigenvector_Integrator_Parameter;

        // The state to be evolved.
        std::array<double,2> state;

        // The field 
        double asymptotic_field;

        /**
            @brief Initialization method to start this object. 
            
            @param configuration    The (numerical) configuration parameters. 
            @param parameter        The (physical) parameters for the model. 
        */
        void initialize (const Configuration& configuration, const Parameter& parameter) {

            this->parameter = parameter;
            this-> configuration = configuration;
            asymptotic_field = configuration.field_perturbation;

            // The eigenvector at start is a free parameter so we fix to 1.0.
            const double eigenvector = 1.0;
            const double prefactor = 
                (1.0 + parameter.trajectory [0].potential_0prime) * 
                (1.0 + parameter.trajectory [0].potential_0prime) /
                parameter.s_constant / parameter.dimension_factor;
            const double eigenvector_derivative = 
                -prefactor * (parameter.dimension + parameter.eigenvalue);
            const double eigenvector_2derivative = 0.5 * prefactor * prefactor * 
                (eigenvector + parameter.dimension) * (eigenvector * parameter.implied_s_factor);
            // Introduces perturbation to the initial condition. 
            const double eigenvector_perturbed = eigenvector + 
                configuration.field_perturbation * eigenvector_derivative + 
                0.5 * std::pow (configuration.field_perturbation, 2.0) * eigenvector_2derivative;
            const double eigenvector_derivative_perturbed = 
                eigenvector_derivative + configuration.field_perturbation * eigenvector_2derivative;

            this->state = {eigenvector_perturbed, eigenvector_derivative_perturbed};
        }


        /**
            @brief The actual differential equation that this model follows. 

            @param state               The state to be evolved.
            @param state_derivative    The derivative of the state. 
            @param field               The field at point of evaluation. 
        */
        void ODE_step (
            const std::array<double,2> &state, 
            std::array<double,2> &state_derivative, 
            const double field) {

            const double potential_contribution = 
                1.0 + get_potential (field, 1) + 2.0 * field * get_potential (field, 2);
            const double prefactor = std::pow (potential_contribution, 2.0) / 
                parameter.dimension_factor / parameter.s_constant;
            const double state_2derivative = - 0.5 / field * (
                state [1] + prefactor * ((parameter.dimension + parameter.eigenvalue) * state [0] -
                (parameter.dimension - parameter.implied_s_factor) * field * state [1]));
            state_derivative [0] = state [1];
            state_derivative [1] = state_2derivative;
        }

        /**
            @brief Checker whether a termination event was triggered. 

            @return    A boolean whether termination event was triggered. 
        */
        bool termination_event() {
            const bool is_termination_event =
                std::abs (state [1]) > configuration.practically_infinity ||
                std::abs (state [0]) > configuration.practically_infinity;
            return is_termination_event; 
        }
    
        /**
            @brief Getter method for the state after the model was integrated. 

            @return    A double-value element corresponding to the terminal state. 
        */
        double get_result () {
            return state [0];
        }

        /**
            @brief Dummy method to satisfy Integrator_Model_Concept. 
        */
        void on_success_step () {}


    private:
    
        // The physical parameters for the model.
        Parameter parameter;

        // The numerical configurations for the integration. 
        Configuration configuration;

        /**
            @brief Getter method for the potential and its derivatives. 
            
            @param x             The field at point of evaluation. 
            @param derivative    The derivative in question. 
            @return              The potential and/or derivatives. 
        */
        double get_potential (const double x, const int derivative) {

            
    
            // First we already deal with boundary issues from above. If field 
            // provided above biggest maximum, use asymptotic continuation.
            if (x > parameter.trajectory.at (parameter.trajectory.size () - 1).field) {
                return parameter_builder::get_asymtotic_continuation (x, derivative, parameter);
            }
            
            // If not, first index closest to the solution. 
            int i = 0;
            while (i < parameter.trajectory.size () && x > parameter.trajectory [i].field) {
                i++;
            }
            // Handles out of boundary from below. 
            if (i==0) i=1;
            // Chooses whether to take the zeroth, first, or second derivative.
            double potential;
            double potential_previous;
            if (derivative == 0) {
                potential = parameter.trajectory [i].potential_0prime;
                potential_previous = parameter.trajectory [i-1].potential_0prime;
            } else if (derivative == 1) {
                potential = parameter.trajectory [i].potential_1prime;
                potential_previous = parameter.trajectory [i-1].potential_1prime;
            } else if (derivative == 2) {
                potential = parameter.trajectory [i].potential_2prime;
                potential_previous = parameter.trajectory [i-1].potential_2prime;
            }
            // Interpolates between (i-1) and (i). 
            const double overstep = (x - parameter.trajectory [i-1].field) / 
                (parameter.trajectory [i].field - parameter.trajectory [i-1].field);
            const double y = potential_previous + overstep * (potential - potential_previous);
            return y;
        }
};

#endif