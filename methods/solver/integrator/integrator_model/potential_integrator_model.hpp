#ifndef POTENTIAL_INTEGRATOR_HPP
#define POTENTIAL_INTEGRATOR_HPP

#include "structure.hpp"
#include <array>

/**
    @brief A class representing an integrator model for the integral of the 
           differential equation for the potential at the fixed point.
        
    @see Integrator_Model_Concept, Solver_Concept, integrator
*/
class Potential_Integrator_Model {

    public:

        // Definitions of parameter and result element. 
        using Result_Element = Potential_Integrator_Result_Element;
        using Parameter = Potential_Integrator_Parameter;

        // The field and the state to be evolved. 
        double asymptotic_field;
        std::array<double,2> state;

        /**
            @brief Method to get the trajectory, that is, the evolution of the potential
                   and its derivatives as a function of the field. 

            @return    A vector containing result elements. 
        */
        std::vector<Potential_Integrator_Result_Element>& get_trajectory () {
            return trajectory;
        }
    
        /**
            @brief Initialization method to start this object. 
            
            @param configuration    The (numerical) configuration parameters. 
            @param parameter        The (physical) parameters for the model. 
        */
        void initialize (const Configuration& configuration, const Parameter& parameter) {

            this->parameter = parameter;
            this->configuration = configuration; 
            trajectory.clear ();
            asymptotic_field = configuration.field_perturbation;

            const double potential = parameter.dimension_factor * parameter.s_constant * 
                parameter.symmetry_factor_N / parameter.dimension / (1.0 + parameter.sigma);
            const double potential_derivative = parameter.sigma;
            const double potential_2derivative = 
                - parameter.implied_s_factor *(1.0 + parameter.sigma) * (1.0 + parameter.sigma) * 
                parameter.sigma / (parameter.s_constant * parameter.dimension_factor);
            // Introduces perturbation to the initial condition.
            const double potential_perturbed = 
                potential + configuration.field_perturbation * potential_derivative + 
                0.5 * std::pow (configuration.field_perturbation, 2.0) * potential_2derivative;
            const double potential_derivative_perturbed = 
                potential_derivative + configuration.field_perturbation * potential_2derivative;
                
            this->state = {potential_perturbed, potential_derivative_perturbed};
        }

        /**
            @brief The actual differential equation that this model follows. 

            @param state               The state to be evolved.
            @param state_derivative    The derivative of the state. 
            @param field               The field at point of evaluation. 
        */
        void ODE_step 
            (const std::array<double,2> &state, 
            std::array<double,2> &state_derivative,
            const double field) {

            const double symmetry_contribution = - (parameter.symmetry_factor_N - 1.0) * 
                parameter.dimension_factor * parameter.s_constant / (1.0 + state [0]);
            const double denominator = symmetry_contribution + parameter.dimension * state [0] - 
                (parameter.dimension - parameter.implied_s_factor) * state [1] * field;
            const double state_2derivative = 0.5 / field * (
                parameter.dimension_factor * parameter.s_constant / denominator - 1.0 - state [1]);
            state_derivative [0] = state [1];
            state_derivative [1] = state_2derivative;
        }

        /**
            @brief Checker whether a termination event was triggered. 

            @return    A boolean whether termination event was triggered. 
        */
        bool termination_event () {

            const double symmetry_contribution = - (parameter.symmetry_factor_N - 1.0) * 
                parameter.dimension_factor * parameter.s_constant / (1.0 + state [0]);
            const double denominator = symmetry_contribution + parameter.dimension * state [0] -
                (parameter.dimension - parameter.implied_s_factor) * state [1] * asymptotic_field;
            const double state_2derivative =
                    0.5 / asymptotic_field * (parameter.dimension_factor * 
                    parameter.s_constant / denominator - 1.0 - state [1]);
            const double the_real_denominator =
                1.0 + state [1] + 2.0 * asymptotic_field * state_2derivative;
            const bool is_termination_event = 
                std::abs (denominator) < configuration.practically_zero ||
                std::abs (the_real_denominator) < configuration.practically_zero ||
                std::abs (1.0 + state [1]) < configuration.practically_zero;
            return is_termination_event;
        }

        /**
            @brief Getter method for the state after the model was integrated. 

            @return    A double-value element corresponding to the terminal state. 
        */
        double get_result ()  {
            return asymptotic_field;
        }

        /**
            @brief Method to define  an action in case of successful step. 
        */
        void on_success_step () {

            const double symmetry_contribution = - (parameter.symmetry_factor_N - 1.0) *
                parameter.dimension_factor * parameter.s_constant / (1.0 + state [0]);
            const double denominator = symmetry_contribution + parameter.dimension * state [0] -
                (parameter.dimension - parameter.implied_s_factor) * state [1] * asymptotic_field;
            const double state_2derivative = 
                0.5 / asymptotic_field * (parameter.dimension_factor * 
                parameter.s_constant / denominator - 1.0 - state [1]);
            const double the_real_denominator = 
                1 + state [1] + 2 * asymptotic_field * state_2derivative; 

            Potential_Integrator_Result_Element trajectory_element;
            trajectory_element.field = asymptotic_field;
            trajectory_element.potential_0prime = state [0];
            trajectory_element.potential_1prime = state [1];
            trajectory_element.potential_2prime = state_2derivative;
            trajectory_element.denominator = denominator;
            trajectory_element.the_real_denominator = the_real_denominator;
            trajectory.push_back (trajectory_element);
        }

    private:

        // The trajectory of the potential as function of the field.  
        std::vector<Result_Element> trajectory;
        // The (physical) model parameters and (numerical) configuration.
        Parameter parameter;
        Configuration configuration;
};

#endif