#ifndef STANDARD_POTENTIAL_INTEGRATOR_HPP
#define STANDARD_POTENTIAL_INTEGRATOR_HPP

#include "structure.hpp"
#include <array>
#include <cmath>

/**
    @brief A class representing an integrator model for the integral of the 
           differential equation for the potential at the fixed point.
        
    @see Integrator_Model_Concept, Solver_Concept, integrator
*/
class Standard_Potential_Integrator_Model {

    public:

        // Definitions of parameter and result element. 
        using Result_Element = Potential_Integrator_Result_Element;
        using Parameter = Potential_Integrator_Parameter;

        // The field and the state to be evolved. 
        double asymptotic_field;
        std::array<double,2> state;
        double threshold;
        bool save_trajectory;

        /**
            @brief Method to get the trajectory, that is, the evolution of the potential
                   and its derivatives as a function of the field. 

            @return    A vector containing result elements. 
        */
        std::vector<Potential_Integrator_Result_Element>& get_trajectory () {

            if (!save_trajectory) {
                return trajectory;
            }

            for (int i = 1; i < trajectory.size () - 1; i++) {
                trajectory [i].potential_2prime =
                    (trajectory[i + 1].potential_1prime - trajectory[i - 1].potential_1prime) /
                    (trajectory[i + 1].field - trajectory[i - 1].field);
            }
            trajectory[0].potential_2prime = parameter.sigma;
            trajectory[trajectory.size () - 1].potential_2prime = 
                (trajectory.back ().potential_1prime - 
                trajectory.at (trajectory.size () - 2).potential_1prime) / 
                (trajectory.back ().field - trajectory.at (trajectory.size () - 2).field);
                
            for (int i = 0; i < trajectory.size (); i++) {
                const double field = trajectory[i].field;
                const double potential_0prime_value = trajectory[i].potential_0prime;
                const double potential_1prime_value = trajectory[i].potential_1prime;
                const double potential_2prime_value = trajectory[i].potential_2prime;

                trajectory [i].denominator = parameter.dimension * potential_0prime_value - 
                    0.5 * (parameter.dimension - parameter.implied_s_factor) * 
                    field * potential_1prime_value;


                const double derivative_threshold = 1e-6;
                double shape_1prime = 0.0;
                if (std::abs (potential_1prime_value) > derivative_threshold) {
                    shape_1prime = potential_2prime_value * field / potential_1prime_value;
                    const double max_exponent = 10.0;
                    const double min_exponent = -10.0;
                    shape_1prime = std::max (min_exponent, std::min (max_exponent, shape_1prime));
                }
                trajectory[i].potential_1prime_shape = shape_1prime;

                double shape_0prime = 0.0;
                if (std::abs (potential_0prime_value) > derivative_threshold) {
                    shape_0prime = potential_1prime_value * field / potential_0prime_value;
                    const double max_exponent = 10.0;
                    const double min_exponent = -10.0;
                    shape_0prime = std::max (min_exponent, std::min (max_exponent, shape_0prime));
                } else {
                    shape_0prime = shape_1prime + 1.0;
                }
                trajectory[i].potential_0prime_shape = shape_0prime;

                trajectory[i].potential_2prime_shape = shape_1prime - 1.0;
            }
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
            threshold = configuration.field_threshold;

            const double potential = parameter.dimension_factor * parameter.s_constant / 
                parameter.dimension / (1.0 + parameter.sigma);
            const double potential_1derivative = 0.0;

            this->state = {potential, potential_1derivative};
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

            const double denominator = parameter.dimension * state [0] - 
                0.5 * (parameter.dimension - parameter.implied_s_factor) * field * state [1];

            // If very big denominator, don't compute the rest, just throw big number. 
            if (denominator <= configuration.practically_zero) {
                state_derivative [0] = configuration.practically_infinity;
                state_derivative [1] = configuration.practically_infinity;
                return; 
            }
            
            const double potential_2prime = 
                parameter.dimension_factor * parameter.s_constant / denominator - 1.0;

            state_derivative [0] = state [1];
            state_derivative [1] = potential_2prime;
        }

        /**
            @brief Checker whether a termination event was triggered. 

            @return    A boolean whether termination event was triggered. 
        */
        bool termination_event () {

            const double denominator = parameter.dimension * state [0] - 
                0.5 * (parameter.dimension - parameter.implied_s_factor) * 
                asymptotic_field * state [1];

            const bool is_termination_event = 
                denominator <= configuration.practically_zero || 
                std::abs (state [0]) >= configuration.practically_infinity || 
                std::abs (state [1]) >= configuration.practically_infinity;
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

            const double denominator = parameter.dimension * state [0] - 
                0.5 * (parameter.dimension - parameter.implied_s_factor) * 
                asymptotic_field * state [1];


            Potential_Integrator_Result_Element trajectory_element;
            trajectory_element.field = asymptotic_field;
            trajectory_element.potential_0prime = state [0];
            trajectory_element.potential_1prime = state [1];
            trajectory_element.the_real_denominator = denominator;
            trajectory_element.denominator = denominator;
            trajectory.push_back (trajectory_element);
        }

        double get_threshold () {
            return threshold;
        }

    private:

        // The trajectory of the potential as function of the field.  
        std::vector<Result_Element> trajectory;
        // The (physical) model parameters and (numerical) configuration.
        Parameter parameter;
        Configuration configuration;
};

#endif