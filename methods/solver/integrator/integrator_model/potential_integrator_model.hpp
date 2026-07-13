#ifndef POTENTIAL_INTEGRATOR_HPP
#define POTENTIAL_INTEGRATOR_HPP

#include "structure.hpp"
#include <array>
#include <cmath>

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
        double threshold;
        bool save_trajectory;
        double latest_3derivative;

        /**
            @brief Method to get the trajectory, that is, the evolution of the potential
                   and its derivatives as a function of the field. 

            @return    A vector containing result elements. 
        */
        std::vector<Potential_Integrator_Result_Element>& get_trajectory () {

            if (!save_trajectory) {
                return trajectory;
            }

            double potential_0prime = 0;
            trajectory[0].potential_0prime = 0;
            for (int i = 1; i < trajectory.size (); i++) {
                potential_0prime +=
                    (trajectory[i].field - trajectory[i - 1].field) *
                    (trajectory[i].potential_1prime + trajectory[i - 1].potential_1prime) / 2;
                trajectory[i].potential_0prime = potential_0prime;
            }

            for (int i = 0; i < trajectory.size (); i++) {
                const double field = trajectory[i].field;
                const double potential_0prime_value = trajectory[i].potential_0prime;
                const double potential_1prime_value = trajectory[i].potential_1prime;
                const double potential_2prime_value = trajectory[i].potential_2prime;

                trajectory[i].denominator =
                    1.0 + potential_1prime_value + 2.0 * field * potential_2prime_value;

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
            latest_3derivative = 0.0;

            const double factor = parameter.s_constant * parameter.dimension_factor;

            const double potential_derivative = parameter.sigma;
            const double potential_2derivative = - parameter.implied_s_factor *  
                parameter.sigma * (1.0 + parameter.sigma) * (1.0 + parameter.sigma) / 
                3.0 / parameter.dimension_factor / parameter.s_constant;
                
            this->state = {potential_derivative, factor * potential_2derivative};
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

            const double term = - parameter.implied_s_factor * state [0] + 
                (parameter.dimension - parameter.implied_s_factor) * field * state [1];
            const double denominator = 1 + state [0] + 2 * field * state [1];

            // If very big denominator, don't compute the rest, just throw big number. 
            if (denominator <= configuration.practically_zero) {
                state_derivative [0] = configuration.practically_infinity;
                state_derivative [1] = configuration.practically_infinity;
                return; 
            }

            //const double prefactor = 1.0 / parameter.dimension_factor / parameter.s_constant;
            const double prefactor= 1.0;
            const double state_3derivative = 0.5 / field * (
                term * denominator * denominator * prefactor - 3 * state [1]);

            latest_3derivative = state_3derivative;
            state_derivative [0] = state [1];
            state_derivative [1] = state_3derivative;
        }

        /**
            @brief Checker whether a termination event was triggered. 

            @return    A boolean whether termination event was triggered. 
        */
        bool termination_event () {
            
            const double denominator =
                1.0 + state [0] + 2.0 * asymptotic_field * state [1];

            // When denominator goes to infinity, the RHS of the equation goes to zero and the 
            // asymptotic limit is reached. This asymptote must NEVER be crossed. Thus, this term 
            // must always stay positive.   
            const double asymptotic_term = parameter.implied_s_factor * state [0] + 
                (parameter.dimension - parameter.implied_s_factor) * asymptotic_field * state[1];
            const bool is_termination_event = 
                //std::abs (asymptotic_term) <= configuration.practically_zero ||
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

            const double the_real_denominator = 
                1 + state [0] + 2 * asymptotic_field * state [1]; 

            Potential_Integrator_Result_Element trajectory_element;


            trajectory_element.potential_0prime =
                trajectory.empty() ? 0.0 : 
                    trajectory.back().potential_0prime + 0.5 * 
                    (state[0] + trajectory.back ().potential_1prime) *
                    (asymptotic_field - trajectory.back().field);


            trajectory_element.field = asymptotic_field;
            trajectory_element.potential_1prime = state [0];
            trajectory_element.potential_2prime = state [1];
            trajectory_element.potential_3prime = latest_3derivative;
            trajectory_element.the_real_denominator = the_real_denominator;


            trajectory_element.denominator = the_real_denominator;

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