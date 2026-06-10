#ifndef INITIAL_CONDITION_SOLVER_HPP
#define INITIAL_CONDITION_SOLVER_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"
#include "solver/solver/utils/helper_utilities.hpp"

#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"

/**
    @brief A solver for the initial condition problem. 
    @see Solver_Concept
*/
class Initial_Condition_Solver {
        
    public:

        // some definitions.
        using Parameter = Initial_Condition_Solver_Parameter;
        using Result_Element = Initial_Condition_Solver_Result_Element;
        using Output_Element = Initial_Condition_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Initial_Condition_Solver (const Configuration& configuration) {
                this->configuration = configuration;
        }

        /**
            @brief getter method for the results of the solver. 

            @return    vector of output elements corresponding to solver results. 
        */
        std::vector<Output_Element>& get_result () {
            return output;
        }

        /**
            @brief Initialization method to reset the solver. 
        */
        void initialize () {
            output.clear ();
        }

        /**
            @brief A method to process an input parameter, that is, given an input to this solver,
                   this method processes the input and inserts results as output elements. 
            
            @param param    The (physical) model parameters for this solver, for which to process.
        */
        void process_parameter (const Parameter& param) {

            std::vector<Result_Element> result;
            const std::vector<double> sigma_range = helper_utilities::populate_range (param);

            for (int s = 0; s < param.s_factor_delta; s++) {
                const double s_factor = param.s_factor_minima + 
                    (param.s_factor_maxima - param.s_factor_minima) * 
                    s / (param.s_factor_delta - 1.0);

                Potential_Integrator_Parameter integrator_parameter = 
                    parameter_builder::make_potential_integrator_parameter (param);
                integrator_parameter.s_factor = s_factor;

                
                std::vector<double> asymptotic_field;
                std::vector<double> recalculated_anomalous_dimension;

                for (double sigma : sigma_range) {

                    integrator_parameter.sigma = sigma;
                    Potential_Integrator_Model integrator_model;
                    integrator::integrate 
                        (configuration, integrator_parameter, integrator_model);
                    // integrator::integrate (configuration, integrator_parameter, integrator_model);
                    asymptotic_field.push_back (integrator_model.get_result ());
                    recalculated_anomalous_dimension.push_back 
                        (integrator_parameter.anomalous_dimension);
                }

                const std::vector<int> indicies = 
                    helper_utilities::find_spike (sigma_range, asymptotic_field);

                // Needs to create the element. 
                for (int index : indicies) {
                    Result_Element result_element;
                    result_element.spike = sigma_range.at (index + 1);
                    result_element.s_factor = s_factor;
                    result_element.anomalous_dimension = 
                        recalculated_anomalous_dimension.at (index + 1);
                    result.push_back (result_element);
                }
            }
            // Inserts output element to results. 
            Output_Element output_element;
            output_element.result = result;
            output_element.parameter = param;
            output.push_back (output_element);
        }

    private:

        // Numerical configurations and the output. 
        Configuration configuration;
        std::vector<Output_Element> output;
};

#endif