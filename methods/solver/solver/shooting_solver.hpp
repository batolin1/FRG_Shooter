#ifndef SHOOTING_SOLVER_HPP
#define SHOOTING_SOLVER_HPP

#include <vector>
#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"
#include "solver/solver/utils/helper_utilities.hpp"
#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "logger.hpp"

/**
    @brief A solver for the shooting problem. 
    @see Solver_Concept
*/
class Shooting_Solver {


    public:

        // Some definitions
        using Parameter = Shooting_Solver_Parameter;
        using Result_Element = Shooting_Solver_Result_Element;
        using Output_Element = Shooting_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Shooting_Solver (const Configuration& configuration) {
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

        void write_log (const Parameter& param) {
            std::ostringstream oss;
            oss << "Initiating shooting solver with the following parameters:\n" 
                << "dimension: " << param.dimension << "\n" 
                << "anomalous_dimension: " << param.anomalous_dimension << "\n" 
                << "s_factor: " << param.s_factor << "\n" 
                << "sigma_minima: " << param.sigma_minima << "\n" 
                << "sigma_maxima: " << param.sigma_maxima << "\n" 
                << "sigma_steps: " << param.sigma_steps << "\n" 
                << "maximum_iterations: " << param.maximum_iterations << "\n" 
                << "tolerance: " << param.tolerance << "\n";
            Logger::instance ().log (oss.str ());
        }

        /**
            @brief A method to process an input parameter, that is, given an input to this solver,
                   this method processes the input and inserts results as output elements. 
            
            @param param    The (physical) model parameters for this solver, for which to process.
        */
        void process_parameter (const Parameter& param) {

            write_log (param);

            std::vector<Result_Element> result;
            std::vector<Result_Element> spikes;

            for (int i = 0; i < param.sigma_steps; ++i) {

                const double sigma = param.sigma_minima +
                (param.sigma_maxima - param.sigma_minima) * i / (param.sigma_steps - 1.0);

                Potential_Integrator_Parameter integrator_parameter = 
                    parameter_builder::make_potential_integrator_parameter (param);

                integrator_parameter.sigma = sigma;
                Potential_Integrator_Model integrator_model;
                // No need to save trajectory in this case. 
                integrator_model.save_trajectory = false;
                integrator::integrate_model (configuration, 
                    integrator_parameter, integrator_model, 
                    param.maximum_iterations, param.tolerance);

                // Actually get asymptotic field at the point where cutoff occurs.
                const double asymptotic_field = integrator_model.get_result ();

                Result_Element result_element;
                result_element.asymptotic_field = asymptotic_field;
                result_element.sigma = sigma;
                // The self-consistently-calculated anomalous dimension.
                result_element.anomalous_dimension = integrator_parameter.anomalous_dimension;
                result.push_back (result_element);
            }

            // Find now the spikes.
            const std::vector<int> indexes = 
                helper_utilities::find_spike (configuration, result);

            // Store spikes to the result array. 
            for (int index : indexes) {
                spikes.push_back (result [index]); 
            }
            
            // The result elements. 
            Output_Element output_element;
            output_element.parameter = param;
            output_element.result = result;
            output_element.spikes = spikes; 
            output.push_back (output_element);
        }

    private:

        // Numerical configurations and the output. 
        Configuration configuration;
        std::vector<Output_Element> output;
};

#endif