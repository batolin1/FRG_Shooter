#ifndef GRID_SEARCH_SOLVER_HPP
#define GRID_SEARCH_SOLVER_HPP

#include <vector>

#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"

#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "logger.hpp"

/**
    @brief A solver for the shooting problem. 
    @see Solver_Concept
*/
class Grid_Search_Solver {


    public:

        // Some definitions
        using Parameter = Grid_Search_Solver_Parameter;
        using Result_Element = Grid_Search_Solver_Result_Element;
        using Output_Element = Grid_Search_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Grid_Search_Solver (const Configuration& configuration) {
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

        std::vector<double> run_for_single_iteration (
            const Parameter& param, const double sigma_minima, const double sigma_maxima, 
            double& anomalous_dimension_minima, double& anomalous_dimension_maxima) {

            double best_sigma = 0;
            double best_anomalous_dimension = 0;
            double best_asymptotic_field = 0;

            // If we do NOT want to grid-search the anomalous dimension, explicitly set it here. 
            int anomalous_dimension_steps = param.anomalous_dimension_steps;
            if (!param.search_anomalous_dimension) {
                anomalous_dimension_minima = param.anomalous_dimension_minima;
                // Minima = maxima, intentional.
                anomalous_dimension_maxima = param.anomalous_dimension_minima;
                // Only a single anomalous dimension value in "grid".
                anomalous_dimension_steps = 1;
            }

            int count = 0;
            const int total_count = (int) (param.sigma_steps * param.anomalous_dimension_steps);
            const int interval = (int) (0.01 * total_count);

            for (int i = 0; i < param.sigma_steps; i++) {
                for (int j = 0; j < param.anomalous_dimension_steps; j++) {

                    const double sigma = sigma_minima + 
                        (sigma_maxima - sigma_minima) * i / (param.sigma_steps - 1.0);
                    const double anomalous_dimension = anomalous_dimension_steps == 1.0 ? 
                        anomalous_dimension_minima :
                        anomalous_dimension_minima + 
                        (anomalous_dimension_maxima - anomalous_dimension_minima) * 
                        j / (anomalous_dimension_steps - 1.0);
                    Potential_Integrator_Parameter integrator_parameter = 
                        parameter_builder::make_potential_integrator_parameter (param);
                    integrator_parameter.sigma = sigma;
                    integrator_parameter.anomalous_dimension = anomalous_dimension;
                    Potential_Integrator_Model integrator_model;
                    // No need to save trajectory in this case. 
                    integrator_model.save_trajectory = false;
                    integrator::integrate (configuration, integrator_parameter, integrator_model);

                    // Actually get asymptotic field at the point where cutoff occurs.
                    const double asymptotic_field = integrator_model.get_result ();

                    Result_Element result_element;
                    result_element.asymptotic_field = asymptotic_field;
                    result_element.sigma = sigma;
                    result_element.anomalous_dimension = anomalous_dimension;
                    result.push_back (result_element);   
                    
                    // Update best estimate for deciding the window again. 
                    if (asymptotic_field > best_asymptotic_field) {
                        best_sigma = sigma;
                        best_anomalous_dimension = anomalous_dimension;
                        best_asymptotic_field = asymptotic_field;
                    }
                }
            }

            std::vector<double> best_estimate;
            best_estimate.push_back (best_sigma);
            best_estimate.push_back (best_anomalous_dimension);
            best_estimate.push_back (best_asymptotic_field);
            return best_estimate;

        }

        void write_log (const Parameter& param) {
            std::ostringstream oss;
            const std::string search_anomalous_dimension = 
                param.search_anomalous_dimension ? "true" : "false"; 
            oss << "Initiating grid search with the following parameters:" << "\n"
                << "dimension: " << param.dimension << "\n"
                << "s_factor: " << param.s_factor << "\n"
                << "sigma_minima: " << param.sigma_minima << "\n"
                << "sigma_maxima: " << param.sigma_maxima << "\n"
                << "sigma_steps: " << param.sigma_steps << "\n"
                << "anomalous_dimension_minima: " << param.anomalous_dimension_minima << "\n"
                << "anomalous_dimension_maxima: " << param.anomalous_dimension_maxima << "\n"
                << "anomalous_dimension_steps: " << param.anomalous_dimension_steps << "\n"
                << "number_of_iterations: " << param.number_of_iterations << "\n"
                << "search anomalous dimension: " << search_anomalous_dimension << "\n"
                << "window_range: " << param.window_range << "\n";
            Logger::instance ().log (oss.str ());
        }

        /**
            @brief A method to process an input parameter, that is, given an input to this solver,
                   this method processes the input and inserts results as output elements. 
            
            @param param    The (physical) model parameters for this solver, for which to process.
        */
        void process_parameter (const Parameter& param) {

            write_log (param);

            // Re-initializes the results. 
            result.clear ();

            double sigma_minima = param.sigma_minima;
            double sigma_maxima = param.sigma_maxima;
            double anomalous_dimension_minima = param.anomalous_dimension_minima;
            double anomalous_dimension_maxima = param.anomalous_dimension_maxima;

            for (int iteration = 1; iteration < param.number_of_iterations + 1; iteration++) {
                const std::vector<double> best_estimate = run_for_single_iteration (
                    param, sigma_minima, sigma_maxima, 
                    anomalous_dimension_minima, anomalous_dimension_maxima);

                // At end of each iteration, update the bounds. 
                const double sigma_range = sigma_maxima - sigma_minima;
                const double anomalous_dimension_range = 
                    anomalous_dimension_maxima - anomalous_dimension_minima;
                const double sigma_delta = 
                    0.5 * param.window_range * sigma_range;
                const double anomalous_dimension_delta = 
                    0.5 * param.window_range * anomalous_dimension_range;
                sigma_minima = best_estimate[0] - sigma_delta;
                sigma_maxima = best_estimate[0] + sigma_delta;
                anomalous_dimension_minima = best_estimate[1] - anomalous_dimension_delta;
                anomalous_dimension_maxima = best_estimate[1] + anomalous_dimension_delta;
                // Do NOT allow it to cross this particular range, as discontinuity exists as -1.
                const double MINIMAL_POSSIBLE_SIGMA = -0.98;
                if (sigma_minima < MINIMAL_POSSIBLE_SIGMA) {
                    sigma_minima = MINIMAL_POSSIBLE_SIGMA;
                }
            }


            // The result elements. 
            Output_Element output_element;
            output_element.parameter = param;
            output_element.result = result;
            output.push_back (output_element);
        }




    private:

        // Numerical configurations and the output. 
        Configuration configuration;
        std::vector<Output_Element> output;
        std::vector<Result_Element> result;
};

#endif