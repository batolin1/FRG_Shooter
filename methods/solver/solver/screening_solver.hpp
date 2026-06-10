#ifndef SCREENING_SOLVER_HPP
#define SCREENING_SOLVER_HPP

#include <vector>

#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"

#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "solver/solver/grid_search_solver.hpp"
#include "logger.hpp"

#include <random>


/**
    @brief A solver for the shooting problem. 
    @see Solver_Concept
*/
class Screening_Solver {


    public:

        // Some definitions
        using Parameter = Screening_Solver_Parameter;
        using Result_Element = Screening_Solver_Result_Element;
        using Output_Element = Screening_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Screening_Solver (const Configuration& configuration) {
                this->configuration = configuration;
                std::random_device rd;
                random_engine.seed (rd ());
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

        void execute_for_single_step (const Parameter& param, double& temperature,
            double& current_sigma, double& current_anomalous_dimension, double& asymptotic_field) {
  

                double new_sigma = current_sigma + 
                    (2.0 * uniform_01_distribution(random_engine) - 1.0) * temperature * 
                    (param.sigma_maxima - param.sigma_minima);
                double new_eta   = param.search_anomalous_dimension ? 
                    current_anomalous_dimension + 
                    (2.0 * uniform_01_distribution(random_engine) - 1.0) * temperature * 
                    (param.anomalous_dimension_maxima - param.anomalous_dimension_minima) : 
                    param.anomalous_dimension_minima;

                // Clamps.
                new_sigma = std::clamp (new_sigma, param.sigma_minima, param.sigma_maxima);
                new_eta = std::clamp (
                    new_eta, param.anomalous_dimension_minima, param.anomalous_dimension_maxima);

                Potential_Integrator_Parameter integrator_parameter = 
                    parameter_builder::make_potential_integrator_parameter (param);

                integrator_parameter.sigma = new_sigma;
                integrator_parameter.anomalous_dimension = new_eta;
                Potential_Integrator_Model  integrator_model;
                // No need to save trajectory in this case.
                integrator_model.save_trajectory = false; 
                integrator::integrate (configuration, integrator_parameter, integrator_model);

                // Actually get asymptotic field at the point where cutoff occurs.
                const double new_asymptotic_field = integrator_model.get_result ();
                const double delta = 
                    (std::abs(new_asymptotic_field) - std::abs(asymptotic_field)) / 
                    std::abs (asymptotic_field);

                // Accept improvement always; accept worse with probability exp(delta/T)
                if (delta > 0 || std::exp(delta / temperature) > 
                    uniform_01_distribution(random_engine)) {
                        current_sigma = new_sigma;
                        current_anomalous_dimension = new_eta;
                        asymptotic_field = new_asymptotic_field;
                        // If accepted, also store results. We won't store bad results as they 
                        // don't really help for anything. 
                        sigma_values.push_back (current_sigma);
                        anomalous_dimension_values.push_back (current_anomalous_dimension);
                        asymptotic_field_values.push_back (asymptotic_field);
                }
                // Cools-down the temperature. 
                temperature *= param.cooling_rate;
        }

        std::vector<double> execute_for_single_iteration (const Parameter& param, 
            const int iteration, const double number_of_steps,
            const double start_temperature, const double cooling_rate, 
            double& current_sigma, double& current_anomalous_dimension) {

            // Clears the vectors for storing sampled results. 
            sigma_values.clear ();
            anomalous_dimension_values.clear ();
            asymptotic_field_values.clear ();
        
            double temperature = start_temperature;
            
            // Executes the very first run separatedly
            Potential_Integrator_Parameter integrator_parameter = 
                parameter_builder::make_potential_integrator_parameter (param);
            integrator_parameter.sigma = current_sigma;
            integrator_parameter.anomalous_dimension = current_anomalous_dimension;
            Potential_Integrator_Model integrator_model;
            integrator_model.save_trajectory = false; // No need to save trajectory in this case.
            integrator::integrate (configuration, integrator_parameter, integrator_model);
            double asymptotic_field = integrator_model.get_result ();

            // stores data from first run. 
            sigma_values.push_back (current_sigma);
            anomalous_dimension_values.push_back (current_anomalous_dimension);
            asymptotic_field_values.push_back (asymptotic_field);

            for (int step = 1; step < number_of_steps; step++) {
                execute_for_single_step (param, temperature, current_sigma, 
                    current_anomalous_dimension, asymptotic_field); 
            }

            // Actually creates the result elements and stores. 
            for (int i = 0; i < sigma_values.size (); i++) {
                Result_Element result_element;
                result_element.asymptotic_field = asymptotic_field_values[i];
                result_element.sigma = sigma_values[i];
                result_element.anomalous_dimension = anomalous_dimension_values[i];
                    result.push_back (result_element);   
            }

            // The best estimates. 
            std::vector<double> best_estimates; 
            best_estimates.push_back (current_sigma);
            best_estimates.push_back (current_anomalous_dimension);
            best_estimates.push_back (asymptotic_field);
            return best_estimates;
        }

        void write_log (const Parameter& param) {
            std::ostringstream oss;
            const std::string run_subprocess = param.run_subprocess ? "true" : "false";
            const std::string run_window_search = param.run_window_search ? "true" : "false";
            const std::string search_anomalous_dimension = 
                param.search_anomalous_dimension ? "true" : "false";

            oss << "Initiating screening solver with the following parameters:" << "\n"
                << "dimension: " <<  param.dimension << "\n"
                << "s_factor: " <<  param.s_factor << "\n"
                << "sigma_minima: " <<  param.sigma_minima << "\n"
                << "sigma_maxima: " <<  param.sigma_maxima << "\n"
                << "anomalous_dimension_minima: " <<  param.anomalous_dimension_minima << "\n"
                << "anomalous_dimension_maxima: " <<  param.anomalous_dimension_maxima << "\n"
                << "temperature: " <<  param.temperature << "\n"
                << "cooling_rate: " <<  param.cooling_rate << "\n"
                << "number_of_iterations: " <<  param.number_of_iterations << "\n"
                << "number_of_steps: " <<  param.number_of_steps << "\n"
                << "temperature_subprocess: " <<  param.temperature_subprocess << "\n"
                << "cooling_rate_subprocess: " <<  param.cooling_rate_subprocess << "\n"
                << "number_of_steps_subprocess: " <<  param.number_of_steps_subprocess << "\n"
                << "number_of_steps_window_search: " <<  param.number_of_steps_window_search << "\n"
                << "run_subprocess: " <<  run_subprocess << "\n"
                << "run_window_search: " <<  run_window_search << "\n" 
                << "search_anomalous_dimension: " << search_anomalous_dimension << "\n";
            Logger::instance ().log (oss.str ());
        }

        /**
            @brief A method to process an input parameter, that is, given an input to this solver,
                   this method processes the input and inserts results as output elements. 
            
            @param param    The (physical) model parameters for this solver, for which to process.
        */
        void process_parameter (const Parameter& param) {

            write_log (param);

            // Create results' container. 
            result.clear ();

            // Initiates the distributions. 
            sigma_distribution = std::uniform_real_distribution<double> 
                (param.sigma_minima, param.sigma_maxima);
            if (param.search_anomalous_dimension) {
                anomalous_dimension_distribution = std::uniform_real_distribution<double> 
                    (param.anomalous_dimension_minima, param.anomalous_dimension_maxima); 
            }
            
            Logger::instance ().log ("Initiating simulated annealing (main process).");
            // Repeat simulated annealing over a number of iterations.
            for (int iteration = 0; iteration < param.number_of_iterations; iteration++) {
                double current_sigma = sigma_distribution (random_engine);
                double current_anomalous_dimension = param.search_anomalous_dimension ? 
                    anomalous_dimension_distribution (random_engine) : 
                    param.anomalous_dimension_minima;
                const double temperature = param.temperature;
                const double cooling_rate = param.cooling_rate;
                const double number_of_steps = param.number_of_steps;
                const std::vector<double> best_estimates = execute_for_single_iteration (
                    param, iteration, number_of_steps, temperature, cooling_rate, 
                    current_sigma, current_anomalous_dimension);
                best_sigma_values.push_back (best_estimates [0]);
                best_anomalous_dimension_values.push_back (best_estimates [1]);
                best_asymptotic_field_values.push_back (best_estimates [2]);
            }

            // Only if required, execute a second annealign to improve on the estimates.  
            if (param.run_subprocess) {
                Logger::instance ().log ("Initiating simulated annealing (subprocess).");
                // Now, we actually want to improve the best estimates themselves. First we run a 
                // secondary annealing at a lower temperature. 
                for (int iteration = 0; iteration < best_sigma_values.size (); iteration++) {
                    double current_sigma = best_sigma_values [iteration];
                    double current_anomalous_dimension = param.search_anomalous_dimension ? 
                        best_anomalous_dimension_values [iteration] : 
                        param.anomalous_dimension_minima;
                    // We have two cooling rates and temperatures, one for the search process,
                    // and one for the improvement process. 
                    const double temperature = param.temperature_subprocess;
                    const double cooling_rate = param.cooling_rate_subprocess;
                    const double number_of_steps = param.number_of_steps_subprocess;
                    const std::vector<double> best_estimates = execute_for_single_iteration 
                        (param, iteration, number_of_steps, temperature, cooling_rate, 
                        current_sigma, current_anomalous_dimension);
                    best_sigma_values [iteration] = best_estimates [0];
                    best_anomalous_dimension_values [iteration] = best_estimates [1];
                    best_asymptotic_field_values [iteration] = best_estimates [2];   
                }
            }

            if (param.run_window_search) {
                Logger::instance ().log ("Initiating window search.");
                // We run the window search using the grid search solver. 
                Grid_Search_Solver grid_search_solver (configuration); 
                for (int i = 0; i < best_sigma_values.size (); i++) {
                    const double sigma = best_sigma_values [i];
                    const double anomalous_dimension = best_anomalous_dimension_values [i];
                    Grid_Search_Solver_Parameter grid_search_parameter = 
                        parameter_builder::build_grid_search_solver_parameter 
                            (param, sigma, anomalous_dimension, configuration);
                    grid_search_solver.process_parameter (grid_search_parameter);
                }

                // For each output from window search, store elements to result elements. 
                for (const Grid_Search_Solver_Output_Element& e : grid_search_solver.get_result ()) {
                    double best_sigma = 0;
                    double best_anomalous_dimension = 0;
                    double best_asymptotic_field = 0;
                    for (const Grid_Search_Solver_Result_Element& r : e.result) {
                        Result_Element result_element;
                        result_element.asymptotic_field = r.asymptotic_field;
                        result_element.sigma = r.sigma;
                        result_element.anomalous_dimension = r.anomalous_dimension;
                        result.push_back (result_element);
                        if (r.asymptotic_field > best_asymptotic_field) {
                            best_sigma = r.sigma;
                            best_anomalous_dimension = r.anomalous_dimension;
                            best_asymptotic_field = r.asymptotic_field; 
                        }
                    }
                    best_sigma_values.push_back (best_sigma);
                    best_anomalous_dimension_values.push_back (best_anomalous_dimension);
                    best_asymptotic_field_values.push_back (best_asymptotic_field);
                }
            
            }

            // Stores the best results. 
            for (int i = 0; i < best_sigma_values.size (); i++) {
                Result_Element result_element;
                result_element.asymptotic_field = best_asymptotic_field_values[i];
                result_element.sigma = best_sigma_values[i];
                result_element.anomalous_dimension = best_anomalous_dimension_values[i];
                best_result.push_back (result_element);   
            }

        
            Output_Element output_element;
            output_element.parameter = param;
            output_element.result = result;
            output_element.best_result = best_result;
            output.push_back (output_element);
        }

    private:

        // Numerical configurations and the output. 
        Configuration configuration;
        std::vector<Output_Element> output;
        std::vector<Result_Element> result;
        std::vector<Result_Element> best_result;
        std::mt19937 random_engine;
        std::uniform_real_distribution<double> anomalous_dimension_distribution;
        std::uniform_real_distribution<double> sigma_distribution;
        std::uniform_real_distribution<double> uniform_01_distribution {0.0, 1.0};
        std::vector<double> sigma_values;
        std::vector<double> anomalous_dimension_values;
        std::vector<double> asymptotic_field_values;
        std::vector<double> best_sigma_values; 
        std::vector<double> best_anomalous_dimension_values; 
        std::vector<double> best_asymptotic_field_values; 
};

#endif