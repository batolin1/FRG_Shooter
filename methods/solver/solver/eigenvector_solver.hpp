#ifndef EIGENVECTOR_SOLVER_HPP
#define EIGENVECTOR_SOLVER_HPP

#include <vector>
#include <algorithm>

#include "solver/integrator/utils/parameter_builder.hpp"
#include "solver/solver/utils/helper_utilities.hpp"
#include "solver/solver/utils/grid_relaxer.hpp"
#include "structure.hpp"
#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "solver/integrator/integrator_model/eigenvector_integrator_model.hpp"
#include "solver/solver/grid_search_solver.hpp"
#include "logger.hpp"

/**
    @brief A solver for the eigenvector problem. 
    @see Solver_Concept
*/
class Eigenvector_Solver {

    public:

        // some definitions.
        using Parameter = Eigenvector_Solver_Parameter;
        using Result_Element = Eigenvector_Solver_Result_Element;
        using Trajectory_Element = Potential_Integrator_Result_Element;
        using Output_Element = Eigenvector_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Eigenvector_Solver (const Configuration& configuration) {
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
            oss << "Initiating eigenvector solver with the following parameters:" << "\n"
                << "dimension: " << param.dimension << "\n"
                << "anomalous_dimension: " << param.anomalous_dimension << "\n"
                << "s_factor: " << param.s_factor << "\n"
                << "sigma: " << param.sigma << "\n"
                << "eigenvalue_maxima: " << param.eigenvalue_maxima << "\n"
                << "eigenvalue_minima: " << param.eigenvalue_minima << "\n"
                << "eigenvalue_steps: " << param.eigenvalue_steps << "\n"; 
            Logger::instance ().log (oss.str ());
        }

        void write_grid_search_log (const Parameter& param, const double asymptotic_field, 
            const double sigma, const double anomalous_dimension) {
                std::ostringstream oss;
                oss << "Following grid search calculation in eigenvector solver, obtained the " 
                    << "following parameter improvement for comparison: \n" 
                    << "sigma input: " << param.sigma << "\n"
                    << "sigma output: " << sigma << "\n"
                    << "anomalous dimension input: " << param.anomalous_dimension << "\n"
                    << "anomalous dimension output: " << anomalous_dimension << "\n"
                    << "asymptotic field: " << asymptotic_field << "\n"
                    << "Proceeding with calculation.\n";
                Logger::instance (). log (oss.str ());
            }

        /**
            @brief A method to process an input parameter, that is, given an input to this solver,
                   this method processes the input and inserts results as output elements. 
            
            @param param    The (physical) model parameters for this solver, for which to process.
        */
        void process_parameter (const Parameter& param) {

            write_log (param);

            // First we want to improve the results (if possible) or, at the very least, make sure
            // we are starting at very good initial conditions. Therefore run grid search solver in 
            // the vicinity of proposed initial condition. 
            Grid_Search_Solver grid_search_solver (configuration);
            Grid_Search_Solver_Parameter grid_search_parameter = 
                parameter_builder::make_grid_search_solver_parameter (param, configuration);
            grid_search_solver.process_parameter (grid_search_parameter);
            double asymptotic_field = 0;
            double sigma = param.sigma;
            double anomalous_dimension = param.anomalous_dimension;
            for (const Grid_Search_Solver_Result_Element& r : 
                grid_search_solver.get_result (). back ().result) {
                    if (r.asymptotic_field > asymptotic_field) {
                        asymptotic_field = r.asymptotic_field;
                        sigma = r.sigma;
                        anomalous_dimension = r.anomalous_dimension;
                    }
            }

            write_grid_search_log (param, asymptotic_field, sigma, anomalous_dimension);
            
            // Integrates the potential again to actually obtain the trajectory. 
            Potential_Integrator_Parameter potential_integrator_parameter = 
                parameter_builder::make_potential_integrator_parameter (param);
            potential_integrator_parameter.sigma = sigma;
            potential_integrator_parameter.anomalous_dimension = anomalous_dimension;
            Potential_Integrator_Model potential_integrator_model;
            potential_integrator_model.save_trajectory = true;
            integrator::integrate 
                (configuration, potential_integrator_parameter, potential_integrator_model);
            std::vector<Trajectory_Element> trajectory = 
                potential_integrator_model.get_trajectory ();
                
            // After the trajectory reaches asymptotic behaviour, we want to actually cut it off. 
            // We want very granular results *before* asymptote and afterwards we just care about 
            // the average behaviour. 
            int cutoff_index = trajectory.size ();
            for (int i = 0; i < trajectory.size (); i++) {
                if (trajectory[i].potential_1prime > configuration.patching_threshold || 
                    trajectory[i].potential_2prime > configuration.patching_threshold) {
                    cutoff_index = i; 
                    break;
                }
            }

            // Next relax it and patch. 
            trajectory = grid_relaxer::relax_trajectory
                (trajectory, configuration.relaxation_grid_size);
            trajectory = helper_utilities::patch_asymptote 
                (trajectory, param, configuration);

            // Actually reached the eigenvalue integrator. 
            std::vector<Result_Element> result;

            // Loops over the proposed eigenvalues. 
            for (int i = 0; i < param.eigenvalue_steps; ++i) {

                double eigenvalue = param.eigenvalue_minima +
                    (param.eigenvalue_maxima - param.eigenvalue_minima) *
                    i / (param.eigenvalue_steps - 1.0);

                // prepares the integrator parameter.
                Eigenvector_Integrator_Parameter eigenvector_integrator_parameter = 
                    parameter_builder::make_eigenvector_integrator_parameter (param);
                eigenvector_integrator_parameter.eigenvalue = eigenvalue;
                eigenvector_integrator_parameter.sigma = sigma;
                eigenvector_integrator_parameter.anomalous_dimension = anomalous_dimension;
                    
                Eigenvector_Integrator_Model eigenvector_integrator_model;
                eigenvector_integrator_parameter.trajectory = trajectory;

                integrator::integrate (
                    configuration, eigenvector_integrator_parameter, eigenvector_integrator_model);
                const double asymptotic_eigenvector = eigenvector_integrator_model.get_result ();
                
                Result_Element result_element;
                result_element.asymptotic_eigenvector = asymptotic_eigenvector;
                result_element.eigenvalue = eigenvalue;
                result.push_back (result_element);
            }

            // Inserts as output element
            Output_Element output_element;
            output_element.parameter = param;
            output_element.trajectory = trajectory;
            output_element.result = result;
            // Insert (re-calculated) anomalous dimension.
            output_element.anomalous_dimension = 
                potential_integrator_parameter.anomalous_dimension;
            // Inserts the crossings-at-zero values. 
            helper_utilities::find_zeroes (output_element);
            output.push_back (output_element);
        }

    private:

        // The configurations and the output elements. 
        Configuration configuration;
        std::vector<Output_Element> output;
};

#endif