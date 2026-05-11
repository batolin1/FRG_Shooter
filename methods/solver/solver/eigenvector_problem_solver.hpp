#ifndef EIGENVECTOR_PROBLEM_SOLVER_HPP
#define EIGENVECTOR_PROBLEM_SOLVER_HPP

#include <vector>
#include <algorithm>

#include "solver/integrator/utils/parameter_builder.hpp"
#include "solver/solver/utils/helper_utilities.hpp"
#include "structure.hpp"
#include "solver/integrator/integrator.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "solver/integrator/integrator_model/eigenvector_integrator_model.hpp"

/**
    @brief A solver for the eigenvector problem. 
    @see Solver_Concept
*/
class Eigenvector_Problem_Solver {

    public:

        // some definitions.
        using Parameter = Eigenvector_Problem_Solver_Parameter;
        using Result_Element = Eigenvector_Problem_Solver_Result_Element;
        using Trajectory_Element = Potential_Integrator_Result_Element;
        using Output_Element = Eigenvector_Problem_Solver_Output_Element;

        /**
            @brief Instantiation method. 

            @param configuration    The numerical configurations for the solver.
        */
        Eigenvector_Problem_Solver (const Configuration& configuration) {
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

            Potential_Integrator_Parameter potential_integrator_parameter = 
                parameter_builder::make_potential_integrator_parameter (param);
            // O(N) models not implemented for this solver. 
            potential_integrator_parameter.symmetry_factor_N = 1;

            // Here we add our integrator model. 
            Potential_Integrator_Model potential_integrator_model;
            integrator::integrate_model 
                (configuration, potential_integrator_parameter, potential_integrator_model);
            // integrator::integrate (
            //     configuration, potential_integrator_parameter,potential_integrator_model);
            std::vector<Trajectory_Element> trajectory = 
                potential_integrator_model.get_trajectory ();

            // Removes the asymptote.
            helper_utilities::eliminate_nonphysical_asymptote (trajectory);

            std::vector<Result_Element> result;

            // Loops over the proposed eigenvalues. 
            for (int i = 0; i < param.eigenvalue_delta; ++i) {

                double eigenvalue = param.eigenvalue_minima +
                    (param.eigenvalue_maxima - param.eigenvalue_minima) *
                    i / (param.eigenvalue_delta - 1.0);

                // prepares the integrator parameter.
                Eigenvector_Integrator_Parameter eigenvector_integrator_parameter = 
                    parameter_builder::make_eigenvector_integrator_parameter (param);
                eigenvector_integrator_parameter.eigenvalue = eigenvalue;
                // Not implemented for O(N) models. 
                eigenvector_integrator_parameter.symmetry_factor_N = 1;
                // Use anomalous dimension from re-calculation.
                eigenvector_integrator_parameter.anomalous_dimension = 
                    potential_integrator_parameter.anomalous_dimension;

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