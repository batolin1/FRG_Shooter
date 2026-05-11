#ifndef PARAMETER_BUILDER_HPP
#define PARAMETER_BUILDER_HPP

#include <cmath>

#include "structure.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"
#include "solver/solver/utils/helper_utilities.hpp"

namespace parameter_builder {

    /**
        @brief Helper method to calculate additional, implied parameters, for integrator parameters. 
        
        @see Integrator, Integrator_Parameter_Concept

        @param parameter    The parameter in question. 
    */
    template <typename Parameter>
    inline void calculate_additional_parameters (Parameter& parameter) {

        constexpr double PI = 3.14159265358979323846;

        parameter.dimension_factor = 
            4.0 / parameter.dimension / std::pow (2.0, parameter.dimension + 1.0) /
            std::pow(PI, parameter.dimension * 0.5) / std::tgamma(0.5 * parameter.dimension);
        parameter.anomalous_constant =
            1.0 - parameter.anomalous_dimension / (2.0 + parameter.dimension);
        parameter.implied_s_factor = parameter.s_factor - parameter.anomalous_dimension;
        parameter.s_constant = parameter.s_factor / 2.0 * parameter.anomalous_constant;
    }

    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Shooting_Problem_Solver_parameter object.

        @see Shooting_Problem_Solver, Potential_Integrator_Model

        @param param    The shooting problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Shooting_Problem_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;
        integrator_parameter.symmetry_factor_N = param.symmetry_factor_N;

        return integrator_parameter;
    }

    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Initial_Condition_Problem_Solver_parameter object.

        @see Initial_Condition_Problem_Solver, Potential_Integrator_Model

        @param param    The initial condition problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Initial_Condition_Problem_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.symmetry_factor_N = param.symmetry_factor_N;

        return integrator_parameter;
    }

    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Eigenvector_Problem_Solver_parameter object.

        @see Eigenvector_Problem_Solver, Eigenvector_Integrator_Model

        @param param    The eigenvector problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Eigenvector_Problem_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;
        integrator_parameter.sigma = param.sigma;

        return integrator_parameter;
    }

    /**
        @brief Helper method to make a Eigenvector_Integrator_Parameter object, given a 
               Eigenvector_Problem_Solver_parameter object.

        @see Eigenvector_Problem_Solver, Eigenvector_Integrator_Model

        @param param    The eigenvector problem solver parameter. 
        @return         An eigenvector integrator parameter. 
    */
    inline Eigenvector_Integrator_Parameter make_eigenvector_integrator_parameter (
        const Eigenvector_Problem_Solver_Parameter& param) {

        Eigenvector_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;
        integrator_parameter.sigma = param.sigma;

        return integrator_parameter;
    }

    inline double recalculate_anomalous_dimension (
        const Potential_Integrator_Parameter& parameter, Potential_Integrator_Model& model) {

        constexpr double PI = 3.14159265358979323846;
    
        // Get the trajectories and eliminate the "fake" parts. 
        std::vector<Potential_Integrator_Result_Element> trajectory = model.get_trajectory ();

        // Unless this only works for actual, physical trajectories?
        helper_utilities::eliminate_nonphysical_asymptote (trajectory);

        // Finds the potential minima, and the respective rho and U'' (rho) at the minima. 
        // (IS THIS THE SAME RULE FOR TRI-CRITICAL AND MULTI-CRITICAL MODELS??)
        double potential_minima = trajectory.at (0).potential_0prime;
        double field_at_minima = trajectory.at (0). field;
        double potential_2prime =  trajectory.at (0).potential_2prime;
        for (Potential_Integrator_Result_Element& result_element : trajectory) {
            if (result_element.potential_0prime < potential_minima) {
                potential_minima = result_element.potential_0prime;
                field_at_minima = result_element.field;
                potential_2prime = result_element.potential_2prime;
            }
        }

        // Actually calculates anomalous dimension. 
        const double anomalous_dimension = 
            16.0 / parameter.dimension / std::pow (2.0, parameter.dimension + 1.0) /
            std::pow(PI, parameter.dimension * 0.5) / std::tgamma(0.5 * parameter.dimension) * 
            field_at_minima * potential_2prime * potential_2prime / 
            std::pow ((1.0 + 2.0 * field_at_minima * potential_2prime), 2.0) * 
            (1 - parameter.anomalous_dimension / (parameter.s_factor + parameter.dimension));

        return anomalous_dimension;
    }

        /**
            @brief Method to get the asymptotic continuation of the potential and its derivatives, 
                   given some field input x, with asymptotic formula. 
                
            @param x             The field. 
            @param derivative    zeroth, first, or second derivative.
            @return              The potential at the nth derivative, given by asymptote. 
        */
        inline double get_asymtotic_continuation (const double x, const int derivative,
            const Eigenvector_Integrator_Parameter& parameter) {

            const double factor = 
                parameter.dimension / (parameter.dimension - parameter.implied_s_factor);

            // The last element available. 
            Potential_Integrator_Result_Element last_element = 
                parameter.trajectory.at (parameter.trajectory.size () - 2);

            const double ratio = x / last_element.field;

            if (derivative == 0) {
                return last_element.potential_0prime * std::pow (ratio, factor);
            } else if (derivative == 1) {
                return last_element.potential_1prime * std::pow (ratio, factor - 1.0);
            } else if (derivative == 2) {
                return last_element.potential_2prime * std::pow (ratio, factor - 2.0);
            }
            return 0;
        }
}

#endif