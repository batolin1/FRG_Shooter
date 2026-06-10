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
            1.0 - parameter.anomalous_dimension / (parameter.s_factor + parameter.dimension);
        parameter.implied_s_factor = parameter.s_factor - parameter.anomalous_dimension;
        parameter.s_constant = parameter.s_factor / 2.0 * parameter.anomalous_constant;
    }

    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Shooting_Solver_Parameter object.

        @see Shooting_Solver_Parameter, Potential_Integrator_Model

        @param param    The shooting problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Shooting_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;

        return integrator_parameter;
    }

    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Initial_Condition_Solver_parameter object.

        @see Initial_Condition_Solver, Potential_Integrator_Model

        @param param    The initial condition problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Initial_Condition_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;

        return integrator_parameter;
    }

    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Screening_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.s_factor= param.s_factor;

        return integrator_parameter;
    }

    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Grid_Search_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.s_factor= param.s_factor;

        return integrator_parameter;
    }
    
    /**
        @brief Helper method to make a Potential_Integrator_Parameter object, given a 
               Eigenvector_Solver_parameter object.

        @see Eigenvector_Solver, Eigenvector_Integrator_Model

        @param param    The eigenvector problem solver parameter. 
        @return         A potential integrator parameter. 
    */
    inline Potential_Integrator_Parameter make_potential_integrator_parameter (
        const Eigenvector_Solver_Parameter& param) {
        
        Potential_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;
        integrator_parameter.sigma = param.sigma;

        return integrator_parameter;
    }

    /**
        @brief Helper method to make a Eigenvector_Integrator_Parameter object, given a 
               Eigenvector_Solver_parameter object.

        @see Eigenvector_Solver, Eigenvector_Integrator_Model

        @param param    The eigenvector problem solver parameter. 
        @return         An eigenvector integrator parameter. 
    */
    inline Eigenvector_Integrator_Parameter make_eigenvector_integrator_parameter (
        const Eigenvector_Solver_Parameter& param) {

        Eigenvector_Integrator_Parameter integrator_parameter;

        integrator_parameter.dimension = param.dimension;
        integrator_parameter.anomalous_dimension = param.anomalous_dimension;
        integrator_parameter.s_factor = param.s_factor;
        integrator_parameter.sigma = param.sigma;

        return integrator_parameter;
    }

    inline Grid_Search_Solver_Parameter make_grid_search_solver_parameter (
        const Eigenvector_Solver_Parameter& param, const Configuration& config) {

            Grid_Search_Solver_Parameter grid_search_parameter;

            const double sigma = param.sigma;
            const double anomalous_dimension = param.anomalous_dimension;
            const double TEN_PERCENT = 0.1;

            // Only search for anomalous dimension when it is explicitly non-zero. 
            bool search_anomalous_dimension = true;
            if (std::abs (anomalous_dimension) < config.practically_zero) {
                search_anomalous_dimension = false;
            }
            grid_search_parameter.search_anomalous_dimension = search_anomalous_dimension;

            const double sigma_delta = std::abs (sigma) * TEN_PERCENT;
            const double anomalous_dimension_delta = std::abs (anomalous_dimension) * TEN_PERCENT;

            // Make sure window size is odd, this way the annealing estimations lie inside grid.
            const int window_size = config.window_size % 2 == 0 ? 
                config.window_size + 1.0 : config.window_size;

            grid_search_parameter.dimension = param.dimension;
            grid_search_parameter.s_factor = param.s_factor;
            grid_search_parameter.sigma_minima = sigma - sigma_delta;
            grid_search_parameter.sigma_maxima = sigma + sigma_delta;
            grid_search_parameter.anomalous_dimension_minima = 
                anomalous_dimension - anomalous_dimension_delta;
            grid_search_parameter.anomalous_dimension_maxima = 
                anomalous_dimension + anomalous_dimension_delta;
            grid_search_parameter.sigma_steps = window_size; 
            grid_search_parameter.anomalous_dimension_steps = window_size; 
            grid_search_parameter.window_range = config.window_range;

            // This one is not available !! I think in config?
            grid_search_parameter.number_of_iterations = config.number_recalculations;

        return grid_search_parameter;
    }


    inline Grid_Search_Solver_Parameter build_grid_search_solver_parameter (
        const Screening_Solver_Parameter& param, const double sigma, 
        const double anomalous_dimension, const Configuration& config) {

            Grid_Search_Solver_Parameter grid_search_parameter; 

            const double TEN_PERCENT = 0.1;
            const double sigma_delta = std::abs (sigma) * TEN_PERCENT;
            const double anomalous_dimension_delta = std::abs (anomalous_dimension) * TEN_PERCENT;

            // Make sure window size is odd, this way the annealing estimations lie inside grid.
            const int window_size = config.window_size % 2 == 0 ? 
                config.window_size + 1.0 : config.window_size;

            grid_search_parameter.dimension = param.dimension;
            grid_search_parameter.s_factor = param.s_factor;
            grid_search_parameter.sigma_minima = sigma - sigma_delta;
            grid_search_parameter.sigma_maxima = sigma + sigma_delta;
            grid_search_parameter.anomalous_dimension_minima = 
                anomalous_dimension - anomalous_dimension_delta;
            grid_search_parameter.anomalous_dimension_maxima = 
                anomalous_dimension + anomalous_dimension_delta;
            grid_search_parameter.sigma_steps = window_size; 
            grid_search_parameter.anomalous_dimension_steps = window_size; 
            grid_search_parameter.window_range = config.window_range;
            grid_search_parameter.number_of_iterations = param.number_of_steps_window_search;
            grid_search_parameter.search_anomalous_dimension = param.search_anomalous_dimension;

            return grid_search_parameter;
        } 



};

#endif