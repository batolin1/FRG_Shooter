#ifndef PARAMETER_BUILDER_HPP
#define PARAMETER_BUILDER_HPP

#include <cmath>
#include <stdexcept>
#include <limits>
#include "structure.hpp"
#include "solver/integrator/integrator_model/potential_integrator_model.hpp"

namespace parameter_builder {

    /**
        @brief Helper method to calculate additional, implied parameters, for integrator parameters. 
        
        @see Integrator, Integrator_Parameter_Concept

        @param parameter    The parameter in question. 
    */
    template <typename Parameter>
    inline void calculate_additional_parameters (Parameter& parameter) {

        constexpr double PI = 3.14159265358979323846;

        parameter.dimension_factor = 4.0 / parameter.dimension / 
            std::pow((4.0 * PI), parameter.dimension * 0.5) / 
            std::tgamma(0.5 * parameter.dimension);
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

            grid_search_parameter.search_anomalous_dimension = param.search_anomalous_dimension;

            const double sigma_delta = std::abs (sigma) * config.grid_extension;
            const double anomalous_dimension_delta = 
                std::abs (anomalous_dimension) * config.grid_extension;

            // Make sure window size is odd, this way the annealing estimations lie inside grid.
            const int window_size = config.window_size % 2 == 0 ? 
                config.window_size + 1.0 : config.window_size;

            grid_search_parameter.dimension = param.dimension;
            grid_search_parameter.s_factor = param.s_factor;
            grid_search_parameter.sigma_minima = 
                // minima smaller than -0.95 ? truncate. Otherwise keep. 
                (sigma - sigma_delta) < -0.95 ? -0.95 : sigma - sigma_delta;
            // No bounds on the upper side. 
            grid_search_parameter.sigma_maxima = sigma + sigma_delta;
            grid_search_parameter.anomalous_dimension_minima = 
                // Anomalous dimension less than zero? truncate at zero. 
                (anomalous_dimension - anomalous_dimension_delta)  < 0.0 ? 0.0 : 
                anomalous_dimension - anomalous_dimension_delta;
            grid_search_parameter.anomalous_dimension_maxima = 
                // Maxima greater than one? truncate. 
                (anomalous_dimension + anomalous_dimension_delta) > 1.0 ? 1.0 : 
                anomalous_dimension + anomalous_dimension_delta;
            grid_search_parameter.sigma_steps = window_size; 
            grid_search_parameter.anomalous_dimension_steps = window_size; 
            grid_search_parameter.window_range = config.window_range;

            // Are we NOT searching anomalous dimension? In that case do NOT apply the deltas. 
            grid_search_parameter.anomalous_dimension_minima = 
                !param.search_anomalous_dimension ? anomalous_dimension : 
                grid_search_parameter.anomalous_dimension_minima;

            // This one is not available !! I think in config?
            grid_search_parameter.number_of_iterations = config.number_recalculations;

        return grid_search_parameter;
    }


    inline Grid_Search_Solver_Parameter build_grid_search_solver_parameter (
        const Screening_Solver_Parameter& param, const double sigma, 
        const double anomalous_dimension, const Configuration& config) {

            Grid_Search_Solver_Parameter grid_search_parameter; 

            const double sigma_delta = std::abs (sigma) * config.grid_extension;
            const double anomalous_dimension_delta = 
                std::abs (anomalous_dimension) * config.grid_extension;

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

    /**
        @brief Helper structure to allocate points where the potential U'(rho) = 0. 
    */
    struct CrossingPoint {
        double field;
        double potential_0prime;
        double potential_1prime;
        double potential_2prime;
        double potential_3prime;
        bool is_minimum;
        // Where it occurs 
        int index;
    };

    /**
        @brief Helper methods to find points where U' (rho) = 0.

        @param parameter    The potential integrator model data used to calculate trajectory.
        @param model        The integration model from which the trajectory is extracted.
        @return             The points where zero crossings occur.   
    */
    inline std::vector<CrossingPoint> find_crossing_points (
        const std::vector<Potential_Integrator_Result_Element>& trajectory) {

        std::vector<CrossingPoint> crossing_points;

        if (trajectory.empty ()) {
            return crossing_points;
        }

        // Consider the origin as a candidate (symmetric phase / unbroken minimum).
        
        CrossingPoint cp;
        cp.field            = trajectory.front().field;
        cp.potential_0prime = trajectory.front().potential_0prime;
        cp.potential_1prime = trajectory.front().potential_1prime;
        cp.potential_2prime = trajectory.front().potential_2prime;
        cp.potential_3prime = trajectory.front().potential_3prime;
        cp.is_minimum       = (trajectory.front().potential_2prime > 0.0);
        cp.index = -1;
        crossing_points.push_back(cp);

        // All the other candidates
        for (int i = 1; i < trajectory.size(); i++) {
            const bool sign_flip = 
                trajectory[i].potential_1prime * trajectory[i-1].potential_1prime < 0;
            if (!sign_flip) continue;

            // Linear interpolation to find more precise zero location.
            const double f0 = trajectory[i-1].field;
            const double f1 = trajectory[i  ].field;
            const double u0 = trajectory[i-1].potential_1prime;
            const double u1 = trajectory[i  ].potential_1prime;
            const double t  = -u0 / (u1 - u0); // fraction in [0,1]
            const double field_zero = f0 + t * (f1 - f0);

            // Interpolate all quantities at the zero.
            const double p0 = trajectory[i-1].potential_0prime + 
                t * (trajectory[i].potential_0prime - trajectory[i-1].potential_0prime);
            const double p2 = trajectory[i-1].potential_2prime + 
                t * (trajectory[i].potential_2prime - trajectory[i-1].potential_2prime);
            const double p3 = trajectory[i-1].potential_3prime + 
                t * (trajectory[i].potential_3prime - trajectory[i-1].potential_3prime);

            CrossingPoint cp;
            cp.field            = field_zero;
            cp.potential_0prime = p0;
            cp.potential_1prime = 0.0;
            cp.potential_2prime = p2;
            cp.potential_3prime = p3;
            cp.is_minimum       = (p2 > 0.0);
            cp.index = i;
            crossing_points.push_back(cp);
        }

        return crossing_points;
    }
     
    /**
        @brief Helper method to re-calculate the anomalous dimension self-consistently.

        @param parameter    The potential integration parameters.
        @param model        The integration model which we are integrating. 
        @return             The new self-consistently calculated anomalous dimension. 
    */
    inline double recalculate_anomalous_dimension (
        const Potential_Integrator_Parameter& parameter, 
        Potential_Integrator_Model& model) {

        const std::vector<Potential_Integrator_Result_Element>& trajectory = model.get_trajectory();

        const std::vector<CrossingPoint>& crossing_points = 
            find_crossing_points (trajectory);

        if (crossing_points.empty()) {
            return parameter.anomalous_dimension;
        }

        // Among all minima (V'' > 0), find the one with lowest potential value V.
        // This is the global minimum — the physically meaningful evaluation point for eta.
        const CrossingPoint* best = nullptr;
        double lowest_potential   = std::numeric_limits<double>::max();

        for (const CrossingPoint& cp : crossing_points) {
            if (!cp.is_minimum) continue;
            if (cp.potential_0prime < lowest_potential) {
                lowest_potential = cp.potential_0prime;
                best = &cp;
            }
        }

        // Fallback: if no proper minimum found, use the point with smallest |V'|.
        if (best == nullptr) {
            double smallest_uprime = std::numeric_limits<double>::max();
            for (const CrossingPoint& cp : crossing_points) {
                if (std::abs(cp.potential_1prime) < smallest_uprime) {
                    smallest_uprime = std::abs(cp.potential_1prime);
                    best = &cp;
                }
            }
        }

        // Still nothing — return unchanged.
        if (best == nullptr) {
            return parameter.anomalous_dimension;
        }

        const double rho0 = best->field;
        const double upp  = best->potential_2prime;
        const double uppp = best->potential_3prime;

        const double num = 3.0 * upp + 2.0 * rho0 * uppp;

        const double mod_num = 3.0 + 6.0 * rho0 * upp + 2.0 * std::pow (rho0 * upp, 2.0);


        // Eq. (2) from Codello / eq. (48) from Defenu et al.:
        // eta = 4 * rho0 * (V'')^2 / (1 + 2*rho0*V'')^2
        const double denom = 1.0 + 2.0 * rho0 * upp;
        if (std::abs(denom) < 1e-10) {
            return parameter.anomalous_dimension;
        }

        // const double anomalous_dimension = 
        //     4.0 * rho0 * upp * upp / std::pow (denom, 2.0);


        // THis one is a LEGITIMATE solution. 
        const double anomalous_dimension = 2.0 * rho0 * upp * upp * mod_num /
            denom / denom / (1.0 + rho0 * upp) / (1.0 + rho0 * upp);



        // const double anomalous_dimension = 2.0 * rho0 * num * num / std::pow (denom, 4.0);
        
        // Simple safe-guard. -> calculated anomalous dimension too big or too small?
        // In this case don't even bother recomputing. This is hard-coded. 
        if (anomalous_dimension > 1.0 || anomalous_dimension < 0.00001) {
            return parameter.anomalous_dimension;
        }

        return anomalous_dimension;
    }
};

#endif