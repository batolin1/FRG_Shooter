#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <boost/numeric/odeint.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>

#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"
#include "concept.hpp"

using namespace boost::numeric::odeint;

namespace integrator {

    /**
        @brief An integrator method to realize integration of integrator models. 

        @see Integrator_Model_Concept, Solver_Concept, Configuration. 

        @tparam Parameter           A template for parameters for the integrator model. 
        @tparam Integrator_Model    A template for integrator models. 
        @param  configuration       The (numerical) configurations for the integration.
        @param  parameter           The (physical) parameters for the integration. 
        @param  model               The integratormodel in question. 
    */
    template <typename Parameter, typename Integrator_Model>
    requires Integrator_Model_Concept<Integrator_Model>
    void integrate (
        const Configuration& configuration, 
        Parameter& parameter, 
        Integrator_Model& model) {

        const double PI = 3.14159265358979323846;

        // Creates instance of the solver.
        typedef runge_kutta_fehlberg78<std::array<double,2>> stepper_type;
        auto stepper = make_controlled (
            configuration.tolerance, 
            configuration.absolute_tolerance, 
            stepper_type ());

        // Fixes integration step and starting parameters.
        double integration_time_step = configuration.integration_time_default;

        // Calculates auxiliary parameters and initializes integrator model. 
        parameter_builder::calculate_additional_parameters (parameter);
        model.initialize (configuration, parameter);

        // The ODE.
        const auto ODE = [&model] (
            const std::array<double, 2>& x, std::array<double, 2>& dxdt, double t) {
                
            model.ODE_step (x, dxdt, t);
        };
        
        const int max_steps = 100000;
        int step = 0;

        // Loops over until the terminal value is reached or until 
        // a termination event is triggered. 
        while (model.asymptotic_field < model.get_threshold () && step < max_steps) {

            step++;

            if (model.termination_event ()) break;

            if (std::abs (model.asymptotic_field) > configuration.field_threshold) {
                std::cout << "Field entering integrator: " << model.asymptotic_field << std::endl;
            }

            // Adaptative step must be controlled near boundary. 
            const double remaining = model.get_threshold () - model.asymptotic_field;
            integration_time_step = std::min (integration_time_step, 0.9 * remaining);

            auto result = stepper.try_step (
                ODE, model.state, model.asymptotic_field, integration_time_step);


            if (result == success) {
                model.on_success_step ();
            }
            // Reduces timestep on failure. 
            else {
                integration_time_step *= 0.5;
                // Too small timestep? - Stop. 
                if (integration_time_step < 1e-8) {
                    return; 
                }
            }
        }
    }

    /**
        @brief An integrator model to realize integrations, particularly when the anomalous 
               dimension is meant to be re-calculated adaptitatively. Note that this method 
               integrates only potentials, but not eigenvectors. 
        @param config                The configurations for this integrator. 
        @param param                 The parameters for this integration. 
        @param model                 The model to be integrated. 
        @param maximum_iterations    The maximum number of times the anomalous dimension will be 
                                     adaptitatively re-calculated. 
        @param tolerance             The tolerance below which the anomalous dimension is treated 
                                     as if it had converged. 

    */
    void integrate_model (
        const Configuration& config, 
        Potential_Integrator_Parameter& param, 
        Potential_Integrator_Model& model, 
        const int maximum_iterations, 
        const double tolerance) {

        const double damping = config.damping;

        double old_anomalous_dimension = param.anomalous_dimension;
        int counter = 0;
        while (counter < maximum_iterations) {
            counter++;
            // Realize the integration
            integrate (config, param, model);
            // (re)Calculate anomalous dimension
            const double anomalous_dimension = 
                parameter_builder::recalculate_anomalous_dimension (param, model);

            param.anomalous_dimension =  
                damping * anomalous_dimension + (1.0 - damping) * old_anomalous_dimension;

            // If newly-calculated parameter good already, no need to keep calculating...
            const bool threshold_reached = 
                std::abs (old_anomalous_dimension - param.anomalous_dimension) < tolerance;  
            // If estimate already fairly bad, stop running. 
            // const bool bad_estimate = model.asymptotic_field < 1.0;
            if (threshold_reached) {
                break;
            }
            old_anomalous_dimension = param.anomalous_dimension;
        }
    }

};

#endif