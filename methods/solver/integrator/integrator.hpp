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
        typedef runge_kutta_cash_karp54<std::array<double,2>> stepper_type;
        auto stepper = make_controlled (1e-6, 1e-8, stepper_type ());

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
        
        // Loops over until the terminal value is reached or until 
        // a termination event is triggered. 
        while (model.asymptotic_field < configuration.field_threshold) {

            if (model.termination_event ()) break;

            auto result = stepper.try_step (
                ODE, model.state, model.asymptotic_field, integration_time_step);

            if (result == success) {
                model.on_success_step ();
            }
            // Reduces timestep on failure. 
            else {
                integration_time_step *= 0.5;
            }
        }
    }

    void integrate_model (
        const Configuration& config, 
        Potential_Integrator_Parameter& param, 
        Potential_Integrator_Model& model) {

        const double damping = 0.5;

        double old_anomalous_dimension = param.anomalous_dimension;
        int counter = 0;
        while (counter < config.maximum_iterations) {
            // Realize the integration
            integrate (config, param, model);
            // (re)Calculate anomalous dimension
            const double anomalous_dimension = 
                parameter_builder::recalculate_anomalous_dimension (param, model);
                
            param.anomalous_dimension =  
                damping * anomalous_dimension + (1.0 - damping) * old_anomalous_dimension;
            if (std::abs (old_anomalous_dimension - anomalous_dimension) < config.tolerance) {
                old_anomalous_dimension = anomalous_dimension;
                break;
            }
            counter++;
            old_anomalous_dimension = anomalous_dimension;
        }
        if (counter >= config.maximum_iterations) {
            std::cout << "Triggered maximum iterations at sigma = " << param.sigma << std::endl;
        }
    }

};

#endif