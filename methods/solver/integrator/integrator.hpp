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
            
            // if ((int) (step * 1000) % max_steps == 0) {
            //     std::cout << step << " out of " << max_steps << std::endl;
            //     std::cout << model.asymptotic_field << " of max " << model.get_threshold () << std::endl;
            // }

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
                if (integration_time_step < 1e-8) {
                    // std::cout << "Reached few too tiny steps!" << std::endl;
                    // std::cout << "Field" << model.asymptotic_field << std::endl;
                    // std::cout << "Potential" << model.state[0] << std::endl;
                    // std::cout << "Primed Potential" << model.state[1] << std::endl;
                    return; 
                }
            }
        }
    }
};

#endif