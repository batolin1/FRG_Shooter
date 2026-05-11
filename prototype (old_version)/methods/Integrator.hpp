#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <boost/numeric/odeint.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>

// Labels. 
using namespace boost::numeric::odeint;
using integrable_element = std::array<double, 2>;
using trajectory = std::vector<double>;

/**
    An abstract class for realizing ODE integrations. 
*/
class Integrator {

public:

    virtual ~Integrator() = default;

    /**
        Computes the asymptotic value of the ODE state.
    */
    double compute_asymptotic_value () {

        // Creates instance of the solver.
        typedef runge_kutta_cash_karp54<integrable_element> stepper_type;
        auto stepper = make_controlled(1e-6, 1e-9, stepper_type());

        // Fixes integration step and starting parameters.
        double integration_time_step = INTEGRATION_TIME_DEFAULT;
        asymptotic_field = FIELD_PERTURBATION;
        update_derived_parameters();
        set_initial_condition();

        // Integrates up to the threshold
        while (asymptotic_field < FIELD_THRESHOLD) {
            // Terminates execution if a termination event is triggered
            // and returns variable of interest. 
            if (termination_event()) {
                return result();
            }

            // Realizes one "infinitesimal" step. 

            const auto ODE = [this](
                const integrable_element &x,
                integrable_element &dxdt,
                double t) {
                ODE_step (x, dxdt, t);
            };

            const auto step_result =
                stepper.try_step (
                    ODE, 
                    state, 
                    asymptotic_field, 
                    integration_time_step);

            // Realizes updates depending on whether step was successful. 
            if (step_result == success) {
                on_success_step();
            } else {
                integration_time_step *= 0.5;
            }
        }
        return result();
    }

        /**
            Sets a particular configuration concerning numerical stability, 
            that is, those are "spurious" parameters that only exist due to
            computational limitations. 
            @param practically_zero             The definition of zero. 
            @param practically_infinity         The definition of infinity. 
            @param integration_time_default     The definition of an 
                                                infinitesimal timestep.
            @param field_perturbation           The definition of zero in the 
                                                context of the field's 
                                                initial condition.
            @param field_threshold              The definition of infinity in 
                                                the context of the field
                                                's asymptotic value.
        */
        void set_configuration (
            const double practically_zero,
            const double practically_infinity,
            const double integration_time_default,
            const double field_perturbation, 
            const double field_threshold) {
                PRACTICALLY_ZERO = practically_zero;
                PRACTICALLY_INFINITY = practically_infinity;
                INTEGRATION_TIME_DEFAULT = integration_time_default;
                FIELD_PERTURBATION = field_perturbation;
                FIELD_THRESHOLD = field_threshold;
            }

protected:

    // the integrable element for the ODE solver, and the variable to keep 
    // track of the asymptotic field value. 
    integrable_element state;
    double asymptotic_field;

    // "Spurious" configuration parameters due to numerical limitations, 
    // in order to ensure numerical stability. 
    double PRACTICALLY_ZERO = 1e-6;
    double PRACTICALLY_INFINITY = 1e+6;
    double INTEGRATION_TIME_DEFAULT = 1e-6;
    double FIELD_PERTURBATION = 1e-8;
    double FIELD_THRESHOLD = 100;
    double PI = 3.14159265358979323846;

    // Physical parameters.
    double dimension;
    double anomalous_dimension;
    double s_factor;
    double sigma;
    double symmetry_factor_N;

    // Calculated parameters 
    double dimension_factor;
    double anomalous_constant;
    double implied_s_factor;
    double s_constant;

    // Output trajectories
    trajectory field;
    trajectory potential_0prime;
    trajectory potential_1prime;
    trajectory potential_2prime;

    /**
        Given input, this method calculates additional parameters to complete
        the integration step. 
    */
    void update_derived_parameters() {
        dimension_factor = 4.0 / dimension /
            std::pow(2.0, dimension + 1.0) /
            std::pow(PI, dimension * 0.5) /
            std::tgamma(0.5 * dimension);
        anomalous_constant =
            1.0 - anomalous_dimension / (2.0 + dimension);
        implied_s_factor = s_factor - anomalous_dimension;
        s_constant = s_factor / 2.0 * anomalous_constant;
    }

    /**
        Virtual method to set the initial condition given input. 
    */
    virtual void set_initial_condition() = 0;

    /**
        Virtual method, the actual differential equation to be integrated. 
        @param state               The state to be integrated.
        @param state_derivative    The derivative of the state integrated.
        @param field               The field value at the particular step. 
    */
    virtual void ODE_step (
        const integrable_element&, 
        integrable_element&, 
        double) = 0;

    /**
        Virtual method to set the definition of a termination event. 
    */
    virtual bool termination_event() = 0;

    /**
        Virtual method to set the definition of the result.
    */
    virtual double result() = 0;

    /**
        Virtual method to set the action after a success. 
    */
    virtual void on_success_step() {}
};

#endif