#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <boost/numeric/odeint.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>

using namespace boost::numeric::odeint;
using integrable_element = std::array<double, 2>;
using trajectory = std::vector<double>;

class Integrator {

public:

    virtual ~Integrator() = default;

    double compute_asymptotic_value () {

        typedef runge_kutta_cash_karp54<integrable_element> stepper_type;
        auto stepper = make_controlled(1e-6, 1e-9, stepper_type());

        double integration_time_step = INTEGRATION_TIME_DEFAULT;

        update_derived_parameters();
        set_initial_condition();
        asymptotic_wavefunction = WAVEFUNCTION_PERTURBATION;

        while (asymptotic_wavefunction < WAVEFUNCTION_THRESHOLD) {

            if (termination_event()) {
                return result();
            }

            const auto ODE = [this](
                const integrable_element &x,
                integrable_element &dxdt,
                double t) {
                ODE_step (x, dxdt, t);
            };

            const auto step_result =
                stepper.try_step(ODE, state, asymptotic_wavefunction, integration_time_step);

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
            @param wavefunction_perturbation    The definition of zero in the 
                                                context of the wavefunction's 
                                                initial condition.
            @param wavefunction_threshold       The definition of infinity in 
                                                the context of the wavefunction
                                                's asymptotic value.
        */
        void set_configuration (
            const double practically_zero,
            const double practically_infinity,
            const double integration_time_default,
            const double wavefunction_perturbation, 
            const double wavefunction_threshold) {
                PRACTICALLY_ZERO = practically_zero;
                PRACTICALLY_INFINITY = practically_infinity;
                INTEGRATION_TIME_DEFAULT = integration_time_default;
                WAVEFUNCTION_PERTURBATION = wavefunction_perturbation;
                WAVEFUNCTION_THRESHOLD = wavefunction_threshold;
            }

protected:

    // ===== shared state =====
    integrable_element state;
    double asymptotic_wavefunction;

    // ===== configuration =====
    double PRACTICALLY_ZERO = 1e-6;
    double PRACTICALLY_INFINITY = 1e+6;
    double INTEGRATION_TIME_DEFAULT = 1e-6;
    double WAVEFUNCTION_PERTURBATION = 1e-8;
    double WAVEFUNCTION_THRESHOLD = 100;

    // ===== shared physics params =====
    double dimension;
    double anomalous_dimension;
    double s_factor;
    double sigma;
    double symmetry_factor_N;

    // ===== derived =====
    double dimension_factor;
    double anomalous_constant;
    double implied_s_factor;
    double s_constant;
    trajectory wavefunction;
    trajectory potential_0prime;
    trajectory potential_1prime;
    trajectory potential_2prime;

    double PI = 3.14159265358979323846;

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

    // ===== hooks (MUST be implemented) =====
    virtual void set_initial_condition() = 0;
    virtual void ODE_step (const integrable_element&, integrable_element&, double) = 0;
    virtual bool termination_event() = 0;
    virtual double result() = 0;

    // ===== optional hook =====
    virtual void on_success_step() {}
};

#endif