#ifndef INTEGRATOR_EIGENVECTOR_HPP
#define INTEGRATOR_EIGENVECTOR_HPP

#include "Integrator.hpp"

/**
    This class realizes the integration for the ODE corresponding to the 
    linearized eigenvector equation, given a "guess" for the eigenvalue. 
*/
class Integrator_Eigenvector : public Integrator {

    public:

        /**
            Constructs an instance of this class. For more detail regarding the
            input parameters, see comments on "Eigenperturbation_Method.hpp".
            Crucially, sigma is the index parameter for the particular family  
            of initial conditions that integrates the ODE corresponding to the 
            RG flow of the potential at a fixed point. It thus characterizes
            the particular fixed point being investigated. The wavefunction and
            the potential_prime trajectories correspond to the ODE evolution of
            the potential and derivatives at the fixed point of the RG flow. 
            Those trajectories are required as the ODE for the eigenvector 
            depends on those trajectories. 
            @param dimension              The dimension for the system in question. 
            @param anomalous_dimension    The anomalous dimension. 
            @param s_factor               The s-factor for the model.
            @param symmetry_factor_N      The N factor for an O(N) model.
            @param sigma                  The parameter sigma corresponding 
                                          to the solution of the ODE for the 
                                          potential RG flow at a critical pt.
            @param eigenvalue             The "guess" eigenvalue corresponding
                                          to a potential solution of this ODE.
            @param wavefunction           The trajectory for wavefunctions. 
            @param potential_0prime       The trajectory for the potential
                                          along the solution
                                          to the potential RG flow ODE. 
            @param potential_1prime       The trajectory for the potential's 
                                          first derivative along the solution
                                          to the potential RG flow ODE. 
            @param potential_2prime       The trajectory for the potential's 
                                          second derivative along the solution
                                          to the potential RG flow ODE. 
        */
        void initialize (
            // Simulation parameters (from input file).
            const double dimension,
            const double anomalous_dimension,
            const double s_factor,
            const double sigma, 
            const double eigenvalue, 
            const trajectory& wavefunction, 
            const trajectory& potential_0prime, 
            const trajectory& potential_1prime, 
            const trajectory& potential_2prime) {
                this->dimension = dimension;
                this->anomalous_dimension = anomalous_dimension;
                this->s_factor = s_factor;
                this->sigma = sigma;
                this->eigenvalue = eigenvalue;
                this->wavefunction = wavefunction;
                this->potential_0prime = potential_0prime;
                this->potential_1prime = potential_1prime;
                this->potential_2prime = potential_2prime;
        }

    private:

        double eigenvalue;


        /**
            Given a family labelled by an eigenvalue, calculates the initial
            conditions for the eigenvectors corresponding to this family.
            @return               an array containing the first and second 
                                  derivatives of the potential at the initial
                                  condition and for the family sigma. 
        */
        void set_initial_condition () override {
            // The eigenvector at start is a free parameter so we fix to 1.0.
            const double eigenvector = 1.0;
            const double prefactor = 
                (1.0 + potential_0prime[0]) * 
                (1.0 + potential_0prime[0]) /
                 s_constant / dimension_factor;
            const double eigenvector_derivative = 
                -prefactor * (dimension + eigenvalue);
            const double eigenvector_2derivative = 
                0.5 * prefactor * prefactor * 
                (eigenvector + dimension) * 
                (eigenvector * implied_s_factor);
            // Introduces perturbation to the initial condition. 
            const double eigenvector_perturbed = eigenvector + 
                WAVEFUNCTION_PERTURBATION * eigenvector_derivative + 0.5 *
                std::pow (WAVEFUNCTION_PERTURBATION, 2.0) * eigenvector_2derivative;
            const double eigenvector_derivative_perturbed = eigenvector_derivative + 
                WAVEFUNCTION_PERTURBATION * eigenvector_2derivative;

            this->state = {eigenvector_perturbed, eigenvector_derivative_perturbed};
        }

        /**
            The ordinary differential equation for the eigenvectors of the RG.
            @param eigenvector               Array containing the eigenvector 
                                             and its first derivative. 
            @param eigenvector_derivative    Array for the first and second
                                             derivatives of the eigenvector.
            @param wavefunction              The wavefunction at this step. 
        */
        void ODE_step 
            (const integrable_element &state, 
            integrable_element &state_derivative,
            const double field) override {
                const double potential_contribution = 
                    1.0 + get_potential (field, 1) + 
                    2.0 * field * get_potential (field, 2);
                const double prefactor = 
                    std::pow (potential_contribution, 2.0) / 
                    dimension_factor / 
                    s_constant;
                const double state_2derivative = - 0.5 / field * (
                    state[1] + prefactor * (
                        (dimension + eigenvalue) * state[0] -
                        (dimension - implied_s_factor) * field * state[1]
                    )
                );
                state_derivative[0] = state [1];
                state_derivative[1] = state_2derivative;
        }


        bool termination_event() override {
            const bool is_termination_event =
            std::abs (state[1]) > PRACTICALLY_INFINITY ||
            std::abs (state[0]) > PRACTICALLY_INFINITY;
            return is_termination_event; 
        }
    
        double result() override {
            return state[0];
        }

        /** 
            An interpolation method to, given a wavefunction, extract the
            potentials and its derivatives. 
            @param x             The wavefunction at which the interpolation
                                 is being evaluated. 
            @param derivative    Whether to extract zeroth, first or 
                                 second derivative. 
            @return              Returns the interpolated value for the 
                                 potential or its derivatives at x.
        */
        double get_potential (const double x, const int derivative) {
            // Firs find index closest to the solution. 
            int i = 0;
            while (i < wavefunction.size () && x > wavefunction[i]) {
                i++;
            }
            // Handles out of boundary issues.
            if (i >= wavefunction.size ()) i = wavefunction.size () - 1;
            if (i==0) i=1;
            // Chooses whether to take the zeroth, first, or second derivative.
            const trajectory* potential = nullptr;
            if (derivative == 0) {
                potential = &potential_0prime;
            } else if (derivative == 1) {
                potential = &potential_1prime;
            } else if (derivative == 2) {
                potential = &potential_2prime;
            }
            // Interpolates between (i-1) and (i). 
            const double overstep = 
                (x - wavefunction[i-1]) / 
                (wavefunction[i] - wavefunction[i-1]);
            const double y = (*potential)[i-1] + 
                overstep * ((*potential)[i] - (*potential)[i-1]);
            return y;
        }
};

#endif