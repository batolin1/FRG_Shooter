#ifndef HELPER_UTILITIES_HPP
#define HELPER_UTILITIES_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include "logger.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"

namespace helper_utilities {

    /**
        Method to eliminate non-physical asymptotic contributions 
        to the potential. 
        @param field               The field trajectory.
        @param potential_0prime    The potential trajectory.
        @param potential_0prime    The potential derivative trajectory.
        @param potential_0prime    The potential 2-derivative trajectory.
        @param denominator         The denominator.
    */
    template <typename Configuration, typename Trajectory_Element, typename Parameter>
    std::vector<Trajectory_Element> patch_asymptote (std::vector<Trajectory_Element>& trajectory,
        Parameter& param, Configuration& configuration) {

        std::vector<Trajectory_Element> patched_trajectory = trajectory;

        const int traj_size = static_cast<int> (patched_trajectory.size ());

        // Compute shapes using central differences. 
        for (int i = 1; i < traj_size - 1; i++) {
            const double up_next = patched_trajectory [i+1].potential_1prime;
            const double up_prev = patched_trajectory [i-1].potential_1prime;
            const double f_next  = patched_trajectory [i+1].field;
            const double f_prev  = patched_trajectory [i-1].field;
            const double up      = patched_trajectory [i].potential_1prime;
            if (std::abs (up) > 1e-10) {
                const double dy = (up_next - up_prev) / (f_next - f_prev);
                patched_trajectory [i].potential_1prime_shape =
                    dy * patched_trajectory [i].field / up;
            } else {
                patched_trajectory [i].potential_1prime_shape = 0.0;
            }
        }

        // Find zero crossings. 
        const std::vector<parameter_builder::CrossingPoint>& crossing_points = 
            parameter_builder::find_crossing_points (patched_trajectory);

        // Trajectory empty? - logs a warning
        if (crossing_points.empty ()) {
            std::ostringstream oss;
            oss << "Warning: trajectory with s_factor = " << param.s_factor << " anomalous "
                << "dimension= " << param.anomalous_dimension << " dimension= " << param.dimension 
                << "could not find minimas or maximas for this trajectory.";
            Logger::instance (). log (oss.str ());
        }

        // Trajectory has only one element and is a maxima? Logs a warning. 
        if (crossing_points.size () == 1 && crossing_points.back ().potential_2prime < 0) {
            std::ostringstream oss;
            oss << "Warning: trajectory with s_factor = " << param.s_factor << " anomalous "
                << "dimension= " << param.anomalous_dimension << " dimension= " << param.dimension 
                << "could not find a miima.";
            Logger::instance (). log (oss.str ());
        }

        // Always pick the last minima to start patching from. 
        // Loop over all points. Is minima and it the index further located? -> accept.
        int zero_cross_idx = -1;
        for (const parameter_builder::CrossingPoint& cp : crossing_points) {
            if (cp.potential_2prime > 0 && cp.index > zero_cross_idx) {
                zero_cross_idx = cp.index;
            }
        }
        // Could not find one? -> Then it is at the origin! Also throws warning. 
        if (zero_cross_idx == -1) {
            zero_cross_idx = 0;
            std::cout << "WARNING: No fixed point found." << std::endl; 
        }

        // Start search slightly beyond zero crossing.
        const int search_start = zero_cross_idx + 3;
        // End search: is there still a maxima after the found minima? Use it. 
        // Otherwise, just extend trajectory from minima. 
        int search_end = 
            (crossing_points.back ().index > zero_cross_idx && 
            crossing_points.back ().potential_2prime < 0.0) ? 
            crossing_points.back ().index : zero_cross_idx * 3.0;

        // Block if shape drops below zero. 
        for (int i = search_start; i < traj_size; i++) {
            if (patched_trajectory[i].potential_1prime_shape < 0.0) {
                search_end = i;
                break;
            }
        }

        // Find best patch point with respect to the expected shape. 

        const double expected_shape =
            (param.s_factor - param.anomalous_dimension) /
            (param.dimension - (param.s_factor - param.anomalous_dimension));

        int patch_idx = search_start;
        for (int i = search_start + 1; i < search_end - 1; i++) {
            const double shape_prev =
                patched_trajectory [i-1].potential_1prime_shape;
            const double shape_curr =
                patched_trajectory [i].potential_1prime_shape;
            const double shape_next =
                patched_trajectory [i+1].potential_1prime_shape;

            const double err_prev = std::abs (shape_prev - expected_shape);
            const double err_curr = std::abs (shape_curr - expected_shape);
            const double err_next = std::abs (shape_next - expected_shape);

            if (err_prev < err_curr && err_curr < err_next) {
                patch_idx = i - 1;
                break;
            }
        }

        std::ostringstream oss;

        // Safety warning - Details the quality of the patching at asymptote. 
        const double match_quality = std::abs(
            patched_trajectory [patch_idx].potential_1prime_shape - expected_shape);
        oss << "Patching at index: " << patch_idx
            << " field: " << patched_trajectory[patch_idx].field
            << " expected_shape: " << expected_shape
            << " actual_shape: " << patched_trajectory[patch_idx].potential_1prime_shape
            << " delta: " << match_quality << "\n";
        if (match_quality > 0.1) {
            oss << "WARNING: poor asymptote match (delta = "
                << match_quality << ")" << "\n";
        }
        // Logs.
        Logger::instance ().log (oss.str ());

        // Now we actually build the power-law tail for U', U''. Intentioanlly leave U = 0.0 since 
        // This is rebuilt at the trajectory relaxation step anyway. 
        const double rho_patch = patched_trajectory [patch_idx].field;
        const double uprime_patch = patched_trajectory [patch_idx].potential_1prime;
        const double C = uprime_patch / std::pow (rho_patch, expected_shape);

        const double rho_ext_limit = patched_trajectory[zero_cross_idx].field * 3.0;
        const double step = patched_trajectory[1].field - patched_trajectory [0].field;

        std::vector<Trajectory_Element> patched_tail;
        for (double rho = rho_patch; rho <= rho_ext_limit; rho += step) {
            Trajectory_Element e;
            e.field = rho;
            e.potential_1prime = C * std::pow(rho, expected_shape);
            e.potential_2prime = expected_shape * C * std::pow(rho, expected_shape - 1.0);
            e.potential_0prime = 0.0;
            // Beyond the threshold? just stop patching. 
            if (std::abs (e.potential_1prime) > configuration.patching_threshold ||
                std::abs (e.potential_2prime) > configuration.patching_threshold) {
                break;
            }
            patched_tail.push_back (e);
        }

        // Join the tail to the original trajectory. 
        patched_trajectory.resize (patch_idx);
        patched_trajectory.insert (
            patched_trajectory.end (), patched_tail.begin (), patched_tail.end ());

        const int total = static_cast<int> (patched_trajectory.size ());

        // Anchor: trust the ODE value only at the very start, before any drift has accumulated.
        patched_trajectory [0].potential_0prime = trajectory [0].potential_0prime;

        for (int i = 1; i < total; i++) {
            const double up_prev = patched_trajectory [i-1].potential_1prime;
            const double up_curr = patched_trajectory [i].potential_1prime;
            const double dx = patched_trajectory [i].field - patched_trajectory [i-1].field;
            patched_trajectory [i].potential_0prime =
                patched_trajectory [i-1].potential_0prime + 0.5 * dx * (up_prev + up_curr);
        }
        return patched_trajectory;
    }

    /**
        Helper method to find the spikes, that is, the points where potentially a transition occurs. 
        @param x_values    The x_values for the dataset. 
        @param y_values    The y_values for the dataset.
        @return            The lits of (potential) spikes (as x-values).
    */
    template<typename Result_Element> 
    std::vector<int> find_spike (
        const Configuration& config, const std::vector<Result_Element>& result) {

        const int n = result.size ();
        if (n < 2) {
            throw std::invalid_argument ("Invalid input sizes");
        }

        // Simultaneous-sort by sigma (x), keeping track of original indices.
        std::vector<int> indices (n);
        for (int i = 0; i < n; ++i) indices [i] = i;
        std::sort(indices.begin (), indices.end (),
            [&](int a, int b) { return result [a].sigma < result [b].sigma; });
        std::vector<double> x(n), y(n);
        for (int i = 0; i < n; ++i) {
            x [i] = result [indices [i]].sigma;
            y [i] = result [indices [i]].asymptotic_field;
        }

        // Finds gradients.
        std::vector<double> gradient(n);
        gradient [0] = (y [1] - y [0]) / (x [1] - x [0]);
        for (int i = 1; i < n - 1; ++i) {
            const double h0 = x [i] - x [i-1];
            const double h1 = x [i+1] - x [i];
            gradient [i] = 
                (h0 * h0 * y [i+1] + (h1*h1 - h0*h0) * y [i] - h1 * h1 * y [i-1]) / 
                (h0 * h1 * (h0 + h1));
        }
        gradient [n-1] = (y [n-1] - y [n-2]) / (x [n-1] - x [n-2]);

        // Percentile threshold (linear interpolation).
        std::vector<double> grad_copy = gradient;
        std::sort(grad_copy.begin (), grad_copy.end ());
        const double pos = (config.value_percentile / 100.0) * (n - 1);
        const int lower = static_cast<int> (std::floor (pos));
        const int upper = static_cast<int> (std::ceil (pos));
        const double weight = pos - lower;
        const double value_threshold =
            grad_copy [lower] * (1.0 - weight) + grad_copy [upper] * weight;

        // Find spikes — full range, with sign-based index shift.

        std::vector<int> spikes;
        bool in_spike = false;
        for (int i = 0; i < n; ++i) {
            bool condition =
                (std::abs (gradient [i]) > config.gradient_threshold); //&&
                //(std::abs(gradient[i]) > value_threshold) &&
                //(y[i] > value_threshold);

            if (condition && !in_spike &&
                std::abs (x [i]) < config.upper_threshold &&
                std::abs (x [i]) > config.lower_threshold) {

                int spike_pos = i; // position in sorted x/y
                if (x [i] < 0) {
                    if (i+1 < n) {
                        spike_pos = i+1;
                    } else {
                        throw std::out_of_range (
                            "Spike detected at last point with negative x; "
                            "no i+1 to shift to (matches Python's IndexError case)");
                    }
                }
                spikes.push_back (indices [spike_pos]);  // map back to original index
                in_spike = true;

            } else if (!condition) {
                in_spike = false;
            }
        }
        return spikes;
    }

    /**
        @brief Helper method to find zeroes to an eigenvector problem output element.
        @see Eigenvector_Solver

        @param output_element    The output element in question
    */
    void find_zeroes (Eigenvector_Solver_Output_Element& output_element) {

        std::vector<double> zero_points; 

        double previous_eigenvector = output_element.result [0].asymptotic_eigenvector;
        
        for (const Eigenvector_Solver_Result_Element& r : output_element.result) {
            bool is_sign_flip = r.asymptotic_eigenvector * previous_eigenvector < 0.0;
            if (is_sign_flip) {
                zero_points.push_back (r.eigenvalue);
            }
            previous_eigenvector = r.asymptotic_eigenvector;
        }
        output_element.zero_points = zero_points;
    }
}

#endif