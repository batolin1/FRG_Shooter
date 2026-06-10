#ifndef HELPER_UTILITIES_HPP
#define HELPER_UTILITIES_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace helper_utilities {


    // Fixed parameters for defining a spike.
    const double GRADIENT_THRESHOLD = 1000.0;
    const double VALUE_PERCENTILE = 68.0;
    const double LOWER_THRESHOLD = 0.025;
    const double UPPER_THRESHOLD = 0.975;

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
    std::vector<Trajectory_Element> patch_asymptote (
        std::vector<Trajectory_Element>& trajectory, 
        Parameter& param,
        Configuration& configuration) {

        std::vector<Trajectory_Element> patched_trajectory = trajectory;

        const int traj_size = static_cast<int>(patched_trajectory.size());

        // =========================================================
        // STEP 1: Recompute shape using central differences
        // =========================================================
        for (int i = 1; i < traj_size - 1; i++) {
            const double up_next = patched_trajectory[i+1].potential_1prime;
            const double up_prev = patched_trajectory[i-1].potential_1prime;
            const double f_next  = patched_trajectory[i+1].field;
            const double f_prev  = patched_trajectory[i-1].field;
            const double up      = patched_trajectory[i].potential_1prime;
            if (std::abs(up) > 1e-10) {
                const double dy = (up_next - up_prev) / (f_next - f_prev);
                patched_trajectory[i].potential_1prime_shape = dy * patched_trajectory[i].field / up;
            } else {
                patched_trajectory[i].potential_1prime_shape = 0.0;
            }
        }

        // =========================================================
        // STEP 2: Find zero crossing where U' BECOMES POSITIVE
        // (mirrors Julia: findfirst(u_prime_vals .> 0))
        // =========================================================
        int zero_cross_idx = traj_size - 1;
        for (int i = 0; i < traj_size; i++) {
            if (patched_trajectory[i].potential_1prime > 0) {
                zero_cross_idx = i;
                break;
            }
        }

        // =========================================================
        // STEP 3: Search window — 3 points after zero crossing
        // =========================================================
        const int search_start = std::min(zero_cross_idx + 3, traj_size - 10);

        // The search end is the point where, having the maxima already been reached, the shape coeffs. 
        // start going to negative values again. 
        // Guess search end region.
        int search_end   = (traj_size > zero_cross_idx * 3.0) ?
            zero_cross_idx * 3.0 : traj_size * 0.9;
        // We already crossed the zero-point here, so we are safely out of divergences.
        for (int i = search_start; i < traj_size; i++) {
            if (patched_trajectory[i].potential_1prime_shape < 0) {
                search_end = i;
                break;
            }
        }

        const double expected_shape = (param.s_factor - param.anomalous_dimension) / 
            (param.dimension - (param.s_factor - param.anomalous_dimension));

    // =========================================================
    // STEP 4: Find best patch point — stable -- This IS the only criteria. 
    // =========================================================
    int patch_idx = search_start;
    double best_score = std::numeric_limits<double>::max();

    /**
    Instead of finding BEST stability score, can find first index where stability is 
    shape prev <= shape curr <= shape next. 
    */
    for (int i = search_start + 1; i < search_end - 1; i++) {

        // The shapes
        const double shape_prev = patched_trajectory[i-1].potential_1prime_shape;
        const double shape_curr = patched_trajectory[i  ].potential_1prime_shape;
        const double shape_next = patched_trajectory[i+1].potential_1prime_shape;

        // I think this portion doesnt work too well. 
        // The errors
        const double err_prev = std::abs(shape_prev - expected_shape);
        const double err_curr = std::abs(shape_curr - expected_shape);
        const double err_next = std::abs(shape_next - expected_shape);

        // Local minima of error reached? stop ! 
        if (err_prev <= err_curr && err_curr <= err_next) {
            patch_idx = i - 1;
            break;
        }

    }

    // Safety warning
    const double match_quality = std::abs(
        patched_trajectory[patch_idx].potential_1prime_shape - expected_shape);
    std::cout << "Patching at index: " << patch_idx
            << " field: "            << patched_trajectory[patch_idx].field
            << " expected_shape: "   << expected_shape
            << " actual_shape: "     << patched_trajectory[patch_idx].potential_1prime_shape
            << " delta: "            << match_quality
            << std::endl;
    if (match_quality > 0.1) {
        std::cout << "WARNING: poor asymptote match (delta = " 
                << match_quality << ")" << std::endl;
    }

        const double actual_shape = patched_trajectory[patch_idx].potential_1prime_shape;
        // =========================================================
        // STEP 5: Throw away tail, append clean power-law grid
        // (mirrors Julia: vcat(rho_vals[1:patch_idx], rho_asymp))
        // =========================================================
        const double rho_patch    = patched_trajectory[patch_idx].field;
        const double uprime_patch = patched_trajectory[patch_idx].potential_1prime;
        const double C            = uprime_patch / std::pow(rho_patch, expected_shape);

        // Build clean power-law tail with fixed step
        const double rho_ext_limit = patched_trajectory [zero_cross_idx].field * 3.0;  // Let's say at most 3x the zero-crossing field.
        const double step          = patched_trajectory[1].field - patched_trajectory[0].field; // Same step as original trajectory.  

        std::vector<Trajectory_Element> patched_tail;
        for (double rho = rho_patch; rho <= rho_ext_limit; rho += step) {
            Trajectory_Element e;
            e.field             = rho;
            e.potential_1prime  = C * std::pow(rho, expected_shape);
            // Truncate patched tail when potential_1prime becomes too big. 
            if (std::abs(e.potential_1prime) > configuration.patching_threshold || 
                std::abs (e.potential_2prime) > configuration.patching_threshold) {
                break;
            }
            /**
            Here I think this is calculating the derivatives different from Vally's !!
            */
            e.potential_2prime  = expected_shape * C * std::pow(rho, expected_shape - 1.0);
            e.potential_0prime  = patched_trajectory[patch_idx].potential_0prime +
                C * rho_patch / (expected_shape + 1.0) *
                (std::pow(rho / rho_patch, expected_shape + 1.0) - 1.0);
            patched_tail.push_back(e);
        }

        // Truncate original patched_trajectory at patch point and append clean tail
        patched_trajectory.resize(patch_idx);
        patched_trajectory.insert(patched_trajectory.end(), patched_tail.begin(), patched_tail.end());

        return patched_trajectory;
    }

    /**
        Helper method to populate the values for sigma. 
    */
    template<typename Parameter>
    static std::vector<double> populate_range (const Parameter& param) {
        std::vector<double> sigma_range;
        for (int i = 0; i < param.sigma_steps; ++i) {
            double sigma = param.sigma_minima +
                (param.sigma_maxima - param.sigma_minima) * i / (param.sigma_steps - 1.0);
            sigma_range.push_back (sigma);
        }
        return sigma_range;
    }


    /**
        Helper method to find the spikes, that is, the potential points 
        where a phase transition occurs. 
        @param x_values    The x_values for the dataset. 
        @param y_values    The y_values for the dataset.
        @return            The lits of (potential) spikes (as x-values).
    */
    std::vector<int> find_spike (
        const std::vector<double>& x_values, const std::vector<double>& y_values) {
        
            const int n = x_values.size ();
            if (n != y_values.size () || n < 2) {
                throw std::invalid_argument ("Invalid input sizes");
            }
        
            // 1) Simultaneous-sort x and y.
            std::vector<int> indices (n);
            for (int i = 0; i < n; ++i) indices [i] = i;
        
            std::sort(indices.begin (), indices.end (),
                [&] (int a, int b) { return x_values [a] < x_values [b]; });
        
            std::vector<double> x (n), y (n);
            for (int i = 0; i < n; ++i) {
                x [i] = x_values [indices [i]];
                y [i] = y_values [indices [i]];
            }
        
            // 2) Find the gradients
            std::vector<double> gradient(n);
        
            gradient[0] = (y[1] - y[0]) / (x[1] - x[0]);

            for (int i = 1; i < n - 1; ++i) {

                const double dx1 = x[i]   - x[i - 1];
                const double dx2 = x[i+1] - x[i];

                const double dy1 = (y[i]   - y[i - 1]) / dx1;
                const double dy2 = (y[i+1] - y[i])     / dx2;

                gradient[i] =
                    (dx2 * dy1 + dx1 * dy2) / (dx1 + dx2);
            }

            gradient[n-1] =
                (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
        
            // 3) Compute percentiles
            std::vector<double> grad_copy = gradient;
            std::sort(grad_copy.begin (), grad_copy.end ());
        
            const double pos =
                (VALUE_PERCENTILE / 100.0) * (n - 1);

            const int lower =
                static_cast<int>(std::floor(pos));

            const int upper =
                static_cast<int>(std::ceil(pos));

            const double weight = pos - lower;

            const double value_threshold =
                grad_copy[lower] * (1.0 - weight) +
                grad_copy[upper] * weight;
        
            // 4) Find (potential) spikes.
            std::vector<int> spikes;
            bool in_spike = false;
        
            for (int i = 0; i < n - 1; ++i) {
                bool condition =
                    (std::abs(gradient [i]) > GRADIENT_THRESHOLD) &&
                    (y [i] > value_threshold);
        
                if (condition && !in_spike &&
                    std::abs(x [i]) < UPPER_THRESHOLD && std::abs(x [i]) > LOWER_THRESHOLD) 
                {
                    spikes.push_back (i);
                    in_spike = true;
                }
                else if (!condition) {
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

        double previous_eigenvector = output_element.result[0].asymptotic_eigenvector;
        
        for (const Eigenvector_Solver_Result_Element& r : output_element.result) {
            bool is_sign_flip = r.asymptotic_eigenvector * previous_eigenvector < 0;
            if (is_sign_flip) {
                zero_points.push_back (r.eigenvalue);
            }
            previous_eigenvector = r.asymptotic_eigenvector;
        }
        output_element.zero_points = zero_points;
    }
}

#endif