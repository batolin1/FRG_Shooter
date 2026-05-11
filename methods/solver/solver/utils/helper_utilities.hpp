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
        template <typename Trajectory_Element>
        void eliminate_nonphysical_asymptote (std::vector<Trajectory_Element>& trajectory) {

            int termination_index = trajectory.size ();

            std::vector<double> denominator;

            // Pick out the ("fake") denominator elements. 
            for (const Trajectory_Element& element : trajectory) {
                denominator.push_back (element.denominator);
            }

            // We find the global minima before the global maxima. This is the true minima. 
            auto maxima = std::max_element (
                denominator.begin (), denominator.begin () + denominator.size ());
            int maximum_element_index = std::distance (denominator.begin (), maxima);
            // If the maxima is at zero, don't actually resize anything. 
            if (maximum_element_index == 0) return;
            auto minima = std::min_element (
                denominator.begin (), denominator.begin () + maximum_element_index);
            termination_index = std::distance (denominator.begin (), minima);
            if (termination_index == 0) return;
            // Resize to truncate trajectory to this minima.
            trajectory.resize (termination_index);
        }

        /**
            Helper method to populate the values for sigma. 
        */
        template<typename Parameter>
        static std::vector<double> populate_range (const Parameter& param) {
            std::vector<double> sigma_range;
            for (int i = 0; i < param.sigma_delta; ++i) {
                double sigma = param.sigma_minima +
                    (param.sigma_maxima - param.sigma_minima) * i / (param.sigma_delta - 1.0);
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
            @see Eigenvector_Problem_Solver

            @param output_element    The output element in question
        */
        void find_zeroes (Eigenvector_Problem_Solver_Output_Element& output_element) {

            std::vector<double> zero_points; 

            double previous_eigenvector = output_element.result[0].asymptotic_eigenvector;
            
            for (const Eigenvector_Problem_Solver_Result_Element& r : output_element.result) {
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