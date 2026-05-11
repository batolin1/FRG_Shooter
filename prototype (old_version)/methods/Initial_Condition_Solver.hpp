#ifndef INITIAL_CONDITION_SOLVER_HPP
#define INITIAL_CONDITION_SOLVER_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "Shooting_Solver.hpp"

class Initial_Condition_Solver : public Shooting_Solver {

    private:

        // Fixed parameters for defining a spike.
        const double GRADIENT_THRESHOLD = 1000.0;
        const double VALUE_PERCENTILE = 50.0;
        const double LOWER_THRESHOLD = 0.025;
        const double UPPER_THRESHOLD = 0.975;

        // Running parameters for solver, user input. 
        double s_factor_minima;
        double s_factor_maxima;
        int s_factor_delta;

    public:

        /**
            Override implementation to execute eigenperturbation solver.
            @see Solver for more details. 
        */
        int execute (
            const std::string& input_filename,
            const std::string& output_filename,
            const std::string& configuration_filename) override {

            // Sets configuration up.
            CONFIGURATION_FILENAME = configuration_filename;
            read_configuration_from_file 
                (CONFIGURATION_FILENAME, integrator_potential);
            
            // Opens input and opens or creates output file.
            std::ifstream infile (input_filename);
            std::ofstream outfile (output_filename);

            // Throws an error if input file could not be openend. 
            if (!infile.is_open ()) {
                std::cerr << "Error opening input file\n";
                return 1;
            }

            // Reads file line-by-line.
            std::string line;
            while (std::getline (infile, line)) {
                // Tries to assign rows to the respective variables. 
                try {
                    read_parameters_from_line (line, 9);
                // If for a particular row any of the token provided could not be
                // converted, ignore row and throws a warning to the user. 
                } catch (...) {
                    std::cerr << "Conversion error in line: " << line << std::endl;
                    continue;
                }
                // Actually runs the shooting algorithm and store results to  
                // file, as comma-separated values; by ranging over the chosen
                // values of sigma, and evaluating the asymptotic field 
                // for each of them. 
                for (int s = 0; s < s_factor_delta; s++) {
                    const double s_factor = s_factor_minima + 
                        (s_factor_maxima - s_factor_minima) * 
                        s / (s_factor_delta - 1.0);
                    const std::vector<double> asymptotic_field = 
                        find_asymptotic_field (s_factor);
                    const std::vector<double> sigma_range = populate_range ();
                    const std::vector<double> initial_condition = 
                        find_spike (sigma_range, asymptotic_field);
                    for (double spike : initial_condition) {
                        outfile << dimension << ","
                        << anomalous_dimension << ","
                        << s_factor << ","
                        << spike << "\n";
                    }
                }
            }
            return 0;
        }

        /**
            Helper method to populate the values for sigma. 
        */
        std::vector<double> populate_range () {
            std::vector<double> sigma_range;
            for (int i = 0; i < sigma_delta; ++i) {
                double sigma = sigma_minima +
                    (sigma_maxima - sigma_minima) * 
                    i / (sigma_delta - 1.0);
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
        std::vector<double> find_spike (
            const std::vector<double>& x_values,
            const std::vector<double>& y_values) {
            
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
            
                for (int i = 1; i < n - 1; ++i) {
                    gradient [i] = (y [i+1] - y [i-1]) / (x [i+1] - x [i-1]);
                }
                // (endpoints)
                gradient [0] = (y [1] - y [0]) / (x [1] - x [0]);
                gradient [n-1] = (y [n-1] - y [n-2]) / (x [n-1] - x [n-2]);
            
                // 3) Compute percentiles
                std::vector<double> grad_copy = gradient;
                std::sort(grad_copy.begin (), grad_copy.end ());
            
                int idx = static_cast<int>(VALUE_PERCENTILE / 100.0 * (n - 1));
                double value_threshold = grad_copy [idx];
            
                // 4) Find (potential) spikes.
                std::vector<double> spikes;
                bool in_spike = false;
            
                for (int i = 0; i < n - 1; ++i) {
                    bool condition =
                        (std::abs(gradient [i]) > GRADIENT_THRESHOLD) &&
                        (y [i] > value_threshold);
            
                    if (condition && !in_spike &&
                        std::abs(x [i]) < UPPER_THRESHOLD /* && std::abs(x [i]) > LOWER_THRESHOLD */) 
                    {
                        spikes.push_back(x [i+1]);
                        in_spike = true;
                    }
                    else if (!condition) {
                        in_spike = false;
                    }
                }
                return spikes;
        }

        /**
            Override of virtual method in abstract class. 
            @see Solver for more details. 
        */
        void set_parameters (const std::vector<std::string>& token) override {
            dimension = std::stod (token[0]);
            anomalous_dimension = std::stod (token[1]);
            s_factor_minima = std::stod (token[2]);
            s_factor_maxima = std::stod (token[3]);
            s_factor_delta= std::stoi (token[4]);
            sigma_minima = std::stod (token[5]);
            sigma_maxima = std::stod (token[6]);
            sigma_delta = std::stoi (token[7]);
            symmetry_factor_N = std::stod (token[8]);
        }

};

#endif