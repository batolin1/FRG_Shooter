#ifndef INITIAL_CONDITION_SOLVER_HPP
#define INITIAL_CONDITION_SOLVER_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "Shooting_Solver.hpp"


/**
    The class executes the shooting method executes as follows:
    1) The input data is read from "input_identifier.txt" and the results are 
       spilled to "output_identifier.txt".
    2) Input: Only a single row is provided as input and a single file is 
       provided as output. The input is the following:
        dimension
            The dimension of the model being analysed; example, a 
            two-dimensional or three-dimensional model. Fractional 
            dimensions are also accepted.
        anomalous_dimension *  
            If the anomalous dimension of the model is already known in 
            advance, user can include it. In general recommended to set to 
            zero. 
        s_factor_minima             
            For long-range ising models, corresponds to "s" or "sigma" in 
            literature. Set to 2.0 to solve for a short-range ising model.
            This is the minimum value for the run. 
        s_factor_maxima
            For long-range ising models, corresponds to "s" or "sigma" in 
            literature. Set to 2.0 to solve for a short-range ising model.
            This is the maximum value for the run. 
        s_factor_delta
            For long-range ising models, corresponds to "s" or "sigma" in 
            literature. Set to 2.0 to solve for a short-range ising model.
            This is the interval for the run.
        nth_critical_point
            The critical point that this algorithm is analysing (example, 
            bi-critical, tri-critical, etc. n=1 is the gaussian fixed point).
        sigma_minima           
            Sigma is the index for the families of initial-conditions in the
            solution-space. This parameter is the minimum sigma in the range
            of families that the algorithm will try to shoot for.
        sigma_maxima          
            The maximum sigma in the range of families that the algorithm will
            try to shoot for.
        number_of_steps        
            The number of steps between the minima and maxima of sigma the 
            algorithm will try for.
        symmetry_factor_N**  
            For O(N) models, this is N. Set to 1.0 for the long-range model.
    * Note that this particular implementation does not solve for the anomalous 
    dimension self-consistently. Therefore, it is generally not recommended for 
    models with anomalous dimension; instead, this parameter should generally 
    be set to 0.0 and include an anomalous dimension only for running 
    validations when the anomalous dimension for the particular model known. 
    **The algorithm is unstable for large N, as this leads to suppression of 
    the potential's second derivative. This is most suitable for ranges such 
    as 1 <= N <= 20. On the other hand, one observes that for large N, the 
    model quickly converges to a mean-field universality class; that one
    should not concern him/herself too deeply with this.
    For each row in the input text file, the algorithm executes one run of the 
    shooting algorithm, within the range from sigma minima and sigma maxima. 
    The results are stored to the output file. The first four columns in the  
    output file are simply labels (dimension, anomalous_dimension, s_factor, 
    and symmetry_factor_N) that allows one to map a particular row to a 
    particular solution, whereas the last two columns correpsond to the values
    of sigma and the respective wavefunction the algorithm integrated to before
    diverging. By plotting the asymptotic wavefunction as function of sigma, 
    peaks on the graph give the values of sigma, that is, the initial 
    condition, for which the wavefunction does not diverge, and therefore those
    represent real physical solutions to the RG-flow at fixed point. Each of 
    the peaks in the graph therefore represent (potentially) a critical point.
    Given the input, the algorithm will run the shooting method over the range
    [s_factor_minima, s_factor_maxima, s_factor_delta] and identify the initial
    condition sigma for which the nth_critical_point for this particular value 
    of the s_factor is a solution to the RG equation. 
    The algorithm will output a text file in the style of 
    "input_eigenperturbation.txt", that can be directly inputted for the 
    identification of RG eigenvalues.
*/
class Initial_Condition_Solver : public Shooting_Solver {


    public:

        /**
            Simple static method to execute the shooting. 
        */
        int execute (
            const std::string& input_filename,
            const std::string& output_filename,
            const std::string& configuration_filename) override {

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
                    read_parameters_from_line (line, 10);
                // If for a particular row any of the token provided could not be
                // converted, ignore row and throws a warning to the user. 
                } catch (...) {
                    std::cerr << "Conversion error in line: " << line << std::endl;
                    continue;
                }

                // Actually runs the shooting algorithm and store results to  
                // file, as comma-separated values; by ranging over the chosen
                // values of sigma, and evaluating the asymptotic wavefunction 
                // for each of them. 
                for (int s = 0; s < s_factor_delta; s++) {
                    const double s_factor = s_factor_minima + 
                        (s_factor_maxima - s_factor_minima) * 
                        s / (s_factor_delta - 1.0);
                    const std::vector<double> asymptotic_wavefunction = 
                        find_asymptotic_wavefunction (s_factor);
                    const std::vector<double> sigma_range = populate_range
                        (sigma_minima, sigma_maxima, sigma_delta);
                    const std::vector<double> initial_condition = 
                        find_spike (sigma_range, asymptotic_wavefunction);
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

        static std::vector<double> populate_range (
            const double sigma_minima,
            const double sigma_maxima, 
            const int sigma_delta) {

            std::vector<double> sigma_range;

            for (int i = 0; i < sigma_delta; ++i) {
                double sigma = sigma_minima +
                    (sigma_maxima - sigma_minima) * 
                    i / (sigma_delta - 1.0);
                sigma_range.push_back (sigma);
            }
                
            return sigma_range;
        }

        static std::vector<double> find_spike (
            const std::vector<double>& x_values,
            const std::vector<double>& y_values) {

            
                const double gradient_threshold = 1000.0;
                const double value_percentile = 50.0;
                const double lower_threshold = 0.025;
                const double upper_threshold = 0.68;
            
                const int n = x_values.size();
                if (n != y_values.size() || n < 2) {
                    throw std::invalid_argument("Invalid input sizes");
                }
            
                // ---- 1. Sort x and y together ----
                std::vector<int> indices(n);
                for (int i = 0; i < n; ++i) indices[i] = i;
            
                std::sort(indices.begin(), indices.end(),
                    [&](int a, int b) { return x_values[a] < x_values[b]; });
            
                std::vector<double> x(n), y(n);
                for (int i = 0; i < n; ++i) {
                    x[i] = x_values[indices[i]];
                    y[i] = y_values[indices[i]];
                }
            
                // ---- 2. Compute gradient (finite differences) ----
                std::vector<double> gradient(n);
            
                for (int i = 1; i < n - 1; ++i) {
                    gradient[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]);
                }
                // endpoints (simple forward/backward)
                gradient[0] = (y[1] - y[0]) / (x[1] - x[0]);
                gradient[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
            
                // ---- 3. Compute percentile (median here = 50th) ----
                std::vector<double> grad_copy = gradient;
                std::sort(grad_copy.begin(), grad_copy.end());
            
                int idx = static_cast<int>(value_percentile / 100.0 * (n - 1));
                double value_threshold = grad_copy[idx];
            
                // ---- 4. Detect spikes ----
                std::vector<double> spikes;
                bool in_spike = false;
            
                for (int i = 0; i < n - 1; ++i) {
                    bool condition =
                        (std::abs(gradient[i]) > gradient_threshold) &&
                        (y[i] > value_threshold);
            
                    if (condition && !in_spike &&
                        std::abs(x[i]) < upper_threshold /* && std::abs(x[i]) > lower_threshold */) 
                    {
                        spikes.push_back(x[i+1]);  // store y-value of spike
                        in_spike = true;
                    }
                    else if (!condition) {
                        in_spike = false;
                    }
                }
                return spikes;
        }

        void set_parameters (const std::vector<std::string>& token) override {
            dimension = std::stod (token[0]);
            anomalous_dimension = std::stod (token[1]);
            s_factor_minima = std::stod (token[2]);
            s_factor_maxima = std::stod (token[3]);
            s_factor_delta= std::stoi (token[4]);
            nth_critical_point = std::stoi (token[5]);
            sigma_minima = std::stod (token[6]);
            sigma_maxima = std::stod (token[7]);
            sigma_delta = std::stoi (token[8]);
            symmetry_factor_N = std::stod (token[9]);
        }

};

#endif