#ifndef SHOOTING_METHOD_HPP
#define SHOOTING_METHOD_HPP

#include "Solver.hpp"

/**
    The class executes the shooting method as follows:
    1) The input data is read from "shooting_input.txt" and the results are 
       spilled to "shooting_output.txt". The actual numerical analysis of the
       output file is carried out via python file "shooting_plotter.py".
    2) Input: each row will run one execution of the shooting algorithm. The 
       user must provide, as comma-separated values, the following input:
        dimension              
            The dimension of the model being analysed; example, a 
            two-dimensional or three-dimensional model. Fractional 
            dimensions are also accepted.
        anomalous_dimension *  
            If the anomalous dimension of the model is already known in 
            advance, user can include it. In general recommended to set to 
            zero. 
        s_factor               
            For long-range ising models, corresponds to "s" or "sigma" in 
            literature. Set to 2.0 to solve for a short-range ising model.
        symmetry_factor_N**  
            For O(N) models, this is N. Set to 1.0 for the long-range model.
        sigma_minima           
            Sigma is the index for the families of initial-conditions in the
            solution-space. This parameter is the minimum sigma in the range
            of families that the algorithm will try to shoot for.
        sigma_maxima          
            The maximum sigma in the range of families that the algorithm will
            try to shoot for.
        sigma_delta        
            The number of steps between the minima and maxima of sigma the 
            algorithm will try for.
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
    of sigma and the respective field the algorithm integrated to before
    diverging. By plotting the asymptotic field as function of sigma, 
    peaks on the graph give the values of sigma, that is, the initial 
    condition, for which the field does not diverge, and therefore those
    represent real physical solutions to the RG-flow at fixed point. Each of 
    the peaks in the graph therefore represent (potentially) a critical point.
*/
class Shooting_Solver : public Solver {

    protected:

        // Running variables for the algorithm, user input.
        double sigma_minima;
        double sigma_maxima;
        int sigma_delta;

    public:

        /**
            Override implementation to execute shooting solver.
            @see Solver for more details. 
        */
        int execute (
            const std::string& input_filename,
            const std::string& output_filename,
            const std::string& configuration_filename) override {

            // Assigns the configuration file. 
            CONFIGURATION_FILENAME = configuration_filename;
     
            // Opens input and opens or creates output file.
            std::ifstream infile (input_filename);
            std::ofstream outfile (output_filename);

            // Throws an error if input file could not be openend. 
            if (!infile.is_open ()) {
                std::cerr << "Error opening input file\n";
                return 1;
            }

            // Reads configuration. 
            const bool success = read_configuration_from_file 
                (CONFIGURATION_FILENAME, integrator_potential);
            if (!success) {
                std::cerr << "Error reading configuration file\n";
                return 1;
            }

            // Reads file line-by-line. and assigns to class variables.
            std::string line;
            while (std::getline (infile, line)) {
                try {
                    read_parameters_from_line (line, 7);
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
                const std::vector<double> asymptotic_field = 
                    find_asymptotic_field (s_factor);
                // Saves results.
                for (int i = 0; i < sigma_delta; ++i) {
                    const double sigma = sigma_minima +
                        (sigma_maxima - sigma_minima) * 
                        i / (sigma_delta - 1.0);
                    outfile << dimension << ","
                            << anomalous_dimension << ","
                            << s_factor << ","
                            << symmetry_factor_N << ","
                            << sigma << ","
                            << asymptotic_field[i] << "\n";
                }
                outfile << "\n";
            }
            return 0;
        }

        /**
            Method to loop over sigma and evaluate the field for each case.
            @param s_factor    The s-factor corresponding to this run. 
            @return            A vector containing the asymptotic fields 
                               for the range of sigmas chosen. 
        */
        std::vector<double> find_asymptotic_field (const double s_factor) {

            std::vector<double> asymptotic_field;

            // Loops over the range of sigmas.
            for (int i = 0; i < sigma_delta; ++i) {
                double sigma = sigma_minima +
                    (sigma_maxima - sigma_minima) * 
                    i / (sigma_delta - 1.0);
                
                // Creates an instance of the integrator_potential and computes
                // asymptotic field. 
                integrator_potential.initialize (
                    dimension, anomalous_dimension, s_factor, 
                    symmetry_factor_N, sigma);

                // Adds to vector.
                asymptotic_field.push_back (
                    integrator_potential.compute_asymptotic_value ()
                );
            }
            return asymptotic_field;
        }

        /**
            Override of virtual method in abstract class. 
            @see Solver for more details. 
        */
        void set_parameters (const std::vector<std::string>& token) override {
            dimension = std::stod (token[0]);
            anomalous_dimension = std::stod (token[1]);
            s_factor = std::stod (token[2]);
            symmetry_factor_N = std::stod (token[3]);
            sigma_minima = std::stod (token[4]);
            sigma_maxima = std::stod (token[5]);
            sigma_delta = std::stoi (token[6]);
        }
};

#endif