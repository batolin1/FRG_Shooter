#ifndef SHOOTING_METHOD_HPP
#define SHOOTING_METHOD_HPP

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include "Integrator_Potential_Flow.hpp"


/**
    The class executes the shooting method executes as follows:
    1) The input data is read from "input_shooting.txt" and the results are 
       spilled to "output_shooting.txt". The actual numerical analysis of the
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
        sigma_minima           Sigma is the index for the families of initial-
                            conditions in the solution-space. This parameter
                            is the minimum sigma in the range of families 
                            that the algorithm will try to shoot for.
        sigma_maxima           The maximum sigma in the range of families 
                            that the algorithm will try to shoot for.
        number_of_steps        The number of steps between the minima and
                            maxima of sigma the algorithm will try for.
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
*/
class Shooting_Method {

    public:

        /**
            Simple static method to execute the shooting. 
        */
        static int execute () {

            // Input parameters.
            double dimension;
            double anomalous_dimension;
            double s_factor;
            double symmetry_factor_N;
            double sigma_minima;
            double sigma_maxima;
            int number_of_steps;
            
            // Opens input and opens or creates output file.
            std::ifstream infile ("input_shooting.txt");
            std::ofstream outfile ("output_shooting.txt");
            // Throws an error if input file could not be openend. 
            if (!infile.is_open ()) {
                std::cerr << "Error opening input file\n";
                return 1;
            }

            // Reads file line-by-line.

            std::string line;

            while (std::getline (infile, line)) {
    
                // Ignores empty lines.
                if (line.empty ()) continue;

                std::stringstream ss (line);
                std::string value;
                std::vector<std::string> token;

                // Comma-separated inputs.
                while (std::getline(ss, value, ',')) {
                    token.push_back(value);
                }
                // Ir row given in incorrect format, warns user. 
                if (token.size() != 7) {
                    std::cerr << "Invalid line: " << line << std::endl;
                    continue;
                }
                // Tries to assign rows to the respective variables. 
                try {
                    dimension = std::stod (token[0]);
                    anomalous_dimension = std::stod (token[1]);
                    s_factor = std::stod (token[2]);
                    symmetry_factor_N = std::stod (token[3]);
                    sigma_minima = std::stod (token[4]);
                    sigma_maxima = std::stod (token[5]);
                    number_of_steps = std::stoi (token[6]);
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
                for (int i = 0; i < number_of_steps; ++i) {
                    double sigma = sigma_minima +
                        (sigma_maxima - sigma_minima) * 
                        i / (number_of_steps - 1.0);
                    
                    // Creates an instance of the integrator and computes
                    // asymptotic wavefunction. 
                    Integrator_Potential_Flow integrator (
                        dimension, anomalous_dimension, s_factor, 
                        symmetry_factor_N, sigma);
                    double asymptotic_wavefunction = 
                        integrator.compute_asymptotic_wavefunction ();
                    // Store result as a row to the output file. 
                    outfile << dimension << ","
                            << anomalous_dimension << ","
                            << s_factor << ","
                            << symmetry_factor_N << ","
                            << sigma << ","
                            << asymptotic_wavefunction << "\n";
                }
                outfile << "\n";
            }
            return 0;
        }


};

#endif