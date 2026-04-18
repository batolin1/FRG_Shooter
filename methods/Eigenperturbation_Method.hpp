#ifndef EIGENPERTURBATION_METHOD_HPP
#define EIGENPERTURBATION_METHOD_HPP

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include "Integrator_Eigenvector.hpp"
#include "Integrator_Potential_Flow.hpp"

/**
    This class executes the shooting method for the eigenperturbations, 
    that is to say, for calculating the eigenvalues of the RG flow. In 
    "input_eigenperturbation.txt", the user must supply the input 
    parameters for the calculation. Each row will represent one execution 
    and/or calculation of this shooting algorithm, and each input parameter 
    is provided by columns, in comma-separated variables. The input parameters
    are the following, that is to say:
        1) dimension
            The dimensionality of the model, example a one-dimensional,
            two-dimensional, or three-dimensional scalar theory. 
        2) anomalous dimension
            Specifically for short range models, the anomalous dimension.
            For long range models, this parameter should be set to 0.0. 
        3) s-factor 
            The s-factor specifying the power-law dependency of the 
            interaction to the distance. For the short-range model, 
            set this parameter to 2.0.
        4) sigma
            The initial condition which triggers a physical fixed-point
            solution to the RG-flow ODE. This parameter is obtained, 
            specifically by first running the Shooting_Method. See more 
            details to follow. 
    The user must first solve the shooting problem in order to find physical
    fixed-point solutions to the RG-flow ODE. Each of those sigmas will 
    (potentially) correspond to a critical point of the theory. The problem is
    solved self-consistently in the sense that first the values of sigma must 
    be found and then provided as input to this class. Only one value of sigma 
    is chosen for each execution. The user must therefore choose one of the 
    critical points at a time, and this class will output the results for the 
    shooting problem for this particular critical point specified by sigma. 
    The output is provided in "eigenperturbation_output.txt". For each row, the 
    first four columns are the labels "dimension", "anomalous_dimension", 
    "s_factor" and "sigma" corresponding to the particular simulation, and the 
    last two rows are the eigenvalues and the and the asymptotic eigenvalues 
    respectively corresponding to the same. This data is processed in more 
    details on the python script "eigenperturbation_plotter.py" and 
    specifications are contained there.
    By plotting the asymptotic eigenvectors as function of the eigenvalues, the
    asymptotic condition implies that the eigenvalues for which the 
    eigenvectors equal to zero correspond to real, physical solutions to the 
    system; that is to say, those are the real RG eigenvalues of the system.
*/
class Eigenperturbation_Method {

    public:

        /**
            Static method to realize the execution. 
        */
        static int execute () {
            
            // Input parameters to be read from the text file. 
            double dimension;
            double anomalous_dimension;
            double s_factor;
            double sigma;
            // Those trajectories are calculated from the Shooting_Method.
            // Those correspond to the evolution of the potential and its 
            // derivatives, as the wavefunction is integrated outwards.
            trajectory wavefunction;
            trajectory potential_0prime;
            trajectory potential_1prime;
            trajectory potential_2prime;
            
            // Opens input and opens or creates output file.
            std::ifstream infile ("input_eigenperturbation.txt");
            std::ofstream outfile ("output_eigenperturbation.txt");
            // Throws an error if input file could not be openend. 
            if (!infile.is_open ()) {
                std::cerr << "Error opening input file\n";
                return 1;
            }
        
            // Reads the input file by line. 

            std::string line;
        
            while (std::getline (infile, line)) {
        
                if (line.empty ()) continue;
        
                std::stringstream ss (line);
                std::string value;
                std::vector<std::string> token;
        
                // Comma-separated entries.
                while (std::getline(ss, value, ',')) {
                    token.push_back(value);
                }
                // Ir row given in incorrect format, igore row and warns user.
                if (token.size() != 4) {
                    std::cerr << "Invalid line: " << line << std::endl;
                    continue;
                }
                // Tries to assign rows to the respective variables. 
                try {
                    dimension = std::stod (token[0]);
                    anomalous_dimension = std::stod (token[1]);
                    s_factor = std::stod (token[2]);
                    sigma = std::stod (token[3]);
                // If for a particular row any of the tokens provided could not be
                // converted, ignore row and throws a warning to the user. 
                } catch (...) {
                    std::cerr << "Conversion error in line: " << line << std::endl;
                    continue;
                }
        
                // Creates an instance of the integrator for the RG flow at the
                // and extracts the trajectories for the potential and its 
                // derivatives, as the ODE is integrated from the wavefunction 
                // towards outwards.
                Integrator_Potential_Flow integrator (
                    dimension, anomalous_dimension, s_factor, 
                    1.0, sigma);
                double result = integrator.compute_asymptotic_wavefunction ();
                wavefunction = integrator.get_wavefunction ();
                potential_0prime = integrator.get_potential_0prime ();
                potential_1prime = integrator.get_potential_1prime ();
                potential_2prime = integrator.get_potential_2prime ();

                // Most eigenvalues are expected to stay within a small range. 
                // Therefore, we constrain to [-5, +5]. In reality, one 
                // remembers that exclusively that eigenvalues smaller than 
                // zero are relevant eigendirections; and the positive 
                // solutions are irrelevant eigendirections that provide higher
                // order corrections.
                const double eigenvalue_minima = -5.0;
                const double eigenvalue_maxima = +5.0;
                const int number_of_steps = 4000;
                // Loops over the proposed eigenvalues. 
                for (int i = 0; i < number_of_steps; ++i) {
                    double eigenvalue = eigenvalue_minima +
                        (eigenvalue_maxima - eigenvalue_minima) * 
                        i / (number_of_steps - 1.0);
        
                    // Solves the shooting problem for the eigenvectors.
                    Integrator_Eigenvector integrator (
                        dimension, anomalous_dimension, s_factor,
                        sigma, eigenvalue, wavefunction, potential_0prime,
                        potential_1prime, potential_2prime);
                    
                    // Adds the result as a row to the file. 
                    double asymptotic_eigenvector = 
                        integrator.compute_asymptotic_eigenvector ();
                    outfile << dimension << ","
                            << anomalous_dimension << ","
                            << s_factor << ","
                            << sigma << ","
                            << eigenvalue << ","
                            << asymptotic_eigenvector << "\n";
                }
            }
            return 0;       
        }
        

};

#endif