#ifndef EIGENPERTURBATION_SOLVER_HPP
#define EIGENPERTURBATION_SOLVER_HPP

#include "Solver.hpp"

/**
    This class executes the eigen-perturbation method, that is to say,  
    for calculating the eigenvalues of the RG flow. In 
    "input_eigenperturbation.txt", the user must supply the input parameters
    for the calculation. A row will represent one execution of the algorithm.
    Each input parameter is provided by columns, in comma-separated variables. 
    The input parameters are the following, that is to say:
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
class Eigenperturbation_Solver : public Solver {

    private:

        // Running parameters for this solver, user input. 
        double eigenvalue_minima;
        double eigenvalue_maxima;
        int eigenvalue_delta;

    public:

        /**
            Override implementation to execute eigenperturbation solver.
            @see Solver for more details. 
        */
        int execute (
            const std::string& input_filename,
            const std::string& output_filename,
            const std::string& configuration_filename) override {
            
            // Sets the configuration file. 
            CONFIGURATION_FILENAME = configuration_filename;

            // Opens input and opens or creates output file.
            std::ifstream infile (input_filename);
            std::ofstream outfile (output_filename);

            // Throws an error if input file could not be openend. 
            if (!infile.is_open ()) {
                std::cerr << "Error opening input file\n";
                return 1;
            }
        
            // Reads the input file by line. 

            std::string line;
        
            while (std::getline (infile, line)) {
        
                // Tries to assign rows to the respective variables. 
                try {
                    read_parameters_from_line (line, 7);
                // If for a particular row any of the tokens provided could not be
                // converted, ignore row and throws a warning to the user. 
                } catch (...) {
                    std::cerr << "Conversion error in line: " << line << std::endl;
                    continue;
                }
        
                // Sets configuration to the integrator of the RG potential.
                success = read_configuration_from_file 
                    (CONFIGURATION_FILENAME, integrator_potential);
                if (!success) {
                    std::cerr << "Error reading configuration (shooting) file\n";
                    return 1;
                }
                // Sets configuration to the integrator of the eigenvector.
                integrator_eigenvector.set_configuration(
                    practically_zero,
                    practically_infinity,
                    integration_time_default,
                    field_perturbation,
                    field_threshold
                );

                // Computes the trajectories for the potential and 
                // its derivatives.
                integrator_potential.initialize (
                    dimension, 
                    anomalous_dimension, 
                    s_factor, 
                    1.0, 
                    sigma);
                integrator_potential.compute_asymptotic_value ();

                // Saves the trajectory for the potential if requested.
                if (save_trajectories) {
                    save_trajectory_to_file ();
                }

                // Loops over the proposed eigenvalues. 
                for (int i = 0; i < eigenvalue_delta; ++i) {
                    double eigenvalue = eigenvalue_minima +
                        (eigenvalue_maxima - eigenvalue_minima) * 
                        i / (eigenvalue_delta - 1.0);
        
                    // Solves the shooting problem for the eigenvectors.
                    integrator_eigenvector.initialize (
                        dimension,
                        anomalous_dimension, 
                        s_factor,
                        sigma, 
                        eigenvalue, 
                        integrator_potential.get_field (),
                        integrator_potential.get_potential_0prime (),
                        integrator_potential.get_potential_1prime (),
                        integrator_potential.get_potential_2prime ());
                    
                    // Adds the result as a row to the file. 
                    double asymptotic_eigenvector = 
                        integrator_eigenvector.compute_asymptotic_value ();
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

        /**
            Override of virtual method in abstract class.
            @see Solver for more details. 
        */
        void set_parameters (const std::vector<std::string>& token) override {
            dimension = std::stod (token[0]);
            anomalous_dimension = std::stod (token[1]);
            s_factor = std::stod (token[2]);
            sigma = std::stod (token[3]);
            eigenvalue_minima = std::stod (token[4]);
            eigenvalue_maxima = std::stod (token[5]);
            eigenvalue_delta = std::stoi (token[6]);
        }

        /**
            Helper method to save trajectory to file. 
        */
        void save_trajectory_to_file () {

            // For each trajectory, creates a specific filename and creates it.
            std::stringstream string_stream;
            string_stream
                    << "output-files/trajectories/"
                    << "dimension" << "=" << dimension << "_"
                    << "anomalous-dimension" << "=" << anomalous_dimension << "_"
                    << "s-factor" << "=" << s_factor << "_"
                    << "sigma" << "=" << sigma
                    << ".txt";
            const std::string file_path_and_name = string_stream.str();
            std::ofstream outfile (file_path_and_name);
            // Warns user if path not found (file not created).
            if (!outfile.is_open()) {
                std::cerr << "Error opening file: " << file_path_and_name << "\n";
                return;
            }
            // Loops over the elements, saving to file. 
            const int number_of_elements = integrator_potential.get_field ().size ();
            for (int i = 0; i < number_of_elements; i++) {
                outfile << integrator_potential.get_field() [i] << ",";
            }
            outfile << "\n";
            for (int i = 0; i < number_of_elements; i++) {
                outfile << integrator_potential.get_potential_0prime () [i] << ",";
            }
            outfile << "\n";
            for (int i = 0; i < number_of_elements; i++) {
                outfile << integrator_potential.get_potential_1prime () [i] << ",";
            }
            outfile << "\n";
            for (int i = 0; i < number_of_elements; i++) {
                outfile << integrator_potential.get_potential_2prime () [i] << ",";
            }
            outfile << "\n";
        }
};

#endif