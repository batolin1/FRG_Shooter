#ifndef Solver_HPP
#define Solver_HPP

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include "Integrator_Eigenvector.hpp"
#include "Integrator_Potential.hpp"

/**
    An abstract solver class to read inputs from files, realize particular 
    operations, and save outputs. 
*/
class Solver {

    public:

        virtual ~Solver() = default;

        /**
            Virtual method to execute a solver action.
            @param input_filename            The filename to access inputs.
            @param output_filename           The filename to save outputs.
            @param configuration_filename    The filename to access configs.
            @return                          Whether operation succeeded.
        */
        virtual int execute (
            const std::string&, 
            const std::string&, 
            const std::string&) = 0;
        
    protected:

        // Physical parameters
        double dimension;
        double anomalous_dimension;
        double s_factor;
        double symmetry_factor_N;
        double sigma;

        // Configuration parameters
        std::string CONFIGURATION_FILENAME;
        double practically_zero;
        double practically_infinity;
        double integration_time_default;
        double field_perturbation;
        double field_threshold;
        bool save_trajectories;
        bool success;
        
        // Integrators
        Integrator_Eigenvector integrator_eigenvector;
        Integrator_Potential integrator_potential;

        /**
            A method to, given a line, read parameters and assign them to the 
            respective class variables. 
            @param line                         The line.
            @param expected_number_of_tokens    The number of tokens the line 
                                                is expected to have. 
        */
        void read_parameters_from_line (
            const std::string& line, 
            const int expected_number_of_tokens) {

            // Ignores empty line.
            if (line.empty ()) return;

            std::stringstream ss (line);
            std::string value;
            std::vector<std::string> token;

            // Comma-separated entries.
            while (std::getline (ss, value, ',')) {
                token.push_back (value);
            }
            // Ir row given in incorrect format, igore row and warns user.
            if (token.size () != expected_number_of_tokens) {
                std::cerr << "Invalid line: " << line << std::endl;
                return;
            }
            // Sets the parameters to class variables.
            set_parameters (token);
        }
       
        /**
            A method to, given the configuration filename, read the 
            configuration file and assign the corresponding 
            parameters to an @Integrator object. 
            @param filename      The filename for the configurations.
            @param integrator    The @Integrator object to be assigned.
            @return              Boolean whether the operation was successful.
        */
        bool read_configuration_from_file (
            const std::string& filename,
            Integrator& integrator) {

            // Reads file
            std::ifstream configuration_file (filename);
            // Throws warning if reading fails.
            if (!configuration_file.is_open ()) {
                std::cerr << "Error opening configuration file\n";
                return false;
            }

            // Reads the first line only. 
            std::string line;
            std::getline (configuration_file, line);
            std::stringstream ss (line);
            std::string value;
            std::vector<std::string> token;

            // Comma-separated inputs.
            while (std::getline (ss, value, ',')) {
                token.push_back (value);
            }

            // Ir row given in incorrect format, warns user. 
            if (token.size () != 6) {
                std::cerr << "Invalid configuration file: " << line << std::endl;
                return false;
            }

            // Assign to class variables and to integrator config.
            try {
                practically_zero = std::stod (token[0]);
                practically_infinity = std::stod (token[1]);
                integration_time_default = std::stod (token[2]);
                field_perturbation = std::stod (token[3]);
                field_threshold= std::stod (token[4]);
                save_trajectories = std::stoi (token[5]);
                integrator.set_configuration (
                    practically_zero, 
                    practically_infinity,
                    integration_time_default, 
                    field_perturbation, 
                    field_threshold);
            // Warns user if fails. 
            } catch (...) {
                std::cerr << "Error in configuration file: " << line << std::endl;
                return false;
            }
            return true;
        }

        /**
            Virtual method to set global parameters given a provided token
            @param token    The token containing the parameters, read from a 
                            line. 
        */
        virtual void set_parameters (const std::vector<std::string>&) = 0;
};

#endif