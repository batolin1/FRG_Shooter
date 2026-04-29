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

class Solver {

    public:

        virtual ~Solver() = default;

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

        // Running parameters
        double sigma_minima;
        double sigma_maxima;
        int sigma_delta;
        double eigenvalue_minima;
        double eigenvalue_maxima;
        int eigenvalue_delta;
        double s_factor_minima;
        double s_factor_maxima;
        int s_factor_delta;
        int nth_critical_point;

        // Configuration parameters
        double practically_zero;
        double practically_infinity;
        double integration_time_default;
        double wavefunction_perturbation;
        double wavefunction_threshold;
        bool save_trajectories;
        bool success;
        std::string CONFIGURATION_FILENAME;

        // Integrators
        Integrator_Eigenvector integrator_eigenvector;
        Integrator_Potential integrator_potential;

        void read_parameters_from_line (
            const std::string& line, 
            const int expected_number_of_tokens) {

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

            set_parameters (token);
        }
       
        bool read_configuration_from_file (
            const std::string& filename,
            Integrator& integrator) {

            std::ifstream configuration_file (filename);

            if (!configuration_file.is_open ()) {
                std::cerr << "Error opening configuration file\n";
                return false;
            }

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

            try {
                practically_zero = std::stod (token[0]);
                practically_infinity = std::stod (token[1]);
                integration_time_default = std::stod (token[2]);
                wavefunction_perturbation = std::stod (token[3]);
                wavefunction_threshold= std::stod (token[4]);
                save_trajectories = std::stoi (token[5]);
                integrator.set_configuration (
                    practically_zero, 
                    practically_infinity,
                    integration_time_default, 
                    wavefunction_perturbation, 
                    wavefunction_threshold);
            } catch (...) {
                std::cerr << "Error in configuration file: " << line << std::endl;
                return false;
            }
            return true;
        }

        virtual void set_parameters (const std::vector<std::string>&) = 0;
};

#endif