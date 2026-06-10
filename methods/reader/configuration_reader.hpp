#ifndef CONFIGURATION_READER_HPP
#define CONFIGURATION_READER_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <exception>

#include "structure.hpp"
#include "reader/file_parsing_exception.hpp"

namespace configuration_reader {

    /**
        @brief
            a method to read configuration parameters from a configuration file
            and return a structure containing those parameters. 

        @param filename    The string representing the filename. 
        @return            The configuration. 
    */
    Configuration read (const std::string& filename) {

        Configuration configuration;
        std::ifstream configuration_file (filename);

        if (!configuration_file.is_open ()) {
            throw File_Parsing_Exception (
                "Could not find configuration file with provided filename.");
        }

        // Reads the first line only. 
        std::string line;
        std::getline (configuration_file, line);
        // If first line is comment, skip to second line actually. 
        if (line.front () == '/') {
            std::getline (configuration_file, line);
        }
        std::stringstream ss (line);
        std::string value;
        std::vector<std::string> token;

        // Comma-separated inputs.
        while (std::getline (ss, value, ',')) {
            token.push_back (value);
        }

        // Ir row given in incorrect format, throws error. 
        if (token.size () != 12) {
            throw File_Parsing_Exception 
                ("Incorrect number of tokens in configuration file");
        }

        // Assign to class variables and to integrator config.
        try {
            configuration.practically_zero = std::stod (token [0]);
            configuration.practically_infinity = std::stod (token [1]);
            configuration.integration_time_default = std::stod (token [2]);
            configuration.field_perturbation = std::stod (token [3]);
            configuration.field_threshold = std::stod (token [4]);
            configuration.patching_threshold = std::stod (token [5]);
            configuration.relaxation_grid_size = std::stoi (token [6]);
            configuration.number_recalculations = std::stoi (token [7]);
            configuration.tolerance = std::stod (token [8]);
            configuration.absolute_tolerance = std::stod (token [9]);
            configuration.window_size = std::stoi (token [10]);
            configuration.window_range = std::stod (token [11]);
        // Throws error if fails.
        } catch (const std::exception& e) {
            throw File_Parsing_Exception (
                "Could not parse input in configuration file.");
        }
        return configuration;
    }
};

#endif

