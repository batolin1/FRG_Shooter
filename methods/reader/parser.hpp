#ifndef PARSER_HPP
#define PARSER_HPP

#include <vector>
#include <string>

#include "structure.hpp"
#include "serialization.hpp"



/**
    @brief A parser for input parameters for the shooting problem. 
    @see Parser_Concept
*/
struct Shooting_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Shooting_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 6) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.anomalous_dimension = std::stod (token [1]);
            param.s_factor = std::stod (token [2]);
            param.sigma_minima = std::stod (token [3]);
            param.sigma_maxima = std::stod (token [4]);
            param.sigma_steps = std::stoi (token [5]);
        } catch (...) {
            return false;
        }
        return true;
    }
};


/**
    @brief A parser for input parameters for the screening problem. 
    @see Parser_Concept
*/
struct Screening_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Screening_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 17) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.s_factor = std::stod (token [1]);
            param.sigma_minima = std::stod (token [2]);
            param.sigma_maxima = std::stod (token [3]);
            param.anomalous_dimension_minima = std::stod (token [4]);
            param.anomalous_dimension_maxima = std::stod (token [5]);
            param.temperature = std::stod (token [6]);
            param.cooling_rate = std::stod (token [7]);
            param.number_of_iterations = std::stoi (token [8]);
            param.number_of_steps = std::stoi (token [9]);
            param.temperature_subprocess = std::stod (token [10]);
            param.cooling_rate_subprocess = std::stod (token [11]);
            param.number_of_steps_subprocess = std::stoi (token [12]);
            param.number_of_steps_window_search = std::stoi (token [13]);
            param.run_subprocess = token [14] == "true" ? true : false;
            param.run_window_search = token [15] == "true" ? true : false; 
            param.search_anomalous_dimension = token [10] == "true" ? true : false;
        } catch (...) {
            return false;
        }
        return true;
    }
};

/**
    @brief A parser for input parameters for the grid search problem. 
    @see Parser_Concept
*/
struct Grid_Search_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Grid_Search_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 11) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.s_factor = std::stod (token [1]);
            param.sigma_minima = std::stod (token [2]);
            param.sigma_maxima = std::stod (token [3]);
            param.sigma_steps = std::stoi (token [4]);
            param.anomalous_dimension_minima = std::stod (token [5]);
            param.anomalous_dimension_maxima = std::stod (token [6]);
            param.anomalous_dimension_steps = std::stoi (token [7]);
            param.number_of_iterations = std::stoi (token [8]);
            param.window_range = std::stod (token [9]);
            param.search_anomalous_dimension = token [10] == "true" ? true : false;
        } catch (...) {
            return false;
        }
        return true;
    }
};


/**
    @brief A parser for input parameters for the initial condition problem. 
    @see Parser_Concept
*/
struct Initial_Condition_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Initial_Condition_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 8) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.anomalous_dimension = std::stod (token [1]);
            param.s_factor_minima = std::stod (token [2]);
            param.s_factor_maxima = std::stod (token [3]);
            param.s_factor_delta= std::stoi (token [4]);
            param.sigma_minima = std::stod (token [5]);
            param.sigma_maxima = std::stod (token [6]);
            param.sigma_steps = std::stoi (token [7]);
        } catch (...) {
            return false;
        }
        return true;
    }
};


/**
    @brief A parser for input parameters for the eigenvector problem. 
    @see Parser_Concept
*/
struct Eigenvector_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Eigenvector_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row.
        if (token.size () != 7) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.anomalous_dimension = std::stod (token [1]);
            param.s_factor = std::stod (token [2]);
            param.sigma = std::stod (token [3]);
            param.eigenvalue_minima = std::stod (token [4]);
            param.eigenvalue_maxima = std::stod (token [5]);
            param.eigenvalue_steps = std::stoi (token [6]);
        } catch (...) {
            return false;
        }

        return true;
    }
};

struct Instruction_Parameter_Parser {

    static bool parse_token (

        const std::vector<std::string>& token, Instruction_Parameter& param) {
            
        // Ir row given in incorrect format, igore row.
        if (token.size () != 4) {
            return false;
        }
        
        try {
            param.input_filename = token [0];
            param.output_filename = token [1];
            param.configuration_filename = token [2];
            param.solver_name = token [3];
            param.id = serialization::generate_uuid ();
        } catch (...) {
            return false;
        }

        return true;
    }

};

#endif