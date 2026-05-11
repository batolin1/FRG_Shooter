#ifndef PARSER_HPP
#define PARSER_HPP

#include <vector>
#include <string>

#include "structure.hpp"



/**
    @brief A parser for input parameters for the shooting problem. 
    @see Parser_Concept
*/
struct Shooting_Problem_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Shooting_Problem_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 7) {
            return false;
        }
        
        try {
            param.dimension = std::stod (token [0]);
            param.anomalous_dimension = std::stod (token [1]);
            param.s_factor = std::stod (token [2]);
            param.symmetry_factor_N = std::stod (token [3]);
            param.sigma_minima = std::stod (token [4]);
            param.sigma_maxima = std::stod (token [5]);
            param.sigma_delta = std::stoi (token [6]);
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
struct Initial_Condition_Problem_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Initial_Condition_Problem_Solver_Parameter& param) {
            
        // Ir row given in incorrect format, igore row and warns user.
        if (token.size () != 9) {
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
            param.sigma_delta = std::stoi (token [7]);
            param.symmetry_factor_N = std::stod (token [8]);
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
struct Eigenvector_Problem_Parameter_Parser {

    static bool parse_token (
        const std::vector<std::string>& token, Eigenvector_Problem_Solver_Parameter& param) {
            
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
            param.eigenvalue_delta = std::stoi (token [6]);
        } catch (...) {
            return false;
        }

        return true;
    }
};

#endif