#ifndef PARAMETER_READER_HPP
#define PARAMETER_READER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "concept.hpp"
#include "reader/file_parsing_exception.hpp"
#include "logger.hpp"

namespace parameter_reader {

    /**
        @brief 
            A parameter reader method to read parameters from a file and return those as 
            parameter objects to be used in problem solvers. 
            
        @see Parser_Concept, Solver_Concept, Pipeline

        @tparam Parameter    The template for parameters. 
        @tparam Parser       The template for a parser. 
        @param  filename     The filename to access and read the data from 
        @return              A vector of parameters, where each parameter correspond
                             to a single line in the input file. 
    */
    template <typename Parameter, typename Parser>
    requires Parser_Concept<Parser, Parameter>
    static std::vector<Parameter> read (const std::string& filename) {

        std::ostringstream oss;
        oss << "Reading and parsing from filename " << filename << "\n";
        Logger::instance().log (oss.str()); 

        std::vector<Parameter> result;

        std::ifstream file (filename);

        if (!file) {
            throw File_Parsing_Exception ("Cannot open file: " + filename);
        }

        std::string line;

        while (std::getline (file, line)) {

            if (line.empty ()) continue;

            if (line.front () == '/') continue; // Skip comments.

            std::stringstream ss (line);
            std::string value;
            std::vector<std::string> token;

            // Comma-separated entries.
            while (std::getline (ss, value, ',')) {
                token.push_back (value);
            }

            Parameter param;

            if (Parser::parse_token (token, param)) {
                result.push_back (param);
            }
        }
        return result;  
    }
};

#endif