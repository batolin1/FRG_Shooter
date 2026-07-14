#ifndef RESULT_WRITER_HPP
#define RESULT_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

#include "concept.hpp"
#include "structure.hpp"
#include "solver/integrator/utils/parameter_builder.hpp"
#include "logger.hpp"

namespace result_writer {

    /**
        @brief A writer method to write output elements from problem solvers to a text file.

        @see Solver_Concept, Formatter_Concept, Pipeline

        @tparam Formatter         The template for formatters. 
        @tparam Output_Element    The template for an output element. 
        @param  filename          The string corresponding to the filename. 
        @param  output            The vector containing output elements. 
    */
    template <typename Formatter, typename Output_Element>
    requires Formatter_Concept<Formatter, Output_Element>
    void write (const std::string& filename, const std::vector<Output_Element>& output) {

        // Logs to logger that this action is taking place. 
        std::ostringstream oss;
        oss << "Formatting and writting from filename " << filename << "\n";
        Logger::instance().log (oss.str()); 

        std::ofstream file (filename);

        if (!file) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        for (const Output_Element& output_element : output) {
            file << Formatter::format (output_element) << "\n";
        }
    }
};

#endif