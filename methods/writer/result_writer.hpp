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

    /**
        @brief A helper method to write an analytical continuation for the trajectory of the 
               potential and its derivatives using the asymptotic formula. This method is used here
               for visualization and valdiation purposes only and it does not affect calculation.

        @see Pipeline, parameter_builder, Eigenvector_Solver, Eigenvector_Integrator

        @param output    The vector containing output elements from eigenvector problem. 
    */
    template <typename Output_Element>
    void make_extended_trajectory (std::vector<Output_Element>& output) {

        // for (Output_Element& o : output) {
        //     Eigenvector_Integrator_Parameter param = 
        //         parameter_builder::make_eigenvector_integrator_parameter (o.parameter);
        //     param.anomalous_dimension = o.anomalous_dimension;
        //     parameter_builder::calculate_additional_parameters (param);
        //     param.trajectory = o.trajectory;
        //     param.eigenvalue = 0;
        //     double start_field = o.trajectory.back ().field;
        //     // Extend by a third. 
        //     const double end_field = start_field * 1.5; 
        //     while (start_field < end_field) {
        //         Potential_Integrator_Result_Element trajectory_element;
        //         start_field += 10e-3;
        //         trajectory_element.field = start_field;
        //         trajectory_element.potential_0prime = 
        //             parameter_builder::get_asymtotic_continuation (start_field, 0, param);
        //         trajectory_element.potential_1prime = 
        //             parameter_builder::get_asymtotic_continuation (start_field, 1, param);
        //         trajectory_element.potential_2prime = 
        //             parameter_builder::get_asymtotic_continuation (start_field, 2, param);
        //         o.trajectory.push_back (trajectory_element);
        //     }
        // }
    }

};

#endif