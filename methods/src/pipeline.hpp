#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include <string>
#include <vector>

#include "structure.hpp"
#include "reader/parameter_reader.hpp"
#include "writer/result_writer.hpp"
#include "concept.hpp"


/**
    @brief A pipeline structure to execute the pipeline, that is, read the data from file,
           calculate numerical results, and store them back onto output files. 
    
    @see parameter_reader, result_writer, configuration_reader, Solver_Concept
*/
template <
    Solver_Concept Solver, 
    typename Parser, 
    typename Formatter, 
    typename Trajectory_Formatter = void,
    typename Zero_Crossing_Formatter = void> 
requires 
    Parser_Concept<Parser, typename Solver::Parameter> && 
    Formatter_Concept<Formatter, typename Solver::Output_Element> 
struct Pipeline {

    using Parameter = typename Solver::Parameter;
    using Output_Element = typename Solver::Output_Element;

    /**
        @brief Method to run the pipieline 

        @param input_file    The string to the directory containing the input file.
        @param output_file   The string to the directory containing the output file. 
        @param config        The (numerical) configurations for the simulation. 
    */
    static int run(
        const std::string& input_file, 
        const std::string& output_file,
        const Configuration& config) {

        // reads input file and parses onto parameter objects for the solver.
        std::vector<Parameter> params = parameter_reader::read<Parameter, Parser> (input_file);

        // Initializes solver and process parameters.

        Solver solver (config);
        solver.initialize ();

        for (const Parameter& p : params) {
            solver.process_parameter (p);
        }

        // Writes results to file.
        result_writer::write<Formatter, Output_Element> (output_file, solver.get_result ());

        // Saves additional data if required -- The potential and derivatives as function of field. 
        if constexpr (!std::is_same_v<Trajectory_Formatter, void>) {
            if (config.save_trajectories) {
                result_writer::make_extended_trajectory (solver.get_result ());
                result_writer::write<Trajectory_Formatter, Output_Element>(
                    output_file + "_trajectory", solver.get_result ());
            }
        }

        // Saves additional data if required -- The points of zero crossing.
        if constexpr (!std::is_same_v<Zero_Crossing_Formatter, void>) {
            result_writer::write<Zero_Crossing_Formatter, Output_Element>(
                output_file + "_zero_crossings", solver.get_result ());
        }

        return 0;
    }
};

#endif