#include <iostream>
#include <stdexcept>
#include <future>
#include <fstream>

#include "pipeline.hpp"
#include "structure.hpp"
#include "reader/configuration_reader.hpp"
#include "reader/parser.hpp"
#include "reader/parameter_reader.hpp"
#include "writer/formatter.hpp"
#include "solver/solver/shooting_solver.hpp"
#include "solver/solver/eigenvector_solver.hpp"
#include "solver/solver/initial_condition_solver.hpp"
#include "solver/solver/screening_solver.hpp"
#include "solver/solver/grid_search_solver.hpp"
#include "logger.hpp"


int process_request (const Instruction_Parameter& instruction, const int thread_number) {

    const std::string eigenvector_solver = "eigenvector_solver";
    const std::string grid_search_solver = "grid_search_solver";
    const std::string screening_solver = "screening_solver";
    const std::string shooting_solver = "shooting_solver";
    const std::string initial_condition_solver = "initial_condition_solver";

    const Configuration config = configuration_reader::read (instruction.configuration_filename);
    const std::string input_filename = instruction.input_filename;
    const std::string output_filename = instruction.output_filename;

    std::ostringstream oss;

    oss << "Processing the instruction for thread number " << thread_number << ": " << "\n"
        << "id: " << instruction.id << "\n"
        << "label: " << instruction.solver_name << "\n"
        << "input_filename: " << instruction.input_filename << "\n"
        << "output_filename: " << instruction.output_filename << "\n"
        << "configuration_filename: " << instruction.configuration_filename;

    Logger::instance().log (oss.str()); 

    if (instruction.solver_name == eigenvector_solver) {
        return 
            Pipeline<
                Eigenvector_Solver,
                Eigenvector_Parameter_Parser,
                Eigenvector_Result_Formatter ,
                Eigenvector_Problem_Zeroes_Formatter,
                Eigenvector_Problem_Trajectory_Formatter
                >::run (
                    input_filename,
                    output_filename,
                    config);
    } else if (instruction.solver_name == grid_search_solver) {
        return 
            Pipeline< 
                Grid_Search_Solver,
                Grid_Search_Parameter_Parser,
                Grid_Search_Result_Formatter
                >::run(
                    input_filename,
                    output_filename,
                    config);
    } else if (instruction.solver_name == screening_solver) {
        return 
            Pipeline< 
                Screening_Solver,
                Screening_Parameter_Parser,
                Screening_Result_Formatter,
                Screening_Best_Result_Formatter
                >::run(
                    input_filename,
                    output_filename,
                    config);
    } else if (instruction.solver_name == shooting_solver) {
        return
            Pipeline<
                Shooting_Solver,
                Shooting_Parameter_Parser,
                Shooting_Result_Formatter 
                >::run(
                    "input_files/input_shooting.txt",
                    "output_files/output_shooting.txt",
                    config);
    } else if (instruction.solver_name == initial_condition_solver) {
        return
            Pipeline<
                Initial_Condition_Solver,
                Initial_Condition_Parameter_Parser,
                Initial_Condition_Result_Formatter 
                >::run(
                    input_filename,
                    output_filename,
                    config);
    } else {
        std::cout << "Could not understand " << instruction.solver_name << std::endl;
    }
    return 1;
}

int main () {

    // Reads the instructions
    std::string input_file = "instruction.txt";
    std::vector<Instruction_Parameter> instructions =
        parameter_reader::read<Instruction_Parameter, Instruction_Parameter_Parser>(input_file);

    std::ostringstream oss;
    oss << "Instructions provided for " << instructions.size () << " total requests."; 
    Logger::instance().log (oss.str()); 

    // Use multithreading. 
    std::vector <std::future<void>> futures;

    int thread_number = 1;

    // For each instruction, processes request given by user.
    for (const Instruction_Parameter& instruction : instructions) {
        // Adds to asynch tasks. 
        futures.push_back(
            std::async(std::launch::async,
                [instruction, thread_number]() {
                    process_request(instruction, thread_number);
                }
            )
        );
        thread_number++;
    }

    // Send. 
    for (auto& future : futures) {
        future.get();
    }

    std::ofstream log_file("logs.txt");
    for (const std::string& s : Logger::instance().get_messages()) {
        log_file << s << "\n";
    }

    return 0;
}