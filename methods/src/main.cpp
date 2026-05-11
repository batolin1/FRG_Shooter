#include <iostream>
#include <stdexcept>

#include "pipeline.hpp"
#include "structure.hpp"
#include "reader/configuration_reader.hpp"
#include "reader/parser.hpp"
#include "writer/formatter.hpp"
#include "solver/solver/shooting_problem_solver.hpp"
#include "solver/solver/eigenvector_problem_solver.hpp"
#include "solver/solver/initial_condition_problem_solver.hpp"


int main() {

    try {
        Configuration config =
            configuration_reader::read("configurations/configuration.txt");
        
         std::cout << "Select solver:\n";
        std::cout << "1. Shooting method\n";
        std::cout << "2. Eigenvector method\n";
        std::cout << "3. Initial Condition method\n";
        std::cout << "0. Exit\n";
        std::cout << "Choice: ";

        int choice;
        std::cin >> choice;

        switch (choice) {

            case 1:
                return
                    Pipeline<
                        Shooting_Problem_Solver,
                        Shooting_Problem_Parameter_Parser,
                        Shooting_Problem_Result_Formatter
                        >::run(
                            "input_files/input_shooting.txt",
                            "output_files/output_shooting.txt",
                            config);

            case 2:
            return 
                Pipeline<
                    Eigenvector_Problem_Solver,
                    Eigenvector_Problem_Parameter_Parser,
                    Eigenvector_Problem_Result_Formatter,
                    Eigenvector_Problem_Trajectory_Formatter,
                    Eigenvector_Problem_Zeroes_Formatter
                    >::run (
                        "input_files/input_eigenperturbation.txt",
                        "output_files/output_eigenperturbation.txt",
                        config);

            case 3:
                return
                    Pipeline<
                        Initial_Condition_Problem_Solver,
                        Initial_Condition_Problem_Parameter_Parser,
                        Initial_Condition_Problem_Result_Formatter
                        >::run(
                            "input_files/input_initial_condition.txt",
                            "output_files/output_initial_condition.txt",
                            config);

            case 0:
                std::cout << "Exiting.\n";
                return 0;

            default:
                std::cerr << "Invalid choice.\n";
                return 1;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}