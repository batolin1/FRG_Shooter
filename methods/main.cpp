#include "Shooting_Solver.hpp"
#include "Eigenperturbation_Solver.hpp"
#include "Initial_Condition_Solver.hpp"

/**
    The main executor for this programme. User selects whether to run the 
    shooting method or the eigenperturbation method from the terminal. Notice
    that the "input_shooting.txt" and "input_eigenperturbation.txt" files must
    be aleady set up, otherwise an error is thrown. For more details on what 
    each of the methods execute, read comments from the paper and furthermore 
    from "Shooting_Method.hpp" and "Eigenperturbation_method.hpp".
*/
int main() {

    // Simple selection menu for choosing between running the computation for
    // the shooting method of for eigenperturbation method. 
    
    std::string input_file;
    std::string output_file;
    const std::string output_directory = "output-files/";
    const std::string input_directory = "input-files/";
    const std::string  configuration_filename = 
        "configurations/configuration.txt";

    int choice;
    
    std::cout << "=========================\n";
    std::cout << " Choose method:\n";
    std::cout << "1 - Shooting method\n";
    std::cout << "2 - Eigenperturbation method\n";
    std::cout << "3 - Initial condition algorithm\n";
    std::cout << "=========================\n";
    std::cout << "Enter choice: ";

    std::cin >> choice;

    switch (choice) {
        case 1: {
            std::cout << "Running Shooting method...\n" << std::endl;
            Shooting_Solver fixed_point_solver;
            input_file = input_directory + "input_shooting.txt";
            output_file = output_directory + "output_shooting.txt";
            fixed_point_solver.execute 
                (input_file, output_file, configuration_filename);
            std::cout << "Execution completed.";
            return 0;
        }
        case 2: {
            std::cout << "Running Eigenvector method...\n";
            Eigenperturbation_Solver eigenvalue_solver;
            input_file = input_directory + "input_eigenperturbation.txt";
            output_file = output_directory + "output_eigenperturbation.txt";
            eigenvalue_solver.execute
                (input_file, output_file, configuration_filename);
            std::cout << "Execution completed.";
            return 0;
        }
        case 3: {
        std::cout << "Running Initial Condition method...\n";
        Initial_Condition_Solver initial_condition_solver;
        input_file = input_directory + "input_identifier.txt";
        output_file = output_directory + "output_identifier.txt";
        initial_condition_solver.execute
            (input_file, output_file, configuration_filename);
        std::cout << "Execution completed.";
        return 0;
        }
        default: {
            std::cout << "Invalid choice.\n";
            return 1;
        }
    }
    return 0;
}