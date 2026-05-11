#include "Shooting_Solver.hpp"
#include "Eigenperturbation_Solver.hpp"
#include "Initial_Condition_Solver.hpp"

/**
    Instructions:
    1) User selects which of the solvers to be executed. Inputs are configured
       from the "input-files" directory. Output is saved on the "output-files"
       directory. Configurations can be set up from the "configurations" 
       directory. These refer to numerical-stability configurations that affect
       only the simulation speed and accuracy but not the actual physics. 
    2) "Shooting solver", given inputs, identifies the points of multicritical 
       phase transitions. "Shooting_plotter.py" then generates the spike plots
       corresponding to those phase transitions. 
    3) "Eigenperturbation solver", given inputs, finds the behaviour of the 
        asymptotic eigenvector, which is then used to find the RG eigenvectors.
        This is later plotted on "Shooting_plotter.py" too. 
    4) "Initial Condition Solver", given inputs, specifically lists out all of
        the (potential) fixed points over a range of s_factors. This can then 
        be provided as input to the "Eigenperturbation solver" to find the 
        corresponding RG eigenvalues. 
*/
int main () {

    // The labels for the directories 
    std::string input_file;
    std::string output_file;
    const std::string output_directory = "output-files/";
    const std::string input_directory = "input-files/";
    const std::string  configuration_filename = 
        "configurations/configuration.txt";


    // Simple selection menu.
    int choice;
    std::cout << "=========================\n";
    std::cout << " Choose solver:\n";
    std::cout << "1 - Fixed point solver\n";
    std::cout << "2 - Eigenperturbation solver\n";
    std::cout << "3 - Initial condition solver\n";
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
        input_file = input_directory + "input_initial_condition.txt";
        output_file = output_directory + "output_initial_condition.txt";
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