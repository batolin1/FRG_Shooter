#include "Shooting_Method.hpp"
#include "Eigenperturbation_Method.hpp"

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
    
    int choice;
    
    std::cout << "=========================\n";
    std::cout << " Choose method:\n";
    std::cout << "1 - Shooting method\n";
    std::cout << "2 - Eigenperturbation method\n";
    std::cout << "=========================\n";
    std::cout << "Enter choice: ";

    std::cin >> choice;

    switch (choice) {
        case 1:
            std::cout << "Running Shooting method...\n" << std::endl;
            Shooting_Method::execute ();
            std::cout << "Execution completed.";
            return 0;
        case 2:
            std::cout << "Running Eigenvector method...\n";
            Eigenperturbation_Method::execute (); 
            std::cout << "Execution completed.";
            return 0;
        default:
            std::cout << "Invalid choice.\n";
            return 1;
    }
    return 0;
}