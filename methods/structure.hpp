#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <vector>

// #####################################################
// Begin of problem solver parameters.
// #####################################################

/**
    @brief Input parameters for the eigenvalue problem solver, i.e., user inputs.
    @see Eigenvector_Solver, Eigenvector_Parameter_Parser
*/
struct Eigenvector_Solver_Parameter {
    double dimension;
    double anomalous_dimension;
    double s_factor;
    double sigma;
    double eigenvalue_maxima;
    double eigenvalue_minima;
    int eigenvalue_steps;
};

/**
    @brief Input parameters for the initial condition problem solver, i.e., user inputs.
    @see Initial_Condition_Solver, Initial_Condition_Parameter_Parser
*/
struct Initial_Condition_Solver_Parameter  {
    double dimension;
    double anomalous_dimension;
    double s_factor_minima;
    double s_factor_maxima;
    int s_factor_delta;
    double sigma_minima;
    double sigma_maxima;
    int sigma_steps;
};

/**
    @brief Input parameters for the shooting problem solver, i.e., user inputs.
    @see Shooting_Solver, Shooting_Parameter_Parser
*/
struct Shooting_Solver_Parameter {
    double dimension;
    double anomalous_dimension;
    double s_factor;
    double sigma_minima;
    double sigma_maxima;
    int sigma_steps;
};


/**
    @brief Input parameters for the screening problem solver, i.e., user inputs.
    @see Screening_Solver, Screening_Parameter_Parser
*/
struct Screening_Solver_Parameter {
    double dimension;
    double s_factor;
    double sigma_minima;
    double sigma_maxima;
    double anomalous_dimension_minima;
    double anomalous_dimension_maxima;
    double temperature;
    double cooling_rate;
    int number_of_iterations;
    int number_of_steps;
    double temperature_subprocess;
    double cooling_rate_subprocess;
    int number_of_steps_subprocess;
    int number_of_steps_window_search;
    bool run_subprocess;
    bool run_window_search; 
    bool search_anomalous_dimension;
};

/**
    @brief Input parameters for the grid search problem solver, i.e., user inputs.
    @see Grid_search_Solver, Grid_search_Parameter_Parser
*/
struct Grid_Search_Solver_Parameter {
    double dimension;
    double s_factor;
    double sigma_minima;
    double sigma_maxima;
    int sigma_steps;
    double anomalous_dimension_minima;
    double anomalous_dimension_maxima;
    int anomalous_dimension_steps;
    int number_of_iterations;
    double window_range;
    bool search_anomalous_dimension;
};

// #####################################################
// End of problem solver parameters.
// #####################################################



// #####################################################
// Begin of integrator parameters and result elements.
// #####################################################

/**
    @brief A structure for physical integration parameters for the 
           potential integrator model. 

    @see Integrator_Parameter_Concept, Potential_Integrator_Model
*/
struct Potential_Integrator_Parameter {
    double dimension;
    double anomalous_dimension;
    double s_factor;
    double sigma;
    double dimension_factor;
    double anomalous_constant;
    double implied_s_factor; 
    double s_constant;
};

/**
    @brief A structure for result elements for the potential integrator model. 
    @see Potential_Integrator_Model
*/
struct Potential_Integrator_Result_Element {
    double field;
    double potential_0prime;
    double potential_1prime;
    double potential_2prime;
    double potential_0prime_shape;
    double potential_1prime_shape;
    double potential_2prime_shape;
    double denominator;
    double the_real_denominator;
};

/**
    @brief A structure for physical integration parameters for the 
           eigenvector integrator model. 
    @see Integrator_Parameter_Concept, Potential_Integrator_Parameter, 
         Eigenvector_Integrator_Model
*/
struct Eigenvector_Integrator_Parameter : Potential_Integrator_Parameter {
    std::vector<Potential_Integrator_Result_Element> trajectory;
    double eigenvalue;
};

// #####################################################
// End of integrator parameters and result elements.
// #####################################################



// #####################################################
// Begin of problem solver result elements.
// #####################################################

/**
    @brief Result Element for the eigenvector problem. 
    @see Eigenvector_Problem_Output_Element, Eigenvector_Solver
*/
struct Eigenvector_Solver_Result_Element {
    double eigenvalue;
    double asymptotic_eigenvector;
};

/**
    @brief Result Element for the initial condition problem. 
    @see Initial_Condition_Problem_Output_Element, Initial_Condition_Solver
*/
struct Initial_Condition_Solver_Result_Element {
    double spike;
    double s_factor;
    double anomalous_dimension;
};

/**
    @brief Result Element for the shooting problem. 
    @see Shooting_Problem_Output_Element, Shooting_Solver
*/
struct Shooting_Solver_Result_Element {
    double asymptotic_field;
    double sigma;
    double anomalous_dimension;
};

/**
    @brief Result Element for the screening problem. 
    @see Screening_Problem_Output_Element, Screening_Solver
*/
struct Screening_Solver_Result_Element {
    double asymptotic_field;
    double sigma;
    double anomalous_dimension;
};

/**
    @brief Result Element for the screening problem. 
    @see Grid_Search_Problem_Output_Element, Grid_Search_Solver
*/
struct Grid_Search_Solver_Result_Element {
    double asymptotic_field;
    double sigma;
    double anomalous_dimension;
};

// #####################################################
// Begin of problem solver output elements.
// #####################################################

/**
    @brief The output from a eigenvector problem solver calculation. 
    @see Eigenvector_Solver, Eigenvector_Result_Formatter 
*/
struct Eigenvector_Solver_Output_Element {
    Eigenvector_Solver_Parameter parameter;
    double anomalous_dimension;
    std::vector<Eigenvector_Solver_Result_Element> result;
    std::vector<Potential_Integrator_Result_Element> trajectory;
    std::vector<double> zero_points;
};

/**
    @brief The output from a initial condition problem solver calculation. 
    @see Initial_Condition_Solver, Initial_Condition_Result_Formatter 
*/
struct Initial_Condition_Solver_Output_Element {
    Initial_Condition_Solver_Parameter parameter;
    std::vector<Initial_Condition_Solver_Result_Element> result;
};

/**
    @brief The output from a shooting problem solver calculation. 
    @see Shooting_Solver, Shooting_Result_Formatter 
*/
struct Shooting_Solver_Output_Element {
    Shooting_Solver_Parameter parameter;
    std::vector<Shooting_Solver_Result_Element> result;
};

/**
    @brief The output from a screening problem solver calculation. 
    @see Screening_Solver, Screening_Result_Formatter 
*/
struct Screening_Solver_Output_Element {
    Screening_Solver_Parameter parameter;
    std::vector<Screening_Solver_Result_Element> result;
    std::vector<Screening_Solver_Result_Element> best_result;
};

/**
    @brief The output from a screening problem solver calculation. 
    @see Grid_Search_Solver, Grid_Search_Result_Formatter 
*/
struct Grid_Search_Solver_Output_Element {
    Grid_Search_Solver_Parameter parameter;
    std::vector<Grid_Search_Solver_Result_Element> result;
};


// #####################################################
// End of problem solver output elements.
// #####################################################



// #####################################################
// Begin of miscellaneous.
// #####################################################

/**
    @brief A configuration structure with numerical configurations. 
    @see configuration_reader, Pipeline
*/
struct Configuration {
    double practically_zero;
    double practically_infinity;
    double integration_time_default;
    double field_perturbation;
    double field_threshold;
    double patching_threshold;
    int relaxation_grid_size;
    int number_recalculations;
    double tolerance;
    double absolute_tolerance;
    int window_size;
    double window_range; 
};

struct Instruction_Parameter {
    std::string id;
    std::string input_filename;
    std::string output_filename;
    std::string configuration_filename;
    std::string solver_name;
};

// #####################################################
// End of miscellaneous.
// #####################################################

#endif