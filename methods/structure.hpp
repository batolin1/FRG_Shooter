#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <vector>

// #####################################################
// Begin of problem solver parameters.
// #####################################################

/**
    @brief Input parameters for the eigenvalue problem solver, i.e., user inputs.
    @see Eigenvector_Problem_Solver, Eigenvector_Problem_Parameter_Parser
*/
struct Eigenvector_Problem_Solver_Parameter {
    double dimension;
    double anomalous_dimension;
    double symmetry_factor_N;
    double s_factor;
    double sigma;
    double eigenvalue_maxima;
    double eigenvalue_minima;
    int eigenvalue_delta;
};

/**
    @brief Input parameters for the initial condition problem solver, i.e., user inputs.
    @see Initial_Condition_Problem_Solver, Initial_Condition_Problem_Parameter_Parser
*/
struct Initial_Condition_Problem_Solver_Parameter  {
    double dimension;
    double anomalous_dimension;
    double symmetry_factor_N;
    double s_factor_minima;
    double s_factor_maxima;
    int s_factor_delta;
    double sigma_minima;
    double sigma_maxima;
    int sigma_delta;
};

/**
    @brief Input parameters for the shooting problem solver, i.e., user inputs.
    @see Shooting_Problem_Solver, Shooting_Problem_Parameter_Parser
*/
struct Shooting_Problem_Solver_Parameter {
    double dimension;
    double anomalous_dimension;
    double symmetry_factor_N;
    double s_factor;
    double sigma_minima;
    double sigma_maxima;
    int sigma_delta;
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
    double symmetry_factor_N;
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
    @see Eigenvector_Problem_Output_Element, Eigenvector_Problem_Solver
*/
struct Eigenvector_Problem_Solver_Result_Element {
    double eigenvalue;
    double asymptotic_eigenvector;
};

/**
    @brief Result Element for the initial condition problem. 
    @see Initial_Condition_Problem_Output_Element, Initial_Condition_Problem_Solver
*/
struct Initial_Condition_Problem_Solver_Result_Element {
    double spike;
    double s_factor;
    double anomalous_dimension;
};

/**
    @brief Result Element for the shooting problem. 
    @see Shooting_Problem_Output_Element, Shooting_Problem_Solver
*/
struct Shooting_Problem_Solver_Result_Element {
    double asymptotic_field;
    double sigma;
    double anomalous_dimension;
};

// #####################################################
// Begin of problem solver output elements.
// #####################################################

/**
    @brief The output from a eigenvector problem solver calculation. 
    @see Eigenvector_Problem_Solver, Eigenvector_Problem_Result_Formatter
*/
struct Eigenvector_Problem_Solver_Output_Element {
    Eigenvector_Problem_Solver_Parameter parameter;
    double anomalous_dimension;
    std::vector<Eigenvector_Problem_Solver_Result_Element> result;
    std::vector<Potential_Integrator_Result_Element> trajectory;
    std::vector<double> zero_points;
};

/**
    @brief The output from a initial condition problem solver calculation. 
    @see Initial_Condition_Problem_Solver, Initial_Condition_Problem_Result_Formatter
*/
struct Initial_Condition_Problem_Solver_Output_Element {
    Initial_Condition_Problem_Solver_Parameter parameter;
    std::vector<Initial_Condition_Problem_Solver_Result_Element> result;
};

/**
    @brief The output from a shooting problem solver calculation. 
    @see Shooting_Problem_Solver, Shooting_Problem_Result_Formatter
*/
struct Shooting_Problem_Solver_Output_Element {
    Shooting_Problem_Solver_Parameter parameter;
    std::vector<Shooting_Problem_Solver_Result_Element> result;
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
    double tolerance;
    int maximum_iterations;
    bool save_trajectories;
};

// #####################################################
// End of miscellaneous.
// #####################################################

#endif