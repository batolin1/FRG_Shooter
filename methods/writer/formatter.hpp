#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include <string>
#include <sstream>

#include "structure.hpp"



/**
    @brief A formatter for output elements related to the shooting problem solver.
    @see Formatter_Concept
*/
struct Shooting_Problem_Result_Formatter {

    static std::string format (const Shooting_Problem_Solver_Output_Element& output_element) {

        std::ostringstream oss;

        for (const Shooting_Problem_Solver_Result_Element& r : output_element.result) {
            oss << output_element.parameter.dimension << ","
            << output_element.parameter.anomalous_dimension << ","
            << output_element.parameter.s_factor << ","
            << output_element.parameter.symmetry_factor_N << ","
            << r.anomalous_dimension << ","
            << r.sigma << ","
            << r.asymptotic_field
            << "\n";
        }

        return oss.str ();
    }
};



/**
    @brief A formatter for output elements related to the initial condition problem solver.
    @see Formatter_Concept
*/
struct Initial_Condition_Problem_Result_Formatter {

    static std::string format (
        const Initial_Condition_Problem_Solver_Output_Element& output_element) {

        std::ostringstream oss;

        for (const Initial_Condition_Problem_Solver_Result_Element& r : output_element.result) {
            oss << output_element.parameter.dimension << ","
            << output_element.parameter.anomalous_dimension << ","
            << r.anomalous_dimension << ","
            << r.s_factor << ","
            << r.spike
            << "\n";
        }

        return oss.str ();
    }
};



/**
    @brief A formatter for output elements related to the eigenvector problem solver.
    @see Formatter_Concept
*/
struct Eigenvector_Problem_Result_Formatter {

    static std::string format (const Eigenvector_Problem_Solver_Output_Element& output_element) {

        std::ostringstream oss;

        for (const Eigenvector_Problem_Solver_Result_Element& r : output_element.result) {
            oss << output_element.parameter.dimension << ","
            << output_element.parameter.anomalous_dimension << ","
            << output_element.parameter.s_factor << ","
            << output_element.parameter.sigma << ","
            << output_element.anomalous_dimension << ","
            << r.eigenvalue << ","
            << r.asymptotic_eigenvector
            << "\n";
        }

        return oss.str ();
    }
};



/**
    @brief A formatter for trajectory elements related to the eigenvector problem solver.
    @see Formatter_Concept
*/
struct Eigenvector_Problem_Trajectory_Formatter {

    static std::string format (const Eigenvector_Problem_Solver_Output_Element& output_element) {

        std::ostringstream oss;

        for (const Potential_Integrator_Result_Element& t : output_element.trajectory) {
            oss << output_element.parameter.dimension << ","
            << output_element.parameter.anomalous_dimension << ","
            << output_element.parameter.s_factor << ","
            << output_element.parameter.sigma << ","
            << output_element.anomalous_dimension << ","
            << t.field << ","
            << t.potential_0prime << ","
            << t.potential_1prime << ","
            << t.potential_2prime << ","
            << t.denominator << ","
            << t.the_real_denominator
            << "\n";
        }

        return oss.str ();
    }
};

/**
    @brief A formatter for the zeroes related to the eigenvector problem solver.
    @see Formatter_Concept
*/
struct Eigenvector_Problem_Zeroes_Formatter {

    static std::string format (const Eigenvector_Problem_Solver_Output_Element& output_element) {

        std::ostringstream oss;

        for (const double value : output_element.zero_points) {
            oss << output_element.parameter.dimension << ","
            << output_element.parameter.anomalous_dimension << ","
            << output_element.parameter.s_factor << ","
            << output_element.parameter.sigma << ","
            << output_element.anomalous_dimension << ","
            << value
            << "\n";
        }

        return oss.str ();
    }
};

#endif