#ifndef CONCEPTS_HPP
#define CONCEPTS_HPP

#include <concepts>
#include <string>
#include <vector>

/**
    @brief The concept of a formatter to convert output elements to strings.     
    @see result_writer, Solver_Concept
*/
template <typename T, typename Output_Element>
concept Formatter_Concept = requires (const Output_Element& o) {

    { T::format (o) } -> std::same_as<std::string>;
};

/**
    @brief The concept of a parser to convert input elements into parameters for problem solvers.
    @see parameter_reader, Solver_Concept
*/
template <typename T, typename Parameter>
concept Parser_Concept = requires (const std::vector<std::string>& token, Parameter& param) {

    { T::parse_token (token, param) } -> std::same_as<bool>;
};

/**
    @brief The concept of a solver to solve different types of problems. 
    @see parameter_reader, prameter_writer, integrator
*/
template <typename S>
concept Solver_Concept = requires(S solver, typename S::Parameter param) {

    typename S::Parameter;
    typename S::Result_Element;
    typename S::Output_Element;

    { solver.initialize () } -> std::same_as<void>;
    { solver.process_parameter (param) } -> std::same_as<void>;
    { solver.get_result () } -> std::same_as<std::vector<typename S::Output_Element>&>;
};

/**
    @brief The concept of an integrator model for realizing integrations. 
    @see integrator, Solver_Concept
*/
template <typename I>
concept Integrator_Model_Concept = 
requires (
    I i,
    typename I::Parameter p,
    const Configuration& c,
    std::array<double,2> s,
    std::array<double,2> ds,
    double t) {

    typename I::Parameter;

    { i.asymptotic_field } -> std::convertible_to<double>;
    { i.state } -> std::same_as<std::array<double,2>&>;
    { i.initialize (c, p) } -> std::same_as<void>;
    { i.termination_event () } -> std::same_as<bool>;
    { i.on_success_step () } -> std::same_as<void>;
    { i.get_result () } -> std::convertible_to<double>;
    { i.ODE_step (s, ds, t) } -> std::same_as<void>;
};

/**
    @brief The concept of an integrator parameter. 
    @see Integrator_Model_Concept
*/
template<typename P> 
concept Integrator_Parameter_Concept = requires (P p) {
    { p.dimension } -> std::convertible_to<double>;
    { p.anomalous_dimension } -> std::convertible_to<double>;
    { p.symmetry_factor_N } -> std::convertible_to<double>;
    { p.s_factor } -> std::convertible_to<double>;
    { p.sigma } -> std::convertible_to<double>;
    { p.dimension_factor } -> std::convertible_to<double>;
    { p.anomalous_constant } -> std::convertible_to<double>;
    { p.implied_s_factor } -> std::convertible_to<double>;
    { p.s_constant } -> std::convertible_to<double>;
};

#endif