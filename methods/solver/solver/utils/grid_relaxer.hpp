#ifndef GRID_RELAXER_HPP
#define GRID_RELAXER_HPP

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>

#include "structure.hpp"

/**
    @brief For the purposes of the eigenvalue problem, after the trajectory for the potential is 
           obtained, a "grid relaxation" process must be executed to make the steps even. This class 
           is largely a copy of whatever Python's ODE solver does for grid relaxation.
*/
namespace grid_relaxer {

    /**
        @brief A structure to record a single step. 
    */
    struct HermiteStep {
        double x0 = 0.0;
        double x1 = 0.0;
        double h  = 0.0;

        std::array<double,2> y0 {};
        std::array<double,2> y1 {};

        std::array<double,2> dy0 {};
        std::array<double,2> dy1 {};

        /**
            @brief This method evaluates the trajectory at a particular value, and returns an array with 
                the respective output values. 
            @param x    The value at which the trajectory is evaluated. 
            @return     The output array. 
                
        */
        std::array<double,2> eval (double x) const {
            if (x<=x0) return y0;
            if (x>=x1) return y1;

            double t = (x - x0) / h;

            double h00 = (1 + 2*t) * (1 - t) * (1 - t);
            double h10 = t * (1 - t) * (1 - t);
            double h01 = t * t * (3 - 2*t);
            double h11 = t * t * (t - 1);

            std::array<double,2> out;

            for (int i = 0; i < 2; ++i) {
                out [i] =
                    h00 * y0 [i] +
                    h10 * h * dy0 [i] +
                    h01 * y1 [i] +
                    h11 * h * dy1 [i];
            }

            return out;
        }
    };

    /**
        @brief A structure to represents the trajectory. 
    */
    struct DenseTrajectory {
        std::vector<HermiteStep> steps;

        /**
            @brief Given input, returns the hermite step corresponding to the particular input. 
            @param x    The x-value in the trajectory. 
            @return     The respective hermite step. 
        */
        const HermiteStep& find_step (double x) const {
            const double eps = 1e-12;

            // Clamp slightly out-of-range values caused by floating-point rounding
            const double x_lo = steps.front ().x0;
            const double x_hi = steps.back ().x1;

            if (x < x_lo - eps || x > x_hi + eps) {
                throw std::out_of_range (
                    "DenseTrajectory::find_step: x out of trajectory range");
            }

            for (const auto& s : steps) {
                if (x >= s.x0 && x <= s.x1 + eps)
                    return s;
            }

            // Should never reach here after the range check above
            return steps.back ();
        }

        std::array<double,2> eval (double x) const {
            return find_step (x).eval (x);
        }
    };

    /**
        @brief Given the input trajectory, builds the dense trajectory with the hermite steps. 
        @param traj    The input trajectory. 
        @return        The dense trajectory. 
    */
    inline DenseTrajectory build_dense (
        const std::vector<Potential_Integrator_Result_Element>& traj) {

        if (traj.size () < 2) {
            throw std::invalid_argument("Trajectory too small");
        }

        DenseTrajectory dense;
        dense.steps.reserve(traj.size () - 1);

        for (int i = 0; i+1 < traj.size (); ++i) {

            const auto& a = traj [i];
            const auto& b = traj [i+1];

            HermiteStep s;
            s.x0 = a.field;
            s.x1 = b.field;
            s.h  = s.x1 - s.x0;

            // state mapping:
            // [0] = U'   (potential_1prime)
            // [1] = U    (potential_0prime)
            s.y0 = { a.potential_1prime, a.potential_0prime };
            s.y1 = { b.potential_1prime, b.potential_0prime };

            // derivatives of the above w.r.t. field:
            // d(U')/dx = U''
            // d(U)/dx  = U'
            s.dy0 = { a.potential_2prime, a.potential_1prime };
            s.dy1 = { b.potential_2prime, b.potential_1prime };

            dense.steps.push_back (s);
        }

        return dense;
    }


    /**
        @brief The actual relaxer for the trajecotry, takes as input a trajectory, and relaxes it 
               onto a homogenously-spaced grif, given inputs.
        @param trajectory    The trajectory to be relaxed. 
        @param N_grid        The size of the relaxed trajectory. 
        @return              Returns the relaxed trajectory. 
    */
    inline std::vector<Potential_Integrator_Result_Element> relax_trajectory ( 
        const std::vector<Potential_Integrator_Result_Element>& trajectory, int N_grid = 500) {

        if (trajectory.size () < 2)
            return trajectory;

        DenseTrajectory dense = build_dense (trajectory);

        const double x_min = trajectory.front ().field;
        const double x_max = trajectory.back ().field;
        const double dx = (x_max - x_min) / (N_grid - 1);

        std::vector<Potential_Integrator_Result_Element> result;
        result.resize (N_grid);

        // Interpolates U' onto uniform grid.
        for (int i = 0; i < N_grid; ++i) {

            double x = x_min + i * dx;
            auto y = dense.eval(x);

            result [i].field = x;
            result [i].potential_1prime = y [0];   
            result [i].potential_0prime = 0.0;    
            result [i].potential_2prime = 0.0;    
        }

        // Reconstructs U by integrating U' with trapezoid rule.
        result [0].potential_0prime = trajectory.front().potential_0prime;

        for (int i = 1; i < N_grid; ++i) {
            double up_prev = result [i-1].potential_1prime;
            double up_curr = result [i].potential_1prime;
            result [i].potential_0prime =
                result [i-1].potential_0prime + 0.5 * dx * (up_prev + up_curr);
        }

        // Reconstructs U'' using central differences. 

        // Interior: central difference
        for (int i = 1; i < N_grid - 1; ++i) {
            result [i].potential_2prime =
                (result [i+1].potential_1prime - result [i-1].potential_1prime) / (2.0 * dx);
        }

        // Boundaries: first-order one-sided differences
        result [0].potential_2prime =
            (result [1].potential_1prime - result [0].potential_1prime) / dx;

        result [N_grid-1].potential_2prime =
            (result [N_grid-1].potential_1prime - result [N_grid-2].potential_1prime) / dx;

        // Computes other parameters. Use epsilon to avoid possible division by zero. 
        const double eps = 1e-8;

        for (int i = 0; i < N_grid; ++i) {

            double x = result [i].field;
            double u = result [i].potential_0prime;
            double up = result [i].potential_1prime;
            double upp = result [i].potential_2prime;

            result [i].denominator = 1.0 + up + 2.0 * x * upp;
            result [i].potential_1prime_shape = (std::abs (up) > eps) ? (x * upp / up) : 0.0;
            result [i].potential_0prime_shape = (std::abs (u) > eps) ? (up * x / u) : 0.0;
            result [i].potential_2prime_shape = result [i].potential_1prime_shape - 1.0;
        }

        return result;
    }

}

#endif