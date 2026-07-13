#ifndef GRID_RELAXER_HPP
#define GRID_RELAXER_HPP

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>

#include "structure.hpp"

namespace grid_relaxer {

struct HermiteStep {
    double x0 = 0.0;
    double x1 = 0.0;
    double h  = 0.0;

    std::array<double,2> y0{};
    std::array<double,2> y1{};

    std::array<double,2> dy0{};
    std::array<double,2> dy1{};

    std::array<double,2> eval(double x) const {
        if (x <= x0) return y0;
        if (x >= x1) return y1;

        double t = (x - x0) / h;

        double h00 = (1 + 2*t) * (1 - t) * (1 - t);
        double h10 = t * (1 - t) * (1 - t);
        double h01 = t * t * (3 - 2*t);
        double h11 = t * t * (t - 1);

        std::array<double,2> out;

        for (int i = 0; i < 2; ++i) {
            out[i] =
                h00 * y0[i] +
                h10 * h * dy0[i] +
                h01 * y1[i] +
                h11 * h * dy1[i];
        }

        return out;
    }
};

// ============================================================
// 2. Internal dense structure
// ============================================================
struct DenseTrajectory {
    std::vector<HermiteStep> steps;

    // Returns the step containing x, or throws if x is out of range.
    const HermiteStep& find_step(double x) const {
        const double eps = 1e-12;

        // Clamp slightly out-of-range values caused by floating-point rounding
        const double x_lo = steps.front().x0;
        const double x_hi = steps.back().x1;

        if (x < x_lo - eps || x > x_hi + eps) {
            throw std::out_of_range(
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

// ============================================================
// 3. Build dense trajectory from input vector
// ============================================================
inline DenseTrajectory build_dense (const std::vector<Potential_Integrator_Result_Element>& traj) {
    if (traj.size () < 2)
        throw std::invalid_argument("Trajectory too small");

    DenseTrajectory dense;
    dense.steps.reserve(traj.size () - 1);

    for (size_t i = 0; i + 1 < traj.size (); ++i) {

        const auto& a = traj [i];
        const auto& b = traj [i + 1];

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

        dense.steps.push_back(s);
    }

    return dense;
}


inline std::vector<Potential_Integrator_Result_Element> relax_trajectory( 
    const std::vector<Potential_Integrator_Result_Element>& trajectory, int N_grid = 500) {
    if (trajectory.size () < 2)
        return trajectory;

    DenseTrajectory dense = build_dense(trajectory);

    const double x_min = trajectory.front().field;
    const double x_max = trajectory.back().field;
    const double dx = (x_max - x_min) / (N_grid - 1);

    std::vector<Potential_Integrator_Result_Element> result;
    result.resize (N_grid);

    // --------------------------------------------------------
    // Step 1: interpolate potential_1prime (U') onto the uniform
    //         grid. This is the only quantity we trust directly
    //         from the Hermite interpolation.
    // --------------------------------------------------------
    for (int i = 0; i < N_grid; ++i) {
        double x = x_min + i * dx;
        auto y = dense.eval(x);

        result [i].field            = x;
        result [i].potential_1prime = y [0];   // U'  — interpolated
        result [i].potential_0prime = 0.0;    // filled in Step 2
        result [i].potential_2prime = 0.0;    // filled in Step 3
    }

    // --------------------------------------------------------
    // Step 2: reconstruct potential_0prime (U) by integrating
    //         potential_1prime using the trapezoid rule.
    //         Anchor the starting value from the original
    //         trajectory at x_min.
    // --------------------------------------------------------
    result [0].potential_0prime = trajectory.front().potential_0prime;

    for (int i = 1; i < N_grid; ++i) {
        double up_prev = result [i-1].potential_1prime;
        double up_curr = result [i  ].potential_1prime;
        result [i].potential_0prime =
            result [i-1].potential_0prime + 0.5 * dx * (up_prev + up_curr);
    }

    // --------------------------------------------------------
    // Step 3: reconstruct potential_2prime (U'') from the
    //         uniform-grid U' values using central differences,
    //         with first-order one-sided differences at the
    //         boundaries (avoids the flat-step artefact).
    // --------------------------------------------------------

    // Interior: central difference
    for (int i = 1; i < N_grid - 1; ++i) {
        result [i].potential_2prime =
            (result [i+1].potential_1prime -
             result [i-1].potential_1prime) / (2.0 * dx);
    }

    // Boundaries: first-order one-sided differences
    result [0].potential_2prime =
        (result [1].potential_1prime -
         result [0].potential_1prime) / dx;

    result [N_grid-1].potential_2prime =
        (result [N_grid-1].potential_1prime -
         result [N_grid-2].potential_1prime) / dx;

    // --------------------------------------------------------
    // Step 4: shape diagnostics (numerically stable)
    // --------------------------------------------------------
    const double eps = 1e-8;

    for (int i = 0; i < N_grid; ++i) {

        double x   = result [i].field;
        double u   = result [i].potential_0prime;
        double up  = result [i].potential_1prime;
        double upp = result [i].potential_2prime;

        result [i].denominator =
            1.0 + up + 2.0 * x * upp;

        result [i].potential_1prime_shape =
            (std::abs (up) > eps) ? (x * upp / up) : 0.0;

        result [i].potential_0prime_shape =
            (std::abs (u) > eps) ? (up * x / u) : 0.0;

        result [i].potential_2prime_shape =
            result [i].potential_1prime_shape - 1.0;
    }

    return result;
}

} // namespace grid_relaxer

#endif