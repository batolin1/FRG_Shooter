#ifndef GRID_RELAXER_HPP
#define GRID_RELAXER_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "structure.hpp"

namespace grid_relaxer {

    /**
        @brief Simple cubic spline interpolator.
    */
    class CubicSpline {
    private:
        std::vector<double> x, y;
        std::vector<double> m;

    public:
        CubicSpline(const std::vector<double>& x_vals, const std::vector<double>& y_vals) 
            : x(x_vals), y(y_vals) {
            if (x.size() != y.size() || x.size() < 2) {
                throw std::invalid_argument("CubicSpline: x and y must have same size >= 2");
            }
            compute_slopes();
        }

        double evaluate(double x_val) const {
            if (x_val <= x.front()) return y.front();
            if (x_val >= x.back()) return y.back();

            auto it = std::lower_bound(x.begin(), x.end(), x_val);
            int i = std::distance(x.begin(), it);
            if (i == 0) i = 1;
            if (i >= static_cast<int>(x.size())) i = x.size() - 1;

            int i0 = i - 1, i1 = i;
            double dx = x[i1] - x[i0];
            double t = (x_val - x[i0]) / dx;

            double h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
            double h10 = t * (1.0 - t) * (1.0 - t);
            double h01 = t * t * (3.0 - 2.0 * t);
            double h11 = t * t * (t - 1.0);

            return h00 * y[i0] + h10 * dx * m[i0] + h01 * y[i1] + h11 * dx * m[i1];
        }

    private:
        void compute_slopes() {
            int n = x.size();
            m.assign(n, 0.0);
            for (int i = 1; i < n - 1; ++i) {
                double dx0 = x[i] - x[i-1];
                double dx1 = x[i+1] - x[i];
                double dy0 = (y[i] - y[i-1]) / dx0;
                double dy1 = (y[i+1] - y[i]) / dx1;
                m[i] = (dy0 * dx1 + dy1 * dx0) / (dx0 + dx1);
            }
            m[0] = (y[1] - y[0]) / (x[1] - x[0]);
            m[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
        }
    };

std::vector<Potential_Integrator_Result_Element> relax_trajectory(
    const std::vector<Potential_Integrator_Result_Element>& trajectory,
    int N_grid = 500) {

    if (trajectory.empty() || trajectory.size() < 2) {
        return trajectory;
    }

    // =========================================================
    // STEP 1: Extract raw data from trajectory
    // =========================================================
    std::vector<double> field_vals, up_vals, upp_vals;
    for (const auto& elem : trajectory) {
        field_vals.push_back(elem.field);
        up_vals.push_back(elem.potential_1prime);
        upp_vals.push_back(elem.potential_2prime);
    }

    // =========================================================
    // STEP 2: Build uniform grid
    // =========================================================
    const double rho_min = field_vals.front();
    const double rho_max = field_vals.back();
    std::vector<double> phi_grid(N_grid);
    for (int i = 0; i < N_grid; ++i) {
        phi_grid[i] = rho_min + (rho_max - rho_min) * i / (N_grid - 1.0);
    }

    // =========================================================
    // STEP 3: Interpolate U' and U'' onto uniform grid
    // =========================================================
    CubicSpline interp_up (field_vals, up_vals);
    CubicSpline interp_upp(field_vals, upp_vals);

    std::vector<Potential_Integrator_Result_Element> result(N_grid);
    for (int i = 0; i < N_grid; ++i) {
        result[i].field            = phi_grid[i];
        result[i].potential_1prime = interp_up.evaluate (phi_grid[i]);
        result[i].potential_2prime = interp_upp.evaluate(phi_grid[i]);
    }

    // =========================================================
    // STEP 4: Integrate U from U' using trapezoidal rule
    // Preserve the original starting value of U
    // =========================================================
    result[0].potential_0prime = trajectory.front().potential_0prime;
    for (int i = 1; i < N_grid; ++i) {
        result[i].potential_0prime = result[i-1].potential_0prime +
            (result[i].field - result[i-1].field) *
            (result[i].potential_1prime + result[i-1].potential_1prime) / 2.0;
    }

    // =========================================================
    // STEP 5: Recompute shape exponents using central differences
    // Much more stable than U'' * field / U' near zero crossings
    // =========================================================
    const double derivative_threshold = 1e-6;

    // Boundary points — forward/backward difference
    auto compute_shape = [&](int i) -> double {
        double dy;
        if (i == 0) {
            dy = (result[1].potential_1prime - result[0].potential_1prime) /
                 (result[1].field - result[0].field);
        } else if (i == N_grid - 1) {
            dy = (result[N_grid-1].potential_1prime - result[N_grid-2].potential_1prime) /
                 (result[N_grid-1].field - result[N_grid-2].field);
        } else {
            dy = (result[i+1].potential_1prime - result[i-1].potential_1prime) /
                 (result[i+1].field - result[i-1].field);
        }
        const double up = result[i].potential_1prime;
        if (std::abs(up) > derivative_threshold) {
            return std::max(-10.0, std::min(10.0, dy * result[i].field / up));
        }
        return 0.0;
    };

    auto compute_shape_0prime = [&](int i) -> double {
        const double u0p = result[i].potential_0prime;
        const double up  = result[i].potential_1prime;
        if (std::abs(u0p) > derivative_threshold) {
            return std::max(-10.0, std::min(10.0, up * result[i].field / u0p));
        }
        return compute_shape(i) + 1.0;
    };

    for (int i = 0; i < N_grid; ++i) {
        const double field = result[i].field;
        const double up    = result[i].potential_1prime;
        const double upp   = result[i].potential_2prime;

        result[i].denominator          = 1.0 + up + 2.0 * field * upp;
        result[i].the_real_denominator = result[i].denominator;
        result[i].potential_1prime_shape = compute_shape(i);
        result[i].potential_0prime_shape = compute_shape_0prime(i);
        result[i].potential_2prime_shape = compute_shape(i) - 1.0;
    }

    return result;
}
}

#endif