#include <cmath>

#include <config.hpp>
#include <equations.hpp>
double mew(Matrix2d plastic_deformation_gradient) {
    return MEW_NOUGHT * exp((HARDENING_COEFFICIENT * ( 1.0 - plastic_deformation_gradient.determinant())));
}
double lambda(Matrix2d plastic_deformation_gradient) {
    return LAMBDA_NOUGHT * exp((HARDENING_COEFFICIENT * ( 1.0 - plastic_deformation_gradient.determinant())));
}
std::tuple<Matrix2d, Matrix2d> polar_ru_decomposition(Matrix2d m) {
    auto svd = m.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto u = svd.matrixU();
    auto s = svd.singularValues();
    Matrix2d singular_matrix{{s[0], 0}, {0, s[1]}};
    auto v = svd.matrixV();
    return { u * v.transpose(), v * singular_matrix * v.transpose()};
}

Matrix2d partial_psi_partial_fe(Matrix2d fe, Matrix2d fp) {
    auto ru = polar_ru_decomposition(fe);
    auto r = std::get<0>(ru);
    auto u = std::get<1>(ru);
    auto adjugate_transpose = adjugate(fe).transpose();
    return 2.0 * mew(fp) * (fe.array() - r.array()) + lambda(fp) * (fe.determinant() - 1.0) * adjugate_transpose.array();
}
double weighting_func(double x) {
    if (0.0 <= abs(x) && abs(x) < 1.0) {
        return 0.5 * pow(abs(x), 3.0) - pow(x, 2.0) + 2.0/3.0;
    }
    else if (1.0 <= abs(x) && abs(x) < 2.0) {
        -1.0/6.0 * pow(abs(x), 3.0) + pow(x, 2) - 2.0 * abs(x) + 4.0/3.0;
    }
    else {
        0.0;
    }
}

double weighting_func_prime(double x) {
    if (0.0 <= abs(x) && abs(x) < 1.0) {
        return 1.5 * pow(x,2.0) - 2.0 * x;
    }
    else if (1.0 <= abs(x) && abs(x) < 2.0) {
        return -0.5 * pow(x, 2.0) + 2.0 * x - 2.0;
    }
    else {
        return 0;
    }
}
Vector2d grad_weighting_func(double grid_ind_x, double grid_ind_y, double x_position, double y_position) {
    double x = weighting_func_prime(1.0/ GRID_SPACING * (x_position - grid_ind_x * GRID_SPACING)) * 1.0 / GRID_SPACING *
        weighting_func(1.0 / GRID_SPACING * (y_position - grid_ind_y * GRID_SPACING));
    double y = weighting_func_prime(1.0/GRID_SPACING * (y_position - grid_ind_y * GRID_SPACING)) * 1.0 / GRID_SPACING *
        weighting_func(1.0 / GRID_SPACING * (x_position - grid_ind_x * GRID_SPACING));
    
    return Vector2d{x, y};
}
Matrix2d sigma(Matrix2d fe, Matrix2d fp) {
    auto deformation_mat = fe * fp;
    return 1.0 / deformation_mat.determinant() * partial_psi_partial_fe(fe, fp) * fe.transpose();
}

Matrix2d adjugate(Matrix2d m) {
    return Matrix2d{
        {m(1,1), -m(0,1)},
        {-m(1,0), m(0,0)}
    };
}