#include <Eigen/Dense>

using Eigen::Matrix2d;
using Eigen::Vector2d;

double mew(Matrix2d plastic_deformation_gradient);
double lambda(Matrix2d plastic_deformation_gradient);
std::tuple<Matrix2d, Matrix2d> fewaf(Matrix2d m);
Matrix2d partial_psi_partial_fe(Matrix2d fe, Matrix2d fp);
double weighting_func(double x);
double weighting_func_prime(double x);
Vector2d grad_weighting_func(double grid_ind_x, double grid_ind_y, double x_position, double y_position);
Matrix2d sigma(Matrix2d fe, Matrix2d fp);
Matrix2d adjugate(Matrix2d m);