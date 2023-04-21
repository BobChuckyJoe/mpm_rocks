#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void test_mat(Matrix3d mat) {
    JacobiSVD<Matrix3d> svd(mat, ComputeFullU | ComputeFullV);
    std::cout << svd.singularValues() << std::endl;
    std::cout << "matrix u: " << std::endl;
    std::cout << svd.matrixU() << std::endl;
    std::cout << "matrix v: " << std::endl;
    std::cout << svd.matrixV() << std::endl;

    std::cout << "U * S * V^T: " << std::endl;
    Matrix3d usv_t = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
    std::cout << svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose() << std::endl;
    std::cout << "diff" << std::endl;
    Matrix3d diff = mat - usv_t;
    std::cout << diff << std::endl;
    std::cout << "Norm" << std::endl;
    std::cout << diff.norm() << std::endl;
}

int main() {
    Matrix3d a {
        {7271926.302062575, -45083290.57492536, -20809342.52386164},
        {18982634.79295168, -117685413.71732813, -54320704.96914578},
        {4833102.3873939235, -29963471.74739077, -13830405.2539423}
    };
    Matrix3d b {
        { 1095564031748140.5, -2003436232048194.3, -1533397669604290.8},
        { 91662846290414.47, -167617553868139.3, -128295268444698.08},
        { 12762433838222480, -23338473423876540, -17862920632250384}
    };
    
    test_mat(b);

    Matrix3d r {
        { 0.9999965895157381, 0.000000004624588069912158,     -0.0026116961716095653 },
        { 0.000027607009689941742, 0.9999441116803974, 0.010572263407384649 },
        { 0.00261155025719159,      -0.010572299451968872,         0.9999407013866136}
    };
    Matrix3d u {
        {1.0001183247237462, -0.000016901223321090454,    0.0019267415428965975},
    {-0.000016901223321034943,       0.9999853959472098,    -0.005935553860422698},
    {0.0019267415428965975,   -0.0059355538604227015,       0.7884486284173051}, 
    };
    Matrix3d fe {
        { 1.0001098126155128, -0.0000011814514224149678,   -0.00012484735462328506},
        {0.00003101823900741668,        0.9998669462130203,     0.0024073120335517126},
        {0.004538666818820031,     -0.016507395941659366,        0.7884694854029726},
    };

    Matrix3d reconstructed = r * u;
    Matrix3d diff = fe - reconstructed;
    std::cout << "diff" << std::endl;
    std::cout << diff << std::endl;
    std::cout << "Norm" << std::endl;
    std::cout << diff.norm() << std::endl;

    Matrix3d large_dg {
        { 27775335.729523204, 4928793.958431432, -4582412.15875561},
        {-1729686.4701394828, -306935.89343111374, 285366.6181921465},
        {1713430.6180917455, 304052.0468986154, -282684.0463539241},
    };
    std::cout << "large_dg" << std::endl;
    std::cout << large_dg.determinant() << std::endl;
}