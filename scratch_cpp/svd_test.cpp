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

}