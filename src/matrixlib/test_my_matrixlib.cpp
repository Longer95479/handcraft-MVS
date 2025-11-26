#include "my_matrix_lib.h"

using namespace Loong;

void test(void)
{
    Matrix<double, 4, 3> A({{1, 2, 3},
                            {4, 5, 6},
                            {7, 8, 9},
                            {10, 11, 12}});
    
    Matrix<double, 4, 3> B;
    B = A;

    Matrix<double, 4, 3> C(B);

    A(1, 1) = 13;

    std::cout << "A:" << std::endl << A << std::endl;
    std::cout << "B:" << std::endl << B << std::endl;
    std::cout << "C:" << std::endl << C << std::endl;
    std::cout << "B.block<2, 1>(1, 1)" << std::endl << B.block<2, 1>(1, 1) << std::endl;

    B.block<2, 1>(1, 1) = Matrix<double, 2, 1>({{954}, {79}});
    std::cout << B << std::endl;

    std::cout << "A + B:" << std::endl << A + B << std::endl;
    std::cout << "A - B:" << std::endl << A - B << std::endl;
    Matrix<double, 3, 3> D;
    D = A.transpose() * A;
    std::cout << "D = A.transpose() * A:" << std::endl << D << std::endl;
    std::cout << "D.row(1): " << std::endl << D.row(1) << std::endl;
    std::cout << "D.col(1): " << std::endl << D.col(1) << std::endl;
    std::cout << "D.col(1).norm(): " << std::endl << D.col(1).norm() << std::endl << std::endl;
    std::cout << "D.col(1) / D.col(1).norm(): " << std::endl << D.col(1) / D.col(1).norm() << std::endl;
            
    auto QR = A.qrDecompose();
    auto Q = QR.first;
    auto R = QR.second;
    auto b = Matrix<double, 4, 1>({{3}, {6}, {9}, {12}});
    auto x = A.solveUsingHouseholderQR(b);

    std::cout << "Q:" << std::endl << Q << std::endl;
    std::cout << "R:" << std::endl << R << std::endl;
    std::cout << "QR:" << std::endl << Q * R << std::endl;
    std::cout << "QTQ:" << std::endl << Q.transpose() * Q << std::endl;

    std::cout << "x:" << std::endl << x << std::endl;
    std::cout << "Ax:" << std::endl << A * x << std::endl;
    std::cout << "b:" << std::endl << b << std::endl;

    Matrix<double, 4, 4> U;
    Matrix<double, 4, 3> C0;
    Matrix<double, 3, 3> V;
    A.svdDecompose(10, &U, &C0, &V);
    std::cout << "U:" << std::endl << U << std::endl;
    std::cout << "C:" << std::endl << C0 << std::endl;
    std::cout << "V:" << std::endl << V << std::endl;
    std::cout << "UCVT:" << std::endl << U * C0 * V.transpose() << std::endl;
    std::cout << "UTU:" << std::endl << U.transpose() * U << std::endl;
    std::cout << "VTV:" << std::endl << V.transpose() * V << std::endl;

}


int main(void)
{
    try {
        test();
    }
    catch (const char* e) {
        std::cout << e << std::endl;
    }

    return 0;
}

