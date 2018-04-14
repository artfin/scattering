#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <complex>

#include <Eigen/Dense>

const double V0 = 8e-4;
const double E = 4.1e-4;
const double a = 3.0;
const double m = 1822.0;

const double k = std::sqrt(2 * m * E);

void fill_ML_matrix( Eigen::Matrix2cd & ML )
{
    ML(0, 0) = 1;
    ML(0, 1) = 1;
    ML(1, 0) = std::complex<double>(0, 1) * k;
    ML(1, 1) = - std::complex<double>(0, 1) * k;
}

void fill_MR_matrix( Eigen::Matrix2cd & MR )
{
    if ( E > V0 )
    {
        const double kprime = std::sqrt(2 * m * (E - V0));
        MR(0, 0) = 1.0;
        MR(0, 1) = 1.0;
        MR(1, 0) = std::complex<double>(0, 1) * kprime;
        MR(0, 1) = - std::complex<double>(0, 1) * kprime;
    }
    else
    {
        const double kappa = std::sqrt(2 * m * (V0 - E));
        MR(0, 0) = 1.0;
        MR(0, 1) = 1.0;
        MR(1, 0) = kappa;
        MR(0, 1) = -kappa;
    }

}

void fill_Mtilde_matrix( Eigen::Matrix2cd & Mtilde )
{
    if ( E > V0 )
    {
        const double kprime = std::sqrt(2 * m * (E - V0));
        Mtilde(0, 0) = std::exp( std::complex<double>(0, 1) * kprime * a);
        Mtilde(0, 1) = std::exp(-std::complex<double>(0, 1) * kprime * a);
        Mtilde(1, 0) = std::complex<double>(0, 1) * kprime * \
                std::exp( std::complex<double>(0, 1) * kprime * a);
        Mtilde(1, 1) = -std::complex<double>(0, 1) * kprime * \
                std::exp(-std::complex<double>(0, 1) * kprime * a);
    }
    else
    {
        const double kappa = std::sqrt(2 * m * (V0 - E));
        Mtilde(0, 0) = std::exp( kappa * a );
        Mtilde(0, 1) = std::exp(-kappa * a );
        Mtilde(1, 0) = kappa * std::exp( kappa * a );
        Mtilde(1, 1) =-kappa * std::exp(-kappa * a );
    }
}

Eigen::Vector2cd find_s_r( Eigen::Matrix2cd const & M )
{
    Eigen::Matrix2cd A;
    A(0, 0) = -std::exp( std::complex<double>(0, 1) * k * a);
    A(0, 1) = M(0, 1);
    A(1, 0) = -std::complex<double>(0, 1) * k * std::exp( std::complex<double>(0, 1) * k * a );
    A(1, 1) = M(1, 1);

    Eigen::Vector2cd b;
    b(0) = -M(0, 0);
    b(1) = -M(1, 0);

    Eigen::Vector2cd sol = A.colPivHouseholderQr().solve(b);
    return sol;
}

int main()
{

    Eigen::Matrix2cd ML; // комплексная матрица 2x2
    fill_ML_matrix( ML );

    Eigen::Matrix2cd MR;
    fill_MR_matrix( MR );

    Eigen::Matrix2cd Mtilde;
    fill_Mtilde_matrix( Mtilde );

    Eigen::Matrix2cd M = Mtilde * MR.inverse() * ML;
    std::cout << M << std::endl;

    Eigen::Vector2cd sol = find_s_r( M );
    std::cout << "|S|^2: " << std::abs(sol(0)) << std::endl;
    std::cout << "|R|^2: " << std::abs(sol(1)) << std::endl;

    return 0;
}

