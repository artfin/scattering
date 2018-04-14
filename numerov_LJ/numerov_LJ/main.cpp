#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <algorithm>

const double lj_eps = 5.9; // meV
const double lj_sigma = 3.57; // A
// returns potential in meV
double potential( const double r )
{
    double ratio = lj_sigma / r;
    return lj_eps * ( std::pow(ratio, 12) - 2 * std::pow(ratio, 6) );
}

// hbar^2 / 2m для H-Kr в единицах meV^-1 * sigma^-2
const double hbar2m = 2.08; // ??

double F( const int l, const double r, const double E )
{
    return potential(r) + hbar2m * l * (l + 1) / r / r - E;
}

// assumed that psi already contains needed [0] and [1] values
void integrate( std::vector<double> & psi, const double E, const int l, const double rmin, const double h )
{
    assert( psi.size() >= 2 );

    const double h2 = h * h;
    const double hh12 = h2 / 12;

    std::vector<double> grid(psi.size());
    for ( size_t k = 0; k < grid.size(); ++k )
        grid[k] = rmin + k * h;

    std::vector<double> w(psi.size());
    w[0] = 1.0 - hh12 / hbar2m * F(l, grid[0], E);
    w[1] = 1.0 - hh12 / hbar2m * F(l, grid[1], E);

    for ( size_t k = 2; k < w.size(); ++k )
    {
        w[k] = 1.0 - hh12 / hbar2m * F(l, grid[k], E);
        psi[k] = ((12.0 - 10.0 * w[k - 1]) * psi[k - 1] - w[k - 2] * psi[k - 2]) / w[k];
    }
}

// spherical Bessel function jl [Stegun, Abramowitz, p.457]
double jl(const int l, const double x)
{
    std::vector<double> v(l + 1);

    v[0] = std::sin(x) / x;
    if ( l == 0 ) return v[0];

    v[1] = std::sin(x) / x /x - std::cos(x) / x;
    if ( l == 1 ) return v[1];

    for ( int k = 2; k <= l; ++k )
        v[k] = (2 * k - 1) / x * v[k - 1] - v[k - 2];

    return v[l];
}

// spherical Bessel function nl [...]
double nl(const int l, const double x)
{
    std::vector<double> v(l + 1);

    v[0] = - std::cos(x) / x;
    if ( l == 0 ) return v[0];

    v[1] = -std::cos(x) / x / x - std::sin(x) / x;
    if ( l == 1 ) return v[1];

    for ( int k = 2; k <= l; ++k )
        v[k] = (2 * k - 1) / x * v[k - 1] - v[k - 2];

    return v[l];
}


int main()
{
    const double xmin = 1.8;
    const double xmax = 18.0;
    const int mesh = 3000;

    std::vector<double> psi(mesh);

    double elb = 0.3; // energy lower bound
    double eup = 3.5; // energy upper bound
    double ne = 100; // количество разных энергетических точек
    double de = (eup - elb) / ne; // delta e

    int lmax = 6;

    double h;
    double rmax, r1, r2;
    int count1, count2;

    // полное сечение рассеяния
    std::vector<double> total_cross(ne);
    // сечение рассеяния при фиксированной энергии, разных l
    std::vector<double> cross_l(lmax + 1);

    // энергии рассеяния
    std::vector<double> energies(ne);

    std::clock_t start = std::clock();

    std::string filename = "../LJ/cross";

    // итерируемся по энергиям
    for ( int k = 0; k < ne; ++k )
    {
        // текущая энергия
        double e = elb + de * k;
        //std::cout << "Energy: " << e << std::endl;
        energies[k] = e;

        double q = std::sqrt(e / hbar2m);

        // половина длины волны де-бройля
        double lmbd2 = M_PI / std::sqrt(e) * std::sqrt(hbar2m);
        //std::cout << "De-broigle wavelength: " << lmbd2 << std::endl;

        // обнуляем сечения рассеяния при разных l
        std::fill(cross_l.begin(), cross_l.end(), 0.0);

        // итерируемся по значениям углового момента
        for ( int l = 0; l <= lmax; ++l ) // !!!!
        {
            // перед каждой итерацией обнуляем волновую функцию
            std::fill( psi.begin(), psi.end(), 0.0 );

            rmax = xmax + lmbd2;
            h = (rmax - xmin) / mesh;

            r2 = rmax;
            // номер r2 на сетке
            count2 = mesh - 1;

            r1 = rmax - h * (int)(lmbd2 / h); // dx * (int) (lmbd2/h) даст нам отступ на целое количество единиц h
            // поэтому r1 также лежит на нашей сетке
            // номер r1 на сетке
            count1 = mesh - (int)(lmbd2 / h);

            psi[0] = std::exp ( -std::sqrt( lj_eps / hbar2m * std::pow(lj_sigma, 12) / 25.0 ) * std::pow(xmin, -5.0) );
            psi[1] = std::exp ( -std::sqrt( lj_eps / hbar2m * std::pow(lj_sigma, 12) / 25.0 ) * std::pow(xmin + h, -5.0) );

            integrate(psi, e, l, xmin, h);

            double integral = 0.0;
            for ( size_t kk = 0; kk < psi.size(); ++kk )
                integral += psi[kk] * psi[kk] * h;

            for ( size_t kk = 0; kk < psi.size(); ++kk )
                psi[kk] /= integral;

            /*
            std::ofstream outWF("../LJ/wf" + std::to_string(e) + ".txt");
            for ( size_t kk = 0; kk < psi.size(); ++kk )
                outWF << xmin + kk * h << " " << psi[kk] / (xmin + kk * h) << std::endl;
            outWF.close();
            */

            double kappa = r1 * psi[count2] / (r2 * psi[count1]);
            //std::cout << "r1 = " << r1 << "; r2 = " << r2 << std::endl;
            //std::cout << "psi[count2] = " << psi[count2] << "; psi[count1] = " << psi[count1] << std::endl;
            //std::cout << "kappa: " << kappa << std::endl;

            double tan_delta = (kappa * jl(l, q*r1) - jl(l, q*r2)) / (kappa * nl(l, q*r1) - nl(l, q*r2));
            double delta = std::atan(tan_delta);

            cross_l[l] = 4.0 * M_PI / (q * q) * (2 * l + 1) * std::sin(delta) * std::sin(delta);

            std::ofstream out;
            out.open(filename + std::to_string(l) + ".txt", std::ofstream::out | std::ofstream::app );
            out << energies[k] << " " << cross_l[l] << std::endl;
            out.close();
        }

        total_cross[k] = std::accumulate(cross_l.begin(), cross_l.end(), 0.0);
    }

    std::ofstream out("../LJ/total_crossection.txt");
    for ( size_t k = 0; k < total_cross.size(); ++k )
        out << energies[k] << " " << total_cross[k] << std::endl;

    std::cout << "---------------------" << std::endl;
    std::cout << "Time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    std::string effpot_filename = "../effective_potentials/";
    for ( int l = 0; l < 10; ++l )
    {
        std::string out_ = effpot_filename + std::to_string(l) + ".txt";
        std::ofstream out(out_);
        for ( double x = 1.8; x < 18.0; x += 0.05 )
            out << x << " " << F(l, x, 0.0 ) << std::endl;
        out.close();
    }

    return 0;
}
