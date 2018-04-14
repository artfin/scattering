#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <cassert>
#include <ctime>
#include <iomanip>

const double XMAX =  1.0;
const int nx = 1e4;
const double h = 2 * XMAX / (nx - 1);
const double h2 = h * h;
const double h12 = h / 12.0;

// границы поиска по энергии
const double elb = 0.01;
const double eub = 80.0;

// поиск точного решения с точностью до
const double EPS = 1e-12;

double potential( const double x )
{
    if ( (x >= -XMAX) && (x <= XMAX) )
        return 0.0;
    else
        return std::numeric_limits<double>::infinity();
}


// нормируем квадрат (!) волновой функции
void normalize( std::vector<double> & psi, const double h )
{
    double norm = psi[0]*psi[0] + psi.back() * psi.back();

    assert( psi.size() % 2 == 0 );
    for ( size_t k = 1; k < psi.size() - 3; k += 2 )
        norm = norm + 4.0 * psi[k]*psi[k] + 2.0 * psi[k + 1]*psi[k + 1];

    norm = 1.0 / std::sqrt(norm * h / 3.0);

    for ( size_t k = 0; k < psi.size(); ++k )
        psi[k] *= norm;
}

// input:
//     std::vector<double> потенциал
// 	   std::vector<double> получающаяся волновая функция
// 	   double значение энергии
void integrate( std::vector<double> const & pot, std::vector<double> & psi, const double ee )
{
    // начальные значения волновой функции
    psi[0] = 0.0;
    psi[1] = 1.0;

    // подготовили значения ф[0], ф[1]
    double f0 = 2.0 * (pot[0] - ee); // правая часть

    double q0 = psi[0] * (1.0 - h12 * f0); // вспомогательная функция phi
    double f1 = 2.0 * (pot[1] - ee);

    double q1 = psi[1] * (1.0 - h12 * f1);

    // ф[k + 1] = 2ф[k] - ф[k - 1] + h^2 * f[k] * psi[k]
    double q2;
    for ( size_t k = 2; k < nx; ++k )
    {
        q2 = 2 * q1 - q0 + h2 * f1 * psi[k - 1];

        q0 = q1;
        q1 = q2;
        f1 = 2.0 * (pot[k] - ee);

        psi[k] = q2 / (1.0 - h12 * f1);
    }

    // нормируем волновую функцию
    normalize(psi, h);
}

// ищем границы по энергии, между которыми находится решение уравнения Ш.
std::pair<double, double> search_for_bracketing_energies( std::vector<double> const & pot, const double lower_bound,
                                     const double upper_bound )
{
    const int ne = 50;
    double de = (upper_bound - lower_bound) / ne;
    double e1, e2;
    std::vector<double> psi1(nx);
    std::vector<double> psi2(nx);

    std::pair<double, double> res;

    for ( e1 = lower_bound, e2 = lower_bound + de; e2 <= upper_bound; e1 += de, e2 += de )
    {
        // находим psi1 с энергией e1
        integrate(pot, psi1, e1);
        // находим psi2 с энергией e2
        integrate(pot, psi2, e2);

        //std::cout << "psi1: " << psi1.back() << std::endl;
        //std::cout << "psi2: " << psi2.back() << std::endl;

        if ( psi1.back() * psi2.back() < 0.0 )
        {
            res.first = e1;
            res.second = e2;
            return res;
        }
    }

    std::cout << "Bracketing energies don't found. Increase energy range." << std::endl;
    exit(1);
}

void bisection( std::vector<double> const & pot, std::vector<double> & res,
                std::pair<double, double> const & enrg )
{
    std::vector<double> psi(nx);

    double e1 = enrg.first;
    integrate(pot, psi, e1);
    double b1 = psi.back();

    double e2 = enrg.second;
    integrate(pot, psi, e2);
    double b2 = psi.back();

    assert(b1 * b2 < 0.0);

    double e0;
    double b0 = b2 - b1;

    std::string filename = "../particle_in_box/wf_iter";

    int iter = 0;
    for ( ; (std::abs(b0) > EPS) && (std::abs(e2 - e1) > EPS); ++iter )
    {
        e0 = 0.5 * (e1 + e2);
        std::cout << "e1: " << e1 << "; e2: " << e2 << "; e0: " << e0 << std::endl;

        integrate(pot, psi, e0);
        b0 = psi.back();

        if ( b0 * b1 <= 0.0 )
            e2 = e0;
        else
            e1 = e0;

        std::ofstream out(filename + std::to_string(iter) + ".txt");
        out << "# e0 = " << e0 << std::endl;
        for ( size_t k = 0; k < psi.size(); ++k )
            out << -XMAX + h * k << " " << psi[k] << std::endl;
        out.close();
    }

    std::cout << std::fixed << std::setprecision(9);
    std::cout << "e0 = " << e0 << std::endl;
    std::copy(psi.begin(), psi.end(), res.begin());
}

int main()
{
    std::vector<double> grid(nx);
    for ( size_t k = 0; k < nx; ++k )
        grid[k] = -XMAX + k * h;

    std::vector<double> pot(nx);
    for ( size_t k = 0; k < nx; ++k )
        pot[k] = potential(grid[k]);

    std::pair<double, double> bracketing_energies = search_for_bracketing_energies(pot, elb, eub);

    /*
    bracketing_energies = search_for_bracketing_energies(pot, bracketing_energies.second, eub);
    std::cout << "enrgs: " << bracketing_energies.first << " " << bracketing_energies.second << std::endl;
    bracketing_energies = search_for_bracketing_energies(pot, bracketing_energies.second, eub);
    std::cout << "enrgs: " << bracketing_energies.first << " " << bracketing_energies.second << std::endl;
    bracketing_energies = search_for_bracketing_energies(pot, bracketing_energies.second, eub);
    std::cout << "enrgs: " << bracketing_energies.first << " " << bracketing_energies.second << std::endl;
    bracketing_energies = search_for_bracketing_energies(pot, bracketing_energies.second, eub);
    std::cout << "enrgs: " << bracketing_energies.first << " " << bracketing_energies.second << std::endl;
    */

    std::clock_t start = std::clock();

    bisection(pot, psi, bracketing_energies);

    std::cout << "Time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    /*
    std::ofstream out("../particle_in_box/wf0.txt");
    for ( size_t k = 0; k < nx; ++k )
        out << grid[k] << " " << psi[k] << std::endl;
    out.close();
    */

    return 0;
}

#if 0
    std::vector<double> psi(nx); // получающийся вектор волновой функции

    double start = 0.61;
    double end = 1.81;
    int n = 100;
    double de = (end - start) / n;

    std::vector<double> enrgs(n);
    std::vector<double> aim(n);

    for ( int k = 0; k < n; ++k )
    {
        double e = start + k * de;
        integrate(pot, psi, e);

        enrgs[k] = e;
        aim[k] = psi.back();
    }

    std::ofstream out("../particle_in_box/bisection_method.txt");
    for ( size_t k = 0; k < enrgs.size(); ++k )
        out << enrgs[k] << " " << aim[k] << std::endl;
    out.close();

#endif
