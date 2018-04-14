#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

const double XMAX = 10.0;
const int MESH = 500;

const int NODES = 1;

// количество итераций Нумерова
const int ITERMAX = 1000;
// критерий сходимости shooting method'a Нумерова
const double EPSILON = 1e-10;

int main()
{
    std::vector<double> x(MESH + 1); // ось X
    std::vector<double> y(MESH + 1);
    std::vector<double> p(MESH + 1);
    std::vector<double> f(MESH + 1); // волновая функция
    std::vector<double> vpot(MESH + 1); // потенциал

    double dx = XMAX / MESH;
    double ddx12 = dx * dx / 12.0;

    // setting up X-axis and potential values
    for ( int i = 0; i <= MESH; ++i )
    {
        x[i] = (double) i * dx;
        vpot[i] = 0.5 * x[i] * x[i];
    }

    // setting initial upper/lower bounds for eigenvalue
    double energy_upper_bound = vpot.end()[-1];
    double energy_lower_bound = energy_upper_bound;
    for ( double & e : vpot )
    {
        if ( e > energy_upper_bound ) energy_upper_bound = e;
        if ( e < energy_lower_bound ) energy_lower_bound = e;
    }

    // setting trial energy
    double e = 0.5 * (energy_lower_bound + energy_upper_bound);

    // eigenvalue search
    for ( int iter = 0; (iter < ITERMAX) && (energy_upper_bound - energy_lower_bound > EPSILON); ++iter )
    {
        f[0] = ddx12 * (2.0 * (vpot[0] - e));
        int icl = -1;

        for ( int i = 1; i <= MESH; ++i )
        {

            f[i] = ddx12 * 2.0 * (vpot[i] - e);

            // if f[i] is exactly zero the change of sign is not observed
            // the following line is a trick to prevent missing a change of sign ???
            if ( f[i] == 0.0 )
                f[i] = 1e-20;

            // store the index icl, where the last change of sign has been found
            if ( f[i] != std::copysign(f[i], f[i-1]) )
                icl = i;
        }

        if ( icl >= MESH - 2 )
        {
            std::cout << "Last change of sign too far." << std::endl;
            exit(1);
        }

        for ( int i = 0; i <= MESH; ++i )
            f[i] = 1.0 - f[i];

        std::fill(y.begin(), y.end(), 0.0);

        /*  beware the integer division: 1/2 = 0 !
        if nodes is even, there are 2*hnodes nodes
        if nodes is odd,  there are 2*hnodes+1 nodes (one is in x=0)
        hnodes is thus the number of nodes in the x>0 semi-axis (x=0 excepted) */
        int hnodes = NODES / 2;

        if ( 2 * hnodes == NODES )
        {
            // even number of nodes: wavefunction is even
            y[0] = 1.0;
            // assume f(-1) = f(1)
            y[1] = 0.5 * (12.0 - f[0] * 10.0) * y[0] / f[1];
        }
        else
        {
            // wavefunction is odd
            y[0] = 0.0;
            y[1] = dx;
        }

        // outward integration and count number of crossings ???
        int ncross = 0;
        for ( int i = 1; i <= MESH - 1; ++i )
        {
            y[i + 1] = ((12.0 - f[i] * 10.0) * y[i] - f[i - 1] * y[i - 1]) / f[i + 1];
            if ( y[i] != std::copysign(y[i], y[i + 1]) )
                ++ncross;
        }

        std::cout << std::fixed << std::setprecision(10);
        std::cout << "iter: " << iter << "; ncross: " << ncross << "; e: " << e << std::endl;

        // too many crossings: current energy is too high -> lower the upper bound
        if ( ncross > hnodes )
            energy_upper_bound = e;
        // too few or correct number of crossings: current energy is too low
        else
            energy_lower_bound = e;

        // new trial value:
        e = 0.5 * (energy_lower_bound + energy_upper_bound);
    }
}

