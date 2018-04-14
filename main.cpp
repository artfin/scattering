#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

/* epsilon and sigma are the parameters of Lennard-Jones potential */
static double epsilon = 5.9, sigma=3.57;
/* hbar2m = hbar^2/2m where m = mass of the proton, in meV*A^2 */
static double hbar2m = 2.08;

/* jl - beware! unstable recurrence for large l */

double jl ( int l, double x )
{
  double cos(), sin();              /* Builtin functions */
  double jm1, j0, jp1;
  int lm;

  jm1= std::cos (x) / x;
  j0 = std::sin (x) / x;

  for ( lm = 0; lm < l; lm++ ) {
    jp1= (2*lm+1)/x*j0-jm1;
    jm1= j0;
    j0 = jp1;
  }
  return j0;
}

/* nl - beware! unstable recurrence for large l */

double nl ( int l, double x )
{
  double cos(), sin();              /* Builtin functions */
  double nm1, n0, np1;
  int lm;

  nm1= std::sin (x) / x;
  n0 =-std::cos (x) / x;

  for ( lm = 0; lm < l; lm++ ) {
    np1= (2*lm+1)/x*n0-nm1;
    nm1= n0;
    n0 = np1;
  }
  return n0;
}

int main()
{
    const double xmin = 1.8;
    const double xmax = 18.0;
    const int mesh = 2000;

    const double energy_lower_bound = 0.1;
    const double energy_upper_bound = 3.5;
    const double energy_delta = 0.01;
    const int lmax = 6;

    std::vector<double> r(mesh + 1);
    std::vector<double> u(mesh + 1);
    std::vector<double> p(mesh + 1);
    std::vector<double> f(mesh + 1);
    std::vector<double> vpot(mesh + 1);

    // ??
    int ne = (int) ((energy_upper_bound - energy_lower_bound) / energy_delta) + 1.5;
    std::cout << "ne = " << ne << std::endl;

    double rmin, rmax;
    double lambda2 = M_PI / std::sqrt( (energy_upper_bound + energy_lower_bound)/2.0/hbar2m );

    for ( int n  = 1; n <= ne; ++n )
    {
        double e = energy_lower_bound + (n - 1) * energy_delta;
        double q = std::sqrt(e/hbar2m);

        rmin = xmin;
        rmax = xmax + lambda2;

        dx = (rmax - rmin) / mesh;
        ddx12 = dx*dx / 12.0;

    }
}


#if 0
int main ( )
{
  /* average half-wavelength on the energy range
     recalculating it at each energy produces some numerical noise */
  lambda2 = pi/sqrt( (eup+elw)/2.0/hbar2m );

  for ( n = 1; n <= ne; n++ ) {
    /* initialization */
    e = elw + (n-1)*de;
    q = sqrt(e/hbar2m);

    for ( l = 0; l <= lmax; l++ ) {
      rmin  = xmin;
      cross[l] = 0.0 ;
      /*  also add half a wavelength to allow matching to spherical waves */
      rmax  = xmax + lambda2;
      dx = (rmax-rmin)/mesh;
      ddx12 = dx*dx/12.0;
      r2 = rmax;
      i2 = mesh;
      /* r1 is a second point, half a wavelength from r2=rmax, needed to
     calculate matching to spherical waves */
      r1 = rmax - dx*(int)( lambda2/dx + 0.5);
      i1 = mesh - (int) ( lambda2/dx + 0.5 );

      /* set up the Lennard-Jones potential */

      for ( i = 0; i <= mesh; i++ ) {
    r[i] = rmin + (double)(i) * dx ;
    vpot[i] = epsilon * ( -2.0*pow(sigma/r[i],6) + pow(sigma/r[i],12) ) ;
      }

      /*  set up the f-function used by the Numerov algorithm */

      for ( i = 0; i <= mesh ; i++ ) {
    f[i] = 1.0 - ddx12/hbar2m * ( hbar2m*l*(l+1)/pow(r[i],2) + vpot[i] - e ) ;
    u[i] = 0.0 ;
      }

      starting_points ( r, u );

      /* outward integration */

      for ( i = 1; i < mesh; i++ ) {
    u[i+1] = ((12.0-10.0*f[i])*u[i]-f[i-1]*u[i-1])/f[i+1];
      }

      /*  normalization (simple integral) */

      norm = 0.0 ;
      for ( i = 0; i <= mesh; i++ ) { norm += u[i]*u[i]*dx ;}
      for ( i = 0; i <= mesh; i++ ) { u[i] = u[i]/sqrt(norm);}

      /* matching */

      kappa = r1*u[i2]/(r2*u[i1]) ;
      tandelta = (kappa*jl(l,q*r1)-jl(l,q*r2)) / (kappa*nl(l,q*r1)-nl(l,q*r2));
      delta = atan(tandelta);
      cross[l] = cross[l] + 4*pi/(q*q) * (2*l+1)*pow(sin(delta),2);

      /*  calculation of asymptotic wavefunction p  */

      for ( i = 0; i <= mesh; i++ ) {
    p[i]= sin (q*r[i]-l*pi/2.0 + delta);
      }
      /* normalize p so that it is equal to u for r=r2 */
      for ( i = 0; i <= mesh; i++ ) {
    p[i] = p[i] / (p[i2]/u[i2] );
      }
    }
    /* uncomment to write wavefunction, asymptotic wavefunctions,
       effective potential */

    /* for ( i = 0; i <= mesh; i++ ) {

        effpot = vpot[i] + hbar2m*l*(l+1)/pow( r[i], 2 );
        fprintf (outl[l], "%f %f %f %f\n", r[i], u[i], p[i], effpot);
      }
      fprintf (outl[l], "\n\n", r[i], u[i], p[i], effpot);
    */

    sumcross =0.0;
    for ( l = 0; l <= lmax; l++ ) { sumcross += cross[l] ; }
    fprintf (out, "%f %f", e, sumcross);
    for ( l = 0; l <= lmax; l++ ) { fprintf (out, " %lf", cross[l]); }
    fprintf (out, "\n");
  }
}


/* starting_points - starting wavefunction for Lennard-Jones potential
   assuming dominant r^12 term - not accurate, see Thjissen notes    */

void starting_points ( double* r, double* u)
{
  /*  wave-function in the first two points (r^12 term assumed dominant) */

  u[0] = exp ( -sqrt( epsilon/hbar2m*pow(sigma,12)/25.0 ) * pow(r[0],-5) );
  u[1] = exp ( -sqrt( epsilon/hbar2m*pow(sigma,12)/25.0 ) * pow(r[1],-5) );

  return ;
}
#endif
