// LHAPDF
#include "LHAPDF/LHAPDF.h"

// APFEL++
#include "APFEL/integrator.h"

using namespace std;
#include <unistd.h>
#define GetCurrentDir getcwd

string GetCurrentWorkingDir()
{
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  string current_working_dir(buff);
  return current_working_dir;
}

int main() {
  // Total energy squared
  const double stot = pow(140.7, 2);
  const double Qmin = 1;
  const double Wmin = sqrt(10);
  const int    sgn  = 1;

  // Include new search path in LHAPDF
  const string path = GetCurrentWorkingDir() + "/../";
  LHAPDF::setPaths(path);

  // Open LHAPDF set with the structure functions
  LHAPDF::PDF* dist = LHAPDF::mkPDF("JAM4EIC");

  // Relevant EW parameters and constants
  const double alpha = 1. / 137.036;
  const double MZ    = 91.1876;
  const double MW    = 80.398;
  const double GF    = 1.1663787e-5;
  const double MN    = 0.93891897;

  // Derived parameters.
  const double MZ2 = MZ * MZ;
  const double MW2 = MW * MW;
  const double MN2 = MN * MN;
  const double GF2 = GF * GF;

  // Conversion factor from natural unists to pb.
  const double conv = 0.3893793721e12;
  
  // Integration accuracy
  const double eps = 1e-7;

  // Integration bounds in ln(Q2)
  const double lnQ2min = 2 * log(Qmin);
  const double lnQ2max = log(stot - MN2);

  // Define integrand in ln(Q2)
  const auto dsigmadlnQ2 = [=] (double const& lnQ2) -> double
    {
      // Helpers
      const double Q2 = exp(lnQ2);
      const double Q  = sqrt(Q2);

      // Integration bounds in ln(x). The 1e-9 term in
      // lnmax is required to avoid running out of
      // boudaries (x > 1).
      const double lnxmin = lnQ2 - log(stot - MN2);
      const double lnxmax = log(1 / ( 1 + pow(Wmin / Q, 2) ));

      // Define integrand in ln(x)
      const auto dsigmadlnxdlnQ2 = [=] (double const& lnx) -> double
      {
	const double x  = exp(lnx);
	const double y  = Q2 / x / ( stot - MN2 );
	const double F2 = dist->xfxQ(908, x, Q);
	const double FL = dist->xfxQ(909, x, Q);
	const double F3 = dist->xfxQ(910, x, Q);
        const double F1 = ( ( 1 + 4 * MN2 / Q2 ) * F2 - FL ) / 2 / x;
        const double factor = conv * 4 * M_PI * pow(alpha, 2) / x / y / Q2 / ( stot - MN2 ) / x;
	// The factor x is the jacobian of the x
	// integration due to the fact that dx = x *
	// d(log(x)).
        return x * factor * ( ( 1 - y - pow(x * y, 2) * MN2 / Q2 ) * F2 + pow(y, 2) * x * F1 +  sgn * ( y - pow(y, 2) / 2 ) * x * F3 );
      };

      // Define integral in ln(x) object and integrate
      // it. The factor Q2 is the jacobian of the Q2
      // integration due to the fact that dQ2 = Q2 *
      // d(log(Q2)).
      const apfel::Integrator Ixq{dsigmadlnxdlnQ2};
      return Q2 * Ixq.integrate(lnxmin, lnxmax, eps);
    };

  // Define integral in ln(Q2) object and integrate it
  const apfel::Integrator Iq{dsigmadlnQ2};
  const double xsec = Iq.integrate(lnQ2min, lnQ2max, eps);

  cout << scientific;
  cout << xsec << endl;

  return 0;
}
