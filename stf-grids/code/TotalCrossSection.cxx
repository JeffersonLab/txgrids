/*
 * Code to compute the BGR18 total UHE neutrino-nucleon cross-sections
 * on isoscalar target at NNLO+NLLx in the nf = 6 FONLL scheme using
 * the DIS structure functions tabulated in
 * "NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF".
 *
 * If you use this code in a scientific publication, please cite:
 *
 * "Neutrino Telescopes as QCD Microscopes" Valerio Bertone, Rhorry
 * Gauld, Juan Rojo arXiv:1808.aaaaa
 */

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// APFEL++
#include "APFEL/integrator.h"

using namespace std;

enum process  { NC, CC };
enum particle { NEUTRINO = 1, ANTINEUTRINO = - 1};

int main() {
  // Cross section switches
  const vector<process>  vproc{CC, NC};
  const vector<particle> vnsgn{NEUTRINO, ANTINEUTRINO};

  // Vector of neutrino energies
  const vector<double> EnuTab{5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4,
      1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8,
      5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10, 1e11, 2e11, 5e11,
      1e12, 2e12, 5e12};

  // Open LHAPDF set with the structure functions (isoscalar target)
  vector<LHAPDF::PDF*> dist = LHAPDF::mkPDFs("NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF");

  // Retrieve parameters from the LHAPDF set
  const int    Nrep = dist.size() - 1;
  const double Ql   = dist[0]->qMin();
  const double Qu   = dist[0]->qMax();
  const double xl   = dist[0]->xMin();
  const double xu   = dist[0]->xMax();

  // Relevant EW parameters and constants
  const double MZ   = 91.1876;
  const double MW   = 80.385;
  const double GF   = 1.1663787e-5;
  const double MN   = 0.938272046;

  // Derived parameters.See Eqs. (8.22) and (8.25) of
  // https://arxiv.org/pdf/0709.1075.pdfhttps://arxiv.org/pdf/0709.1075.pdf
  const double MZ2  = MZ * MZ;
  const double MW2  = MW * MW;
  const double MN2  = MN * MN;
  const double GF2  = GF * GF;
  const double coef = GF * pow(dist[0]->quarkMass(6), 2) / 8 / sqrt(2) / pow(M_PI, 2);
  const double rho  = 1 + 3 * coef * ( 1 + coef * ( 19 - 2 * pow(M_PI, 2) ) );

  // Conversion factor from natural unists to pb.
  const double conv = 0.3894e9;
  
  // Integration accuracy
  const double eps = 1e-3;

  // Loop over processes
  for(auto const& proc : vproc)
    // Loop over neutrino and antineutrino
    for(auto const& nsgn : vnsgn)
      {
	// offset: 1000 = NC, 2000 = CC neutrino, 3000 = CC anti-neutrino
	int offset;
	if(proc == 0)
	  {
	    cout << "\nNeutral-current ";
	    offset = 1000;
	    if(nsgn == 1)
	      cout << "neutrino total cross sections" << endl;
	    else
	      cout << "antineutrino total cross sections" << endl;
	  }
	else
	  {
	    cout << "\nCharged-current ";
	    if(nsgn == 1)
	      {
		offset = 2000;
		cout << "neutrino total cross sections" << endl;
	      }
	    else
	      {
		cout << "antineutrino total cross sections" << endl;
		offset = 3000;
	      }
	  }

	cout << scientific;
	cout << "\n  Enu [GeV]  "
	     << "     cross section [pb]     "
	     << "uncertainty[%]"
	     << endl;
	cout << "             "
	     << "    member 0  "
	     << "    average   "
	     << "              "
	     << endl;
	// Loop over the neutrino energies
	for (auto const& Enu : EnuTab)
	  {
	    int counter = 0;
	    double sum  = 0;
	    double sum2 = 0;
	    double mem0;

	    // Loop over set members
	    for (auto const& mem : dist)
	      {
		// Total energy squared
		const double stot = MN2 + 2 * MN * Enu;

		// Integration bounds in ln(Q2)
		const double lnQ2min = 2 * log(Ql);
		const double lnQ2max = min(log(stot - MN2), 2 * log(500 * MW));

		// Define integrand in ln(Q2)
		const auto dsigmadlnQ2 = [=] (double const& lnQ2) -> double
		  {
		    // Integration bounds in ln(x). The 1e-9 term in
		    // lnmax is required to avoid running out of
		    // boudaries (x > 1).
		    const double lnxmin = lnQ2 - log(stot - MN2);
		    const double lnxmax = 0 - 1e-9;

		    // Helpers
		    const double Q2 = exp(lnQ2);
		    const double Q  = sqrt(Q2);

		    // Define integrand in ln(x)
		    const auto dsigmadlnxdlnQ2 = [=] (double const& lnx) -> double
		    {
		      const double x      = max(exp(lnx), xl);
		      const double y      = Q2 / x / ( stot - MN2 );
		      const double omy2   = pow(1 - y, 2);
		      const double Yplus  = 1 + omy2;
		      const double Yminus = 1 - omy2;
		      double fact;

		      if(proc == 0)
			fact = conv * 8 * pow( ( MZ2 * rho  - MW2) * MW2 / MZ2 / Q2 / rho / rho / ( 2 - rho ), 2) * GF2 / M_PI / x;
		      else
			fact = conv * GF2 / 4 / M_PI / x * pow(MW2 / ( MW2 + Q2 ), 2);

		      // The factor x is the jacobian of the x
		      // integration due to the fact that dx = x *
		      // d(log(x)).
		      return x * fact * ( Yplus           * mem->xfxQ(offset + 1, x, Q)
					  - y * y         * mem->xfxQ(offset + 2, x, Q)
					  + nsgn * Yminus * mem->xfxQ(offset + 3, x, Q) );
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

		if(counter == 0)
		  mem0 = xsec;
		else
		  {
		    sum  += xsec;
		    sum2 += pow(xsec, 2);
		  }
		counter++;
	      }
	    const double mean = sum / Nrep;
	    const double unc  = sqrt( sum2 / Nrep - mean * mean );
	    cout << Enu << "  " << mem0 << "  " << mean << "  " << 100 * unc / mean << endl;
	  }
	cout << "\n";
      }

  return 0;
}
