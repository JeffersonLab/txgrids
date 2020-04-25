/*
  Authors: Valerio Bertone and Rabah Abdul Khalek
 */

#include "LHAPDF/LHAPDF.h"
#include "apfel/APFEL.h"
#include <sstream>

// Open LHAPDF set.
const std::string set = "NNPDF31_nnlo_pch_as_0118";
const std::vector<LHAPDF::PDF*> dist = LHAPDF::mkPDFs(set);

// Value of rs
const double rs = 0.5;

template <typename T>
std::string to_string_with_precision(T const& a_value, int const& n = 1)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

extern "C" void externalsetapfelrep_(double* x, double* Q, int* irep, double* xf);

int main()
{
  // Retrieve evolution parameters from the LHAPDF grid to avoid
  // mismatches.
  const int    pto   = dist[0]->orderQCD();
  const double Qref  = 91.2;
  const double asref = dist[0]->alphasQ(Qref);
  const double mc    = dist[0]->quarkThreshold(4);
  const double mb    = dist[0]->quarkThreshold(5);
  const double mt    = dist[0]->quarkThreshold(6);
  const double Qmin  = dist[0]->qMin();
  const double Qmax  = dist[0]->qMax();
  const double xmin  = dist[0]->xMin();
  const int    nrep  = dist.size();

  // Name of the output set
  const std::string outname = set + "_rs_" + to_string_with_precision(rs);

  // Initialize APFEL and write the LHAPDF grid
  APFEL::EnableWelcomeMessage(false);
  APFEL::LockGrids(true);
  APFEL::SetAlphaEvolution("expanded");
  APFEL::SetPDFEvolution("truncated");
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetPDFSet("repexternal");
  APFEL::SetPerturbativeOrder(pto);
  APFEL::SetAlphaQCDRef(asref, Qref);
  APFEL::SetPoleMasses(mc, mb, mt);
  APFEL::SetLHgridParameters(150, 75, xmin, 0.1, 1, 50, Qmin*Qmin, Qmax*Qmax);
  APFEL::LHAPDFgrid(nrep-1, Qmin, outname);
  APFEL::CleanUp();

  return 0;
}

// function to be called by APFEL
void externalsetapfelrep_(double* x, double* Q, int* irep, double* xf)
{
  // Get original PDFs
  std::map<int, double> pdfs = dist[*irep] -> xfxQ(*x, *Q);

  // Define kappa function
  const double kappa = ( rs * ( pdfs[-2] + pdfs[-1] ) - ( pdfs[3] + pdfs[-3] ) ) / ( 2 * ( 1 + rs ) );

  // Fill in array of PDFs
  xf[6-6] = pdfs[-6];
  xf[6-5] = pdfs[-5];
  xf[6-4] = pdfs[-4];
  xf[6-3] = pdfs[-3] + kappa;
  xf[6-2] = pdfs[-2] - kappa;
  xf[6-1] = pdfs[-1] - kappa;
  xf[6-0] = pdfs[21];
  xf[6+1] = pdfs[1];
  xf[6+2] = pdfs[2];
  xf[6+3] = pdfs[3] + kappa;
  xf[6+4] = pdfs[4];
  xf[6+5] = pdfs[5];
  xf[6+6] = pdfs[6];
  xf[6+7] = 0;
}
