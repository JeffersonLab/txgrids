#include <iostream>
#include <cmath>
#include <iomanip>
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
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

int main()
{
  // Include new search path in LHAPDF
  const string path = GetCurrentWorkingDir() + "/";
  LHAPDF::setPaths(path);

  APFEL::SetPDFSet("NNPDF31_nnlo_pch_as_0118_rs_1.0");
  APFEL::SetPerturbativeOrder(1);
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,130,3,1e-9);
  APFEL::SetGridParameters(2,60,5,1e-1);
  APFEL::SetGridParameters(3,20,5,8e-1);

  // Initializes integrals on the grids
  APFEL::InitializeAPFEL();

  double Q;
  cout << "Enter the scale in GeV" << endl;
  cin >> Q;

  APFEL::EvolveAPFEL(Q, Q);

  double momsr = 0.0;
  for (int i = -6; i < 7; i++)
    momsr += APFEL::NPDF(i,2);
  
  momsr += APFEL::Ngamma(2);

  double uvsr = APFEL::NPDF(2,1) - APFEL::NPDF(-2,1);
  double dvsr = APFEL::NPDF(1,1) - APFEL::NPDF(-1,1);
  double svsr = APFEL::NPDF(3,1) - APFEL::NPDF(-3,1);

  cout << setprecision(15);
  cout << "Sum rules at Q = " << Q << " GeV:\n" << endl;
  cout << "- Momentum sum rule         = " << momsr << endl;
  cout << "- Up valence sum rule       = " << uvsr << endl;
  cout << "- Down valence sum rule     = " << dvsr << endl;
  cout << "- Strange valence sum rule  = " << svsr << endl;
  cout << "\n";

  return 0;
}
