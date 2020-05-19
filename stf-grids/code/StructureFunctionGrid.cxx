/*
  Authors: Valerio Bertone and Rabah Abdul Khalek
 */

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL libs
#include "apfel/APFEL.h"

#include <unistd.h>
#define GetCurrentDir getcwd

using namespace std;

string GetCurrentWorkingDir()
{
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  string current_working_dir(buff);
  return current_working_dir;
} 

int main() {
  // Include new search path in LHAPDF
  const string path = GetCurrentWorkingDir() + "/";
  LHAPDF::setPaths(path);

  // Open LHAPDF set
  const string set = "NNPDF31_nnlo_pch_as_0118";
  vector<LHAPDF::PDF*> dist = LHAPDF::mkPDFs(set);

  // Retrieve relevant parameters
  const int    Nrep  = dist.size() - 1;
  const int    pto   = dist[0]->orderQCD();
  const double Qref  = 91.2;
  const double asref = dist[0]->alphasQ(Qref);
  const double mc    = dist[0]->quarkThreshold(4);
  const double mb    = dist[0]->quarkThreshold(5);
  const double mt    = dist[0]->quarkThreshold(6);
  const double Q2min = dist[0]->q2Min();
  const double Qin   = sqrt(Q2min);
  const double Q2max = dist[0]->q2Max();
  const double xmin  = dist[0]->xMin();

  // EW parameters
  /*
  const double s2tw = 0.23126;
  const double aem  = 1. / 127.955;
  const double MZ   = 91.1876;
  const double MW   = MZ * sqrt( 1 - s2tw );
  const double MZ2  = MZ * MZ;
  const double MW2  = MW * MW;
  const double GF   = aem * M_PI / sqrt(2) / MW2 /  ( 1 - MW2 / MZ2 );
  */
  // Prescription is from Denner, eq. (8.25) of
  // https://arxiv.org/pdf/0709.1075.pdf The correction (second term)
  // is: Nc * Sqrt[2] GF mtop^2 / (16 pi^2) It is the universal
  // correction involving tops, relating mw and mz
  const double MZ   = 91.1876;
  const double MW   = 80.385;
  const double MZ2  = MZ * MZ;
  const double MW2  = MW * MW;
  const double GF   = 1.1663787e-5;
  const double drho = 0.00940161;
  const double rho  = 1 + drho;
  const double s2tw = 1 - MW2 / MZ2 / rho;

  APFEL::SetSin2ThetaW(s2tw);
  APFEL::SetZMass(MZ);
  APFEL::SetWMass(MW);
  APFEL::SetGFermi(GF);
  if(pto == 1)
    APFEL::SetMassScheme("FONLL-B");
  else if(pto == 2)
    APFEL::SetMassScheme("FONLL-C");
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetQLimits(Qin,1e5);
  APFEL::LockGrids(true);
  APFEL::SetAlphaEvolution("expanded");
  APFEL::SetPDFEvolution("truncated");
  APFEL::SetPDFSet(set);
  APFEL::SetPerturbativeOrder(pto);
  APFEL::SetAlphaQCDRef(asref,Qref);
  APFEL::SetPoleMasses(mc,mb,mt);
  APFEL::SetLHgridParameters(150, 75, xmin, 0.1, 1, 100, Q2min, Q2max);
  APFEL::EnableTargetMassCorrections(true);
  APFEL::LHAPDFgridStructureFunctions(dist.size()-1, Qin, set + "_SF");

  // Test produced grid
  vector<LHAPDF::PDF*> disttest = LHAPDF::mkPDFs(set + "_SF");
  const vector<double> xtest{0.00001, 0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9};
  const vector<double> Qtest{2, 5, 10, 20, 100, 1000, 10000};
  APFEL::SetReplica(0);
  APFEL::SetProcessDIS("NC");
  //APFEL::SetProjectileDIS("neutrino");

  cout << scientific;
  for (auto const& q: Qtest)
    {
      APFEL::ComputeStructureFunctionsAPFEL(Qin, q);
      for (auto const& x: xtest)
	cout << q << "  " << x << "  " << disttest[0]->xfxQ(1001, x, q) << "  " << APFEL::F2total(x) << "  " << disttest[0]->xfxQ(1001, x, q) / APFEL::F2total(x) << endl;
      cout << "\n";
    }
  APFEL::CleanUp();

  return 0;
}
