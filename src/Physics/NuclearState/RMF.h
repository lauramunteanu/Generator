//____________________________________________________________________________
/*!

\class    genie::RMF

 RMF nuclear model

 Authors: G. Megias, J. Patiño, J. Rosa, S. Dolan, L. Munteanu
          
*/
//____________________________________________________________________________

#ifndef _RMF_H_
#define _RMF_H_

#include <map>
#include <fstream>

#include <TH1D.h>
#include <TFile.h>
#include "TSystem.h"

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class RMFNucleus {

  public:

  RMFNucleus(int pdg);
  ~RMFNucleus();

  int fPDG; // nucleus PDG
  int fNneut;//number of neutrons
  int fNprot;//number of protons
  int fA; //number of nucleons
  int fNneut_shells; //number of neutron shells
  int fNprot_shells; //number of proton shells

  std::string neut_dir;
  std::string prot_dir;
  std::string data_dir;

  std::string neut_ermv_file;
  std::string prot_ermv_file;

  int* neut_shell_occ;// array of shell occupancies for neutrons
  int* prot_shell_occ;// array of shell occupancies for protons
  double* Ermv_neut; // removal energies for neutrons
  double* Ermv_prot; // removal energies for protons
  
  int* neut_shell_n;
  int* neut_shell_l;
  int* neut_shell_2j;

  int* prot_shell_n;
  int* prot_shell_l;
  int* prot_shell_2j;

  int n_mom_bins_neut;
  int n_mom_bins_prot;

  std::string neut_pnucl_file;
  std::string prot_pnucl_file;

  std::vector<double> pnucl_neut;
  std::vector<double> pnucl_prot;

  std::vector< std::vector<double> > pnucl_neut_prob;
  std::vector< std::vector<double> > pnucl_prot_prob;

  TH1D** hist_prob_neut;
  TH1D** hist_prob_prot;

  //double GetShellErmv(int nucleon_pdgc, int index);
  //double GetShellQuantumNumber(char Qnumber, int index);
  //double SampleErmv();
  void FillNucleusInfo();

};

class RMF : public NuclearModelI {

public:
  RMF();
  RMF(string config);
  virtual ~RMF();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & tgt) const;
  //bool           GenerateNucleon (const Target & tgt,
  //                                        double q, double w) const;
  double         Prob            (double p, double w, const Target & t) const;


  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmRMF; 
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

 protected:
  void   LoadConfig (void);
  
private:

  //TH1D * ProbDistro (const Target & t) const;

  mutable map<string, TH1D *> fProbDistroMap;

  map<int, double> fNucRmvE;
  map<int, RMFNucleus*> fNuclearMap;

  double fPMax;
  double fPCutOff;
  bool fUseParametrization;

};

}         // genie namespace
#endif    // _RMF_H_

