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

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class RMF : public NuclearModelI {

public:
  RMF();
  RMF(string config);
  virtual ~RMF();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & tgt,
                                          double q, double w) const;
  //double         Prob            (double mom, double E, const Target & t) const;
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

  double fPMax;
  double fPCutOff;
  bool fUseParametrization;

};

class RMFNucleus {

  RMFNucleus(int pdg);
  virtual ~RMFNucleus();

  const int fPDG; // nucleus PDG
  const int fNneut;//number of neutrons
  const int fNprot;//number of protons
  const int fA; //number of nucleons
  const int fNneut_shells; //number of neutron shells
  const int fNprot_shells; //number of proton shells

  std::string neut_dir;
  std::string prot_dir;
  std::string data_dir;

  std::string neut_file;
  std::string prot_file;

  int neut_shell_occ[fNneut_shells];// array of shell occupancies for neutrons
  int prot_shell_occ[fNprot_shells];// array of shell occupancies for protons
  double Ermv_neut[fNneut_shells]; // removal energies for neutrons
  double Ermv_prot[fNprot_shells]; // removal energies for protons
  
  int neut_shell_n[fNneut_shells];
  int neut_shell_l[fNneut_shells];
  int neut_shell_2j[fNneut_shells];

  int prot_shell_n[fNprot_shells];
  int prot_shell_l[fNprot_shells];
  int prot_shell_2j[fNprot_shells];

  vector<double> pnucl_neut_prob[fNneut_shells];
  vector<double> pnucl_prot_prob[fNprot_shells];

  vector<double> pnucl_neut;
  vector<double> pnucl_prot;

  TH1D* hist_prob_neut[fNneut_shells];
  TH1D* hist_prob_prot[fNprot_shells];

  double GetShellErmv(int nucleon_pdgc, int index);
  double GetShellQuantumNumber(char Qnumber, int index);
  double SampleErmv();


}

}         // genie namespace
#endif    // _RMF_H_

