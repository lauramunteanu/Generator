//____________________________________________________________________________
/*
 RMF nuclear model

 Authors: G. Megias, J. Patiño, J. Rosa, S. Dolan, L. Munteanu
*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/RMF.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
RMF::RMF() :
NuclearModelI("genie::RMF")
{

}
//____________________________________________________________________________
RMF::RMF(string config) :
NuclearModelI("genie::RMF", config)
{

}
//____________________________________________________________________________
RMF::~RMF()
{
  map<string, TH1D*>::iterator iter = fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.end(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//____________________________________________________________________________
bool RMF::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  RMFNucleus* nucleus = fNuclearMap[target.Pdg()];

  //choose a removal energy
  RandomGen * rnd = RandomGen::Instance();
  double nProt = nucleus->fNprot;

  // is this a neutron or a proton?
  bool isNeut = pdg::IsNeutron(target.HitNucPdg());
  bool isProt = pdg::IsProton(target.HitNucPdg());

  assert(isProt||isNeut);

  int shellindex = 0;

  if (isProt){
    int rand_whichshell = rnd->RndGen().Integer(nProt);
    int this_shell_occ = 0;

    while (shellindex <= rand_whichshell ){
      if (this_shell_occ > nucleus->prot_shell_occ[shellindex]) {
        shellindex++;
        this_shell_occ=0;
      }
      else {this_shell_occ ++;}
    }

    fCurrRemovalEnergy = nucleus->Ermv_prot[shellindex];
  }
  else if (isNeut){
    int rand_whichshell = rnd->RndGen().Integer(nNeut);
    int this_shell_occ = 0;

    while (shellindex <= rand_whichshell ){
      if (this_shell_occ > nucleus->neut_shell_occ[shellindex]) {
        shellindex++;
        this_shell_occ=0;
      }
      else {this_shell_occ ++;}
    }

    fCurrRemovalEnergy = nucleus->Ermv_neut[shellindex];
  }

  //-- set momentum vector
  //
  TH1D * prob;
  if (isProt) prob = nucleus->hist_prob_prot[shellindex];
  else if (isNeut) prob = nucleus->hist_prob_neut[shellindex];
  if ( ! prob ) {
    LOG("RMF nuclear model", pNOTICE)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }
  double p = prob->GetRandom();
  LOG("RMF nuclear model", pINFO) << "|p,nucleon| = " << p;

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  fCurrMomentum.SetXYZ(px,py,pz);

  return true;
}
//____________________________________________________________________________
// double RMF::Prob(double mom, double w, const Target & target) const
// {
//   if(w<0) {
//      TH1D * prob_distr = this->ProbDistro(target);
//      int bin = prob_distr->FindBin(mom);
//      double y  = prob_distr->GetBinContent(bin);
//      double dx = prob_distr->GetBinWidth(bin);
//      double prob = y*dx;
//      return prob;
//   }
//   return 1;
// }
//____________________________________________________________________________

//____________________________________________________________________________

//____________________________________________________________________________
void RMF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RMF::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void RMF::LoadConfig(void)
{

  NuclearModelI::LoadConfig() ; 

  LOG("RMF nuclear model", pDEBUG) << "Loading coonfiguration for RMF nuclear model";

  // Load removal energies for specific nuclei  and shells from either the algorithm's
  // configuration file or the UserPhysicsOptions file.

  string data_dir =
      string(gSystem->Getenv("GENIE")) + string("/data/evgen/nucl/RMF/");
  vector<int> PDGlist;
  map<int, RMFNucleus*> fNuclearMap;

  string pdglist_file = data_dir + "ListOfNuclei.dat";

  int thispdg;
  while (std::getline(pdglist_file, thispdg))
  {
    LOG("RMF nuclear model", pDEBUG) << "Adding nucleus " << thispdg << " to list of PDGs";
    PDGlist.push_back(thispdg);
    RMFNucleus* thisnucleus = new RMFNucleus(thispdg);
    fNuclearMap.insert({thispdg, thisnucleus});
  }

  LOG("RMF nuclear model", pDEBUG)
    << "Finished LoadConfig";

}
//____________________________________________________________________________

//____________________________________________________________________________
RMFNucleus::RMFNucleus(int pdg)
{
  fPDG = pdg;
  this->FillNucleusInfo();
}
//____________________________________________________________________________
RMFNucleus::~RMFNucleus()
{

}
//____________________________________________________________________________
RMFNucleus::FillNucleusInfo()
{
  string data_dir =
      string(gSystem->Getenv("GENIE")) + string("/data/evgen/nucl/RMF/") + std::to_string(fPDG) + "/";

  neut_dir = data_dir + "n";
  prot_dir = data_dir + "p";

  neut_ermv_file = neut_dir + "/n.rmf";
  prot_ermv_file = prot_dir + "/p.rmf"

  std::ifstream in_file_ermv_neut( neut_ermv_file.c_str() );
  std::ifstream in_file_ermv_prot( prot_ermv_file.c_str() );

  // Skip the initial comment line
  std::string dummy;
  std::getline(in_file_ermv_neut, dummy);
  std::getline(in_file_ermv_prot, dummy);
  
  in_file_ermv_neut >> fNneut >> fNneut_shells;

  in_file_ermv_prot >> fNprot >> fNprot_shells;

  for (int i_neut_shell=0; i_neut_shell < fNneut_shells; i_neut_shell++ ){
    in_file_ermv_neut >> Ermv_neut[i_neut_shell] >> neut_shell_n[i_neut_shell] >> neut_shell_l[i_neut_shell] >> neut_shell_2j[i_neut_shell] >> neut_shell_occ[i_neut_shell];
  }

  for (int i_prot_shell=0; i_prot_shell < fNprot; i_prot_shell++ ){
    in_file_ermv_prot >> Ermv_prot[i_prot_shell] >> prot_shell_n[i_prot_shell] >> prot_shell_l[i_prot_shell] >> prot_shell_2j[i_prot_shell] >> prot_shell_occ[i_prot_shell];
  }
  
  fA = fNneut + fNprot; 
  
  neut_pnucl_file = neut_dir + "/RMF_n.rmf";
  prot_pnucl_file = prot_dir + "/RMF_p.rmf"

  std::ifstream in_file_pnucl_neut( neut_pnucl_file.c_str() );
  std::ifstream in_file_pnucl_prot( prot_pnucl_file.c_str() );

  double pnucl_neut_holder;
  double pnucl_neut_prob_holder;
  double pnucl_neut_array_holder[fNneut_shells];
  
  while(!in_file_pnucl_neut.eof())
  {
    in_file_pnucl_neut >> pnucl_neut_holder;
    pnucl_neut.push_back(pnucl_neut_holder);
    for(int i_neut_shell=0; i_neut_shell < fNneut_shells; i_neut_shell++ ){
      in_file_pnucl_neut >> pnucl_neut_prob_holder;
      pnucl_neut_array_holder[i_neut_shell] = pnucl_neut_prob_holder;
    }
    pnucl_neut_prob.push_back(pnucl_neut_array_holder);
  }

  double pnucl_prot_holder;
  double pnucl_prot_prob_holder;
  double pnucl_prot_array_holder[fNneut_shells];
  
  while(!in_file_pnucl_prot.eof())
  {
    in_file_pnucl_prot >> pnucl_prot_holder;
    pnucl_prot.push_back(pnucl_prot_holder);
    for(int i_prot_shell=0; i_prot_shell < fNprot_shells; i_prot_shell++ ){
      in_file_pnucl_prot >> pnucl_prot_prob_holder;
      pnucl_prot_array_holder[i_prot_shell] = pnucl_prot_prob_holder;
    }
    pnucl_prot_prob.push_back(pnucl_neut_array_holder);
  }

  //debug file, delete me please
  TFile* fdebug = new TFile(Form("debug_%i.root", fPDG));

  for (int i_neut_shell=0; i_neut_shell < fNneut_shells; i_neut_shell++ ){
    hist_prob_neut[i_neut_shell] = new TH1D(Form("hist_prob_neut_%i", i_neut_shell), Form("hist_prob_neut_%i", i_neut_shell), fNneut_shells, &pnucl_neut[0]);
    for (i_pnucl_bin = 0; i_pnucl_bin < hist_prob_neut[i_neut_shell]->GetNbins(); i_pnucl_bin++){
      hist_prob_neut[i_neut_shell]->SetBinContent(i_pnucl_bin+1, pnucl_neut_prob[i_neut_shell][i_pnucl_bin]);
    }
    hist_prob_neut[i_neut_shell]->Scale(1./hist_prob_neut[i_neut_shell]->Integral());
    fdebug->cd();
    hist_prob_neut[i_neut_shell]->Write(hist_prob_neut[i_neut_shell]->GetName()); 
  }
  
  for (int i_prot_shell=0; i_prot_shell < fNprot_shells; i_prot_shell++ ){
    hist_prob_prot[i_prot_shell] = new TH1D(Form("hist_prob_prot_%i", i_prot_shell), Form("hist_prob_prot_%i", i_prot_shell), fNprot_shells, &pnucl_prot[0]);
    for (i_pnucl_bin = 0; i_pnucl_bin < hist_prob_prot[i_prot_shell]->GetNbins(); i_pnucl_bin++){
      hist_prob_prot[i_prot_shell]->SetBinContent(i_pnucl_bin+1, pnucl_prot_prob[i_prot_shell][i_pnucl_bin]);
    }
    hist_prob_prot[i_prot_shell]->Scale(1./hist_prob_prot[i_prot_shell]->Integral());
    fdebug->cd();
    hist_prob_prot[i_prot_shell]->Write(hist_prob_prot[i_neut_shell]->GetName()); 
  }

  fdebug->Close();

}
//____________________________________________________________________________
