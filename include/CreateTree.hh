#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  
  TTree*  ftree ;
  TString fname ;
  
public:
  
  CreateTree (TString name);
  ~CreateTree () ;
  
  TTree*          GetTree() const { return ftree; };
  TString         GetName() const { return fname; };
  void             AddEnergyDeposit(int index, float deposit);
  void             AddScintillationPhoton(int index);
  void             AddCerenkovPhoton(int index);
  int                Fill();
  bool             Write(TFile *);
  void             Clear() ;
  
  static CreateTree* Instance() { return fInstance; } ;
  static CreateTree* fInstance;
  
  int Event;


  int inputE1Thick;
  int inputE2Thick;
  int inputE1Width;


  
  std::vector<float>* inputMomentum ; // Px Py Pz E
  std::vector<float>* inputInitialPosition ; // x, y, z

  std::vector<float>* primaryMomT1 ; // Px Py Pz E
  std::vector<float>* primaryPosT1 ; // x, y, z

  std::vector<float>* primaryMomE1 ; // Px Py Pz E
  std::vector<float>* primaryPosE1 ; // x, y, z



  //integrated energy in each longitudinal layer
  float depositedEnergyEscapeWorld;

  float depositedEnergyTotal;
  float depositedEnergyECAL_f;
  float depositedEnergyECAL_r;
  float depositedEnergyEcalGap;
  float depositedEnergyEcalDet;
  float depositedEnergyWorld;

  float depositedIonEnergyTotal;
  float depositedIonEnergyTiming_f;
  float depositedIonEnergyTiming_r;
  float depositedIonEnergyECAL_f;
  float depositedIonEnergyECAL_r;
  float depositedIonEnergyEcalGap;
  float depositedIonEnergyEcalDet;
  float depositedIonEnergyWorld;

  float depositedElecEnergyTotal;
  float depositedElecEnergyECAL_f;
  float depositedElecEnergyECAL_r;
  float depositedElecEnergyEcalGap;
  float depositedElecEnergyEcalDet;
  float depositedElecEnergyWorld;
  


  int tot_phot_cer_ECAL_f ;
  int tot_phot_cer_ECAL_r ;







  //energy deposit in each trasnversally segmented channel
  float Edep_ECAL_f_ch[2500];
  float Edep_ECAL_r_ch[2500];






  TH1F* h_phot_cer_lambda_ECAL_f ;
  TH1F* h_phot_cer_lambda_ECAL_r;
  TH1F* h_phot_cer_parentID;


};

#endif
