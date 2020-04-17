#include "CreateTree.hh"
#include <algorithm>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if( fInstance )
  {
    return ;
  }
  
  this -> fInstance = this ;
  this -> fname     = name ;
  this -> ftree     = new TTree (name,name) ;

  
  this -> GetTree ()->Branch ("Event", &this->Event, "Event/I") ;



  this -> GetTree() -> Branch("inputE1Thick",      		&this->inputE1Thick,               	"inputE1Thick/F");
  this -> GetTree() -> Branch("inputE2Thick",      		&this->inputE2Thick,               	"inputE2Thick/F");
  this -> GetTree() -> Branch("inputE1Width",      		&this->inputE1Width,              	"inputE1Width/F");


  
  inputInitialPosition	= new vector<float>(3,0.); 
  inputMomentum		= new vector<float>(4,0.);
  primaryPosT1		= new vector<float>(3,0.); 
  primaryMomT1		= new vector<float>(4,0.); 
  primaryPosE1		= new vector<float>(3,0.); 
  primaryMomE1		= new vector<float>(4,0.);  

  this -> GetTree() -> Branch("inputInitialPosition",    "vector<float>", &inputInitialPosition);
  this -> GetTree() -> Branch("inputMomentum",        "vector<float>", &inputMomentum);
  this -> GetTree() -> Branch("primaryPosT1", 	"vector<float>", &primaryPosT1);
  this -> GetTree() -> Branch("primaryMomT1",	"vector<float>", &primaryMomT1);
  this -> GetTree() -> Branch("primaryPosE1", 	"vector<float>", &primaryPosE1);
  this -> GetTree() -> Branch("primaryMomE1",	"vector<float>", &primaryMomE1);  




  //integrated per longitudinal layer 
  this -> GetTree() -> Branch("depositedEnergyTotal",     	&this->depositedEnergyTotal,        	"depositedEnergyTotal/F");
  this -> GetTree() -> Branch("depositedEnergyEscapeWorld",     	&this->depositedEnergyEscapeWorld,        	"depositedEnergyEscapeWorld/F");

  this -> GetTree() -> Branch("depositedEnergyECAL_f",     &this->depositedEnergyECAL_f,            "depositedEnergyECAL_f/F");


  this -> GetTree() -> Branch("depositedEnergyWorld",     	&this->depositedEnergyWorld,         	"depositedEnergyWorld/F");
  this -> GetTree() -> Branch("depositedEnergyEcalGap",     	&this->depositedEnergyEcalGap,         	"depositedEnergyEcalGap/F");
  this -> GetTree() -> Branch("depositedEnergyEcalDet",     	&this->depositedEnergyEcalDet,         	"depositedEnergyEcalDet/F");


  this -> GetTree() -> Branch("depositedIonEnergyTotal",     	&this->depositedIonEnergyTotal,        	"depositedIonEnergyTotal/F");

  this -> GetTree() -> Branch("depositedIonEnergyECAL_f",     &this->depositedIonEnergyECAL_f,            "depositedIonEnergyECAL_f/F");

  this -> GetTree() -> Branch("depositedIonEnergyWorld",     	&this->depositedIonEnergyWorld,         	"depositedIonEnergyWorld/F");
  this -> GetTree() -> Branch("depositedIonEnergyEcalGap",     	&this->depositedIonEnergyEcalGap,         	"depositedIonEnergyEcalGap/F");
  this -> GetTree() -> Branch("depositedIonEnergyEcalDet",     	&this->depositedIonEnergyEcalDet,         	"depositedIonEnergyEcalDet/F");


  this -> GetTree() -> Branch("depositedElecEnergyTotal",     	&this->depositedElecEnergyTotal,        	"depositedElecEnergyTotal/F");

  this -> GetTree() -> Branch("depositedElecEnergyECAL_f",     &this->depositedElecEnergyECAL_f,            "depositedElecEnergyECAL_f/F");

  this -> GetTree() -> Branch("depositedElecEnergyWorld",     	&this->depositedElecEnergyWorld,         	"depositedElecEnergyWorld/F");
  this -> GetTree() -> Branch("depositedElecEnergyEcalGap",     	&this->depositedElecEnergyEcalGap,         	"depositedElecEnergyEcalGap/F");
  this -> GetTree() -> Branch("depositedElecEnergyEcalDet",     	&this->depositedElecEnergyEcalDet,         	"depositedElecEnergyEcalDet/F");


  //single channels
  this -> GetTree() -> Branch("Edep_ECAL_f_ch",   		&this->Edep_ECAL_f_ch,      			"Edep_ECAL_f_ch[100]/F");

  

  //Cerenkov photons
  this -> GetTree() -> Branch("tot_phot_cer_ECAL_f",        &this->tot_phot_cer_ECAL_f,               "tot_phot_cer_ECAL_f/I");





  h_phot_cer_lambda_ECAL_f    = new TH1F("h_phot_cer_lambda_ECAL_f","",1250, 0.,1250.);

  h_phot_cer_parentID = new TH1F("h_phot_cer_parentID","",600,-300,300);
  h_detected_photon = new TH1F("h_detected_photon","",500,0.,1000.);


  
  this -> Clear() ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



CreateTree::~CreateTree()
{}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



int CreateTree::Fill() 
{ 
  return this -> GetTree() -> Fill(); 
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



bool CreateTree::Write(TFile * outfile)
{
  outfile -> cd();
  ftree -> Write();


  h_phot_cer_lambda_ECAL_f    ->Write();

  h_phot_cer_parentID->Write();
  h_detected_photon->Write();

  return true ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



void CreateTree::Clear()
{
  Event	= 0;
  




  depositedEnergyEscapeWorld=0.;

  depositedEnergyTotal = 0.;
  depositedEnergyECAL_f = 0.;

  depositedEnergyWorld = 0.;
  depositedEnergyEcalGap = 0.;
  depositedEnergyEcalDet = 0.;


  depositedIonEnergyTotal = 0.;
  depositedIonEnergyECAL_f = 0.;

  depositedIonEnergyWorld = 0.;
  depositedIonEnergyEcalGap = 0.;
  depositedIonEnergyEcalDet = 0.;


  depositedElecEnergyTotal = 0.;
  depositedElecEnergyECAL_f = 0.;

  depositedElecEnergyWorld = 0.;
  depositedElecEnergyEcalGap = 0.;
  depositedElecEnergyEcalDet = 0.;


  tot_phot_cer_ECAL_f = 0.;






  for (int iCh = 0; iCh<100; iCh++)
  {
	Edep_ECAL_f_ch[iCh] = 0.;

  }

  

  
  for (int i = 0 ; i < 3 ; ++i) 
  {
    inputInitialPosition -> at(i) = 0.;
    primaryPosT1 -> at(i) = 0.;
    primaryPosE1 -> at(i) = 0.;

  }
  for (int i = 0 ; i < 4 ; ++i) 
  {
    inputMomentum ->at(i) = 0.;
    primaryMomT1 -> at(i) = 0.;
    primaryMomE1 -> at(i) = 0.;
  }




}
