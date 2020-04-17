#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
//#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();



using namespace std;
using namespace CLHEP;



int to_int (string name)
{
  int Result ;             // int which will contain the result
  stringstream convert (name) ;
  string dummy ;           
  convert >> dummy ;       
  convert >> Result ;
  return Result ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


SteppingAction::SteppingAction (DetectorConstruction* detectorConstruction,
                                const G4int& scint, const G4int& cher):
  fDetectorConstruction(detectorConstruction),
  propagateScintillation(scint),
  propagateCerenkov(cher)
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


SteppingAction::~SteppingAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void SteppingAction::UserSteppingAction (const G4Step * theStep)
{
 
  G4Track* theTrack = theStep->GetTrack () ;
  
//  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
//  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();

//  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());

  G4ParticleDefinition* particleType = theTrack->GetDefinition () ;
  //G4int trackID = theTrack->GetTrackID();
  
  G4StepPoint * thePrePoint  = theStep->GetPreStepPoint () ;
  G4StepPoint * thePostPoint = theStep->GetPostStepPoint () ;
  //  const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName  = "" ; if ( thePrePV )  thePrePVName  = thePrePV  -> GetName () ;
  G4String thePostPVName = "" ; if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;
  //  std::cout<<thePrePVName<<std::endl;

//  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();
    
  G4int nStep = theTrack -> GetCurrentStepNumber();

  G4int TrPDGid= theTrack->GetDefinition()->GetPDGEncoding();

  
//        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------
 
  // get position
  

    G4double energy = theStep->GetTotalEnergyDeposit();
    G4double energyIon = theStep->GetNonIonizingEnergyDeposit();
    G4double energyElec=0.;

    if(abs(TrPDGid)==11) {
      energyElec=energy-energyIon;
    }
    //std::cout<<"TrPDGid energy energyIon enegyElec are "<<TrPDGid<<" "<<energy<<" "<<energyIon<<" "<<energyElec<<std::endl;

    CreateTree::Instance() -> depositedEnergyTotal += energy/GeV;
    CreateTree::Instance() -> depositedIonEnergyTotal += energyIon/GeV;
    CreateTree::Instance() -> depositedElecEnergyTotal += energyElec/GeV;

    //    if(thePrePVName.contains("world")) {
      bool haha4=((theStep->GetPostStepPoint())->GetStepStatus())==fWorldBoundary;
      if(haha4) {
	//std::cout<<"leaving "<<std::endl;
	CreateTree::Instance() -> depositedEnergyEscapeWorld += (theStep->GetPostStepPoint())->GetKineticEnergy()/GeV;
      }


  
  // optical photon

  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  {
  
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
    float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV);
      //kill very long wavelengths
    if (photWL> 1000 ||  photWL< 300)   theTrack->SetTrackStatus(fKillTrackAndSecondaries); 

    else{

      if(thePostPVName.contains("ecalDet")) {
	CreateTree::Instance()->h_detected_photon->Fill(photWL);
	  theTrack->SetTrackStatus(fKillTrackAndSecondaries); 
	}
        
    if(
        (nStep == 1) && (processName == "Cerenkov") )
    {
      
      TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());
      //G4ThreeVector haha=theTrackInfo->GetParentMomentum();
      //G4double haha2=theTrackInfo->GetParentEnergy()/GeV;
      //G4double haha3=haha.mag()/GeV;
      //G4double betaa=0.;
      //if(haha2>0) betaa=haha3/haha2;

      G4int aapdgid=theTrackInfo->GetParentPDGid();
      CreateTree::Instance()->h_phot_cer_parentID -> Fill( aapdgid );

      
      //std::cout << " generated Cerenkov photon with parent " << theTrackInfo->GetParentName()<<" "<<aapdgid<<" with beta of "<<betaa<<" and energy "<<haha2<<std::endl;



      if (thePrePVName.contains("ecalCrystalP_f"))
      {      
	CreateTree::Instance()->tot_phot_cer_ECAL_f += 1;
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_f -> Fill( photWL);
      }
      else if(thePrePVName.contains("ecalDet"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else if(thePrePVName.contains("ecalGap"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  //theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else if(thePrePVName.contains("world"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else 
	{
	  std::cout<<"weird PrePVName "<<thePrePVName<<std::endl;
	}
      if( thePostPVName.contains("world") ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      

      if( !propagateCerenkov ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      

    } // end nstep==1 and cerenkov
    else { // later steps

      if(thePrePVName.contains("ecalDet"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else if(thePrePVName.contains("ecalGap"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  //theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else if(thePrePVName.contains("world"))
	{
	  //	  std::cout<<"hit ecal photo detector"<<std::endl;
	  theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
	}
      else 
	{
	  std::cout<<"weird PrePVName "<<thePrePVName<<std::endl;
	}
      if( thePostPVName.contains("world") ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);     
      //      std::cout<<nStep<<" "<<processName<<" "<<thePrePVName<<std::endl;

    }
    }
  }  // end optical photon

    else{




    //ecal
    if( thePrePVName.contains("ecalCrystalP_f") )
    {
      CreateTree::Instance()->depositedEnergyECAL_f += energy/GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_f += energyIon/GeV;
      CreateTree::Instance()->depositedElecEnergyECAL_f += energyElec/GeV;
      for (int iCh = 0; iCh<2500; iCh++)
      {
	if (thePrePVName == Form("ecalCrystalP_f_%d", iCh)) CreateTree::Instance()->Edep_ECAL_f_ch[iCh] += energy/GeV;
      }
    }





    if( thePrePVName.contains("world") )
    {
      CreateTree::Instance() -> depositedEnergyWorld += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyWorld += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyWorld += energyElec/GeV;
    }




    


    if( thePrePVName.contains("ecalGap") )
    {
      CreateTree::Instance() -> depositedEnergyEcalGap += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyEcalGap += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyEcalGap += energyElec/GeV;
    }


    if( thePrePVName.contains("ecalDet") )
    {
      CreateTree::Instance() -> depositedEnergyEcalDet += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyEcalDet += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyEcalDet += energyElec/GeV;
    }



    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon
  
  
  return ;
}

