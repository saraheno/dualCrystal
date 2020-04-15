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
/*
SteppingAction::SteppingAction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;

  config.readInto(core_material, "core_material"); 

  if (core_material == 0)
  {
	  config.readInto(toy_ly,	"toy_ly");
	  config.readInto(toy_decay,	"toy_decay");
	  config.readInto(toy_rise,	"toy_rise");
  }
  

}*/



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
  const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName  = "" ; if ( thePrePV )  thePrePVName  = thePrePV  -> GetName () ;
  G4String thePostPVName = "" ; if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;

//  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();
    
  G4int nStep = theTrack -> GetCurrentStepNumber();

  G4int TrPDGid= theTrack->GetDefinition()->GetPDGEncoding();

  
//        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------
 
  // get position
  G4double global_x = thePrePosition.x()/mm;
  G4double global_y = thePrePosition.y()/mm;
  G4double global_z = thePrePosition.z()/mm;
  

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
      //}

  
  // optical photon

  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  {
  
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();

        
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
      float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV);

      //kill very long wavelengths
      if (photWL> 1000 ||  photWL< 300)  theTrack->SetTrackStatus(fKillTrackAndSecondaries); 



      else if (thePrePVName.contains("corePV_front"))
      {      
	CreateTree::Instance()->tot_phot_cer_Timing_f += 1;
        CreateTree::Instance()->h_phot_cer_lambda_Timing_f -> Fill( photWL );


      }
      else if (thePrePVName.contains("corePV_rear"))
      {      
	CreateTree::Instance()->tot_phot_cer_Timing_r += 1;
        CreateTree::Instance()->h_phot_cer_lambda_Timing_r -> Fill( photWL );
      }

      else if (thePrePVName.contains("ecalCrystalP_f"))
      {      
	CreateTree::Instance()->tot_phot_cer_ECAL_f += 1;
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_f -> Fill( photWL);
      }
      else if (thePrePVName.contains("ecalCrystalP_r"))
      {      
	CreateTree::Instance()->tot_phot_cer_ECAL_r += 1;
        CreateTree::Instance()->h_phot_cer_lambda_ECAL_r -> Fill( photWL );
      }
      else if (thePrePVName.contains("hcalTile_layer"))
      {      
	CreateTree::Instance()->tot_phot_cer_HCAL += 1;
        CreateTree::Instance()->h_phot_cer_lambda_HCAL -> Fill( photWL );
      }


      if( !propagateCerenkov ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      

    }
  }

    else{

   //count tracks before SCEPCAL at the tracker layers
    if (  thePrePVName.contains("world") && thePostPVName.contains("trackerPV_Layer")    )	// interface between T1 and T2
   {
        for (int iLayer = 0; iLayer<6; iLayer++)
	{
		if (thePostPVName == Form("trackerPV_Layer%d", iLayer)    )	// interface between T1 and T2
	        CreateTree::Instance() -> nTracksTRK[iLayer] ++;   //counting tracks crossing the boundary
	}
   }


    //save primary particle position and momentum before timing layer T1 and before ECAL E1
    else if ( thePrePVName.contains("world") && thePostPVName.contains("corePV_front")  ) // interface between world and T1
   {
 
        CreateTree::Instance() -> nTracksT1 ++;   //counting tracks crossing the boundary

	if ( theTrack->GetParentID() == 0 )	// select only the primary particle
	{
	        CreateTree::Instance() -> primaryPosT1->at(0) = global_x;
	        CreateTree::Instance() -> primaryPosT1->at(1) = global_y;
	        CreateTree::Instance() -> primaryPosT1->at(2) = global_z;

	        CreateTree::Instance() -> primaryMomT1->at(0) = thePrePoint->GetMomentum().x()/GeV;
        	CreateTree::Instance() -> primaryMomT1->at(1) = thePrePoint->GetMomentum().y()/GeV;
	        CreateTree::Instance() -> primaryMomT1->at(2) = thePrePoint->GetMomentum().z()/GeV;
	        CreateTree::Instance() -> primaryMomT1->at(3) = thePrePoint->GetTotalEnergy()/GeV;
	}
   }

    else if (  thePrePVName.contains("world") && thePostPVName.contains("corePV_rear")    )	// interface between T1 and T2
   {
        CreateTree::Instance() -> nTracksT2 ++;   //counting tracks crossing the boundary
   }


    else if ( 	( thePrePVName.contains("world")  || thePrePVName.contains("ecalGapP_f") || thePrePVName.contains("ecalDetP_f") ) 
		&& thePostPVName.contains("ecalCrystalP_f") 	// interface between world and E1
       )
   {
        CreateTree::Instance() -> nTracksE1 ++;   //counting tracks crossing the boundary

	if ( theTrack->GetParentID() == 0 )	// select only the primary particle
	{
	        CreateTree::Instance() -> primaryPosE1->at(0) = global_x;
	        CreateTree::Instance() -> primaryPosE1->at(1) = global_y;
	        CreateTree::Instance() -> primaryPosE1->at(2) = global_z;

	        CreateTree::Instance() -> primaryMomE1->at(0) = thePrePoint->GetMomentum().x()/GeV;
        	CreateTree::Instance() -> primaryMomE1->at(1) = thePrePoint->GetMomentum().y()/GeV;
	        CreateTree::Instance() -> primaryMomE1->at(2) = thePrePoint->GetMomentum().z()/GeV;
	        CreateTree::Instance() -> primaryMomE1->at(3) = thePrePoint->GetTotalEnergy()/GeV;
	}
   }

    else if ( 	( thePrePVName.contains("ecalCrystalP_f") || thePrePVName.contains("world") )
		&& thePostPVName.contains("ecalCrystalP_r"  )	
       )// interface between E1 and E2
   {
        CreateTree::Instance() -> nTracksE2 ++;   //counting tracks crossing the boundary
   }



    //tracker
    if( thePrePVName.contains("trackerPV_Layer") )
    {
      for (int iLayer = 0; iLayer<6; iLayer++)
      {
	if (thePrePVName == Form("trackerPV_Layer%d", iLayer)) CreateTree::Instance()->Edep_Tracker_layer[iLayer] += energy/GeV;
      }
    }

    //timing
    if( thePrePVName.contains("corePV_front") )
    {
      CreateTree::Instance()->depositedEnergyTiming_f += energy/GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_f += energyIon/GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_f += energyElec/GeV;
      for (int iBar = 0; iBar<18; iBar++)
      {
	if (thePrePVName == Form("corePV_front_%d", iBar)) CreateTree::Instance()->Edep_Timing_f_ch[iBar] += energy/GeV;
      }
    }
    if( thePrePVName.contains("corePV_rear") )
    {
      CreateTree::Instance()->depositedEnergyTiming_r += energy/GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_r += energyIon/GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_r += energyElec/GeV;
      for (int iBar = 0; iBar<18; iBar++)
      {
	if (thePrePVName == Form("corePV_rear_%d", iBar)) CreateTree::Instance()->Edep_Timing_r_ch[iBar] += energy/GeV;
      }	
    }

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
    if( thePrePVName.contains("ecalCrystalP_r") )
    {
      CreateTree::Instance()->depositedEnergyECAL_r += energy/GeV;
      CreateTree::Instance()->depositedIonEnergyECAL_r += energyIon/GeV;
      CreateTree::Instance()->depositedElecEnergyECAL_r += energyElec/GeV;
      for (int iCh = 0; iCh<2500; iCh++)
      {
	if (thePrePVName == Form("ecalCrystalP_r_%d", iCh)) CreateTree::Instance()->Edep_ECAL_r_ch[iCh] += energy/GeV;
      }
    }    


    //hcal
    if( thePrePVName.contains("hcalTile_layer") )
    {
       CreateTree::Instance()->depositedEnergyHCALAct += energy/GeV;
       CreateTree::Instance()->depositedIonEnergyHCALAct += energyIon/GeV;
       CreateTree::Instance()->depositedElecEnergyHCALAct += energyElec/GeV;
//      for (int iLayer = 0; iLayer<100; iLayer++)
  //    {
//	if (thePrePVName == Form("ecalCrystalP_f_%d", iCh)) CreateTree::Instance()->Edep_ECAL_f_ch[iLayer] += energy/GeV;
//      }
    }
    if( thePrePVName.contains("hcalAbs") )
    {
       CreateTree::Instance()->depositedEnergyHCALPas += energy/GeV;
       CreateTree::Instance()->depositedIonEnergyHCALPas += energyIon/GeV;
       CreateTree::Instance()->depositedElecEnergyHCALPas += energyElec/GeV;
//      for (int iLayer = 0; iLayer<100; iLayer++)
  //    {
//	if (thePrePVName == Form("ecalCrystalP_f_%d", iCh)) CreateTree::Instance()->Edep_ECAL_f_ch[iLayer] += energy/GeV;
//      }
    }


    if( thePrePVName.contains("world") )
    {
      CreateTree::Instance() -> depositedEnergyWorld += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyWorld += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyWorld += energyElec/GeV;
    }



    if( thePrePVName.contains("services") )
    {
      CreateTree::Instance() -> depositedEnergyServices += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyServices += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyServices += energyElec/GeV;
    }

    

    if( thePrePVName.contains("TimingGap") )
    {
      CreateTree::Instance() -> depositedEnergyTimingGap += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergyTimingGap += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergyTimingGap += energyElec/GeV;
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


    if( thePrePVName.contains("solenoid") )
    {
      CreateTree::Instance() -> depositedEnergySolenoid += energy/GeV;
      CreateTree::Instance() -> depositedIonEnergySolenoid += energyIon/GeV;
      CreateTree::Instance() -> depositedElecEnergySolenoid += energyElec/GeV;
    }



    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon
  
  
  return ;
}

