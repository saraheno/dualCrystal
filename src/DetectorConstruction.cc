//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

#include "DetectorConstruction.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"


#include "DetectorConstruction.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>

using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;
  
  config.readInto(checkOverlaps,	"checkOverlaps");
  
  config.readInto(world_material,	"world_material");
  config.readInto(bar_length,  	"bar_length");
  

  // gap between crystal and photodetector
  config.readInto(gap_l,         		"gap_l");
  config.readInto(gap_size_x,    	"gap_size_x");
  config.readInto(gap_size_y,    	"gap_size_y");
  config.readInto(gap_material,  	"gap_material");
  
  config.readInto(det_l,         		"det_l");
  config.readInto(det_material,  	"det_material");
  config.readInto(ecal_det_size,	"ecal_det_size");
  
  config.readInto(depth,			"depth");


  config.readInto(ecal_material, 	"ecal_material");
  config.readInto(ecal_front_length,"ecal_front_length");

  config.readInto(ecal_front_face,	"ecal_front_face");
  config.readInto(ecal_rear_face,	"ecal_rear_face");


  config.readInto(ecal_n_cell,"ecal_n_cell");


  config.readInto(cReffile, "cReffile");
  config.readInto(crystal_reflectivity, "cReflectivity");
  config.readInto(cSurrefind, "cSurrefind");
  config.readInto(cSurtype, "cSurtype");
  config.readInto(cSpecularspike, "cSpecularspike");
  config.readInto(cSpecularlobe, "cSpecularlobe");
  config.readInto(cSigmaalpha, "cSigmaalpha");
  config.readInto(cSpecularspike, "cSpecularspike");
  config.readInto(cSpecularlobe, "cSpecularlobe");
  config.readInto(cSigmaalpha, "cSigmaalpha");
  config.readInto(cBackscatter, "cBackscatter");
  config.readInto(cLambertian, "cLambertian");
  config.readInto(crystalSurfinish, "crystalSurfinish");


  
  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;
  
  //expHall_x = 450.*cm;
  //expHall_y = 450.*cm;
  //expHall_z = 300.*cm;

  expHall_x = sqrt(ecal_n_cell)*ecal_front_length*3;
  expHall_y = sqrt(ecal_n_cell)*ecal_front_length*3;
  expHall_z = 3*(ecal_front_length);
  
  B_field_IsInitialized = false ;
  
  initializeMaterials();


  CreateTree::Instance() -> inputE1Thick 		= ecal_front_length; 

  CreateTree::Instance() -> inputE1Width 		= ecal_front_face; 




}

//---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ---- 



DetectorConstruction::~DetectorConstruction ()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;
  
  
  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------


  //
  //   S U R F A C E S   I N I T I A L I Z A T I O N
  //

  G4LogicalBorderSurface *CrystalSurfaceSide   = NULL;
  G4LogicalBorderSurface *CrystalSurfaceGapf      = NULL;
  G4LogicalBorderSurface *CrystalSurfaceGapr      = NULL;

  G4OpticalSurface *OpWrappingSurface         = NULL;
  G4OpticalSurface *OpCrystalSurface      = NULL;

  
  
  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid * worldS = new G4Box ("worldS", 3 * expHall_x, 3 * expHall_y, 6 * expHall_z) ;
  G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
  G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "worldPV", 0, false, 0, checkOverlaps) ;



  //********************************************
  //  ELECTROMAGNETIC CALORIMETER
  //********************************************





  double pointingAngle = 0.;
  double alveola_thickness = 0.2*mm;

  
  // ECAL solid
  G4Trd* ecalCrystalS_f = new G4Trd("ecalCrystalS_f",0.5*ecal_front_face, 0.5*ecal_rear_face, 0.5*ecal_front_face , 0.5*ecal_rear_face, 0.5*ecal_front_length);

  G4Box* ecalGapS      = new G4Box("ecalGapS",      ecal_det_size*0.5, ecal_det_size*0.5, 0.5*(gap_l-depth) );
  G4Box* ecalDetS      = new G4Box("ecalDetS",      ecal_det_size*0.5, ecal_det_size*0.5, 0.5*(det_l-depth));
  
  
  // ECAL logic
  G4LogicalVolume* ecalCrystalL_f = new G4LogicalVolume(ecalCrystalS_f, EcalMaterial, "ecalCrystalL_f", 0, 0, 0);


  // gaps
  G4LogicalVolume* ecalGapL	 = new G4LogicalVolume(ecalGapS, GaMaterial, "ecalGapL", 0, 0, 0);
  // photodetector
  G4LogicalVolume* ecalDetL      = new G4LogicalVolume(ecalDetS, DeMaterial, "ecalDetL", 0, 0, 0);

  

  // ECAL physical placement
  //const int NECAL_CRYST = 2500;

  const int NECAL_CRYST = ecal_n_cell;

  G4VPhysicalVolume* ecalCrystalP_f[NECAL_CRYST];

  G4VPhysicalVolume* ecalGapP_f[NECAL_CRYST];
  G4VPhysicalVolume* ecalGapP_r[NECAL_CRYST];

  G4VPhysicalVolume* ecalDetP_f[NECAL_CRYST];
  G4VPhysicalVolume* ecalDetP_r[NECAL_CRYST];

  
  char name[60];

  G4double x_pos[NECAL_CRYST];
  G4double y_pos[NECAL_CRYST];
  int nArrayECAL = (int) sqrt(NECAL_CRYST);


  if(surConfig == 0) {
    //Nothing - Crystal completely polished
    
  } else if(surConfig == 1) {
    cout << "Configuring a naked crystal, with only a tiny wrapping" << endl;
    /*-------CRYSTAL SURFACE-------*/
    OpCrystalSurface = new G4OpticalSurface("crystal");
    initializeSurface(OpCrystalSurface, "crystal");
    initializeReflectivitySurface(OpCrystalSurface, "crystal");
  }


  int iCrystal;
  double off;
  for (int iX = 0; iX < nArrayECAL; iX ++)
  {
    for (int iY = 0; iY < nArrayECAL; iY ++)
    {
      
      iCrystal = nArrayECAL*iX + iY;
      
      y_pos[iCrystal] =( iY-nArrayECAL/2)*(ecal_front_face + alveola_thickness);
      off=0;
      if(iY % 2 != 0) off=0.5;
      x_pos[iCrystal] = (iX-off-nArrayECAL/2)*(ecal_front_face + alveola_thickness);    

      G4RotationMatrix* piRotEcal = new G4RotationMatrix;
      piRotEcal->rotateY(-pointingAngle*rad*(iX + 0.5));
      piRotEcal->rotateX(pointingAngle*rad*(iY + 0.5));

      

      cout << " x_pos [" <<setw(3)<<iCrystal << "] = " <<setw(8)<< x_pos[iCrystal] << " :: y_pos[" <<setw(3)<< iCrystal << "] = " <<setw(8)<<y_pos[iCrystal] <<" "<<iX<<" "<<iY<<" "<<off<< endl;
      
      sprintf(name, "ecalCrystalP_f_%d", iCrystal);
      ecalCrystalP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_front_length*0.5), ecalCrystalL_f, name, worldLV, false, 0);


    CrystalSurfaceSide   = new G4LogicalBorderSurface("CrystalSurfaceSide", ecalCrystalP_f[iCrystal], worldPV, OpCrystalSurface);  






      sprintf(name, "ecalGapP_f_%d", iCrystal);
      ecalGapP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], -1.* gap_l*0.5), ecalGapL, name, worldLV, false, 0);

      //    CrystalSurfaceGapf   = new G4LogicalBorderSurface("CrystalSurfaceGapf", ecalCrystalP_f[iCrystal], ecalGapP_f[iCrystal], OpCrystalSurface);  



      ecalGapP_r[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_front_length + gap_l*0.5), ecalGapL, name, worldLV, false, 0);

      //CrystalSurfaceGapr   = new G4LogicalBorderSurface("CrystalSurfaceGapr", ecalCrystalP_f[iCrystal], ecalGapP_r[iCrystal], OpCrystalSurface);  



      sprintf(name, "ecalDetP_f_%d", iCrystal);
      ecalDetP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], -1.* gap_l - det_l*0.5), ecalDetL, name, worldLV, false, 0);

      sprintf(name, "ecalDetP_r_%d", iCrystal);
      ecalDetP_r[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_front_length +gap_l + det_l*0.5), ecalDetL, name, worldLV, false, 0);

    }
  }





  
  
  
  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------
  
  G4Colour white  (1.00, 1.00, 1.00);  // white
  G4Colour gray   (0.50, 0.50, 0.50);  // gray
  G4Colour black  (0.00, 0.00, 0.00);  // black
  G4Colour red    (1.00, 0.00, 0.00);  // red
  G4Colour green  (0.00, 1.00, 0.00);  // green
  G4Colour blue   (0.00, 0.00, 1.00);  // blue
  G4Colour cyan   (0.00, 1.00, 1.00);  // cyan
  G4Colour air    (0.90, 0.94, 1.00);  // cyan
  G4Colour magenta(1.00, 0.00, 1.00);  // magenta 
  G4Colour yellow (1.00, 1.00, 0.00);  // yellow
  G4Colour brass  (0.80, 0.60, 0.40);  // brass
  G4Colour brown  (0.70, 0.40, 0.10);  // brown
  
  G4VisAttributes* VisAttWorld = new G4VisAttributes(black);
  VisAttWorld -> SetVisibility(true) ;
  VisAttWorld -> SetForceWireframe(true) ;
  worldLV -> SetVisAttributes(VisAttWorld) ;
  


  G4VisAttributes* VisCrystalCore = new G4VisAttributes(cyan);
  VisCrystalCore -> SetVisibility(true);
  VisCrystalCore -> SetForceWireframe(true);
  ecalCrystalL_f -> SetVisAttributes(VisCrystalCore);






 
  
  if (B_field_intensity > 0.1 * tesla) ConstructField () ; 
  
  
  
  
  
  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::initializeMaterials ()
{
  //-----------------
  // define materials
  
  WoMaterial = NULL ;
  if      	( world_material == 0 ) WoMaterial = MyMaterials::Vacuum () ;
  else if   ( world_material == 1 ) WoMaterial = MyMaterials::Air () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Wo. material: "<< WoMaterial << G4endl ;
  

  EcalMaterial = NULL ;
  if      ( ecal_material == 1 )  EcalMaterial = MyMaterials::Quartz();
  else if ( ecal_material == 2 )  EcalMaterial = MyMaterials::SiO2();
  else if ( ecal_material == 3 )  EcalMaterial = MyMaterials::SiO2_Ce();
  else if ( ecal_material == 4 )  EcalMaterial = MyMaterials::LuAG_Ce();
  else if ( ecal_material == 5 )  EcalMaterial = MyMaterials::YAG_Ce();
  else if ( ecal_material == 6 )  EcalMaterial = MyMaterials::LSO();
  else if ( ecal_material == 7 )  EcalMaterial = MyMaterials::LYSO();
  else if ( ecal_material == 8 )  EcalMaterial = MyMaterials::LuAG_undoped();
  else if ( ecal_material == 9 )  EcalMaterial = MyMaterials::GAGG_Ce();
  else if ( ecal_material == 10 ) EcalMaterial = MyMaterials::LuAG_Pr();
  else if ( ecal_material == 11 ) EcalMaterial = MyMaterials::PbF2();
  else if ( ecal_material == 12 ) EcalMaterial = MyMaterials::PlasticBC408();
  else if ( ecal_material == 13 ) EcalMaterial = MyMaterials::PlasticBC418();
  else if ( ecal_material == 14 ) EcalMaterial = MyMaterials::PWO();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "ECAL material: "<< EcalMaterial << G4endl ;
  


  GaMaterial = NULL;
  if     ( gap_material == 1 ) GaMaterial = MyMaterials::Air();
  else if( gap_material == 2 ) GaMaterial = MyMaterials::OpticalGrease();
  else if( gap_material == 3 ) GaMaterial = MyMaterials::MeltMount168();
  else if( gap_material == 4 ) GaMaterial = MyMaterials::OpticalGrease155();                            else
    {
      G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
      exit(-1);
    }
  G4cout << "Gap material: " << gap_material << G4endl;

  DeMaterial = NULL;
  if     ( det_material == 1 ) DeMaterial = MyMaterials::Silicon();
  else if( det_material == 2 ) DeMaterial = MyMaterials::Quartz();                                      else if( det_material == 3 ) DeMaterial = MyMaterials::Air();
  else
    {
      G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
      exit(-1);
    }
  G4cout << "Detector material: " << det_material << G4endl;

  
  
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::ConstructField () 
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl ;
  
  static G4TransportationManager * trMgr = G4TransportationManager::GetTransportationManager () ; 
  
  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager * globalFieldMgr = trMgr->GetFieldManager () ;
  
  if( !B_field_IsInitialized )
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522*B_field_intensity,0.0522*B_field_intensity,0.9973*B_field_intensity);   
    
    B_field = new G4UniformMagField (fieldVector) ; 
    globalFieldMgr->SetDetectorField (B_field) ;
    globalFieldMgr->CreateChordFinder (B_field) ;
    globalFieldMgr->GetChordFinder ()->SetDeltaChord (0.005 * mm) ;
    B_field_IsInitialized = true ;
  }
  
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl ;
  return ;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}


// Initialization classes
//
void DetectorConstruction::initializeSurface(G4OpticalSurface *mySurface, string surfaceType)
{
    if(surfaceType == "crystal") {
//         cout << "CRISTALLO " << crystalSurfinish << endl;
        surfinish   = crystalSurfinish;
        RefFile     = cReffile;
        reflectivity    = cReflectivity;
        surrefind   = cSurrefind;
        surtype     = cSurtype;
        specularspike   = cSpecularspike;
        specularlobe    = cSpecularlobe;
        sigmaalpha  = cSigmaalpha;
        backscatter     = cBackscatter;
        lambertian  = cLambertian;
    }



    if(this->surfinish <= 5) {
        G4cout << "Using unified model." << G4endl;
        mySurface -> SetModel(unified);
        switch(this->surtype) {
        case 0:
            mySurface -> SetType(dielectric_metal);
            G4cout << "Surface type: dielectric_metal" << G4endl;
            break;
        case 1:
            mySurface -> SetType(dielectric_dielectric);
            G4cout << "Surface type: dielectric_dielectric" << G4endl;
            break;
        }
    } else if(this->surfinish > 5 && surfaceType == "wrapping") G4cout << "Value not allowed" << G4endl;
    else {
        G4cout << "Using LUT for surface treatment." << G4endl;
        mySurface -> SetModel(LUT);
        mySurface -> SetType(dielectric_LUT);
    }

    switch(this->surfinish) {
    case 0:
        mySurface -> SetFinish(polished);
        G4cout << "Surface finish: polished" << G4endl;
        break;
    case 1:
        mySurface -> SetFinish(polishedfrontpainted);
        G4cout << "Surface finish: polishedfrontpainted" << G4endl;
        break;
    case 2:
        mySurface -> SetFinish(polishedbackpainted);
        G4cout << "Surface finish: polishedbackpainted" << G4endl;
        break;
    case 3:
        mySurface -> SetFinish(ground);
        G4cout << "Surface finish: ground" << G4endl;
        break;
    case 4:
        mySurface -> SetFinish(groundfrontpainted);
        G4cout << "Surface finish: groundfrontpainted" << G4endl;
        break;
    case 5:
        mySurface -> SetFinish(groundbackpainted);
        G4cout << "Surface finish: groundbackpainted" << G4endl;
        break;
    case 17:
        mySurface -> SetFinish(polishedteflonair);
        G4cout << "Surface finish: polishedteflonair" << G4endl;
        break;
    case 18:
        mySurface -> SetFinish(polishedtioair);
        G4cout << "Surface finish: polishedtioair" << G4endl;
        break;
    case 26:
        mySurface -> SetFinish(etchedtioair);
        G4cout << "Surface finish: etchedtioair" << G4endl;
        break;
    case 34:
        mySurface -> SetFinish(groundtioair);
        G4cout << "Surface finish: groundtioair" << G4endl;
        break;
    case 36:
        mySurface -> SetFinish(polishedtyvekair);
        G4cout << "Surface finish: polishedtyvekair" << G4endl;
        break;
    default:
        G4cout << "Surface finish unkown!" << G4endl;
        exit(0);
    }
}


//
// reflectivity
//
void DetectorConstruction::initializeReflectivitySurface(G4OpticalSurface *surface, string surfaceType)
{

    int NumRefl = 0;
    G4double EphotonRefl[1000];
    G4double Refl[1000];
    if(this->RefFile != "none") {
        ifstream myReadFile;
        myReadFile.open(this->RefFile);

        G4cout << "Reflectivities read from file:" << G4endl;
        if(myReadFile.is_open()) {
            while(!myReadFile.eof()) {
                myReadFile >> EphotonRefl[NumRefl];
                if(EphotonRefl[NumRefl] == -1) break;
                myReadFile >> Refl[NumRefl];
                // convert to energy (1eV corresponds to 1239.8 nm: energy [eV]= 1239.8nmeV/lambda[nm])
                EphotonRefl[NumRefl] = 1239.8 / EphotonRefl[NumRefl] * eV;
                NumRefl++;
            }
        } else {
            G4cerr << "<DetectorConstruction> : Could not read file with reflectivities!" << G4endl;
            exit(0);
        }
        myReadFile.close();
    } else {
        EphotonRefl[0] = 0.0001 * eV;
        EphotonRefl[1] = 1.0 * eV;
        EphotonRefl[2] = 4.08 * eV;
        Refl[0] = 0.0; // suppress photons with energy < 1eV (will not be detected)
        Refl[1] = this->crystal_reflectivity;
        Refl[2] = this->crystal_reflectivity;
        NumRefl = 3;


    }
    G4cout << "Reflectivities as a function of the photon energy:" << G4endl;
    for(int i = 0; i < NumRefl; i++) {
        G4cout << i << "   " << EphotonRefl[i] << "   " << Refl[i] << G4endl;
    }


    Ephoton[0] = 0.0001 * eV;
    Ephoton[1] = 1.0 * eV;
    Ephoton[2] = 4.08 * eV;
    G4double RefractiveIndex[3] = {this->surrefind, this->surrefind, this->surrefind};
    G4double SpecularLobe[3]    = {this->specularlobe, this->specularlobe, this->specularlobe};
    G4double SpecularSpike[3]   = {this->specularspike, this->specularspike, this->specularspike};
    G4double Backscatter[3]     = {this->backscatter, this->backscatter, this->backscatter};
    G4double Lambertian[3]      = {this->lambertian, this->lambertian, this->lambertian};
    G4MaterialPropertiesTable *myST = new G4MaterialPropertiesTable();
    G4cout << "Read from config-file: " << G4endl;
    G4cout << "Read SPECULARLOBECONSTANT  : " << SpecularLobe[0] << G4endl;
    G4cout << "Read SPECULARSPIKECONSTANT : " << SpecularSpike[0] << G4endl;
    G4cout << "Read BACKSCATTERCONSTANT   : " << Backscatter[0] << G4endl;
    G4cout << "Read LAMBERTIAN            : " << Lambertian[0] << G4endl;
    G4cout << "Read ref. index            : " << RefractiveIndex[0] << G4endl;

    myST->AddProperty("RINDEX", Ephoton, RefractiveIndex, 3);
    if(this->specularlobe >= 0) {
        G4cout << "Setting SPECULARLOBECONSTANT to : " << SpecularLobe[0] << G4endl;
        myST->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    3);
    }
    if(this->specularspike >= 0) {
        G4cout << "Setting SPECULARSPIKECONSTANT to : " << SpecularSpike[0] << G4endl;
        myST->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   3);
    }
    if(this->backscatter >= 0) {
        G4cout << "Setting BACKSCATTERCONSTANT to : " << Backscatter[0] << G4endl;
        myST->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     3);
    }
    if(this->lambertian >= 0) {
        G4cout << "Setting LAMBERTIAN to : " << Lambertian[0] << G4endl;
        myST->AddProperty("LAMBERTIAN",            Ephoton, Lambertian,      3);
    }


    myST->AddProperty("REFLECTIVITY", EphotonRefl, Refl, NumRefl);
    //try with real and complex index... remove line above and use ones below.
    G4double tyvek_rwavelength[1000] = { 1.24984 * eV, 1.3051 * eV, 1.3776 * eV, 1.45864 * eV, 1.5498 * eV, 1.65312 * eV, 1.7712 * eV, 1.90745 * eV, 2.0664 * eV,
                                         2.25426 * eV, 2.47968 * eV, 2.7552 * eV, 3.0996 * eV, 3.17908 * eV, 3.26274 * eV, 3.35092 * eV, 3.44401 * eV, 3.54241 * eV, 4.13281 * eV
                                       };
    G4double tyvek_rindex[1000] = { 1.37, 1.49, 2.06, 2.61, 2.7, 2.4, 1.8, 1.32, 1.13, 0.907, 0.72, 0.578, 0.456, 0.433, 0.41, 0.382, 0.38, 0.4, 0.276};


    G4double tyvek_cwavelength[1000] = {1.24997 * eV, 1.3051 * eV, 1.34037 * eV, 1.3776 * eV, 1.41696 * eV, 1.45864 * eV, 1.50284 * eV, 1.6 * eV, 1.65312 * eV,
                                        1.74995 * eV, 1.8 * eV, 1.89985 * eV, 1.95005 * eV, 2.05 * eV, 2.1 * eV, 2.19986 * eV, 2.25426 * eV, 2.34997 * eV, 2.4498 * eV,
                                        2.50019 * eV, 2.7 * eV, 2.8 * eV, 2.99986 * eV, 3.19959 * eV, 3.39962 * eV, 3.54241 * eV, 3.69992 * eV, 3.9001 * eV, 4.13281 * eV
                                       };
    G4double tyvek_cindex[1000] = { 9.49, 8.88, 8.49, 8.3, 8.18, 8.22, 8.31, 8.6, 8.62,
                                    8.39, 8.21, 7.82, 7.65, 7.31, 7.15, 6.85, 6.69, 6.42, 6.15,
                                    6.03, 5.58, 5.38, 5.02, 4.71, 4.43, 4.24, 4.06, 3.84, 3.61
                                  };






    //-----------------------------------------------------------------------------------//
    myST->AddProperty("REALRINDEX", tyvek_rwavelength, tyvek_rindex, 19);
    myST->AddProperty("IMAGINARYRINDEX", tyvek_cwavelength, tyvek_cindex, 29);


    G4double air_rwavelength[1000] = {0.1 * eV, 1 * eV, 2 * eV, 3 * eV, 4 * eV};
    G4double air_rindex[1000] = {1, 1, 1, 1, 1};

    G4double air_cwavelength[1000] = {0.1 * eV, 1 * eV, 2 * eV, 3 * eV, 4 * eV};
    G4double air_cindex[1000] = {1, 1, 1, 1, 1};
//     myST->AddProperty("REALRINDEX", air_rwavelength, air_rindex, 5);
//     myST->AddProperty("IMAGINARYRINDEX", air_cwavelength, air_cindex, 5);



    surface->SetMaterialPropertiesTable(myST);
    if(this->sigmaalpha >= 0) surface->SetSigmaAlpha(this->sigmaalpha);


}
