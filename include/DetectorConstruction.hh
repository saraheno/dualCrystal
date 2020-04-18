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
// * institutes,nor the agencies providing financial support for this *
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
// $Id: DetectorConstruction.hh,v 1.5 2006-06-29 17:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <iostream>
#include <string>
#include <fstream>
#include <utility>

#include "ConfigFile.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"



class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  //! ctor
  DetectorConstruction();
  DetectorConstruction(const string& configFileName);
  
  //! dtor
  ~DetectorConstruction();
  
  //! construct method
  G4VPhysicalVolume* Construct();
  
  //! other methods
  void initializeMaterials();
  void ConstructField();
  void SetMaxStep(G4double);
  void initializeSurface(G4OpticalSurface *,string);
  void initializeReflectivitySurface(G4OpticalSurface *, string);
  
  G4OpticalSurface* MakeS_TyvekCrystal();
  G4OpticalSurface* MakeS_ESR();
  G4OpticalSurface* MakeS_IdealTyvekCrystal();
  G4OpticalSurface* MakeS_UnifiedTyvek();
  G4OpticalSurface* MakeS_IdealTyvek();
  G4OpticalSurface* MakeS_Polished();
  G4OpticalSurface* MakeS_IdealPolished();
  G4OpticalSurface* MakeS_Mirror();
  G4OpticalSurface* MakeS_IdealMirror();
  G4OpticalSurface* MakeS_IdealWhiteSurface();
  G4OpticalSurface* MakeS_Absorbing();
  G4OpticalSurface* MakeS_Rough();
  G4OpticalSurface* MakeS_SiPM();
  G4OpticalSurface* MakeS_PCBSurface();
  G4OpticalSurface* MakeS_AirAir();




  Fiber* GetFiber() { return &fib; };
  
  
private:
  G4bool    checkOverlaps;
  
  G4double  expHall_x;
  G4double  expHall_y;
  G4double  expHall_z;
  
  G4int    world_material;
 G4double bar_length;      
  G4int    detector;
  

  
  G4int    gap_material;
  G4double gap_l;

  
  G4int    det_material;
  G4double det_l;

    
  G4double depth;
  G4double cryst_dist;


  G4int ecal_material;
  G4double ecal_front_length;

  G4double ecal_front_face;
  G4double ecal_rear_face;

  G4double ecal_det_size;


  
  G4UniformMagField* B_field;
  G4bool   B_field_IsInitialized; 
  G4double B_field_intensity;     // magnetic field, in units of Tesla
  
  Fiber fib;
  
  //Materials
  G4Material* WoMaterial;
  G4Material* EcalMaterial;
  G4Material* CoMaterial;
  G4Material* GaMaterial;
  G4Material* DeMaterial;



  //surfaces
  int surConfig;



  int crystalSurfinish;
  TString cReffile;
  G4double cReflectivity;
  G4double cSurrefind;
  int cSurtype;
  G4double cSpecularspike;
  G4double cSpecularlobe;
  G4double cSigmaalpha;
  G4double cLambertian;
  G4double cBackscatter;


  int wrappingSurfinish;
  TString wReffile;
  G4double wReflectivity;
  G4double wSurrefind;
  int wSurtype;
  G4double wSpecularspike;
  G4double wSpecularlobe;
  G4double wSigmaalpha;
  G4double wLambertian;
  G4double wBackscatter;
  G4double wrapping_reflectivity;


  TString RefFile;
  G4double reflectivity;
  G4double surrefind;
  int surtype;
  G4double specularspike;
  G4double specularlobe;
  G4double sigmaalpha;
  G4double lambertian;
  G4double backscatter;
  G4double crystal_reflectivity;
  int surfinish;
  G4double Ephoton[3];

  
  G4UserLimits* stepLimit;             // pointer to user step limits

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
