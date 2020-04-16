#ifndef MyMaterials_hh
#define MyMaterials_hh

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "ConfigFile.hh"



class MyMaterials
{
private:
  
public:
  MyMaterials();
  ~MyMaterials();
//  static const G4double kInfinity = 9.0E99; 
  static G4Material* Vacuum();
  static G4Material* Air();
  static G4Material* Water();
  static G4Material* Silicon();
  static G4Material* Aluminium();
  static G4Material* Iron();
  static G4Material* Lead();
  static G4Material* Brass();
  static G4Material* copper();
  static G4Material* Tungsten();
  static G4Material* CopperTungstenAlloy(const G4double& WFrac);
  static G4Material* Quartz();
  static G4Material* SiO2();
  static G4Material* PlasticO2WLS();
  static G4Material* PlasticBC408();
  static G4Material* PlasticBC418();
  static G4Material* OpticalGrease();
  static G4Material* OpticalGrease155();
  static G4Material* MeltMount168();
  static G4Material* LSO();
  static G4Material* LYSO();
  static G4Material* EJ200();
  static G4Material* Acrylic();
  static G4Material* PWO();
  static G4Material* SiO2_Ce();

  static G4Material* LuAG_undoped();
  static G4Material* LuAG_Ce();
  static G4Material* LuAG_Pr();
  static G4Material* YAG_Ce();
  static G4Material* GAGG_Ce();

  static G4Material* PbF2();
  
  static G4double fromNmToEv(G4double wavelength);
  static G4double fromEvToNm(G4double energy);
  static G4double CalculateSellmeier(int size, G4double indexZero, G4double *nVec, G4double *lVec, G4double wavelength);
};

#endif
