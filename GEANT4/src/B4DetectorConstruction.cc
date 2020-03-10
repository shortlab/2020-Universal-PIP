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
// * institutes_Nor the agencies providing financial support for this *
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
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class 

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


#include "G4Region.hh"



// CADMESH //
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4Tet.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4NistManager* nistManager = G4NistManager::Instance();

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 


G4Element* elHe = new  G4Element("Helium", "He", 2, 4.0026*g/mole);
G4Element* elNe = new  G4Element("Neon", "Ne", 10, 20.1797*g/mole);
G4Element* elF = new  G4Element("Fluorine", "F", 9, 18.9984*g/mole);
G4Element* elAr = new  G4Element("Argon", "Ar", 18, 39.948*g/mole);
G4Element* elN = new  G4Element("Nitrogen", "N", 7, 14.0067*g/mole);
G4Element* elKr = new  G4Element("Krypton", "Kr", 36, 83.8*g/mole);
G4Element* elO = new  G4Element("Oxygen", "O", 8, 15.9994*g/mole);
G4Element* elH = new  G4Element("Hydrogen", "H", 1, 1.0079*g/mole);
G4Element* elCl = new  G4Element("Chlorine", "Cl", 17, 35.453*g/mole);
G4Element* elXe = new  G4Element("Xenon", "Xe", 54, 131.293*g/mole);
G4Element* elBr = new  G4Element("Bromine", "Br", 35, 79.904*g/mole);
G4Element* elC = new  G4Element("Carbon", "C", 6, 12.0107*g/mole);
G4Element* elRn = new  G4Element("Radon", "Rn", 86, 222*g/mole);
G4Element* elP = new  G4Element("Phosphorus", "P", 15, 30.9738*g/mole);
G4Element* elI = new  G4Element("Iodine", "I", 53, 126.9045*g/mole);
G4Element* elHg = new  G4Element("Mercury", "Hg", 80, 200.59*g/mole);
G4Element* elS = new  G4Element("Sulfur", "S", 16, 32.065*g/mole);
G4Element* elAs = new  G4Element("Arsenic", "As", 33, 74.9216*g/mole);
G4Element* elSe = new  G4Element("Selenium", "Se", 34, 78.96*g/mole);
G4Element* elZn = new  G4Element("Zinc", "Zn", 30, 65.39*g/mole);
G4Element* elBe = new  G4Element("Beryllium", "Be", 4, 9.0122*g/mole);
G4Element* elAt = new  G4Element("Astatine", "At", 85, 210*g/mole);
G4Element* elAu = new  G4Element("Gold", "Au", 79, 196.9665*g/mole);
G4Element* elTe = new  G4Element("Tellurium", "Te", 52, 127.6*g/mole);
G4Element* elCd = new  G4Element("Cadmium", "Cd", 48, 112.411*g/mole);
G4Element* elIr = new  G4Element("Iridium", "Ir", 77, 192.217*g/mole);
G4Element* elPt = new  G4Element("Platinum", "Pt", 78, 195.078*g/mole);
G4Element* elSb = new  G4Element("Antimony", "Sb", 51, 121.76*g/mole);
G4Element* elOs = new  G4Element("Osmium", "Os", 76, 190.23*g/mole);
G4Element* elPo = new  G4Element("Polonium", "Po", 84, 209*g/mole);
G4Element* elPd = new  G4Element("Palladium", "Pd", 46, 106.42*g/mole);
G4Element* elB = new  G4Element("Boron", "B", 5, 10.811*g/mole);
G4Element* elSi = new  G4Element("Silicon", "Si", 14, 28.0855*g/mole);
G4Element* elFe = new  G4Element("Iron", "Fe", 26, 55.845*g/mole);
G4Element* elGe = new  G4Element("Germanium", "Ge", 32, 72.64*g/mole);
G4Element* elCo = new  G4Element("Cobalt", "Co", 27, 58.9332*g/mole);
G4Element* elW = new  G4Element("Tungsten", "W", 74, 183.84*g/mole);
G4Element* elRe = new  G4Element("Rhenium", "Re", 75, 186.207*g/mole);
G4Element* elCu = new  G4Element("Copper", "Cu", 29, 63.546*g/mole);
G4Element* elMg = new  G4Element("Magnesium", "Mg", 12, 24.305*g/mole);
G4Element* elNi = new  G4Element("Nickel", "Ni", 28, 58.6934*g/mole);
G4Element* elAg = new  G4Element("Silver", "Ag", 47, 107.8682*g/mole);
G4Element* elTa = new  G4Element("Tantalum", "Ta", 73, 180.9479*g/mole);
G4Element* elRh = new  G4Element("Rhodium", "Rh", 45, 102.9055*g/mole);
G4Element* elMn = new  G4Element("Manganese", "Mn", 25, 54.938*g/mole);
G4Element* elPb = new  G4Element("Lead", "Pb", 82, 207.2*g/mole);
G4Element* elRu = new  G4Element("Ruthenium", "Ru", 44, 101.07*g/mole);
G4Element* elSn = new  G4Element("Tin", "Sn", 50, 118.71*g/mole);
G4Element* elBi = new  G4Element("Bismuth", "Bi", 83, 208.9804*g/mole);
G4Element* elTc = new  G4Element("Technetium", "Tc", 43, 98*g/mole);
G4Element* elMo = new  G4Element("Molybdenum", "Mo", 42, 95.94*g/mole);
G4Element* elTi = new  G4Element("Titanium", "Ti", 22, 47.867*g/mole);
G4Element* elHf = new  G4Element("Hafnium", "Hf", 72, 178.49*g/mole);
G4Element* elCr = new  G4Element("Chromium", "Cr", 24, 51.9961*g/mole);
G4Element* elNb = new  G4Element("Niobium", "Nb", 41, 92.9064*g/mole);
G4Element* elV = new  G4Element("Vanadium", "V", 23, 50.9415*g/mole);
//G4Element* elNo = new  G4Element("Nobelium", "No", 102, 259*g/mole);
G4Element* elZr = new  G4Element("Zirconium", "Zr", 40, 91.224*g/mole);
//G4Element* elMd = new  G4Element("Mendelevium", "Md", 101, 258*g/mole);
G4Element* elSc = new  G4Element("Scandium", "Sc", 21, 44.9559*g/mole);
//G4Element* elFm = new  G4Element("Fermium", "Fm", 100, 257*g/mole);
//G4Element* elEs = new  G4Element("Einsteinium", "Es", 99, 252*g/mole);
G4Element* elTh = new  G4Element("Thorium", "Th", 90, 232.0381*g/mole);
//G4Element* elCf = new  G4Element("Californium", "Cf", 98, 251*g/mole);
G4Element* elNp = new  G4Element("Neptunium", "Np", 93, 237*g/mole);
G4Element* elYb = new  G4Element("Ytterbium", "Yb", 70, 173.04*g/mole);
G4Element* elY = new  G4Element("Yttrium", "Y", 39, 88.9059*g/mole);
G4Element* elBk = new  G4Element("Berkelium", "Bk", 97, 247*g/mole);
G4Element* elU = new  G4Element("Uranium", "U", 92, 238.0289*g/mole);
G4Element* elTm = new  G4Element("Thulium", "Tm", 69, 168.9342*g/mole);
G4Element* elGd = new  G4Element("Gadolinium", "Gd", 64, 157.25*g/mole);
G4Element* elCa = new  G4Element("Calcium", "Ca", 20, 40.078*g/mole);
G4Element* elTl = new  G4Element("Thallium", "Tl", 81, 204.3833*g/mole);
G4Element* elEr = new  G4Element("Erbium", "Er", 68, 167.259*g/mole);
G4Element* elPu = new  G4Element("Plutonium", "Pu", 94, 244*g/mole);
G4Element* elHo = new  G4Element("Holmium", "Ho", 67, 164.9303*g/mole);
G4Element* elGa = new  G4Element("Gallium", "Ga", 31, 69.723*g/mole);  
G4Element* elCm = new  G4Element("Curium", "Cm", 96, 247*g/mole);
G4Element* elAl = new  G4Element("Aluminum", "Al", 13, 26.9815*g/mole);
G4Element* elAm = new  G4Element("Americium", "Am", 95, 243*g/mole);
G4Element* elDy = new  G4Element("Dysprosium", "Dy", 66, 162.5*g/mole);
G4Element* elPa = new  G4Element("Protactinium", "Pa", 91, 231.0359*g/mole);
G4Element* elTb = new  G4Element("Terbium", "Tb", 65, 158.9253*g/mole);
G4Element* elIn = new  G4Element("Indium", "In", 49, 114.818*g/mole);
G4Element* elSr = new  G4Element("Strontium", "Sr", 38, 87.62*g/mole);
G4Element* elEu = new  G4Element("Europium", "Eu", 63, 151.964*g/mole);
G4Element* elSm = new  G4Element("Samarium", "Sm", 62, 150.36*g/mole);
G4Element* elPm = new  G4Element("Promethium", "Pm", 61, 145*g/mole);
G4Element* elLa = new  G4Element("Lanthanum", "La", 57, 138.9055*g/mole);
G4Element* elCe = new  G4Element("Cerium", "Ce", 58, 140.116*g/mole);
G4Element* elNd = new  G4Element("Neodymium", "Nd", 60, 144.24*g/mole);
G4Element* elPr = new  G4Element("Praseodymium", "Pr", 59, 140.9077*g/mole);
G4Element* elLu = new  G4Element("Lutetium", "Lu", 71, 174.967*g/mole);
G4Element* elLi = new  G4Element("Lithium", "Li", 3, 6.941*g/mole);
G4Element* elRa = new  G4Element("Radium", "Ra", 88, 226*g/mole);
G4Element* elBa = new  G4Element("Barium", "Ba", 56, 137.327*g/mole);
G4Element* elAc = new  G4Element("Actinium", "Ac", 89, 227*g/mole);
G4Element* elNa = new  G4Element("Sodium", "Na", 11, 22.9897*g/mole);
//G4Element* elLr = new  G4Element("Lawrencium", "Lr", 103, 262*g/mole);
G4Element* elK = new  G4Element("Potassium", "K", 19, 39.0983*g/mole);
G4Element* elRb = new  G4Element("Rubidium", "Rb", 37, 85.4678*g/mole);
G4Element* elFr = new  G4Element("Francium", "Fr", 87, 223*g/mole);
G4Element* elCs = new  G4Element("Cesium", "Cs", 55, 132.9055*g/mole);
//G4Element* elRf = new  G4Element("Rutherfordium", "Rf", 104, 261*g/mole);
//G4Element* elDb = new  G4Element("Dubnium", "Db", 105, 262*g/mole);
//G4Element* elSg = new  G4Element("Seaborgium", "Sg", 106, 266*g/mole);
//G4Element* elBh = new  G4Element("Bohrium", "Bh", 107, 264*g/mole);
//G4Element* elHs = new  G4Element("Hassium", "Hs", 108, 277*g/mole);
//G4Element* elMt = new  G4Element("Meitnerium", "Mt", 109, 268*g/mole);

G4double densityForAll = 4.0;

G4Material* DenseHelium = new G4Material ("DenseHelium", densityForAll*(g/cm3), 1);
G4Material* DenseNeon = new G4Material ("DenseNeon", densityForAll*(g/cm3), 1);
G4Material* DenseFluorine = new G4Material ("DenseFluorine", densityForAll*(g/cm3), 1);
G4Material* DenseArgon = new G4Material ("DenseArgon", densityForAll*(g/cm3), 1);
G4Material* DenseNitrogen = new G4Material ("DenseNitrogen", densityForAll*(g/cm3), 1);
G4Material* DenseKrypton = new G4Material ("DenseKrypton", densityForAll*(g/cm3), 1);
G4Material* DenseOxygen = new G4Material ("DenseOxygen", densityForAll*(g/cm3), 1);
G4Material* DenseHydrogen = new G4Material ("DenseHydrogen", densityForAll*(g/cm3), 1);
G4Material* DenseChlorine = new G4Material ("DenseChlorine", densityForAll*(g/cm3), 1);
G4Material* DenseXenon = new G4Material ("DenseXenon", densityForAll*(g/cm3), 1);
G4Material* DenseBromine = new G4Material ("DenseBromine", densityForAll*(g/cm3), 1);
G4Material* DenseCarbon = new G4Material ("DenseCarbon", densityForAll*(g/cm3), 1);
G4Material* DenseRadon = new G4Material ("DenseRadon", densityForAll*(g/cm3), 1);
G4Material* DensePhosphorus = new G4Material ("DensePhosphorus", densityForAll*(g/cm3), 1);
G4Material* DenseIodine = new G4Material ("DenseIodine", densityForAll*(g/cm3), 1);
G4Material* DenseMercury = new G4Material ("DenseMercury", densityForAll*(g/cm3), 1);
G4Material* DenseSulfur = new G4Material ("DenseSulfur", densityForAll*(g/cm3), 1);
G4Material* DenseArsenic = new G4Material ("DenseArsenic", densityForAll*(g/cm3), 1);
G4Material* DenseSelenium = new G4Material ("DenseSelenium", densityForAll*(g/cm3), 1);
G4Material* DenseZinc = new G4Material ("DenseZinc", densityForAll*(g/cm3), 1);
G4Material* DenseBeryllium = new G4Material ("DenseBeryllium", densityForAll*(g/cm3), 1);
G4Material* DenseAstatine = new G4Material ("DenseAstatine", densityForAll*(g/cm3), 1);
G4Material* DenseGold = new G4Material ("DenseGold", densityForAll*(g/cm3), 1);
G4Material* DenseTellurium = new G4Material ("DenseTellurium", densityForAll*(g/cm3), 1);
G4Material* DenseCadmium = new G4Material ("DenseCadmium", densityForAll*(g/cm3), 1);
G4Material* DenseIridium = new G4Material ("DenseIridium", densityForAll*(g/cm3), 1);
G4Material* DensePlatinum = new G4Material ("DensePlatinum", densityForAll*(g/cm3), 1);
G4Material* DenseAntimony = new G4Material ("DenseAntimony", densityForAll*(g/cm3), 1);
G4Material* DenseOsmium = new G4Material ("DenseOsmium", densityForAll*(g/cm3), 1);
G4Material* DensePolonium = new G4Material ("DensePolonium", densityForAll*(g/cm3), 1);
G4Material* DensePalladium = new G4Material ("DensePalladium", densityForAll*(g/cm3), 1);
G4Material* DenseBoron = new G4Material ("DenseBoron", densityForAll*(g/cm3), 1);
G4Material* DenseSilicon = new G4Material ("DenseSilicon", densityForAll*(g/cm3), 1);
G4Material* DenseIron = new G4Material ("DenseIron", densityForAll*(g/cm3), 1);
G4Material* DenseGermanium = new G4Material ("DenseGermanium", densityForAll*(g/cm3), 1);
G4Material* DenseCobalt = new G4Material ("DenseCobalt", densityForAll*(g/cm3), 1);
G4Material* DenseTungsten = new G4Material ("DenseTungsten", densityForAll*(g/cm3), 1);
G4Material* DenseRhenium = new G4Material ("DenseRhenium", densityForAll*(g/cm3), 1);
G4Material* DenseCopper = new G4Material ("DenseCopper", densityForAll*(g/cm3), 1);
G4Material* DenseMagnesium = new G4Material ("DenseMagnesium", densityForAll*(g/cm3), 1);
G4Material* DenseNickel = new G4Material ("DenseNickel", densityForAll*(g/cm3), 1);
G4Material* DenseSilver = new G4Material ("DenseSilver", densityForAll*(g/cm3), 1);
G4Material* DenseTantalum = new G4Material ("DenseTantalum", densityForAll*(g/cm3), 1);
G4Material* DenseRhodium = new G4Material ("DenseRhodium", densityForAll*(g/cm3), 1);
G4Material* DenseManganese = new G4Material ("DenseManganese", densityForAll*(g/cm3), 1);
G4Material* DenseLead = new G4Material ("DenseLead", densityForAll*(g/cm3), 1);
G4Material* DenseRuthenium = new G4Material ("DenseRuthenium", densityForAll*(g/cm3), 1);
G4Material* DenseTin = new G4Material ("DenseTin", densityForAll*(g/cm3), 1);
G4Material* DenseBismuth = new G4Material ("DenseBismuth", densityForAll*(g/cm3), 1);
G4Material* DenseTechnetium = new G4Material ("DenseTechnetium", densityForAll*(g/cm3), 1);
G4Material* DenseMolybdenum = new G4Material ("DenseMolybdenum", densityForAll*(g/cm3), 1);
G4Material* DenseTitanium = new G4Material ("DenseTitanium", densityForAll*(g/cm3), 1);
G4Material* DenseHafnium = new G4Material ("DenseHafnium", densityForAll*(g/cm3), 1);
G4Material* DenseChromium = new G4Material ("DenseChromium", densityForAll*(g/cm3), 1);
G4Material* DenseNiobium = new G4Material ("DenseNiobium", densityForAll*(g/cm3), 1);
G4Material* DenseVanadium = new G4Material ("DenseVanadium", densityForAll*(g/cm3), 1);
//G4Material* DenseNobelium = new G4Material ("DenseNobelium", densityForAll*(g/cm3), 1);
G4Material* DenseZirconium = new G4Material ("DenseZirconium", densityForAll*(g/cm3), 1);
//G4Material* DenseMendelevium = new G4Material ("DenseMendelevium", densityForAll*(g/cm3), 1);
G4Material* DenseScandium = new G4Material ("DenseScandium", densityForAll*(g/cm3), 1);
//G4Material* DenseFermium = new G4Material ("DenseFermium", densityForAll*(g/cm3), 1);
//G4Material* DenseEinsteinium = new G4Material ("DenseEinsteinium", densityForAll*(g/cm3), 1);
G4Material* DenseThorium = new G4Material ("DenseThorium", densityForAll*(g/cm3), 1);
//G4Material* DenseCalifornium = new G4Material ("DenseCalifornium", densityForAll*(g/cm3), 1);
G4Material* DenseNeptunium = new G4Material ("DenseNeptunium", densityForAll*(g/cm3), 1);
G4Material* DenseYtterbium = new G4Material ("DenseYtterbium", densityForAll*(g/cm3), 1);
G4Material* DenseYttrium = new G4Material ("DenseYttrium", densityForAll*(g/cm3), 1);
G4Material* DenseBerkelium = new G4Material ("DenseBerkelium", densityForAll*(g/cm3), 1);
G4Material* DenseUranium = new G4Material ("DenseUranium", densityForAll*(g/cm3), 1);
G4Material* DenseThulium = new G4Material ("DenseThulium", densityForAll*(g/cm3), 1);
G4Material* DenseGadolinium = new G4Material ("DenseGadolinium", densityForAll*(g/cm3), 1);
G4Material* DenseCalcium = new G4Material ("DenseCalcium", densityForAll*(g/cm3), 1);
G4Material* DenseThallium = new G4Material ("DenseThallium", densityForAll*(g/cm3), 1);
G4Material* DenseErbium = new G4Material ("DenseErbium", densityForAll*(g/cm3), 1);
G4Material* DensePlutonium = new G4Material ("DensePlutonium", densityForAll*(g/cm3), 1);
G4Material* DenseHolmium = new G4Material ("DenseHolmium", densityForAll*(g/cm3), 1);
G4Material* DenseGallium = new G4Material ("DenseGallium", densityForAll*(g/cm3), 1);
G4Material* DenseCurium = new G4Material ("DenseCurium", densityForAll*(g/cm3), 1);
G4Material* DenseAluminum = new G4Material ("DenseAluminum", densityForAll*(g/cm3), 1);
G4Material* DenseAmericium = new G4Material ("DenseAmericium", densityForAll*(g/cm3), 1);
G4Material* DenseDysprosium = new G4Material ("DenseDysprosium", densityForAll*(g/cm3), 1);
G4Material* DenseProtactinium = new G4Material ("DenseProtactinium", densityForAll*(g/cm3), 1);
G4Material* DenseTerbium = new G4Material ("DenseTerbium", densityForAll*(g/cm3), 1);
G4Material* DenseIndium = new G4Material ("DenseIndium", densityForAll*(g/cm3), 1);
G4Material* DenseStrontium = new G4Material ("DenseStrontium", densityForAll*(g/cm3), 1);
G4Material* DenseEuropium = new G4Material ("DenseEuropium", densityForAll*(g/cm3), 1);
G4Material* DenseSamarium = new G4Material ("DenseSamarium", densityForAll*(g/cm3), 1);
G4Material* DensePromethium = new G4Material ("DensePromethium", densityForAll*(g/cm3), 1);
G4Material* DenseLanthanum = new G4Material ("DenseLanthanum", densityForAll*(g/cm3), 1);
G4Material* DenseCerium = new G4Material ("DenseCerium", densityForAll*(g/cm3), 1);
G4Material* DenseNeodymium = new G4Material ("DenseNeodymium", densityForAll*(g/cm3), 1);
G4Material* DensePraseodymium = new G4Material ("DensePraseodymium", densityForAll*(g/cm3), 1);
G4Material* DenseLutetium = new G4Material ("DenseLutetium", densityForAll*(g/cm3), 1);
G4Material* DenseLithium = new G4Material ("DenseLithium", densityForAll*(g/cm3), 1);
G4Material* DenseRadium = new G4Material ("DenseRadium", densityForAll*(g/cm3), 1);
G4Material* DenseBarium = new G4Material ("DenseBarium", densityForAll*(g/cm3), 1);
G4Material* DenseActinium = new G4Material ("DenseActinium", densityForAll*(g/cm3), 1);
G4Material* DenseSodium = new G4Material ("DenseSodium", densityForAll*(g/cm3), 1);
//G4Material* DenseLawrencium = new G4Material ("DenseLawrencium", densityForAll*(g/cm3), 1);
G4Material* DensePotassium = new G4Material ("DensePotassium", densityForAll*(g/cm3), 1);
G4Material* DenseRubidium = new G4Material ("DenseRubidium", densityForAll*(g/cm3), 1);
G4Material* DenseFrancium = new G4Material ("DenseFrancium", densityForAll*(g/cm3), 1);
G4Material* DenseCesium = new G4Material ("DenseCesium", densityForAll*(g/cm3), 1);
//G4Material* DenseRutherfordium = new G4Material ("DenseRutherfordium", densityForAll*(g/cm3), 1);
//G4Material* DenseDubnium = new G4Material ("DenseDubnium", densityForAll*(g/cm3), 1);
//G4Material* DenseSeaborgium = new G4Material ("DenseSeaborgium", densityForAll*(g/cm3), 1);
//G4Material* DenseBohrium = new G4Material ("DenseBohrium", densityForAll*(g/cm3), 1);
//G4Material* DenseHassium = new G4Material ("DenseHassium", densityForAll*(g/cm3), 1);
//G4Material* DenseMeitnerium = new G4Material ("DenseMeitnerium", densityForAll*(g/cm3), 1);

G4int natoms = 1;


DenseHelium->AddElement(elHe, natoms=1);
DenseNeon->AddElement(elNe, natoms=1);
DenseFluorine->AddElement(elF, natoms=1);
DenseArgon->AddElement(elAr, natoms=1);
DenseNitrogen->AddElement(elN, natoms=1);
DenseKrypton->AddElement(elKr, natoms=1);
DenseOxygen->AddElement(elO, natoms=1);
DenseHydrogen->AddElement(elH, natoms=1);
DenseChlorine->AddElement(elCl, natoms=1);
DenseXenon->AddElement(elXe, natoms=1);
DenseBromine->AddElement(elBr, natoms=1);
DenseCarbon->AddElement(elC, natoms=1);
DenseRadon->AddElement(elRn, natoms=1);
DensePhosphorus->AddElement(elP, natoms=1);
DenseIodine->AddElement(elI, natoms=1);
DenseMercury->AddElement(elHg, natoms=1);
DenseSulfur->AddElement(elS, natoms=1);
DenseArsenic->AddElement(elAs, natoms=1);
DenseSelenium->AddElement(elSe, natoms=1);
DenseZinc->AddElement(elZn, natoms=1);
DenseBeryllium->AddElement(elBe, natoms=1);
DenseAstatine->AddElement(elAt, natoms=1);
DenseGold->AddElement(elAu, natoms=1);
DenseTellurium->AddElement(elTe, natoms=1);
DenseCadmium->AddElement(elCd, natoms=1);
DenseIridium->AddElement(elIr, natoms=1);
DensePlatinum->AddElement(elPt, natoms=1);
DenseAntimony->AddElement(elSb, natoms=1);
DenseOsmium->AddElement(elOs, natoms=1);
DensePolonium->AddElement(elPo, natoms=1);
DensePalladium->AddElement(elPd, natoms=1);
DenseBoron->AddElement(elB, natoms=1);
DenseSilicon->AddElement(elSi, natoms=1);
DenseIron->AddElement(elFe, natoms=1);
DenseGermanium->AddElement(elGe, natoms=1);
DenseCobalt->AddElement(elCo, natoms=1);
DenseTungsten->AddElement(elW, natoms=1);
DenseRhenium->AddElement(elRe, natoms=1);
DenseCopper->AddElement(elCu, natoms=1);
DenseMagnesium->AddElement(elMg, natoms=1);
DenseNickel->AddElement(elNi, natoms=1);
DenseSilver->AddElement(elAg, natoms=1);
DenseTantalum->AddElement(elTa, natoms=1);
DenseRhodium->AddElement(elRh, natoms=1);
DenseManganese->AddElement(elMn, natoms=1);
DenseLead->AddElement(elPb, natoms=1);
DenseRuthenium->AddElement(elRu, natoms=1);
DenseTin->AddElement(elSn, natoms=1);
DenseBismuth->AddElement(elBi, natoms=1);
DenseTechnetium->AddElement(elTc, natoms=1);
DenseMolybdenum->AddElement(elMo, natoms=1);
DenseTitanium->AddElement(elTi, natoms=1);
DenseHafnium->AddElement(elHf, natoms=1);
DenseChromium->AddElement(elCr, natoms=1);
DenseNiobium->AddElement(elNb, natoms=1);
DenseVanadium->AddElement(elV, natoms=1);
//DenseNobelium->AddElement(elNo, natoms=1);
DenseZirconium->AddElement(elZr, natoms=1);
//DenseMendelevium->AddElement(elMd, natoms=1);
DenseScandium->AddElement(elSc, natoms=1);
//DenseFermium->AddElement(elFm, natoms=1);
//DenseEinsteinium->AddElement(elEs, natoms=1);
DenseThorium->AddElement(elTh, natoms=1);
//DenseCalifornium->AddElement(elCf, natoms=1);
DenseNeptunium->AddElement(elNp, natoms=1);
DenseYtterbium->AddElement(elYb, natoms=1);
DenseYttrium->AddElement(elY, natoms=1);
DenseBerkelium->AddElement(elBk, natoms=1);
DenseUranium->AddElement(elU, natoms=1);
DenseThulium->AddElement(elTm, natoms=1);
DenseGadolinium->AddElement(elGd, natoms=1);
DenseCalcium->AddElement(elCa, natoms=1);
DenseThallium->AddElement(elTl, natoms=1);
DenseErbium->AddElement(elEr, natoms=1);
DensePlutonium->AddElement(elPu, natoms=1);
DenseHolmium->AddElement(elHo, natoms=1);
DenseGallium->AddElement(elGa, natoms=1);
DenseCurium->AddElement(elCm, natoms=1);
DenseAluminum->AddElement(elAl, natoms=1);
DenseAmericium->AddElement(elAm, natoms=1);
DenseDysprosium->AddElement(elDy, natoms=1);
DenseProtactinium->AddElement(elPa, natoms=1);
DenseTerbium->AddElement(elTb, natoms=1);
DenseIndium->AddElement(elIn, natoms=1);
DenseStrontium->AddElement(elSr, natoms=1);
DenseEuropium->AddElement(elEu, natoms=1);
DenseSamarium->AddElement(elSm, natoms=1);
DensePromethium->AddElement(elPm, natoms=1);
DenseLanthanum->AddElement(elLa, natoms=1);
DenseCerium->AddElement(elCe, natoms=1);
DenseNeodymium->AddElement(elNd, natoms=1);
DensePraseodymium->AddElement(elPr, natoms=1);
DenseLutetium->AddElement(elLu, natoms=1);
DenseLithium->AddElement(elLi, natoms=1);
DenseRadium->AddElement(elRa, natoms=1);
DenseBarium->AddElement(elBa, natoms=1);
DenseActinium->AddElement(elAc, natoms=1);
DenseSodium->AddElement(elNa, natoms=1);
//DenseLawrencium->AddElement(elLr, natoms=1);
DensePotassium->AddElement(elK, natoms=1);
DenseRubidium->AddElement(elRb, natoms=1);
DenseFrancium->AddElement(elFr, natoms=1);
DenseCesium->AddElement(elCs, natoms=1);
//DenseRutherfordium->AddElement(elRf, natoms=1);
//DenseDubnium->AddElement(elDb, natoms=1);
//DenseSeaborgium->AddElement(elSg, natoms=1);
//DenseBohrium->AddElement(elBh, natoms=1);
//DenseHassium->AddElement(elHs, natoms=1);
//DenseMeitnerium->AddElement(elMt, natoms=1);



  
G4String materialNames[270] = {
"G4_Li",
"G4_N-PENTANE",
"G4_N-HEXANE",
"G4_N-HEPTANE",
"G4_OCTANE",
"G4_DIETHYL_ETHER",
"G4_CYCLOHEXANE",
"G4_ETHYL_ALCOHOL",
"G4_ACETONE",
"G4_METHANOL",
"G4_N-PROPYL_ALCOHOL",
"G4_lN2",
"G4_N-BUTYL_ALCOHOL",
"G4_OCTADECANOL",
"G4_LITHIUM_HYDRIDE",
"G4_K",
"G4_TOLUENE",
"G4_XYLENE",
"G4_BENZENE",
"G4_POLYPROPYLENE",
"G4_RUBBER_BUTYL",
"G4_RUBBER_NATURAL",
"G4_PARAFFIN",
"G4_POLYETHYLENE",
"G4_N,N-DIMETHYL_FORMAMIDE",
"G4_ADIPOSE_TISSUE_ICRP",
"G4_FREON-13",
"G4_STILBENE",
"G4_Na",
"G4_PYRIDINE",
"G4_MIX_D_WAX",
"G4_Fr",
"G4_MS20_TISSUE",
"G4_TISSUE_SOFT_ICRU-4",
"G4_WATER",
"G4_DNA_ADENINE",
"G4_DNA_GUANINE",
"G4_DNA_CYTOSINE",
"G4_DNA_THYMINE",
"G4_DNA_URACIL",
"G4_ANILINE",
"G4_FERROUS_SULFATE",
"G4_CERIC_SULFATE",
"G4_TISSUE_SOFT_ICRP",
"G4_PLASTIC_SC_VINYLTOLUENE",
"G4_BRAIN_ICRP",
"G4_LUNG_ICRP",
"G4_MUSCLE_STRIATED_ICRU",
"G4_TESTIS_ICRP",
"G4_M3_WAX",
"G4_MUSCLE_SKELETAL_ICRP",
"G4_BLOOD_ICRP",
"G4_POLYSTYRENE",
"G4_EYE_LENS_ICRP",
"G4_MUSCLE_WITHOUT_SUCROSE",
"G4_TRIETHYL_PHOSPHATE",
"G4_NYLON-8062",
"G4_SKIN_ICRP",
"G4_AMBER",
"G4_DIMETHYL_SULFOXIDE",
"G4_CHLOROBENZENE",
"G4_MUSCLE_WITH_SUCROSE",
"G4_FREON-12",
"G4_POLYVINYL_BUTYRAL",
"G4_A-150_TISSUE",
"G4_ETHYL_CELLULOSE",
"G4_NYLON-6-6",
"G4_NYLON-6-10",
"G4_lO2",
"G4_NAPHTHALENE",
"G4_POLYACRYLONITRILE",
"G4_LITHIUM_AMIDE",
"G4_PLEXIGLASS",
"G4_POLYVINYL_ACETATE",
"G4_LUCITE",
"G4_NITROBENZENE",
"G4_CELLULOSE_BUTYRATE",
"G4_POLYCARBONATE",
"G4_DICHLORODIETHYL_ETHER",
"G4_RUBBER_NEOPRENE",
"G4_VALINE",
"G4_NEOPRENE",
"G4_THYMINE",
"G4_1,2-DICHLOROETHANE",
"G4_TERPHENYL",
"G4_BAKELITE",
"G4_POLYVINYL_PYRROLIDONE",
"G4_GLYCEROL",
"G4_ANTHRACENE",
"G4_GEL_PHOTO_EMULSION",
"G4_POLYCHLOROSTYRENE",
"G4_POLYVINYL_ALCOHOL",
"G4_POLYVINYL_CHLORIDE",
"G4_1,2-DICHLOROBENZENE",
"G4_CR39",
"G4_URACIL",
"G4_UREA",
"G4_lAr",
"G4_MYLAR",
"G4_DACRON",
"G4_ALANINE",
"G4_CELLULOSE_CELLOPHANE",
"G4_KAPTON",
"G4_NYLON-11_RILSAN",
"G4_POLYOXYMETHYLENE",
"G4_KEVLAR",
"G4_B-100_BONE",
"G4_GLUTAMINE",
"G4_TRICHLOROETHYLENE",
"G4_CHLOROFORM",
"G4_CELLULOSE_NITRATE",
"G4_FREON-13B1",
"G4_Rb",
"G4_Ca",
"G4_CYTOSINE",
"G4_SUCROSE",
"G4_CARBON_TETRACHLORIDE",
"G4_ADENINE",
"G4_TETRACHLOROETHYLENE",
"G4_POLYVINYLIDENE_CHLORIDE",
"G4_GRAPHITE_POROUS",
"G4_Mg",
"G4_C-552",
"G4_POLYVINYLIDENE_FLUORIDE",
"G4_FREON-12B2",
"G4_FREON-13I1",
"G4_VITON",
"G4_BORON_OXIDE",
"G4_Be",
"G4_BONE_COMPACT_ICRU",
"G4_Cs",
"G4_BONE_CORTICAL_ICRP",
"G4_C",
"G4_S",
"G4_LITHIUM_OXIDE",
"G4_POLYTRIFLUOROCHLOROETHYLENE",
"G4_LITHIUM_CARBONATE",
"G4_P",
"G4_GUANINE",
"G4_TEFLON",
"G4_GRAPHITE",
"G4_Pyrex_Glass",
"G4_SODIUM_NITRATE",
"G4_SODIUM_MONOXIDE",
"G4_CONCRETE",
"G4_GYPSUM",
"G4_POTASSIUM_OXIDE",
"G4_SILICON_DIOXIDE",
"G4_Si",
"G4_B",
"G4_GLASS_PLATE",
"G4_TUNGSTEN_HEXAFLUORIDE",
"G4_lKr",
"G4_LITHIUM_TETRABORATE",
"G4_BORON_CARBIDE",
"G4_MAGNESIUM_TETRABORATE",
"G4_SODIUM_CARBONATE",
"G4_Sr",
"G4_LITHIUM_FLUORIDE",
"G4_Al",
"G4_CALCIUM_CARBONATE",
"G4_lXe",
"G4_MAGNESIUM_CARBONATE",
"G4_CALCIUM_SULFATE",
"G4_Sc",
"G4_MAGNESIUM_FLUORIDE",
"G4_BERYLLIUM_OXIDE",
"G4_lBr",
"G4_POTASSIUM_IODIDE",
"G4_CALCIUM_FLUORIDE",
"G4_CALCIUM_OXIDE",
"G4_LITHIUM_IODIDE",
"G4_Ba",
"G4_MAGNESIUM_OXIDE",
"G4_SODIUM_IODIDE",
"G4_PHOTO_EMULSION",
"G4_ALUMINUM_OXIDE",
"G4_CESIUM_FLUORIDE",
"G4_TITANIUM_DIOXIDE",
"G4_Y",
"G4_Se",
"G4_BARIUM_SULFATE",
"G4_CESIUM_IODIDE",
"G4_Ti",
"G4_BARIUM_FLUORIDE",
"G4_I",
"G4_Ra",
"G4_FERRIC_OXIDE",
"G4_Eu",
"G4_GALLIUM_ARSENIDE",
"G4_Ge",
"G4_SILVER_CHLORIDE",
"G4_FERROUS_OXIDE",
"G4_As",
"G4_LANTHANUM_OXYSULFIDE",
"G4_Ga",
"G4_SILVER_IODIDE",
"G4_CALCIUM_TUNGSTATE",
"G4_V",
"G4_La",
"G4_CADMIUM_TELLURIDE",
"G4_GLASS_LEAD",
"G4_Te",
"G4_LANTHANUM_OXYBROMIDE",
"G4_MERCURIC_IODIDE",
"G4_SILVER_HALIDES",
"G4_SILVER_BROMIDE",
"G4_Zr",
"G4_Ce",
"G4_Sb",
"G4_Pr",
"G4_Yb",
"G4_Nd",
"G4_THALLIUM_CHLORIDE",
"G4_BGO",
"G4_Zn",
"G4_FERROBORIDE",
"G4_Cr",
"G4_Pm",
"G4_In",
"G4_Sn",
"G4_Mn",
"G4_GADOLINIUM_OXYSULFIDE",
"G4_Sm",
"G4_Fe",
"G4_CADMIUM_TUNGSTATE",
"G4_Gd",
"G4_STAINLESS-STEEL",
"G4_Tb",
"G4_PbWO4",
"G4_BRASS",
"G4_Dy",
"G4_Nb",
"G4_Cd",
"G4_Ho",
"G4_BRONZE",
"G4_Co",
"G4_Ni",
"G4_Cu",
"G4_Er",
"G4_Po",
"G4_At",
"G4_Tm",
"G4_LEAD_OXIDE",
"G4_Bi",
"G4_Lu",
"G4_Ac",
"G4_Mo",
"G4_Ag",
"G4_URANIUM_OXIDE",
"G4_URANIUM_DICARBIDE",
"G4_Pb",
"G4_Tc",
"G4_Tl",
"G4_Th",
"G4_Pd",
"G4_Ru",
"G4_Rh",
"G4_Hf",
"G4_Hg",
"G4_URANIUM_MONOCARBIDE",
"G4_Pa",
"G4_Ta",
"G4_U",
"G4_W",
"G4_Au",
"G4_Re",
"G4_Pt",
"G4_Ir",
"G4_Os"};



G4Material* elementMaterials[97] = {
DenseHelium,
DenseNeon,
DenseFluorine,
DenseArgon,
DenseNitrogen,
DenseKrypton,
DenseOxygen,
DenseHydrogen,
DenseChlorine,
DenseXenon,
DenseBromine,
DenseCarbon,
DenseRadon,
DensePhosphorus,
DenseIodine,
DenseMercury,
DenseSulfur,
DenseArsenic,
DenseSelenium,
DenseZinc,
DenseBeryllium,
DenseAstatine,
DenseGold,
DenseTellurium,
DenseCadmium,
DenseIridium,
DensePlatinum,
DenseAntimony,
DenseOsmium,
DensePolonium,
DensePalladium,
DenseBoron,
DenseSilicon,
DenseIron,
DenseGermanium,
DenseCobalt,
DenseTungsten,
DenseRhenium,
DenseCopper,
DenseMagnesium,
DenseNickel,
DenseSilver,
DenseTantalum,
DenseRhodium,
DenseManganese,
DenseLead,
DenseRuthenium,
DenseTin,
DenseBismuth,
DenseTechnetium,
DenseMolybdenum,
DenseTitanium,
DenseHafnium,
DenseChromium,
DenseNiobium,
DenseVanadium,
DenseZirconium,
DenseScandium,
DenseThorium,
DenseNeptunium,
DenseYtterbium,
DenseYttrium,
DenseBerkelium,
DenseUranium,
DenseThulium,
DenseGadolinium,
DenseCalcium,
DenseThallium,
DenseErbium,
DensePlutonium,
DenseHolmium,
DenseGallium,
DenseCurium,
DenseAluminum,
DenseAmericium,
DenseDysprosium,
DenseProtactinium,
DenseTerbium,
DenseIndium,
DenseStrontium,
DenseEuropium,
DenseSamarium,
DensePromethium,
DenseLanthanum,
DenseCerium,
DenseNeodymium,
DensePraseodymium,
DenseLutetium,
DenseLithium,
DenseRadium,
DenseBarium,
DenseActinium,
DenseSodium,
DensePotassium,
DenseRubidium,
DenseFrancium,
DenseCesium};

G4Material* G4_Galactic = nistManager->FindOrBuildMaterial("G4_Galactic");



G4double A; G4int Z; G4double d;


d = 2.165*g/cm3;
G4Material* NaCl = new G4Material("Sodium Chloride",d,2);
NaCl ->AddElement(elNa,1); // 1 is the number of Sodium atoms in NaCl 
NaCl ->AddElement(elCl,1);

// cm
G4double World_Half_Thickness = (100.0/2.0);

G4Box* world_solid
= new G4Box("World", World_Half_Thickness*cm, World_Half_Thickness*cm, World_Half_Thickness*cm); // its size
                         
G4LogicalVolume* world_logical
= new G4LogicalVolume(
    world_solid,           // its solid
    G4_Galactic,  // its material
    "World");         // its name
                                   
G4VPhysicalVolume* world_physical
= new G4PVPlacement(
    0,                // no rotation
    G4ThreeVector(),  // at (0,0,0)
    world_logical,          // its logical volume                         
    "World",          // its name
    0,                // its mother  volume
    false,            // no boolean operation
    0,                // copy number
    fCheckOverlaps);  // checking overlaps 

// *******
// 22NaCl source information - cylinder with radius of 1 mm and thickness of 0.6249760658373185  nm to get 10 uCi if 100% purity

// nm
G4double NaCl_Source_thickness = 0.624976;
// cm
G4double NaCl_Source_radius = 0.1;

G4RotationMatrix* myRotation = new G4RotationMatrix();
    myRotation->rotateX(0.0*deg);
    myRotation->rotateY(0.0*deg);
    myRotation->rotateZ(90.0*deg);

G4Tubs* NaCl_Source_solid
= new G4Tubs("NaCl_Source", 0.0*mm, NaCl_Source_radius*cm, (NaCl_Source_thickness/2.0)*nm, 0.0*deg,  360.0*deg);
G4LogicalVolume* NaCl_Source_logical
= new G4LogicalVolume(
    NaCl_Source_solid,      // its solid
    NaCl,                  // its material
    "NaCl_Source");         // its name
G4VPhysicalVolume* NaCl_Source_physical
= new G4PVPlacement(
    myRotation,                // no rotation
    G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm),  
    NaCl_Source_logical,          // its logical volume                         
    "NaCl_Source",          // its name
    world_logical,   // its mother  volume
    false,            // no boolean operation
    0,                // copy number
    fCheckOverlaps);  // checking overlaps 

// *******
// Sample Material to be characterized information - just let it be really thick so that they all anihilate inside, hopefully

// If unning for the whole material database

G4String material_name =  materialNames[(material_num)];
G4Material* Sample_Material = nistManager->FindOrBuildMaterial(material_name);


// If running for the 4 g/cc materials
//G4Material* Sample_Material = elementMaterials[material_num];
//G4String material_name =  Sample_Material->GetName();

cout << "NOW RUNNING! " << material_name << endl;



for( unsigned int mn = 0; mn < sizeof(elementMaterials); mn = mn + 1 ) {
G4Material* sm = elementMaterials[mn];

G4String material_name =  sm->GetName();

cout << material_name << endl;

std::ostringstream commandOS;
commandOS << "MaterialInfo_1000000_1mm_diameter_10uCi_22Na_4g_cc_elements.txt";
std::ofstream ofile;
ofile.open (G4String(commandOS.str()), std::ios::out | std::ios::app);     // ascii file   
ofile << " " << material_name << " " << (sm->GetDensity())/(g/cm3)  << " " << sm->GetElectronDensity()/(1.0/cm3) << " " << sm->GetRadlen()/cm << " " << sm->GetIonisation()->GetMeanExcitationEnergy()/eV << " " << sm->GetIonisation()->GetPlasmaEnergy()/eV << " " << sm->GetIonisation()->GetAdjustmentFactor() << " " << sm->GetIonisation()->GetEnergy1fluct()/eV << " " << sm->GetIonisation()->GetEnergy2fluct()/eV << " " << sm->GetIonisation()->GetEnergy0fluct()/eV << " " << sm->GetIonisation()->GetRateionexcfluct()/eV << " " << sm->GetIonisation()->GetZeffective() << " " << sm->GetIonisation()->GetFermiEnergy()/eV  << " " << sm->GetIonisation()->GetBirksConstant() << " " << sm->GetIonisation()->GetMeanEnergyPerIonPair()/eV << "\n";;
ofile.close(); 



} 








G4Box* SampleToBeSubtracted_solid
    = new G4Box("SampleToBeSubstracted", (World_Half_Thickness - 0.005*World_Half_Thickness)*cm, (World_Half_Thickness - 0.005*World_Half_Thickness)*cm, (NaCl_Source_thickness + NaCl_Source_thickness*0.01)*nm);  // half-thickness
G4Box* SampleToBeSubtractedFrom_solid
    = new G4Box("SampleToBeSubtractedFrom", (World_Half_Thickness - 0.01*World_Half_Thickness)*cm, (World_Half_Thickness - 0.01*World_Half_Thickness)*cm, (World_Half_Thickness - 0.01*World_Half_Thickness)*cm);  // half-thickness

G4VSolid* Sample_solid = new G4SubtractionSolid("Sample", SampleToBeSubtractedFrom_solid, SampleToBeSubtracted_solid, 0, G4ThreeVector(0.,0.,0.));

G4LogicalVolume* Sample_logical
= new G4LogicalVolume(
    Sample_solid,      // its solid
    Sample_Material,     // its material
    "Sample");         // its name                                 
G4VPhysicalVolume* Sample_physical
= new G4PVPlacement(
    0,                // no rotation
    G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm),  
    Sample_logical,           // its logical volume                         
    "Sample",       // its name
    world_logical,    // its mother  volume
    false,            // no boolean operation
    0,                // copy number
    fCheckOverlaps);  // checking overlaps 

 
  G4Colour colorRed(G4Colour::Red());
  G4Colour colorGreen(G4Colour::Green());
  G4Colour colorBlue(G4Colour::Blue());
  G4Colour colorYellow(G4Colour::Yellow());


  G4VisAttributes* detVisAtt= new G4VisAttributes(colorYellow);
  //detVisAtt->SetForceSolid(true);
  NaCl_Source_logical->SetVisAttributes(detVisAtt);


  G4VisAttributes* detVisAtt2= new G4VisAttributes(colorBlue);
  //detVisAtt2->SetForceSolid(true);
  Sample_logical->SetVisAttributes(detVisAtt2);





  return world_physical;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}


 void B4DetectorConstruction::SetParameters(G4int mn) 
    {
        material_num = mn;
    }

/*
/gps/particle proton
/gps/energy 100000.0 eV
/gps/pos/type Volume
/gps/pos/shape Cylinder
/gps/pos/radius 1.875 mm
/gps/pos/halfz 0.5 nm
/gps/ang/type iso
/run/beamOn 50

*/


