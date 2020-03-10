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
// $Id: B4RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class 

#include "B4RunAction.hh"


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B4PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);

    std::fill_n(implantation_x_histogram, 4000, 0);
    std::fill_n(implantation_y_histogram, 4000, 0);
    std::fill_n(implantation_z_histogram, 4000, 0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* run)
{ 


    std::fill_n(implantation_x_histogram, 4000, 0);
    std::fill_n(implantation_y_histogram, 4000, 0);
    std::fill_n(implantation_z_histogram, 4000, 0);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* run)
{ 

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

    // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B4PrimaryGeneratorAction* generatorAction
   = static_cast<const B4PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }


  


if (implantation_z_histogram[0] > 0) {

  int runNum = run->GetRunID();
  std::ostringstream commandOS;
  commandOS << "Annihilation_Coords_1000000_1mm_diameter_10uCi_22Na_Hist_PEN_4gccElements.txt";
  std::ofstream ofile;
  ofile.open (G4String(commandOS.str()), std::ios::out | std::ios::app);     // ascii file   

  for (size_t i = 0; i < 4000; i++) {
  ofile << (runNum+84) << " " << implantation_x_histogram[i] << " "  << implantation_y_histogram[i] << " " << implantation_z_histogram[i] << "\n";
   }

  ofile.close(); 


std::ostringstream commandOS2;
commandOS2 << "MaterialsRun_1000000_1mm_diameter_10uCi_22Na_Hist_PEN_4gccElements.txt";
std::ofstream ofile2;
ofile2.open (G4String(commandOS2.str()), std::ios::out | std::ios::app);     // ascii file   
ofile2 << (run->GetRunID()+84) << "\n";
ofile2.close(); 


}








 // G4cout << "Num fissioned ones: " << numFissionNeutrons << G4endl;
 // G4cout << "Num lost ones: " << numLostNeutrons << G4endl;
        
  // Print
  //  
  // if (IsMaster()) {
  //   G4cout
  //    << G4endl
  //    << "--------------------End of Global Run-----------------------";
  // }
  // else {
  //   G4cout
  //    << G4endl
  //    << "--------------------End of Local Run------------------------";
  // }
  
  // G4cout
  //    << G4endl
  //    << " The run consists of " << nofEvents << " "<< runCondition
  //    << G4endl
  //    << "------------------------------------------------------------"
  //    << G4endl
  //    << G4endl;



    std::fill_n(implantation_x_histogram, 4000, 0);
    std::fill_n(implantation_y_histogram, 4000, 0);
    std::fill_n(implantation_z_histogram, 4000, 0);



}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
