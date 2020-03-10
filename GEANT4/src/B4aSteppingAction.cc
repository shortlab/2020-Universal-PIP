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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class 

#include "B4aSteppingAction.hh"
#include "B4RunAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "B4RunAction.hh"
#include "G4Decay.hh"
//#include "G4RadioActiveDecay.hh"


#include "G4ElectronIonPair.hh"

#include "G4Step.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>



using namespace std;

G4Decay* theDecayProcess = new G4Decay();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction, B4RunAction* runAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    fRunAction(runAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* theStep)
{




  //get event #
  G4int eID = 0;
  G4int runNum = 0;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
  if(evt) eID = evt->GetEventID();
  if(run) runNum = run->GetRunID();
  G4Track *theTrack = theStep->GetTrack(); 


// Collect energy and track length step by step




  G4StepPoint* preStepPoint = theStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = theStep->GetPostStepPoint();




if(preStepPoint->GetProcessDefinedStep() != 0) { 
if (postStepPoint->GetProcessDefinedStep()->GetProcessName() == "annihil") {
        
          /*
          std::ostringstream commandOS;
          commandOS << "Annihilation_Coords_1000000_1mm_diameter_10uCi_22Na.txt";
          std::ofstream ofile;
          ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file    
          ofile << runNum << " " << (int)(postStepPoint->GetPosition().getX()/um) << " " << (int)(postStepPoint->GetPosition().getY()/um) << " " << (int)(postStepPoint->GetPosition().getZ()/um) <<  "\n";
          ofile.close(); */


  int xposition_um = abs((int)(postStepPoint->GetPosition().getX()/um));
  int yposition_um = abs((int)(postStepPoint->GetPosition().getY()/um));
  int zposition_um = abs((int)(postStepPoint->GetPosition().getZ()/um));


  fRunAction->implantation_x_histogram[min(4000,xposition_um)] += 1;
  fRunAction->implantation_y_histogram[min(4000,yposition_um)] += 1;
  fRunAction->implantation_z_histogram[min(4000,zposition_um)] += 1;

  theTrack->SetTrackStatus(fKillTrackAndSecondaries);


}} 






  // trackID == 1 means that you are a trimary particle (i.e. you are a beta that I am shooting in)
/*
  if((preStepPoint->GetProcessDefinedStep() == 0) && (theTrack->GetTrackID()  == 1)) {

  // look at all the stuff as soon as you make it

  if(preStepPoint->GetProcessDefinedStep() == 0) {

          std::ostringstream commandOS;
          commandOS << "22Na_Beta_InitialEnergy_1000000.txt";
          std::ofstream ofile;
          ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file    
          ofile << theTrack->GetParticleDefinition()->GetParticleName() << " " << preStepPoint->GetKineticEnergy()/MeV << "\n";
          ofile.close(); 

}}
*/



 




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......