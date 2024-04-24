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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class


#include "SteppingAction.hh"
#include <G4LogicalVolumeStore.hh>
#include "ClusteringAlgo.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction(),RunInitObserver(),event_action_(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::Initialize()
{
  event_action_ = (EventAction*) G4EventManager::GetEventManager()->
      GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* pStep)
{  
    
if (!fScoringVolume) {
    const DetectorConstruction* detConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }



    //G4LogicalVolume *target_volume = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalOrb");
    
    
    G4LogicalVolume *the_volume = pStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    //G4double edepStep = pStep->GetTotalEnergyDeposit();

    //const G4StepPoint* preStepPoint = pStep->GetPreStepPoint();
   // const G4StepPoint* postStepPoint = pStep->GetPostStepPoint();

    // 获取前一步和后一步的位置
   // G4ThreeVector prePosition = preStepPoint->GetPosition();
    //G4ThreeVector postPosition = postStepPoint->GetPosition();
  
     // G4ThreeVector center(0, 0, 0); // 球体中心
      //G4double radius = 5 * um;      // 球体半径
    //G4double distance = (center - prePosition).mag();
    //G4double distance1 = (center - postPosition).mag();

    //G4double stepLength = pStep->GetStepLength();

    
    

   if (the_volume != fScoringVolume) return;
    
    
    G4double edepStep = pStep->GetTotalEnergyDeposit();

     event_action_->clustering_->RegisterPhysicalDamage(
        pStep->GetPreStepPoint()->GetPosition(),edepStep);


    event_action_->AddEdep(edepStep);
  // get the particle name and parent ID
    G4String particle_name = pStep->GetTrack()->GetDefinition()->GetParticleName();
     G4int parentID = pStep->GetTrack()->GetParentID();
  // collect energy deposited in this step for alpha particles with parentID 0
  if (particle_name != "alpha" )return;
  if( parentID != 0) return;
    
   G4double stepLength = pStep->GetStepLength();
    event_action_->AddStepLength(stepLength);
  
  

 /*if (particle_name == G4String("alpha") && parentID == 0)
 {
    
    if (the_volume != fScoringVolume) return;

    G4double stepLength = pStep->GetStepLength();

    fpEventAction->AddLength(stepLength);

    */
    
    /*if(distance>=radius && distance1<=radius )
    {
       fpEventAction->AddLength(stepLength);

    }

    if(distance<radius && distance1<radius )
    {

      fpEventAction->AddLength(stepLength);

   }
 
     if(distance<radius && distance1>=radius )
     {

     fpEventAction->AddLength(stepLength);
     }
  */
   
 
 /*if(target_volume == the_volume)
    {

      fpEventAction->AddEdep(edepStep);
    }
    
     G4String particle_name = pStep->GetTrack()->GetDefinition()->GetParticleName();
     G4int parentID = pStep->GetTrack()->GetParentID();
    
     G4double stepLength = pStep->GetStepLength();
  if (particle_name == G4String("alpha") && parentID == 0) 
    {
         if(target_volume == the_volume)
           {

            fpEventAction->AddLength(stepLength);
           }
    }

     */
     
      
    
  
}






















/*#include "SteppingAction.hh"
#include <G4LogicalVolumeStore.hh>
//#include "ClusteringAlgo.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
#include "ClusteringAlgo.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
    : G4UserSteppingAction(), RunInitObserver(), event_action_(0) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::Initialize() {
  event_action_ = (EventAction *) G4EventManager::GetEventManager()->
      GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* pStep)

{ 
   if (!fScoringVolume) {
    const DetectorConstruction* detConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }


    G4String particle_name = pStep->GetTrack()->GetDefinition()->GetParticleName();
    G4int parentID = pStep->GetTrack()->GetParentID();
  /*if (pStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() !=
      "Transportation") {*/
   //G4LogicalVolume *target_volume = G4LogicalVolumeStore::GetInstance()->GetVolume("targer");
    /*G4LogicalVolume *target_volume1 = G4LogicalVolumeStore::GetInstance()->GetVolume("logic sugar 4");
    G4LogicalVolume *target_volume3 = G4LogicalVolumeStore::GetInstance()->GetVolume("tin");
    G4LogicalVolume *target_volume4 = G4LogicalVolumeStore::GetInstance()->GetVolume("logic blue sphere");
    G4LogicalVolume *target_volume5 = G4LogicalVolumeStore::GetInstance()->GetVolume("logic pink sphere");
   G4LogicalVolume *the_volume = pStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
   //G4double edepStep = pStep->GetTotalEnergyDeposit();
  
    
  if(the_volume == fScoringVolume)
  {   
    G4double edepStep = pStep->GetTotalEnergyDeposit();
       event_action_->AddEdep(edepStep);

       event_action_->clustering_->RegisterPhysicalDamage(
        pStep->GetPreStepPoint()->GetPosition(),edepStep);
  }
 

    if (particle_name == G4String("alpha") && parentID == 0){


      if (the_volume != fScoringVolume) return;

        G4double stepLength = pStep->GetStepLength();
        event_action_->AddStepLength(stepLength);

    }

    
  
 
  
} */



