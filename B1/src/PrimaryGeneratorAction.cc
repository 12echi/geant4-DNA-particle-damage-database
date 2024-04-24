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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include <cmath>
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), particle_gun_(0) {
  G4int n_particle = 1;
  particle_gun_ = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleDefinition *particle
      = G4ParticleTable::GetParticleTable()->FindParticle("alpha");   //("proton");
   particle_gun_ ->SetParticleDefinition(particle);
   particle_gun_ ->SetParticleEnergy(0.5*MeV);
    particle_gun_->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete particle_gun_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
 
 
  //G4double  theta = G4UniformRand() * M_PI;
  //G4double  phi = G4UniformRand() * M_PI * 2.0;

  //G4double momentum_x = sin(theta) * cos(phi);
  //G4double momentum_y = sin(theta) * sin(phi);
  //G4double momentum_z = cos(theta);

  /*G4double pos_x = 5*um * cos(phi) ; 
  G4double pos_y = 5*um * sin(phi) ;
  G4double pos_z = 0.1 * nm;*/
  
  G4double pos_x = 0 ; 
  G4double pos_y = 0 ;
  G4double pos_z = 0.501 * um;


   //G4ThreeVector momentumDirection = -G4ThreeVector(0,0,1);
   particle_gun_->SetParticlePosition(G4ThreeVector(pos_x, pos_y, pos_z));
   //particle_gun_->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
   particle_gun_->GeneratePrimaryVertex(anEvent);
}





