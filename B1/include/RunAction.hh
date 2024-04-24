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
/// \file RunAction.hh
/// \brief Definition of the B1::RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include <map>
#include <utility>
#include <string>

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
//#include "G4ParticleGun.hh"
#include "globals.hh"

#include "SizeDistribution.hh"

//#include "PrimaryGeneratorAction.hh"
class RunActionMessenger;

class RunAction : public G4UserRunAction {
 public:
  RunAction();
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run *);
  virtual void EndOfRunAction(const G4Run *);

  void SetHistoName(G4String &val) {
    fFileName = val;
  };

  inline void AddEnergyDeposition(G4double edep) {
    energy_deposition_ += edep;
  }

  inline void AddSingleStrandBreak(G4int num) {
    ssb_num_ += num;
  }


  inline void AddComplexSingleStrandBreakNew(G4int num) {
    dsb_num_ += num;
  }

  inline void AddDoubleStrandBreakNew(G4int num) {
    dsbn_num_ += num;
  }

  inline void AddDoubleStrandBreakPlus(G4int num) {
    dsbp_num_ += num;
  }

  inline void AddDoubleStrandBreakPlusPlus(G4int num) {
    dsbpp_num_ += num;
  }

  inline void AddCDSB(G4double  num) {
    CDSB+= num;
  }

  inline void AddDistribution(G4int num) {
    distribution_ += std::to_string(num) + " ";
  }

 private:
  std::map<int, int> String_Distribution(const std::string &str);

  G4String fFileName;
  RunActionMessenger *run_messenger_;

  G4Accumulable<G4double> energy_deposition_;
  G4Accumulable<G4int> ssb_num_;
  G4Accumulable<G4int> dsb_num_;
  G4Accumulable<G4int> dsbn_num_;
  G4Accumulable<G4int> dsbp_num_;
  G4Accumulable<G4int> dsbpp_num_;
  G4Accumulable<G4int> CDSB;
  G4Accumulable<G4double> let;
  G4Accumulable<G4String> distribution_;
};

#endif

