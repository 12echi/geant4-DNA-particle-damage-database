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
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include <cstdio>
#include <string>
#include <sstream>
#include <chrono>

#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4Run.hh"

//#include "Analysis.hh"
#include "ClusteringAlgo.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "SizeDistribution.hh"

#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "TrackerHit.hh" // 如果 TrackerHit 类定义在单独的头文件中，请包含该头文件
#include <vector> // 如果使用 vector 容器来存储 TrackerHit 对象，请包含该头文件
#include "TrackerSD.hh"


#include "G4ios.hh"


using namespace std;
extern G4double Energy;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {
  //default parameter values
 

  // Create clustering algorithm
  // These default values have been tuned for the Physics List G4EmDNAPhysics
  // to reproduce data published by:
  // Francis et al. 2011 Comput. Meth. Programs. Biomed. 2011 101(3)
  clustering_ = new ClusteringAlgo(3.2 * nanometer, 2, 1, 0.05, 5 * eV, 37.5 * eV, 0.65);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {
  delete clustering_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *) {
  fEdep = 0.;
  fTrackLength = 0.;
  clustering_->Purge();

  G4cout << "--------Mass: " << (G4LogicalVolumeStore::GetInstance()->
      GetVolume("Logical targer")->GetMass() / kg) << "--------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event) 
{

 //SD探测器获取信息000000000000000000000000000000000000000000000000

   G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
   auto nhit = hc->GetSize();
   
   G4double maxDistance = 0.0;
   G4ThreeVector maxDistancePoints[2];
   
   for (size_t i=0; i<nhit; ++i) {
        auto hit = static_cast<TrackerHit*>(hc->GetHit(i));
        G4ThreeVector pos1 = hit->GetPos();
        // 使用 pos 中的位置信息进行操作

        for (size_t j = i + 1; j < nhit; ++j) {
        auto hit2 = static_cast<TrackerHit*>(hc->GetHit(j));
        G4ThreeVector pos2 = hit2->GetPos();
        
        // 计算两点之间的直线距离
        double distance = (pos2 - pos1).mag();
        if (distance > maxDistance) {
            maxDistance = distance;
            maxDistancePoints[0] = pos1;
            maxDistancePoints[1] = pos2;
        }
        
        }
      
    }
G4cout << "Point 1: (" << maxDistancePoints[0].x()/um << ", " << maxDistancePoints[0].y()/um << ", " << maxDistancePoints[0].z()/um << ")" << G4endl;
G4cout << "Point 2: (" << maxDistancePoints[1].x()/um << ", " << maxDistancePoints[1].y()/um << ", " << maxDistancePoints[1].z()/um << ")" << G4endl;

G4double distance = (maxDistancePoints[1] - maxDistancePoints[0]).mag() / um;

G4cout << "Distance between Point 1 and Point 2: " << distance << " um" << G4endl;
//000000000000000000000000000000000步长计算轨迹长度

G4double totaLength = 0.;
for (size_t i = 0; i < nhit; ++i) {
    auto hit1 = static_cast<TrackerHit*>(hc->GetHit(i));
    G4double StepLength = hit1->GetStepLength(); // 获取当前击中的能量
    totaLength += StepLength; // 累加能量
}
G4cout << "Total Length: " << totaLength/um << G4endl;



//0000000000000000000000000000000000000000000000
G4double totalEdep = 0.; // 总沉积能量

for (size_t i = 0; i < nhit; ++i) {
    auto hit1 = static_cast<TrackerHit*>(hc->GetHit(i));
    G4double edep = hit1->GetEdep(); // 获取当前击中的能量
    totalEdep += edep; // 累加能量
}

// 输出总沉积能量
G4cout << "Total energy deposited: " << totalEdep/keV << G4endl;

G4double let = (totalEdep/keV)/(distance/um) ;



//000000000000000000000000000000000000000000000000000000000000000


  std::map<G4int, G4int> sizeDistribution = clustering_->RunClustering();


  G4int ssb_num = clustering_->GetSSB();
  G4int dsbn = clustering_->GetDSB();
  G4int dsbn_num = clustering_->GetDSBN();
  G4int dsbp_num = clustering_->GetDSBP();
  G4int dsbpp_num = clustering_->GetDSBPP();
  G4int CDSB = dsbp_num + dsbpp_num;

  RunAction *run_action =
      (RunAction *) G4RunManager::GetRunManager()->GetUserRunAction();
  run_action->AddEnergyDeposition(fEdep);
  run_action->AddSingleStrandBreak(ssb_num);
  run_action->AddComplexSingleStrandBreakNew(dsbn);
  run_action->AddDoubleStrandBreakNew(dsbn_num);
  run_action->AddDoubleStrandBreakPlus(dsbp_num);
  run_action->AddDoubleStrandBreakPlusPlus(dsbpp_num);
  run_action->AddCDSB(CDSB);

  for (auto it = sizeDistribution.begin(); it != sizeDistribution.end(); it++) {
    for (int i = 0; i < it->second; i++) {
      run_action->AddDistribution(it->first);
    }
  }
  //G4double let = (fEdep / keV )/(fTrackLength/ um);
 // run_action->Addlet(let);

  G4double total_energy_deposit = fEdep / keV;
  G4double absorbed_dose =
      (fEdep / joule) / (G4LogicalVolumeStore::GetInstance()->GetVolume("Logical targer")->GetMass() / kg);

  G4double ssb_yield = static_cast<double>(ssb_num)/ absorbed_dose / (6);
  G4double cssbn_yield = static_cast<double>(dsbn) / absorbed_dose / (6);
  G4double dsbn_yield = static_cast<double>(dsbn_num) / absorbed_dose / (6);
  G4double dsbp_yield = static_cast<double>(dsbp_num) / absorbed_dose / (6);
  G4double dsbpp_yield = static_cast<double>(dsbpp_num) / absorbed_dose / (6);

  //G4cout << "-----------------RESULT-----G4cout << "| DSB /SSB |" << dsbn_yield /ssb_yield << G4endl;------------" << G4endl;
  G4cout << "| ssb_yield |" << ssb_yield << " Gbp-1Gy-1" << "| ssb_num |" << ssb_num << G4endl;
  G4cout << "| cssbn_yield |" << cssbn_yield << " Gbp-1Gy-1" << "| cssbn_num |" << dsbn << G4endl;
  G4cout << "| dsbn_yield |" <<dsbn_yield << " Gbp-1Gy-1" << "| dsbn_num |" << dsbn_num << G4endl;
  G4cout << "| dsbp_yield |" << dsbp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbp_num << G4endl;
  G4cout << "| dsbpp_yield |" << dsbpp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbpp_num << G4endl;
  if (dsbn_yield != 0) G4cout << "| DSB /SSB |" << dsbn_yield /ssb_yield << G4endl;
  G4cout << "| energy_deposit |" << total_energy_deposit << " keV" << G4endl;
  G4cout << "| absorbed_dose | " << absorbed_dose << " J/Kg" << G4endl;
  G4cout << "----------------------------------------" << G4endl;


  ofstream DataOutOPut;
              DataOutOPut.open("./data/DataOutOPut.txt", ios::app);
              DataOutOPut << "Energy:\t" << Energy << "\t";
              //DataOutOPut << "inter:\t" << inter<< "\t";
              DataOutOPut << "dsb:\t" << dsbn << "\t";
              DataOutOPut << "cdsbn_num:\t" <<CDSB << "\t";
             // DataOutOPut << "energySD:\t" << totalEdep<< " keV"<< "\t";
              //DataOutOPut << "e_deposit:\t" << fEdep/keV << "keV" <<"\t";
             // DataOutOPut << "Length:\t" <<  fTrackLength / um << " um "<<"\t";
              //DataOutOPut << "Length_track:\t" <<  distance / um << " um "<<"\t";
              //DataOutOPut << "Total Length:\t" <<  totaLength/um << " um "<<"\t";
             // DataOutOPut << "let:\t" << let << "kev/um" << "\t";

              DataOutOPut << G4endl;
              DataOutOPut.close();





}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


