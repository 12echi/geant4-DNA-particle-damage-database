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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// The original DBSCAN clustering algorithm is reproduced by Henri Payno
// and Yann Perrot at Blaise Pascal University, which has been released as
// a Geant4-DNA user example, then published in a file sub-folder named
// “extended/medical/dna/clustering”.
//
// Jing Tang,Qinfeng Xiao et al. improve the DBSCAN algorithm by utilizing
// the KD-Tree to find neighbors of each site to calculate clustered DNA damage
// This work is published:
//
// $Id$
//
/// \file ClusteringAlgo.cc
/// \brief Implementation of the ClustreringAlgo class

#include "ClusteringAlgo.hh"
#include "ClusteringAlgoMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//#include "flann/flann.hpp"

#include <map>
#include <chrono>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgo::ClusteringAlgo(G4double pEps,
                               G4int pMinPts,
                               G4double pSPointsProb,
                               G4double pChemicalPointsProb,
                               G4double pEMinDamage,
                               G4double pEMaxDamage,
                               G4double ind_dam_prob)
    : energy_(pEps),
      min_pts_(pMinPts),
      damage_prob_(pSPointsProb),
      chemical_damage_prob_(pChemicalPointsProb),
      min_damage_energy_(pEMinDamage),
      max_damage_energy_(pEMaxDamage),
      ind_dam_prob_(ind_dam_prob),
      runing_time_(0.0) {
  next_point_id_ = 0;
  messenger_ = new ClusteringAlgoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgo::~ClusteringAlgo() {
  delete messenger_;
  Purge();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random sampling in space
G4bool ClusteringAlgo::IsInSensitiveArea() {
  return damage_prob_ > G4UniformRand();
}

G4bool ClusteringAlgo::IsInChemicalSensitiveArea() {
  return chemical_damage_prob_ >= G4UniformRand();
}

G4bool ClusteringAlgo::IsIndDamProb() {
  return ind_dam_prob_ > G4UniformRand();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random sampling in energy
G4bool ClusteringAlgo::IsEdepSufficient(G4double pEdep) {
  if (pEdep < min_damage_energy_) {
    return false;
  } else if (pEdep > max_damage_energy_) {
    return true;
  } else {
    G4double proba = (pEdep / eV - min_damage_energy_ / eV) /
        (max_damage_energy_ / eV - min_damage_energy_ / eV);
    return (proba > G4UniformRand());
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add an event interaction to the unregistered damage in
// the physical process if good conditions (pos and energy) are met
//

void ClusteringAlgo::RegisterPhysicalDamage(G4ThreeVector pPos, G4double pEdep) {
  if (IsEdepSufficient(pEdep)) {
    if (IsInSensitiveArea()) {
      set_of_points_.push_back(new SBPoint(next_point_id_++, pPos, pEdep));
    }
  }
}
// Add an event interaction to the unregistered damage in
// the chemical process if ·OH radical hit DNA
G4bool ClusteringAlgo::RegisterChemicalDamage(G4ThreeVector pPos) {
  if (IsInChemicalSensitiveArea()) {
    if (IsIndDamProb()) {
      set_of_points_.push_back(new SBPoint(next_point_id_++, pPos, max_damage_energy_));
      return true;
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<G4int, G4int> ClusteringAlgo::RunClustering() {
  

  std::vector<SBPoint*>::iterator itVisitorPt, itObservedPt;
  for(itVisitorPt = set_of_points_.begin();
      itVisitorPt != set_of_points_.end();
      ++itVisitorPt  )
  {
    itObservedPt = itVisitorPt;
    itObservedPt ++;
    while(itObservedPt != set_of_points_.end() )
    {
      // if at least one of the two points has not a cluster
      if(!((*itObservedPt)->HasCluster() && (*itVisitorPt)->HasCluster()))
      {
        if(AreOnTheSameCluster( (*itObservedPt)->GetPosition(),
                                (*itVisitorPt)->GetPosition(),energy_))
        {
          // if none has a cluster. Create a new one
          if(!(*itObservedPt)->HasCluster() && !(*itVisitorPt)->HasCluster())
          {
            // create the new cluster
            set<SBPoint*> clusterPoints;
            clusterPoints.insert((*itObservedPt));
            clusterPoints.insert((*itVisitorPt));
            ClusterSBPoints* lCluster = new ClusterSBPoints(clusterPoints);
            assert(lCluster);
            clusters_.push_back(lCluster);
            assert(lCluster);
            // inform SB point that they are part of a cluster now
            assert(lCluster);
            (*itObservedPt)->SetCluster(lCluster);
            assert(lCluster);
            (*itVisitorPt)->SetCluster(lCluster);
          }else
          {
            // add the point to the existing cluster
            if((*itObservedPt)->HasCluster())
            {
              (*itObservedPt)->GetCluster()->AddSBPoint((*itVisitorPt));
              (*itVisitorPt)->SetCluster((*itObservedPt)->GetCluster());
            }

            if((*itVisitorPt)->HasCluster())
            {
              (*itVisitorPt)->GetCluster()->AddSBPoint((*itObservedPt));
              (*itObservedPt)->SetCluster((*itVisitorPt)->GetCluster());
            }
          }
        }
      }
      ++itObservedPt;
    }
  }

  // associate isolated points and merge clusters
  IncludeUnassociatedPoints();
  MergeClusters();

  // return cluster size distribution
  return GetClusterSizeDistribution();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Try to merge cluster between them, based on the distance between barycenters
void ClusteringAlgo::MergeClusters() {
  std::vector<ClusterSBPoints *>::iterator itCluster1, itCluster2;
  for (itCluster1 = clusters_.begin();
       itCluster1 != clusters_.end();
       ++itCluster1) {
    G4ThreeVector baryCenterClust1 = (*itCluster1)->GetBarycenter();
    itCluster2 = itCluster1;
    itCluster2++;
    while (itCluster2 != clusters_.end()) {
      G4ThreeVector baryCenterClust2 = (*itCluster2)->GetBarycenter();
      // if we can merge both cluster
      if (AreOnTheSameCluster(baryCenterClust1, baryCenterClust2, energy_)) {
        (*itCluster1)->MergeWith(*itCluster2);
        delete *itCluster2;
        clusters_.erase(itCluster2);
        return MergeClusters();
      } else {
        itCluster2++;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgo::IncludeUnassociatedPoints() {
  std::vector<SBPoint *>::iterator itVisitorPt;
  int nbPtSansCluster = 0;
  // Associate all point not in a cluster if possible ( to the first found cluster)
  for (itVisitorPt = set_of_points_.begin();
       itVisitorPt != set_of_points_.end();
       ++itVisitorPt) {
    if (!(*itVisitorPt)->HasCluster()) {
      nbPtSansCluster++;
      FindCluster(*itVisitorPt);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgo::FindCluster(SBPoint *pPt) {
  assert(!pPt->HasCluster());
  std::vector<ClusterSBPoints *>::iterator itCluster;
  for (itCluster = clusters_.begin();
       itCluster != clusters_.end();
       ++itCluster) {
    if ((*itCluster)->HasInBarycenter(pPt, energy_)) {
      (*itCluster)->AddSBPoint(pPt);
      return true;
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgo::AreOnTheSameCluster(G4ThreeVector pPt1,
                                         G4ThreeVector pPt2, G4double pMinDist) {
  G4double x1 = pPt1.x() / nm;
  G4double y1 = pPt1.y() / nm;
  G4double z1 = pPt1.z() / nm;

  G4double x2 = pPt2.x() / nm;
  G4double y2 = pPt2.y() / nm;
  G4double z2 = pPt2.z() / nm;

  // If the two points are closed enough
  return ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) <=
      (pMinDist / nm * pMinDist / nm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Determines the types and numbers of clustered damage including SSB, SSB+, DSB, DSB+,
// and DSB++, as denoted by Nikjoo et al. 1997
G4int ClusteringAlgo::GetSSB() const {
  G4int nbSSB = 0;
  std::vector<SBPoint *>::const_iterator itSDSPt;
  for (itSDSPt = set_of_points_.begin();
       itSDSPt != set_of_points_.end();
       ++itSDSPt) {
    if (!(*itSDSPt)->HasCluster()) {
      nbSSB++;
    }
  }
  return nbSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*G4int ClusteringAlgo::GetComplexSSB() const {
  G4int nbSSB = 0;
  std::vector<ClusterSBPoints *>::const_iterator itCluster;
   for (itCluster = clusters_.begin();
      itCluster != clusters_.end();
      ++itCluster) {
    if ((*itCluster)->IsSSB()) {
      nbSSB++;
   }
 }
  return nbSSB;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetDSB() const {
  G4int nbDSB = 0;
  std::vector<ClusterSBPoints *>::const_iterator itCluster;
  for (itCluster = clusters_.begin();
      itCluster != clusters_.end();
      ++itCluster) {
   if ((*itCluster)->IsDSB()) {
     nbDSB++;
   }
 }
  return nbDSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ClusteringAlgo::IsCSSBN(std::vector<ClusterSBPoints *>::const_iterator it) const {
  std::set<SBPoint *> cluster_points = (*it)->GetRegistredSBPoints();
  G4int strand1_damage_num = 0;
  G4int strand2_damage_num = 0;
  for (auto point = cluster_points.begin(); point != cluster_points.end(); point++) {
    if ((*point)->GetTouchedStrand() == 0) {
      strand1_damage_num++;
    } else {
      strand2_damage_num++;
    }
  }
  if (std::min(strand1_damage_num, strand2_damage_num) == 0 && std::max(strand1_damage_num, strand2_damage_num) >= 1) {
    return true;
  } else {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ClusteringAlgo::IsDSBN(std::vector<ClusterSBPoints *>::const_iterator it) const {
  std::set<SBPoint *> cluster_points = (*it)->GetRegistredSBPoints();
  G4int strand1_damage_num = 0;
  G4int strand2_damage_num = 0;
  for (auto point = cluster_points.begin(); point != cluster_points.end(); point++) {
    if ((*point)->GetTouchedStrand() == 0) {
      strand1_damage_num++;
    } else {
      strand2_damage_num++;
    }
  }
  if (strand1_damage_num == 1 && strand2_damage_num == 1) {
    return true;
  } else {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ClusteringAlgo::IsDSBP(std::vector<ClusterSBPoints *>::const_iterator it) const {
  std::set<SBPoint *> cluster_points = (*it)->GetRegistredSBPoints();
  G4int strand1_damage_num = 0;
  G4int strand2_damage_num = 0;
  for (auto point = cluster_points.begin(); point != cluster_points.end(); point++) {
    if ((*point)->GetTouchedStrand() == 0) {
      strand1_damage_num++;
    } else {
      strand2_damage_num++;
    }
  }
  if (std::min(strand1_damage_num, strand2_damage_num) == 1 && std::max(strand1_damage_num, strand2_damage_num) >= 2) {
    return true;
  } else {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ClusteringAlgo::IsDSBPP(std::vector<ClusterSBPoints *>::const_iterator it) const {
  std::set<SBPoint *> cluster_points = (*it)->GetRegistredSBPoints();
  G4int strand1_damage_num = 0;
  G4int strand2_damage_num = 0;
  for (auto point = cluster_points.begin(); point != cluster_points.end(); point++) {
    if ((*point)->GetTouchedStrand() == 0) {
      strand1_damage_num++;
    } else {
      strand2_damage_num++;
    }
  }
  if (std::min(strand1_damage_num, strand2_damage_num) >= 2) {
    return true;
  } else {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetCSSBN() const {
  G4int nbCSSBN = 0;
  std::vector<ClusterSBPoints *>::const_iterator it;
  for (it = clusters_.begin(); it != clusters_.end(); it++) {
    if (IsCSSBN(it)) {
      nbCSSBN++;
    }
  }
  return nbCSSBN;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetDSBN() const {
  G4int nbDSBN = 0;
  std::vector<ClusterSBPoints *>::const_iterator it;
  for (it = clusters_.begin(); it != clusters_.end(); it++) {
    if (IsDSBN(it)) {
      nbDSBN++;
    }
  }
  return nbDSBN;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetDSBP() const {
  G4int nbDSBP = 0;
  std::vector<ClusterSBPoints *>::const_iterator it;
  for (it = clusters_.begin(); it != clusters_.end(); it++) {
    if (IsDSBP(it)) {
      nbDSBP++;
    }
  }
  return nbDSBP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetDSBPP() const {
  G4int nbDSBPP = 0;
  std::vector<ClusterSBPoints *>::const_iterator it;
  for (it = clusters_.begin(); it != clusters_.end(); it++) {
    if (IsDSBPP(it)) {
      nbDSBPP++;
    }
  }
  return nbDSBPP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<G4int, G4int> ClusteringAlgo::GetClusterSizeDistribution() {
  std::map<G4int, G4int> sizeDistribution;
  sizeDistribution[1] = GetSSB();
  std::vector<ClusterSBPoints *>::const_iterator itCluster;
  for (itCluster = clusters_.begin();
       itCluster != clusters_.end();
       itCluster++) {
    sizeDistribution[(*itCluster)->GetSize()]++;
  }
  return sizeDistribution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgo::Purge() {
  next_point_id_ = 0;
  std::vector<ClusterSBPoints *>::iterator itCluster;
  for (itCluster = clusters_.begin();
       itCluster != clusters_.end();
       ++itCluster) {
    delete *itCluster;
    *itCluster = NULL;
  }
  clusters_.clear();
  std::vector<SBPoint *>::iterator itPt;
  for (itPt = set_of_points_.begin();
       itPt != set_of_points_.end();
       ++itPt) {
    delete *itPt;
    *itPt = NULL;
  }
  set_of_points_.clear();
}

