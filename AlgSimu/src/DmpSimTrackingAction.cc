/*
 *  $Id: DmpSimTrackingAction.cc, 2014-05-14 15:35:24 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 04/03/2014
*/

#include "G4TrackingManager.hh"
#include "G4Track.hh"

#include "DmpSimTrackingAction.h"
#include "DmpDataBuffer.h"
#include "DmpLog.h"

//-------------------------------------------------------------------
DmpSimTrackingAction::DmpSimTrackingAction():fTrackInfor(0){
  fTrackInfor = new DmpEvtMCTrack();
//  gDataBuffer->RegisterObject("Event/MCTruth/TrackVertex",fTrackInfor,"DmpEvtMCTrack");
}

//-------------------------------------------------------------------
void DmpSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack){
  if(aTrack->GetDefinition()->GetPDGCharge() == 0){
    fpTrackingManager->SetStoreTrajectory(false);
  }else{
    fpTrackingManager->SetStoreTrajectory(true);
  }

  
  int trackID = aTrack->GetTrackID();
  int nTrack = fTrackInfor->fTrackID.size();
  //std::cout << "Track ID size =" << nTrack << std::endl;
  for(size_t i=0;i<nTrack;++i){
    if(trackID == fTrackInfor->fTrackID[i]){
      return;
    }
  }

  
  fTrackInfor->fTrackID.push_back(trackID);
  fTrackInfor->fParentID.push_back(aTrack->GetParentID());
  fTrackInfor->fPDGCode.push_back(aTrack->GetDynamicParticle()->GetPDGcode());
  G4ThreeVector pos = aTrack->GetVertexPosition();
  fTrackInfor->fPosition.push_back(TVector3(pos.x()/mm,pos.y()/mm,pos.z()/mm));
  G4ThreeVector dir = aTrack->GetVertexMomentumDirection();
  fTrackInfor->fDirection.push_back(TVector3(dir.x(),dir.y(),dir.z()));
  fTrackInfor->fEnergy.push_back(aTrack->GetVertexKineticEnergy()/MeV);

  //debug
  /*
  if ((aTrack->GetDynamicParticle()->GetPDGcode()==22)&&(aTrack->GetParentID()==1)&&(aTrack->GetVolume()->GetName()=="PhotonGenerator")){
    std::cout << "Track ID:" << trackID << ",Paraent ID:" << aTrack->GetParentID() << ",volume:" << aTrack->GetVolume()->GetName() << ",PZ=" << aTrack->GetVertexPosition().z() << ",VKin=" << aTrack->GetVertexKineticEnergy()/MeV << "MeV" << std::endl;
  }
  */
        /*
  if(aTrack->GetParentID()==0) {        // Create trajectory only for primaries
    fpTrackingManager->SetStoreTrajectory(true);        // then will affect what??
  }else if(aTrack->GetParentID()==1 && aTrack->GetDefinition()->GetPDGCharge() != 0.) { // also e+ and e- from primary photon conversion
    fpTrackingManager->SetStoreTrajectory(true); 
  }else {
    //fpTrackingManager->SetStoreTrajectory(false);
  }
  */
}

void DmpSimTrackingAction::ResetTrackingData()const{
  fTrackInfor->Reset();
}

