/*
 *  $Id: DmpEvtMCTrack.cc, 2014-10-20 22:23:41 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 17/10/2014
*/

#include "DmpEvtMCTrack.h"
#include <iostream>
#include "TMath.h"

ClassImp(DmpEvtMCTrack)


//double DmpEvtMCTrack::_ZCutLow = 0;
//double DmpEvtMCTrack::_ZCutHigh = 500;

//-------------------------------------------------------------------
DmpEvtMCTrack::DmpEvtMCTrack(){
}

//-------------------------------------------------------------------
DmpEvtMCTrack::~DmpEvtMCTrack(){
}

//-------------------------------------------------------------------
DmpEvtMCTrack& DmpEvtMCTrack::operator=(const DmpEvtMCTrack &r){
  Reset();
  fTrackID = r.fTrackID;
  fParentID = r.fParentID;
  fPDGCode = r.fPDGCode;
  fPosition = r.fPosition;
  fDirection = r.fDirection;
  fEnergy = r.fEnergy;
}

//-------------------------------------------------------------------
void DmpEvtMCTrack::LoadFrom(DmpEvtMCTrack *r){
  Reset();
  fTrackID = r->fTrackID;
  fParentID = r->fParentID;
  fPDGCode = r->fPDGCode;
  fPosition = r->fPosition;
  fDirection = r->fDirection;
  fEnergy = r->fEnergy;
}

//-------------------------------------------------------------------
void DmpEvtMCTrack::Reset(){
  //std::cout << "Reset MC Track event class~" << std::endl;
  fTrackID.clear();
  fParentID.clear();
  fPDGCode.clear();
  fPosition.clear();
  fDirection.clear();
  fEnergy.clear();
}

void DmpEvtMCTrack::LoadTracks(DmpEvtMCTrack *tmp,double eLow)const
{
  tmp->Reset();
  int n0 = fTrackID.size();
  for(int i=0;i<n0;++i){
    if(fEnergy.at(i) > eLow){
      tmp->fTrackID.push_back(fTrackID.at(i));
      tmp->fParentID.push_back(fParentID.at(i));
      tmp->fPDGCode.push_back(fPDGCode.at(i));
      tmp->fPosition.push_back(fPosition.at(i));
      tmp->fDirection.push_back(fDirection.at(i));
      tmp->fEnergy.push_back(fEnergy.at(i));
    }
  }
}

void DmpEvtMCTrack::LoadTracks(DmpEvtMCTrack *tmp,double z0,double z1,double eLow)const
{
  tmp->Reset();
  int n0 = fTrackID.size();
  for(int i=0;i<n0;++i){
    TVector3 po = fPosition.at(i);
    if(fEnergy.at(i) > eLow && po.z() > z0 && po.z()<z1){
      tmp->fTrackID.push_back(fTrackID.at(i));
      tmp->fParentID.push_back(fParentID.at(i));
      tmp->fPDGCode.push_back(fPDGCode.at(i));
      tmp->fPosition.push_back(fPosition.at(i));
      tmp->fDirection.push_back(fDirection.at(i));
      tmp->fEnergy.push_back(fEnergy.at(i));
    }
  }
}

void DmpEvtMCTrack::Tracking(DmpEvtMCTrack *l,int parentID,bool onlyFirstElement)const
{
  l->Reset();
  int nT = fTrackID.size();
  for(int i=0;i<nT;++i){
    if(fTrackID.at(i) == parentID){
      l->fTrackID.push_back(fTrackID.at(i));
      l->fParentID.push_back(fParentID.at(i));
      l->fPDGCode.push_back(fPDGCode.at(i));
      l->fPosition.push_back(fPosition.at(i));
      l->fDirection.push_back(fDirection.at(i));
      l->fEnergy.push_back(fEnergy.at(i));
    }
  }

  if(onlyFirstElement) return;

  for(int i=0;i<nT;++i){
    if(fParentID.at(i) == parentID){
      l->fTrackID.push_back(fTrackID.at(i));
      l->fParentID.push_back(fParentID.at(i));
      l->fPDGCode.push_back(fPDGCode.at(i));
      l->fPosition.push_back(fPosition.at(i));
      l->fDirection.push_back(fDirection.at(i));
      l->fEnergy.push_back(fEnergy.at(i));
    }
  }
}

void DmpEvtMCTrack::TrackingMaxE(DmpEvtMCTrack *l,int pdgcode)const
{
  l->Reset();
  int id = -1;
  double eM = -1;
  int nT = fEnergy.size();
  for(int i=0;i<nT;++i){
    if(fEnergy.at(i) > eM && (pdgcode == fPDGCode.at(i) || pdgcode == -12345)){
      id = fTrackID.at(i);
    }
  }
  this->Tracking(l,id);
}

void DmpEvtMCTrack::TrackingMaxPath(DmpEvtMCTrack *l,int pdgcode)const
{
  l->Reset();
  double max_dis = 0;
  int id = 0;
  int nT = fTrackID.size();
  DmpEvtMCTrack *tmp = new DmpEvtMCTrack();
  for(int i=0;i<nT;++i){
    if(fPDGCode.at(i) == pdgcode || pdgcode == -12345){
      this->Tracking(tmp,i);
      double dis = TrackDistance(tmp);
      if(max_dis < dis){
        max_dis = dis;
        id = i;
      }
    }
  }
  delete tmp;
  Tracking(l,id);
}

double DmpEvtMCTrack::TrackDistance(DmpEvtMCTrack *me)const
{
  TVector3 p0 = me->fPosition.at(0);
  TVector3 p1;
  double Max_dis=0;
  int np = me->fPosition.size();
  for(int i=0;i<np;++i){
    p1 = fPosition.at(i);
    double dis = TMath::Power((p1.x()-p0.x()),2) + TMath::Power(p1.y()-p0.y(),2) + TMath::Power(p1.z()-p0.z(),2);
    if(dis > Max_dis){Max_dis = dis;}
  }
  return TMath::Sqrt(Max_dis);
}

double DmpEvtMCTrack::GetEnergyOfMaxETrack(int pdgCode,double zCut)const
{
  int nT = fTrackID.size();
  double mE = 0;
  for(int i=0;i<nT;++i){
    if(pdgCode  == fPDGCode.at(i)){
      if(mE < fEnergy.at(i) && fPosition.at(i).z() > zCut){
        mE = fEnergy.at(i);
      }
    }
  }
  return mE;
}

TVector3 DmpEvtMCTrack::GetP0OfMaxETrack(int pdgCode,double zCut)const
{
  int nT = fTrackID.size();
  double mE = 0;
  TVector3 posi(0,0,0);
  for(int i=0;i<nT;++i){
    if(pdgCode  == fPDGCode.at(i)){
      if(mE < fEnergy.at(i) && fPosition.at(i).z()>zCut){
        mE = fEnergy.at(i);
        posi.SetX(fPosition.at(i).x());
        posi.SetY(fPosition.at(i).y());
        posi.SetZ(fPosition.at(i).z());
      }
    }
  }
  return posi;
}

int DmpEvtMCTrack::GetSDParticleNumber(int pdgCode,double z0, double z1,int parentCode,double eLow)const
{
  int count=0;
  int n0 = fTrackID.size();
  for(int i=0;i<n0;++i){
    if(fEnergy.at(i) > eLow && (fParentID.at(i) == parentCode || parentCode == -1)){
      if(fPosition.at(i).z()>z0 && fPosition.at(i).z()<z1){
        ++count;
      }
    }
  }
//std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<")"<<count<<"/"<<fTrackID.size()<<std::endl;
  return count;
}


std::vector<int> DmpEvtMCTrack::GetParentPDGCode(int pdgCode, double z0, double z1,double eLow)const
{
  std::vector<int> rme;
  int n0 = fTrackID.size();
  for(int i=0;i<n0;++i){
    if(fEnergy.at(i) > eLow){
      if(fPosition.at(i).z()>z0 && fPosition.at(i).z()<z1){
        rme.push_back(fParentID.at(i));
      }
    }
  }
  return rme;
}





