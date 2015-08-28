/*
 *  $Id: DmpEvtMCTrack.h, 2014-10-21 20:33:48 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 17/10/2014
*/

#ifndef DmpEvtMCTrack_H
#define DmpEvtMCTrack_H

#include "TVector3.h"
#include <vector>

class DmpEvtMCTrack : public TObject{
/*
 *   This class used to record:
 *    1. vertex information of all secondray particle
 *          fPosition:  Where does the secondary particle generated?
 *          fDirection: the direction of the new secondary particle
 *          fEnergy:    kinetic energy of it
 *    2. information of a G4Step
 *          fPosition:  position of pre-step-point
 *          fDirection: momentum direction of pre-step-point
 *          fEnergy:    deposited energy of this step
 */
public:
  DmpEvtMCTrack();
  ~DmpEvtMCTrack();
  void Reset();
  DmpEvtMCTrack &operator=(const DmpEvtMCTrack &r);
  void LoadFrom(DmpEvtMCTrack *r);

public:
  void LoadTracks(DmpEvtMCTrack *left,double eLow)const;  // load tracks which fEnergy > eLow
  void LoadTracks(DmpEvtMCTrack *left,double zLow,double zHigh,double eLow=2.5)const;  // load tracks which are generated in Z range by E cut
  void Tracking(DmpEvtMCTrack *left,int parentID,bool onlyFirstElement = false)const;   // if parentID match, load it. NOTE:  the frist element is the vertex which TrackID == parentID
  void TrackingMaxE(DmpEvtMCTrack *left,int PDGCode=-12345)const;        // find the track which has the max energy
  void TrackingMaxPath(DmpEvtMCTrack *left,int PDGCode = -12345)const;    // find track has max path at xyz(0,1,2)

public:
  double TrackDistance(DmpEvtMCTrack *me)const;  // calculate absolute distance of this tracking

//public:
  double  GetEnergyOfMaxETrack(int PDGCode,double zCut=0)const;
  TVector3  GetP0OfMaxETrack(int PDGCode,double zCut=0)const;
//  double  GetEnergyOfMaxPathTrack()const;
//  TVector3  GetStartPointOfMaxPathTrack()const;
//  TVector3  GetDirectionOfMaxPathTrack()const;
//  TVector3  GetDirectionOfMaxETrack()const;

//private:
//   static double  _ZCutLow;     //! only study tracks which are created >Zlow < Zhigh
//   static double  _ZCutHigh;    //! unit: mm
//public:
//   static void SetZRange(double z_low,double z_high) {_ZCutLow = z_low; _ZCutHigh = z_high;}
//   static void GetZ0(){return _ZCutLow;};
//   static void GetZ1(){return _ZCutHigh;};

public:
  int GetSDParticleNumber(int pdgCode, double z0=0, double z1 = 600,int parentParticle=-1,double eLow=0)const;
  std::vector<int> GetParentPDGCode(int pdgCode, double z0=0, double z1 = 600,double eLow=0)const;

public:
  std::vector<int>          fTrackID;           // track ID
  std::vector<int>          fParentID;          // ID of Parent track
  std::vector<int>          fPDGCode;           // pdg code
  std::vector<TVector3>     fPosition;          // unit mm, position x,y,z
  std::vector<TVector3>     fDirection;         // momentum direction
  std::vector<double>       fEnergy;            // unit MeV

  ClassDef(DmpEvtMCTrack,1)
};

#endif


