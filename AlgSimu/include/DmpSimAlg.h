/*
 *  $Id: DmpSimAlg.h, 2014-10-03 22:19:33 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 10/06/2014
*/

#ifndef DmpSimAlg_H
#define DmpSimAlg_H

#include "DmpVAlg.h"

class DmpMetadata;
class DmpSimRunManager;
class DmpSimPrimaryGeneratorAction;
class DmpSimDetector;
class DmpSimTrackingAction;
class DmpSimuEventAction;
class G4PhysListFactory;

class DmpSimAlg : public DmpVAlg{
/*
 *  DmpSimAlg
 *
 *
 */
public:
  DmpSimAlg();
  ~DmpSimAlg();

  void Set(const std::string &option,const std::string &value);
  /*
   * Options:
   *    |--Mode
   *    |       value:
   *    |           "batch"     // default
   *    |           "Vis" or "vis" or "VIS" or "visualization" or "Visualization"
   *    +--Physics
   *    |       value:
   *    |           "QGSP_BIC"  // default
   *    |           "QGSP_BERT"
   *    |           "QGSP_BERT_HP"
   *    |           many values...
   *    +--Gdml
   *    |       value:
   *    |          sub-directory       // default: "FM" (Will read $DMPSWSYS/share/Geometry/FM/DAMPE.gdml)
   *    |          "Path/Filename.gdml"
   *    +--Seed             // default: current time
   *    +--Nud
   *    |  +--DeltaTime     // default: "100"
   *    +--BT
   *    |  +--Magnetic      // default: "0 0 0"
   *    |  +--DAMPE
   *    |  |  +--Translation    // deault: "0 0 0"
   *    |  |  +--Rotation       // deault: "0 0 0"
   */
  void SaveTrackVertex() {_SaveTrackVertex = true;}
  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();
  DmpMetadata* GetMetadata(){return fMetadata;}

public:
  void ResetGPS();
  void SetGPSType(int i);
  void ActiveGPS0();
  void ActiveGPS1();    // spherical surface for Acc. simu.
  void ActiveGPS2();    // beam test

private:
  DmpMetadata           *fMetadata;
  DmpSimRunManager      *fSimRunMgr;    // run manager of simulation
  G4PhysListFactory     *fPhyFactory;
  DmpSimPrimaryGeneratorAction  *fSource;   // particle source
  DmpSimDetector        *fDetector;     // DAMPE detector
  DmpSimTrackingAction  *fTracking;     // tracking action
  DmpSimuEventAction    *fEvent;        // event action

private:
  bool  _SaveTrackVertex;

};

#endif

