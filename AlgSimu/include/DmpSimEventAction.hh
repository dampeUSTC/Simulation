/*
 *  $Id: DmpSimEventAction.hh, 2014-10-29 09:48:49 DAMPE/USTC $
 *  Author(s):
 *    Yifeng Wei (weiyf@mail.ustc.edu.cn) 29/10/2014
*/

#ifndef DmpSimEventAction_H
#define DmpSimEventAction_H 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DmpSimAlg.h"

class DmpSimTrackingAction;
class DmpSimuEventAction : public G4UserEventAction
{
public:

  DmpSimuEventAction(DmpSimTrackingAction*);
  virtual ~DmpSimuEventAction();
  
public:
  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  
private:

  DmpSimAlg *fDmpSimAlg;
  DmpSimTrackingAction *fTrackingAction;
};

#endif

    





