/*
 *  $Id: DmpSimEventAction.cc, 2014-10-29 09:48:35 DAMPE/USTC $
 *  Author(s):
 *    Yifeng Wei (weiyf@mail.ustc.edu.cn) 29/10/2014
*/

#include "DmpSimEventAction.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "DmpCore.h"
#include "DmpSimTrackingAction.h"

DmpSimuEventAction::DmpSimuEventAction(DmpSimTrackingAction* fTracking)
{ 
  // Get pointer to the simulation algorithm
  //fDmpSimAlg = (DmpSimAlg*)gCore->AlgorithmManager()->Get("Sim/SimAlg");
  
  fTrackingAction = fTracking;

}

DmpSimuEventAction::~DmpSimuEventAction()
{
}

void DmpSimuEventAction::BeginOfEventAction(const G4Event* evt)
{
}

void DmpSimuEventAction::EndOfEventAction(const G4Event* evt)
{
}













