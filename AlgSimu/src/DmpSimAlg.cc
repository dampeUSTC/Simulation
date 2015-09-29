/*
 *  $Id: DmpSimAlg.cc, 2014-09-30 00:15:41 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 10/06/2014
*/

#include <boost/lexical_cast.hpp>

#include "DmpSimAlg.h"
#include "DmpSimRunManager.h"
#include "G4PhysListFactory.hh"
#include "DmpSimDetector.h"
#include "DmpSimMagneticField.h"
#include "DmpSimPrimaryGeneratorAction.h"
#include "DmpSimTrackingAction.h"
#include "DmpSimEventAction.hh"
#include "DmpCore.h"
#include "DmpDataBuffer.h"
#include "DmpMetadata.h"
#include "G4UImanager.hh"
#ifdef G4UI_USE_QT
#include "G4UIExecutive.hh"
#endif
#ifdef G4VIS_USE_OPENGLQT
#include "G4VisExecutive.hh"
#endif

//-------------------------------------------------------------------
DmpSimAlg::DmpSimAlg()
 :DmpVAlg("Sim/BootAlg"),
  fSimRunMgr(0),
  fPhyFactory(0),
  fSource(0),
  fDetector(0),
  fTracking(0),
  _SaveTrackVertex(false)
{
  fMetadata = new DmpMetadata();
  gDataBuffer->RegisterObject("Metadata/MCTruth/JobOpt",fMetadata,"DmpMetadata");
  fMetadata->SetOption("Mode","batch");
  fMetadata->SetOption("Physics","QGSP_BIC");
  fMetadata->SetOption("Gdml","FM");        // Fly Mode
  fMetadata->SetOption("Seed",boost::lexical_cast<std::string>(fMetadata->JobTime()));
  fMetadata->SetOption("Nud/DeltaTime","100");  // 100 ns
  gRootIOSvc->Set("Output/FileName","DmpSim_"+fMetadata->GetValue("Seed"));
  gRootIOSvc->Set("Output/Key","sim");
}

//-------------------------------------------------------------------
DmpSimAlg::~DmpSimAlg(){
  //delete fSimRunMgr;
}

//-------------------------------------------------------------------
void DmpSimAlg::Set(const std::string &type,const std::string &argv){
  if("gps/centre"==type || "gps/direction" == type){
    DmpLogWarning<<"Reseting "<<type<<": "<<fMetadata->GetValue(type)<<"\t new value = "<<argv<<DmpLogEndl;
  }
  if("Mode" == type && "batch" != argv){
    gRootIOSvc->Set("Output/FileName","DmpSimVis_"+fMetadata->GetValue("Seed"));
  }
  fMetadata->SetOption(type,argv);
}

//-------------------------------------------------------------------
#include <stdlib.h>     // getenv()
bool DmpSimAlg::Initialize(){
// set seed
  DmpLogCout<<"\tRandom seed: "<<fMetadata->GetValue("Seed")<<DmpLogEndl;      // keep this information in any case
  CLHEP::HepRandom::setTheSeed(boost::lexical_cast<long>(fMetadata->GetValue("Seed")));
// set G4 kernel
  fSimRunMgr = new DmpSimRunManager();
  fPhyFactory = new G4PhysListFactory();            fSimRunMgr->SetUserInitialization(fPhyFactory->GetReferencePhysList(fMetadata->GetValue("Physics")));
  fSource = new DmpSimPrimaryGeneratorAction();     fSimRunMgr->SetUserAction(fSource);      // only Primary Generator is mandatory
  fDetector = new DmpSimDetector();                 fSimRunMgr->SetUserInitialization(fDetector);
  fTracking = new DmpSimTrackingAction();
  if(_SaveTrackVertex) fSimRunMgr->SetUserAction(fTracking);
  fEvent    = new DmpSimuEventAction(fTracking);    fSimRunMgr->SetUserAction(fEvent);
  fSimRunMgr->Initialize();
  fSource->ApplyGPSCommand(); // must after fSimRunMgr->Initialize()
// boot simulation
  if(fMetadata->GetValue("Mode") == "batch"){    // batch mode
    if(fSimRunMgr->ConfirmBeamOnCondition()){   // if not vis mode, do some prepare for this run. refer to G4RunManagr::BeamOn()
      fSimRunMgr->SetNumberOfEventsToBeProcessed(gCore->GetMaxEventNumber());
      fSimRunMgr->ConstructScoringWorlds();
      fSimRunMgr->RunInitialization();
      fSimRunMgr->InitializeEventLoop(gCore->GetMaxEventNumber());
    }else{
      DmpLogError<<"G4RunManager::Initialize() failed"<<DmpLogEndl;
      return false;
    }
  }else{    // visualization mode
    G4UImanager *uiMgr = G4UImanager::GetUIpointer();
#ifdef G4UI_USE_QT
    char *dummyargv[20]={"visual"};
    G4UIExecutive *ui = new G4UIExecutive(1,dummyargv);
#ifdef G4VIS_USE_OPENGLQT
    G4VisExecutive *vis = new G4VisExecutive();
    vis->Initialize();
// *
// *  TODO: publish... check prefix
// *
    G4String prefix = (G4String)getenv("DMPSWWORK")+"/share/Simulation/";
    uiMgr->ApplyCommand("/control/execute "+prefix+"DmpSimVis.mac");
#endif
    if (ui->IsGUI()){
      uiMgr->ApplyCommand("/control/execute "+prefix+"DmpSimGUI.mac");
    }
    ui->SessionStart();
    delete ui;
#ifdef G4VIS_USE_OPENGLQT
    delete vis;
#endif
    gRootIOSvc->FillData("Event");
    gCore->TerminateRun();  // just for check GDML, or run one event while debuging
#endif
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpSimAlg::ProcessThisEvent(){
  fTracking->ResetTrackingData();
  if(fSimRunMgr->SimOneEvent(gCore->GetCurrentEventID())){
    return true;
  }
  return false;
}

//-------------------------------------------------------------------
bool DmpSimAlg::Finalize(){
  if(fMetadata->GetValue("Mode") == "batch"){
    fSimRunMgr->TerminateEventLoop();
    fSimRunMgr->RunTermination();
  }
  if(fTracking){
    delete fTracking;
  }
  if(fDetector){
    delete fDetector;
  }
  if(fSource){
    delete fSource;
  }
  if(fPhyFactory){
    delete fPhyFactory;
  }
  return true;
}

void DmpSimAlg::ResetGPS(){
    fMetadata->ClearOptions("gps/");
}

void DmpSimAlg::ActiveGPS0()
{   // point
  fMetadata->SetOption("gps/centre","0 0 -50 cm");
  fMetadata->SetOption("gps/pos/type","Plane");
  fMetadata->SetOption("gps/pos/shape","Circle");
  fMetadata->SetOption("gps/pos/radius","5 mm");
  fMetadata->SetOption("gps/direction","0 0 1");
  fMetadata->SetOption("gps/particle","mu-");
  fMetadata->SetOption("gps/ene/type","Mono");
  fMetadata->SetOption("gps/ene/mono","10 GeV");
}

void DmpSimAlg::ActiveGPS1()    // beam
{
  fMetadata->SetOption("gps/centre","0 0 -2700 cm");//1700 is not enough for SPS
  fMetadata->SetOption("gps/pos/type","Beam");
  fMetadata->SetOption("gps/pos/shape","Circle");
  fMetadata->SetOption("gps/pos/radius","1 mm");
  fMetadata->SetOption("gps/pos/sigma_r","2 mm");
  fMetadata->SetOption("gps/direction","0 0 1");
  fMetadata->SetOption("gps/particle","mu-");
  fMetadata->SetOption("gps/ene/type","Gauss");
  fMetadata->SetOption("gps/ene/mono","100 GeV");
  fMetadata->SetOption("gps/ene/sigma","10 MeV");
}

void DmpSimAlg::ActiveGPS2()    // plane
{
  fMetadata->SetOption("gps/centre","0 0 -50 cm");
  fMetadata->SetOption("gps/pos/type","Plane");
  fMetadata->SetOption("gps/pos/shape","Circle");
  fMetadata->SetOption("gps/pos/radius","50 cm");
  fMetadata->SetOption("gps/direction","0 0 1");
  fMetadata->SetOption("gps/particle","mu-");
  fMetadata->SetOption("gps/ene/type","Mono");
  fMetadata->SetOption("gps/ene/mono","10 GeV");
}

void DmpSimAlg::ActiveGPS3()    // sphere
{
  fMetadata->SetOption("gps/pos/shape","Sphere");
  fMetadata->SetOption("gps/pos/type","Surface");
  fMetadata->SetOption("gps/pos/centre","0 0 0 cm");
  fMetadata->SetOption("gps/pos/radius","1 m");
  fMetadata->SetOption("gps/ang/type","cos");
  fMetadata->SetOption("gps/particle","mu-");
  fMetadata->SetOption("gps/ene/type","Mono");
  fMetadata->SetOption("gps/ene/mono","10 GeV");
}

void DmpSimAlg::SetGPSType(int i){
  if(i == 0){
    this->ActiveGPS0();
  }else if(i == 1){
    this->ActiveGPS1();
  }else if(i == 2){
    this->ActiveGPS2();
  }else if(i == 3){
    this->ActiveGPS3();
  }
}


