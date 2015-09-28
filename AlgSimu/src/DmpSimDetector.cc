/*
 *  $Id: DmpSimDetector.cc, 2014-09-30 14:07:44 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 26/02/2014
 *    Yifeng Wei (weiyf@mail.ustc.edu.cn) 28/09/2014
*/

#include <stdlib.h>     // getenv()

//#include "G4RotationMatrix.hh"
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "DmpLog.h"
#include "DmpDataBuffer.h"
#include "DmpMetadata.h"

#include "DmpSimDetector.h"
//#include "DmpSimPsdSD.h"
//#include "DmpSimStkSD.h"
#include "DmpSimBgoSD.h"
//#include "DmpSimNudSD.h"

//-------------------------------------------------------------------
DmpSimDetector::DmpSimDetector()
 :fParser(0),
  fPhyVolume(0),
//  fPsdSD(0),
//  fStkSD(0),
  fBgoSD(0)
//  fNudSD(0)
{
  fParser = new G4GDMLParser();
//  fPsdSD = new DmpSimPsdSD();
//  fStkSD = new DmpSimStkSD();
  fBgoSD = new DmpSimBgoSD();
//  fNudSD = new DmpSimNudSD();
  fMetadata = dynamic_cast<DmpMetadata*>(gDataBuffer->ReadObject("Metadata/MCTruth/JobOpt"));
}

//-------------------------------------------------------------------
DmpSimDetector::~DmpSimDetector(){
  delete fParser;
//  delete fPsdSD;
//  delete fStkSD;
  delete fBgoSD;
//  delete fNudSD;
}

//-------------------------------------------------------------------
#include <boost/filesystem.hpp>     // path
G4VPhysicalVolume* DmpSimDetector::Construct(){
  boost::filesystem::path   gdmlFile=fMetadata->GetValue("Gdml");
  if(gdmlFile.extension().string() != ".gdml"){    // argv = sub-directory
    gdmlFile = (std::string)getenv("DMPSWSYS")+"/share/Geometry/"+gdmlFile.string()+"/DAMPE.gdml";
  }
  char *dirTmp = getcwd(NULL,NULL);
  DmpLogInfo<<"Reading GDML:\t"<<gdmlFile.string()<<DmpLogEndl;
  chdir(gdmlFile.parent_path().string().c_str());
  fParser->Read(gdmlFile.filename().string().c_str());
  fPhyVolume = fParser->GetWorldVolume();
  fPhyVolume->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  chdir(dirTmp);

  //AddPhotonGenerator(); // must before adjustment

  Adjustment();

  // *
  // *  TODO: set SD of SubDet at here
  // *
  G4SDManager *mgrSD = G4SDManager::GetSDMpointer();
  if(fParser->GetVolume("Psd_Det")){
    DmpLogInfo<<"Setting Sensitive Detector of Psd"<<DmpLogEndl;
    //mgrSD->AddNewDetector(fPsdSD);
    //fParser->GetVolume("Psd_Strip0")->SetSensitiveDetector(fPsdSD);
    //fParser->GetVolume("Psd_Strip1")->SetSensitiveDetector(fPsdSD);
  }
  if(fParser->GetVolume("Stk_Det")){
          /*
    DmpLogInfo<<" Setting Sensitive Detector of Stk"<<DmpLogEndl;
    mgrSD->AddNewDetector(fStkSD);
    fParser->GetVolume("Stk_Strip")->SetSensitiveDetector(fStkSD);
    */
  }
  if(fParser->GetVolume("Bgo_Det")){
    DmpLogInfo<<"Setting Sensitive Detector of Bgo"<<DmpLogEndl;
    mgrSD->AddNewDetector(fBgoSD);
    fParser->GetVolume("Bgo_Bar")->SetSensitiveDetector(fBgoSD);
  }
  if(fParser->GetVolume("Nud_Det")){
    DmpLogInfo<<"Setting Sensitive Detector of Nud"<<DmpLogEndl;
    //mgrSD->AddNewDetector(fNudSD);
    //fParser->GetVolume("Nud_Block0")->SetSensitiveDetector(fNudSD);
    //fParser->GetVolume("Nud_Block1")->SetSensitiveDetector(fNudSD);
    //fParser->GetVolume("Nud_Block2")->SetSensitiveDetector(fNudSD);
    //fParser->GetVolume("Nud_Block3")->SetSensitiveDetector(fNudSD);
  }
  SetAncillarySD();
  return fPhyVolume;
}

//-------------------------------------------------------------------
void DmpSimDetector::Adjustment()const{
  short nCmd = fMetadata->OptionSize();
  G4VPhysicalVolume *PV = PV = G4PhysicalVolumeStore::GetInstance()->GetVolume("Ancillary_Det_PV",false);
  if(PV == 0 ) return;
  PV->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  for(short i =0; i<nCmd;++i){
    std::string command = fMetadata->GetCommand(i);
    if(command == "BT/DAMPE/Rotation"){
      double rad = 0;
      std::istringstream iss(fMetadata->GetValue(command));
      iss>>rad;
      rad = rad / 180 * 3.141592653;
      static G4RotationMatrix rot;
      rot.rotateY(-rad);
      PV->SetRotation(&rot);
    }else if(command == "BT/DAMPE/Translation"){
      double x=0,y=0.,z=0.0;
      std::istringstream iss(fMetadata->GetValue(command));
      iss>>x>>y>>z;
      G4ThreeVector move(x,y,z);
      G4ThreeVector par = PV->GetTranslation();
      par -= move;
      PV->SetTranslation(par);
    }else if(command == "BT/Magnetic"){
      double x=0,y=0.,z=0.0;
      std::istringstream iss(fMetadata->GetValue(command));
      iss>>x>>y>>z;
      ResetMagnetic(x,y,z);
    }else if(command == "BT/PbGlass/Delete"){
      if(G4PhysicalVolumeStore::GetInstance()->GetVolume("Pb_glass_PV",false)){
        G4PhysicalVolumeStore::GetInstance()->DeRegister(G4PhysicalVolumeStore::GetInstance()->GetVolume("Pb_glass_PV",false));
      }
    }
  }

}

//-------------------------------------------------------------------
#include "G4NistManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
void DmpSimDetector::AddPhotonGenerator()const{
  // get logical volume of ancillary detector
  G4LogicalVolume *ancillaryDet_LV = G4LogicalVolumeStore::GetInstance()->GetVolume("Ancillary_Det",false);
  if(not ancillaryDet_LV){
    return;
  }
  std::string command;
  // get size of PhotonGenerator
  double size_x=50,size_y=50,size_z=0.5;    // unit mm.
  command = "BT/PhotonGenerator/Thickness";
  if(fMetadata->HasCommand(command)){
    std::istringstream iss(fMetadata->GetValue(command));
    iss>>size_z;
  }
  command = "BT/PhotonGenerator/Size";
  if(fMetadata->HasCommand(command)){
    std::istringstream iss(fMetadata->GetValue(command));
    iss>>size_x>>size_y>>size_z;
  }
  // get material of PhotonGenerator
  G4NistManager *materialMgr = G4NistManager::Instance();
  G4Material *mat = materialMgr->FindOrBuildMaterial("G4_Galactic");
  command = "BT/PhotonGenerator/Material";
  if(fMetadata->HasCommand(command)){
    mat = materialMgr->FindOrBuildMaterial(fMetadata->GetValue(command).c_str());
  }
  // position of photon generator
  double pos_x=-300.0, pos_y=0.0, pos_z=-8000.0;
  command = "BT/PhotonGenerator/Offset";
  if(fMetadata->HasCommand(command)){
    double z =0.0;
    std::istringstream iss(fMetadata->GetValue(command));
    iss>>z;
    pos_z += z;
  }
  command = "BT/PhotonGenerator/Position";
  if(fMetadata->HasCommand(command)){
    std::istringstream iss(fMetadata->GetValue(command));
    iss>>pos_x>>pos_y>>pos_z;
  }
  std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<") material = "<<mat->GetName()<<std::endl;
  std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<") position = "<<pos_x<<"\t"<<pos_y<<"\t"<<pos_z<<std::endl;
  std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<") size = "<<size_x<<"\t"<<size_y<<"\t"<<size_z<<std::endl;
  //add PhotonGenerator volume
  G4VSolid *photonGeneratorSolid = new G4Box("PhotonGeneratorBox",size_x*mm,size_y*mm,size_z*mm);
  G4LogicalVolume *photonGeneratorLog = new G4LogicalVolume(photonGeneratorSolid,mat,"PhotonGeneratorLog");
  G4VPhysicalVolume *fPhotonGeneratorPhy = new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y,pos_z),photonGeneratorLog,"PhotonGenerator",ancillaryDet_LV,false,0);
}

//-------------------------------------------------------------------
#include "DmpSimMagneticField.h"
void DmpSimDetector::ResetMagnetic(const double &x,const double &y,const double &z)const{
  G4LogicalVolume *LV = G4LogicalVolumeStore::GetInstance()->GetVolume("B_field",false);
  if(LV){
    static DmpSimMagneticField magneticField(x,y,z);
    static G4FieldManager fieldMgr;
    fieldMgr.SetDetectorField(&magneticField);
    fieldMgr.CreateChordFinder(&magneticField);
    LV->SetFieldManager(&fieldMgr,true);
  }
}

//-------------------------------------------------------------------
#include "DmpSimSSD_SD.h"
void DmpSimDetector::SetAncillarySD(){
  G4SDManager *mgrSD = G4SDManager::GetSDMpointer();
  if(G4LogicalVolumeStore::GetInstance()->GetVolume("AMS_LDR_1",false)){
    DmpLogInfo<<"Setting Sensitive Detector of SSD_1"<<DmpLogEndl;
    static DmpSimSSD_SD AMS_LDR_1_SD("1");
    mgrSD->AddNewDetector(&AMS_LDR_1_SD);
    fParser->GetVolume("AMS_LDR_1")->SetSensitiveDetector(&AMS_LDR_1_SD);
  }
  if(G4LogicalVolumeStore::GetInstance()->GetVolume("AMS_LDR_2",false)){
    DmpLogInfo<<"Setting Sensitive Detector of SSD_2"<<DmpLogEndl;
    static DmpSimSSD_SD AMS_LDR_2_SD("2");
    mgrSD->AddNewDetector(&AMS_LDR_2_SD);
    fParser->GetVolume("AMS_LDR_2")->SetSensitiveDetector(&AMS_LDR_2_SD);
  }
  if(G4LogicalVolumeStore::GetInstance()->GetVolume("AMS_LDR_3",false)){
    DmpLogInfo<<"Setting Sensitive Detector of SSD_3"<<DmpLogEndl;
    static DmpSimSSD_SD AMS_LDR_3_SD("3");
    mgrSD->AddNewDetector(&AMS_LDR_3_SD);
    fParser->GetVolume("AMS_LDR_3")->SetSensitiveDetector(&AMS_LDR_3_SD);
  }
  if(G4LogicalVolumeStore::GetInstance()->GetVolume("AMS_LDR_4",false)){
    DmpLogInfo<<"Setting Sensitive Detector of SSD_4"<<DmpLogEndl;
    static DmpSimSSD_SD AMS_LDR_4_SD("4");
    mgrSD->AddNewDetector(&AMS_LDR_4_SD);
    fParser->GetVolume("AMS_LDR_4")->SetSensitiveDetector(&AMS_LDR_4_SD);
  }
}

