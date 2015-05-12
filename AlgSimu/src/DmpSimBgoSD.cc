/*
 *  $Id: DmpSimBgoSD.cc, 2014-10-06 17:13:16 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/03/2014
*/

#include "G4Step.hh"
#include "G4TouchableHistory.hh"

#include <stdlib.h>     // getenv()

#include "DmpSimBgoSD.h"
#include "DmpEvtBgoHits.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "DmpLog.h"
#include "TF1.h"
#include "TMath.h"

//-------------------------------------------------------------------
//define laugaufun
Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}
//-------------------------------------------------------------------
DmpSimBgoSD::DmpSimBgoSD()
 :G4VSensitiveDetector("BgoSD"),
  fEvtMCBgo(0)
{
  _CaliParPath = (std::string)getenv("DMPSWWORK")+"/share/";
  GetMipPar();
  GetAttPar();
  for (int i=0;i<616;i++){
    RanGaus[i] = new TRandom3();
    UInt_t seed=time(0)+i;
    RanGaus[i]->SetSeed(seed);

    
  }
  fEvtMCBgo = new DmpEvtBgoHits();
  gDataBuffer->RegisterObject("Event/MCTruth/Bgo",fEvtMCBgo,"DmpEvtBgoHits");
  fDigitBgo = new DmpEvtBgoHits();
  gDataBuffer->RegisterObject("Event/MCTruth/BgoFDigit",fDigitBgo,"DmpEvtBgoHits");
}

//-------------------------------------------------------------------
DmpSimBgoSD::~DmpSimBgoSD(){
}

//-------------------------------------------------------------------
#include <boost/lexical_cast.hpp>
G4bool DmpSimBgoSD::ProcessHits(G4Step *aStep,G4TouchableHistory*){
  G4TouchableHistory *theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  std::string barName = theTouchable->GetVolume(1)->GetName();
  barName.assign(barName.end()-4,barName.end());        // get ID
  short barID = boost::lexical_cast<short>(barName);
  barID = DmpBgoBase::ConstructGlobalBarID(barID/100,barID%100);
  if (aStep->GetTotalEnergyDeposit()>0){
    AddThisG4Hit(barID,aStep->GetTotalEnergyDeposit()/MeV,aStep->GetPreStepPoint()->GetPosition());
    Eny2ADC(barID,aStep->GetTotalEnergyDeposit()/MeV,(aStep->GetPreStepPoint()->GetPosition()).x()/mm,(aStep->GetPreStepPoint()->GetPosition()).y()/mm);
  }
  return true;
}

//-------------------------------------------------------------------
void DmpSimBgoSD::Initialize(G4HCofThisEvent*){
  fEvtMCBgo->Reset();
  fDigitBgo->Reset();
  Position.Clear();
  memset(TotalE,0,sizeof(TotalE));
}

//-------------------------------------------------------------------
void DmpSimBgoSD::EndOfEvent(G4HCofThisEvent* HCE){
  //sampling & save
  Sampling();
}

//-------------------------------------------------------------------
void DmpSimBgoSD::AddThisG4Hit(const short &id,const double &e,const G4ThreeVector &in){
  TVector3 pos(in.x(),in.y(),in.z());
  for(size_t i=0;i<fEvtMCBgo->fGlobalBarID.size();i++){
    if(fEvtMCBgo->fGlobalBarID.at(i) == id){
      double totE = e + fEvtMCBgo->fEnergy.at(i);
      fEvtMCBgo->fPosition.at(i) = fEvtMCBgo->fPosition.at(i) * (fEvtMCBgo->fEnergy.at(i) / totE);
      fEvtMCBgo->fPosition.at(i) += pos * (e / totE);
      fEvtMCBgo->fEnergy.at(i) = totE;
      return;   // if found gid, update and return
    }
  }
  // if not, creat a new one

  fEvtMCBgo->fGlobalBarID.push_back(id);
  fEvtMCBgo->fEnergy.push_back(e);
  fEvtMCBgo->fES0.push_back(0);  // TODO, to two sides?
  fEvtMCBgo->fES1.push_back(0);
  fEvtMCBgo->fPosition.push_back(pos);
}
//-------------------------------------------------------------------
void DmpSimBgoSD::GetAttPar(){

  //Get Attenuation coefficients 
  ifstream Apar;
  std::string fn = _CaliParPath+"Simulation/AttPar";
  Apar.open(fn.c_str());
   if(!Apar.good()){
    std::cout<<"Can not open Att Par file!"<<std::endl;
    exit(1);
  } 
  int nGbar=14*22;
   for(int i=0; i<nGbar;i++){
     for(int j=0 ;j<2;j++){ 
      Apar>>AttPar[i][j]; 
      //std::cout<<AttPar[i][j]<<"\t";
    }
    //std::cout<<std::endl;
  }
  Apar.close();
}
//-------------------------------------------------------------------
void DmpSimBgoSD::GetMipPar(){

  //Get MIPs parameters
  ifstream Mpar;
  std::string fn = _CaliParPath+"Simulation/MIPsPar";
  Mpar.open(fn.c_str());
  if(!Mpar.good()){
    std::cout<<"Can not open MIPs Par file!"<<std::endl;
    exit(1);
  }
  TF1 *myMIPs=new TF1("myMIPs",langaufun,0.,5800.,4);
  myMIPs->SetParNames("Width","MP","Area","GSigma");

  int nGpmt=14*2*22;
  for(int i=0;i<nGpmt;i++){
    for(int j=0;j<4;j++){
      Mpar>>MipPar[i][j];  
      myMIPs->SetParameters(MipPar[i]);
      MipPar[i][1]=myMIPs->GetMaximumX(0.8*MipPar[i][1],1.5*MipPar[i][1]);
      //std::cout<<MipPar[i][j]<<"\t";
    }
    //std::cout<<std::endl;
  }
  Mpar.close();
}
//-------------------------------------------------------------------
void DmpSimBgoSD::Eny2ADC(const short &id, const double &e, const double &x,const double &y){
  //step Energy to ADC
  //Get iGpmt, iGbar
  short layer=DmpBgoBase::GetLayerID(id);
  short bar=DmpBgoBase::GetBarID(id);
  short iGbar=layer*22+bar;
  short iGpmt[2]={layer*2*22+bar,(layer*2+1)*22+bar};

  //Set Att coe, MIPs coversion ratio;
  double Dis[2]; //distance between energy deposition and BGO end.
  if(layer%2==0){
    Dis[1]=x+300;
    Dis[0]=300-x;
  }
  else{
    Dis[1]=y+300;
    Dis[0]=300-y;
  }
  double AttCoe[2];  //AttPar[][0]=2/lambda; lambda :cm Dis mm;
  AttCoe[0]=1/TMath::Exp(AttPar[iGbar][0]*Dis[0]/2/10);
  AttCoe[1]=1/TMath::Exp(AttPar[iGbar][0]*Dis[1]/2/10);
  double AttHit[2];
  AttHit[0]=e*AttCoe[0];
  AttHit[1]=e*AttCoe[1];
  TotalE[iGpmt[0]]+=AttHit[0];
  TotalE[iGpmt[1]]+=AttHit[1];
}
//-------------------------------------------------------------------
void DmpSimBgoSD::Sampling(){
  //Sampling with calibrated paramneters
  //Get iGpmt, iGbar
  for(short layer=0;layer<14;layer++){
    for(short bar=0;bar<22;bar++){
      short iGbar=layer*22+bar;
      short iGpmt[2]={layer*2*22+bar,(layer*2+1)*22+bar};
	  short gid_bar=DmpBgoBase::ConstructGlobalBarID(layer,bar);
      if(TotalE[iGpmt[0]]>0){
//	  std::cout<<iGpmt[0]<<std::endl;
//	  std::cout<<TotalE[iGpmt[0]]<<std::endl;
//	  std::cout<<iGpmt[1]<<std::endl;
//	  std::cout<<TotalE[iGpmt[1]]<<std::endl;
        double MipCov[2];
//        MipCov[0]=MipPar[iGpmt[0]][1]*TMath::Exp(AttPar[iGbar][0]*30/2)/22.5;   //non-att normalized ADC counts/MeV;
//        MipCov[1]=MipPar[iGpmt[1]][1]*TMath::Exp(AttPar[iGbar][0]*30/2)/22.5;   //non-att normalized ADC counts/MeV;
//        double Mean[2]={TotalE[iGpmt[0]]*MipCov[0],TotalE[iGpmt[1]]*MipCov[1]};
//        double Sigma[2]={MipPar[iGpmt[0]][3]*TMath::Sqrt(Mean[0]/MipPar[iGpmt[0]][1]),MipPar[iGpmt[1]][3]*TMath::Sqrt(Mean[1]/MipPar[iGpmt[1]][1])};
        MipCov[0]=TMath::Exp(AttPar[iGbar][0]*30/2);   //non-att normalized ADC counts/MeV;
        MipCov[1]=TMath::Exp(AttPar[iGbar][0]*30/2);   //non-att normalized ADC counts/MeV;
        double Mean[2]={TotalE[iGpmt[0]]*MipCov[0],TotalE[iGpmt[1]]*MipCov[1]};
        double Sigma[2]={22.5*MipPar[iGpmt[0]][3]/MipPar[iGpmt[0]][1]*TMath::Sqrt(Mean[0]/22.5),22.5*MipPar[iGpmt[1]][3]/MipPar[iGpmt[1]][1]*TMath::Sqrt(Mean[1]/22.5)};
	  
        //  TRandom *s0=new TRandom();
        //  TRandom *s1=new TRandom();
	    double ES0=RanGaus[iGpmt[0]]->Gaus(Mean[0],Sigma[0]);
	    double ES1=RanGaus[iGpmt[1]]->Gaus(Mean[1],Sigma[1]);
        if(ES0>0&&ES1>0){
	  fDigitBgo->fGlobalBarID.push_back(gid_bar);
          fDigitBgo->fES0.push_back(ES0);
          fDigitBgo->fES1.push_back(ES1);
          fDigitBgo->fEnergy.push_back(TMath::Sqrt(ES0*ES1));
	  double pos=(AttPar[iGbar][0]*TMath::Log(ES0/ES1)+AttPar[iGbar][1])*10-300;//unit: mm
	  if(layer%2==0){
	    Position.SetX(pos);	  
	    double yy=DmpBgoBase::Parameter()->BarCenter(gid_bar).y();
            Position.SetY(yy);
	  }
	  else{
	    Position.SetY(pos);
	    double xx=DmpBgoBase::Parameter()->BarCenter(gid_bar).x();
            Position.SetX(xx);
	  }
	  double zz=DmpBgoBase::Parameter()->LayerCenter(gid_bar).z();
          Position.SetZ(zz);
          fDigitBgo->fPosition.push_back(Position);
        }
      }
    }
  }
}
