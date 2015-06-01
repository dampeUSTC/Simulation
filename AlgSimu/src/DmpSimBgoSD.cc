/*
 *  $Id: DmpSimBgoSD.cc, 2014-10-06 17:13:16 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/03/2014
*/

#include "G4Step.hh"
#include "G4TouchableHistory.hh"

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
  GetPedPar();
  GetMipPar();
  GetAttPar();
  GetDyPar();
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
  Apar.open("../CaliParameter/Attenuation/AttPar");
   if(!Apar.good()){
    std::cout<<"Can not open Att Par file!"<<std::endl;
    exit(1);
  } 
  int nGbar=14*22;
   for(int i=0; i<nGbar;i++){
      std::cout<<(int)(i/22)<<"  "<<i%22<<"\t\t";
     for(int j=0 ;j<2;j++){ 
      Apar>>AttPar[i][j]; 
    }
     // std::cout<<2/AttPar[i][0]*<<"\t"<<600/AttPar[i][1]<<"\t";
      std::cout<<2/AttPar[i][0]<<"\t";
    std::cout<<std::endl;
  }
  Apar.close();
}
//-------------------------------------------------------------------
void DmpSimBgoSD::GetPedPar(){

  //Get Pedestal parameters 
  ifstream Ppar;
  Ppar.open("../CaliParameter/Pedestal/PedPar");
   if(!Ppar.good()){
    std::cout<<"Can not open Pedestal Par file!"<<std::endl;
    exit(1);
  } 
   else{
  int nGdy=14*22*3*2;
   for(int i=0; i<nGdy;i++){
     for(int j=0 ;j<2;j++){ 
      Ppar>>PedPar[i][j]; 
  //    std::cout<<PedPar[i][j]<<"\t";
    }
  //  std::cout<<std::endl;
  }
  Ppar.close();
   }
}
//-------------------------------------------------------------------
void DmpSimBgoSD::GetDyPar(){

  //Get Pedestal parameters 
  ifstream Dpar;
  Dpar.open("../CaliParameter/DyCoe/DyPar");
   if(!Dpar.good()){
    std::cout<<"Can not open DynodeRatios Par file!"<<std::endl;
    exit(1);
  }  
  int nGdy=14*22*2;
   for (int i=0; i<nGdy;i++){
     Dpar>>DyPar25[i][0]>>DyPar25[i][1]>>DyPar58[i][0]>>DyPar58[i][1]; 
  }
  Dpar.close();
}
//-------------------------------------------------------------------
void DmpSimBgoSD::GetMipPar(){

  //Get MIPs parameters
  ifstream Mpar;
  Mpar.open("../CaliParameter/MIPs/MIPsPar_beam");
  if(!Mpar.good()){
    std::cout<<"Can not open MIPs Par file!"<<std::endl;
    exit(1);
  }
  TF1 *myMIPs=new TF1("myMIPs",langaufun,0.,5800.,4);
  myMIPs->SetParNames("Width","MP","Area","GSigma");

  int nGpmt=14*3*22;
/*  for(int i=0;i<nGpmt;i++){
      Mpar>>MipPar[i][1]>>MipPar[i][3]>>MipPar[i][0];
      if(MipPar[i][3]/MipPar[i][1]<0.1&&MipPar[i][1]<100){
      MipPar[i][3]=0.2*MipPar[i][1];
      }
      if(MipPar[i][3]<8){
      MipPar[i][3]=16;
      }
   //   myMIPs->SetParameters(MipPar[i]);
   //   MipPar[i][1]=myMIPs->GetMaximumX(0.8*MipPar[i][1],1.5*MipPar[i][1]);
      //std::cout<<MipPar[i][j]<<"\t";
    //std::cout<<std::endl;
  }*/
 //int layer,side,bar;
 // double peak,width,gsigma;
 double spar[6];
  for(int i=0;i<nGpmt;i++){
    for(int j=0;j<6;j++){
	    Mpar>>spar[j];
    }
    if(((int)spar[1])==2)continue;
    int ipmt=(spar[0]*2+spar[1])*22+spar[2];
    MipPar[ipmt][0]=spar[4];
    MipPar[ipmt][1]=spar[3];
    MipPar[ipmt][3]=spar[5];
    int idy=ipmt*3+2;
    std::cout<<spar[0]<<" "<<spar[1]<<" "<<spar[2]<<" ";
    std::cout<<" "<<23*3*PedPar[idy][1]/MipPar[ipmt][1]<<std::endl;
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
      short iGdy8[2]={iGpmt[0]*3+2,iGpmt[1]*3+2};
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
        //double Sigma[2]={22.5*MipPar[iGpmt[0]][3]/MipPar[iGpmt[0]][1]*TMath::Sqrt(Mean[0]/22.5),22.5*MipPar[iGpmt[1]][3]/MipPar[iGpmt[1]][1]*TMath::Sqrt(Mean[1]/22.5)};
        double Sigma[2];
	short usedy[2];
	for(int iside=0;iside<2;iside++){
	usedy[iside]=ChoiseDynode(Mean[iside],iGpmt[iside]);
	Sigma[iside]=SetSigma(iGpmt[iside],usedy[iside],Mean[iside]);
	} 
//	std::cout<<"BarEnergy"<<":s0 "<<Mean[0]<<", s1 "<<Mean[1]<<std::endl;
//	std::cout<<"BarSigma"<<":s0 "<<Sigma[0]<<", s1 "<<Sigma[1]<<std::endl;
//	std::cout<<"Bardynode"<<":s0 "<<usedy[0]<<", s1 "<<usedy[1]<<std::endl;
//	std::cout<<"MIPsPar_MPV"<<":s0 "<<MipPar[iGpmt[0]][1]<<", s1 "<<MipPar[iGpmt[1]][1]<<std::endl;
//	std::cout<<"MIPsPar_Gsigma"<<":s0 "<<MipPar[iGpmt[0]][3]<<", s1 "<<MipPar[iGpmt[1]][3]<<std::endl;

	  
        //  TRandom *s0=new TRandom();
        //  TRandom *s1=new TRandom();
	    double ES0=RanGaus[iGpmt[0]]->Gaus(Mean[0],Sigma[0]);
	    double ES1=RanGaus[iGpmt[1]]->Gaus(Mean[1],Sigma[1]);
            double cut0=22.5*3*PedPar[iGdy8[0]][1]/MipPar[iGpmt[0]][1];
            double cut1=22.5*3*PedPar[iGdy8[1]][1]/MipPar[iGpmt[1]][1];
           // std::cout<<"cut0: "<<cut0<<"MeV"<<std::endl;
           // std::cout<<"cut1: "<<cut1<<"MeV"<<std::endl;

        if(ES0>cut0&&ES1>cut1){
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
short DmpSimBgoSD::ChoiseDynode(double barE,short gid_pmt){//return the selected dynode
  if(barE<=10000./MipPar[gid_pmt][1]*22.5){return 8;}//set ADC cut at 10000
  else {
    double Maxdy5E=(10000.*DyPar58[gid_pmt][0]+DyPar58[gid_pmt][1])/MipPar[gid_pmt][1]*22.5;
    if(barE<=Maxdy5E){return 5;}
    else { 
      double Maxdy2E=((14000.*DyPar25[gid_pmt][0]+DyPar25[gid_pmt][1])*DyPar58[gid_pmt][0]+DyPar58[gid_pmt][1])/MipPar[gid_pmt][1]*22.5;
      if(barE<=Maxdy2E){return 2;}
      else return 0;///dynode 2 became saturated!!!!!!!!
    } 
   }
} 
double DmpSimBgoSD::SetSigma(short gid_pmt, short usedynode,double barE){
  int gid_dy8=gid_pmt*3+2;
  
  double AbsSigma=22.5*TMath::Sqrt(MipPar[gid_pmt][3]*MipPar[gid_pmt][3]-PedPar[gid_dy8][1]*PedPar[gid_dy8][1])/MipPar[gid_pmt][1]/TMath::Sqrt(22.5);
  if(usedynode==8){
    double noisedy8=22.5*PedPar[gid_dy8][1]/MipPar[gid_pmt][1];
    double statcdy8=AbsSigma*TMath::Sqrt(barE);
    double sigmady8=TMath::Sqrt(statcdy8*statcdy8+noisedy8*noisedy8);
    return sigmady8;
    }
  else {
   if(usedynode==5){
     int gid_dy5=gid_pmt*3+1;
     double noisedy5=22.5*PedPar[gid_dy5][1]*DyPar58[gid_pmt][0]/MipPar[gid_pmt][1];
     double statcdy5=AbsSigma*TMath::Sqrt(barE);
     double sigmady5=TMath::Sqrt(statcdy5*statcdy5+noisedy5*noisedy5);
     return sigmady5;
   } 
   else{
     int gid_dy2=gid_pmt*3;
     double noisedy2=22.5*PedPar[gid_dy2][1]*DyPar58[gid_pmt][0]*DyPar25[gid_pmt][0]/MipPar[gid_pmt][1];
     double statcdy2=AbsSigma*TMath::Sqrt(barE);
     double sigmady2=TMath::Sqrt(statcdy2*statcdy2+noisedy2*noisedy2);
     return sigmady2;
   }  
 } 
}  

