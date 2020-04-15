#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

const float sampling_fraction=100.;

void CC_Ana(const char* inputfilename,const char* outputfilename) {

  TH1F *hTotalE = new TH1F("hTotalE","energy total world /true",600,0.,1.1);
  TH1F *hestE = new TH1F("hestE","energy estimated /true",600,0.,1.1);
  TH2F *hEcalEcalem = new TH2F("hEcalEcalem","ecal versus ecal em /true",600,0.,1.1,600,0.,1.1);
  TH2F *hEcalfEcalem = new TH2F("hEcalfEcalem","em frac versus ecal/true ",600,0.,1.1,600,0.,1.1);
  TH1F *hcerTimingf= new TH1F("hcerTimingf","number cerenkov photons timing front",1000,0.,1000.);
  TH1F *hcerTimingr= new TH1F("hcerTimingr","number cerenkov photons timing read",1000,0.,1000.);
  TH1F *hcerECALf= new TH1F("hcerECALf","number cerenkov photons ECAL front",1000,0.,1000.);
  TH1F *hcerECALr= new TH1F("hcerECALr","number cerenkov photons ECAL read",1000,0.,1000.);

  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");

  vector<float> *inputMomentum = new vector<float>;

  float depositedEnergyEscapeWorld;

  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedEnergyECAL_f,depositedEnergyECAL_r;
  float depositedEnergyEcalGap;
  float depositedEnergyEcalDet;


  float depositedElecEnergyTotal,depositedElecEnergyWorld;
  float depositedElecEnergyECAL_f,depositedElecEnergyECAL_r;
  float depositedElecEnergyEcalGap;
  float depositedElecEnergyEcalDet;



  float depositedIonEnergyTotal,depositedIonEnergyWorld;
  float depositedIonEnergyECAL_f,depositedIonEnergyECAL_r;
  float depositedIonEnergyEcalGap;
  float depositedIonEnergyEcalDet;

  int tot_phot_cer_ECAL_f,tot_phot_cer_ECAL_r;

  t1->SetBranchAddress("inputMomentum",&inputMomentum);

  t1->SetBranchAddress("depositedEnergyEscapeWorld",&depositedEnergyEscapeWorld);

  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyWorld",&depositedEnergyWorld);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_r",&depositedEnergyECAL_r);
  t1->SetBranchAddress("depositedEnergyEcalGap",&depositedEnergyEcalGap);
  t1->SetBranchAddress("depositedEnergyEcalDet",&depositedEnergyEcalDet);


  t1->SetBranchAddress("depositedElecEnergyTotal",&depositedElecEnergyTotal);
  t1->SetBranchAddress("depositedElecEnergyWorld",&depositedElecEnergyWorld);
  t1->SetBranchAddress("depositedElecEnergyECAL_f",&depositedElecEnergyECAL_f);
  t1->SetBranchAddress("depositedElecEnergyECAL_r",&depositedElecEnergyECAL_r);
  t1->SetBranchAddress("depositedElecEnergyEcalGap",&depositedElecEnergyEcalGap);
  t1->SetBranchAddress("depositedElecEnergyEcalDet",&depositedElecEnergyEcalDet);



  t1->SetBranchAddress("depositedIonEnergyTotal",&depositedIonEnergyTotal);
  t1->SetBranchAddress("depositedIonEnergyWorld",&depositedIonEnergyWorld);
  t1->SetBranchAddress("depositedIonEnergyECAL_f",&depositedIonEnergyECAL_f);
  t1->SetBranchAddress("depositedIonEnergyECAL_r",&depositedIonEnergyECAL_r);
  t1->SetBranchAddress("depositedIonEnergyEcalGap",&depositedIonEnergyEcalGap);
  t1->SetBranchAddress("depositedIonEnergyEcalDet",&depositedIonEnergyEcalDet);



  t1->SetBranchAddress("tot_phot_cer_ECAL_f",&tot_phot_cer_ECAL_f);
  t1->SetBranchAddress("tot_phot_cer_ECAL_r",&tot_phot_cer_ECAL_r);


  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    float trueE=9999999.;
    if((*inputMomentum)[3]>0) trueE=(*inputMomentum)[3];
    
    std::cout<<endl<<"event number "<<i<<std::endl;
    std::cout<<(*inputMomentum)[0]<<","<<(*inputMomentum)[1]<<","<<(*inputMomentum)[2]<<","<<(*inputMomentum)[3]<<std::endl;
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    /*
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    std::cout<<"world energy deposited is "<<depositedEnergyWorld<<std::endl;
    std::cout<<"ECAL front energy deposited is "<<depositedEnergyECAL_f<<std::endl;
    std::cout<<"ECAL rear energy deposited is "<<depositedEnergyECAL_r<<std::endl;
    std::cout<<"escape energy deposited is "<<depositedEnergyEscapeWorld<<std::endl;
    std::cout<<"Ecal gap energy deposited is "<<depositedEnergyEcalGap<<std::endl;
    std::cout<<"Ecal det energy deposited is "<<depositedEnergyEcalDet<<std::endl;
    */
    float eee=depositedEnergyECAL_f+depositedEnergyECAL_r+depositedEnergyWorld+depositedEnergyEcalGap+depositedEnergyEcalDet;
    float fff=depositedEnergyTotal+depositedEnergyEscapeWorld;
    float ecaltotal=depositedEnergyECAL_f+depositedEnergyECAL_r;
    //std::cout<<" sum in detectors is "<<eee<<std::endl;
    //std::cout<<" deposited plus escaped is "<<fff<<std::endl;

    float estE= 
      (depositedEnergyECAL_f-depositedIonEnergyECAL_f)+
      (depositedEnergyECAL_r-depositedIonEnergyECAL_r);
    std::cout<<"est E is "<<estE<<std::endl;


    hTotalE->Fill(depositedEnergyTotal/trueE);
    hestE->Fill(estE/trueE);
    float yyy=depositedEnergyECAL_f+depositedEnergyECAL_r;
    float yy2=depositedElecEnergyECAL_f+depositedElecEnergyECAL_r;
    hEcalEcalem->Fill(yyy/trueE,yy2/trueE);
    if(yyy>0) hEcalfEcalem->Fill(yyy/trueE,yy2/yyy);

    hcerECALf->Fill(tot_phot_cer_ECAL_f);
    hcerECALr->Fill(tot_phot_cer_ECAL_r);


		    
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hestE->Write();
  hEcalEcalem->Write();
  hEcalfEcalem->Write();
  hcerECALf->Write();
  hcerECALr->Write();

  out->Close();

}


void CC_Analyzer(bool debug) {
  if(debug) {  
    std::cout<<"running on temp.root"<<std::endl;
    CC_Ana("temp.root","temp_hists.root");
  }
  else {
    std::cout<<"running on normal files"<<std::endl;
  CC_Ana("electrons_10GeV.root","electrons_10GeV_hists.root");
  CC_Ana("electrons_20GeV.root","electrons_20GeV_hists.root");
  CC_Ana("electrons_50GeV.root","electrons_50GeV_hists.root");

  }
}
