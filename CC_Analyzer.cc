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
  TH1F *hHcalPasE = new TH1F("hHcalPasE","energy HCAL passive /true ",600,0.,1.1);
  TH1F *hHcalActE = new TH1F("hHcalActE","ioninging HCAL energy active /true",600,0.,0.1);
  TH2F *hEcalHcal = new TH2F("hEcalHcal","ecal versus hcal passive /true ",600,0.,1.1,600,0.,1.1);
  TH2F *hEcalEcalem = new TH2F("hEcalEcalem","ecal versus ecal em /true",600,0.,1.1,600,0.,1.1);
  TH2F *hEcalfEcalem = new TH2F("hEcalfEcalem","em frac versus ecal/true ",600,0.,1.1,600,0.,1.1);
  TH1F *hcerTimingf= new TH1F("hcerTimingf","number cerenkov photons timing front",1000,0.,1000.);
  TH1F *hcerTimingr= new TH1F("hcerTimingr","number cerenkov photons timing read",1000,0.,1000.);
  TH1F *hcerECALf= new TH1F("hcerECALf","number cerenkov photons ECAL front",1000,0.,1000.);
  TH1F *hcerECALr= new TH1F("hcerECALr","number cerenkov photons ECAL read",1000,0.,1000.);
  TH1F *hcerHCAL= new TH1F("hcerHCAL","number cerenkov photons HCAL",1000,0.,1000.);
  TH1F *hHCALActPas = new TH1F("hHCALActPas","ratio of passive to active HCAL energy",1000,0.,500.);

  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");

  vector<float> *inputMomentum = new vector<float>;

  float depositedEnergyEscapeWorld;

  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedEnergyTiming_f,depositedEnergyTiming_r;
  float depositedEnergyECAL_f,depositedEnergyECAL_r;
  float depositedEnergyHCALAct,depositedEnergyHCALPas;
  float depositedEnergyServices;
  float depositedEnergyTimingGap;
  float depositedEnergyEcalGap;
  float depositedEnergyEcalDet;
  float depositedEnergySolenoid;

  float depositedElecEnergyTotal,depositedElecEnergyWorld;
  float depositedElecEnergyTiming_f,depositedElecEnergyTiming_r;
  float depositedElecEnergyECAL_f,depositedElecEnergyECAL_r;
  float depositedElecEnergyHCALAct,depositedElecEnergyHCALPas;
  float depositedElecEnergyServices;
  float depositedElecEnergyTimingGap;
  float depositedElecEnergyEcalGap;
  float depositedElecEnergyEcalDet;
  float depositedElecEnergySolenoid;


  float depositedIonEnergyTotal,depositedIonEnergyWorld;
  float depositedIonEnergyTiming_f,depositedIonEnergyTiming_r;
  float depositedIonEnergyECAL_f,depositedIonEnergyECAL_r;
  float depositedIonEnergyHCALAct,depositedIonEnergyHCALPas;
  float depositedIonEnergyServices;
  float depositedIonEnergyTimingGap;
  float depositedIonEnergyEcalGap;
  float depositedIonEnergyEcalDet;
  float depositedIonEnergySolenoid;

  int tot_phot_cer_Timing_f,tot_phot_cer_Timing_r,tot_phot_cer_ECAL_f,tot_phot_cer_ECAL_r,tot_phot_cer_HCAL;

  t1->SetBranchAddress("inputMomentum",&inputMomentum);

  t1->SetBranchAddress("depositedEnergyEscapeWorld",&depositedEnergyEscapeWorld);

  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyWorld",&depositedEnergyWorld);
  t1->SetBranchAddress("depositedEnergyTiming_f",&depositedEnergyTiming_f);
  t1->SetBranchAddress("depositedEnergyTiming_r",&depositedEnergyTiming_r);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_r",&depositedEnergyECAL_r);
  t1->SetBranchAddress("depositedEnergyHCALAct",&depositedEnergyHCALAct);
  t1->SetBranchAddress("depositedEnergyHCALPas",&depositedEnergyHCALPas);
  t1->SetBranchAddress("depositedEnergyServices",&depositedEnergyServices);
  t1->SetBranchAddress("depositedEnergyTimingGap",&depositedEnergyTimingGap);
  t1->SetBranchAddress("depositedEnergyEcalGap",&depositedEnergyEcalGap);
  t1->SetBranchAddress("depositedEnergyEcalDet",&depositedEnergyEcalDet);
  t1->SetBranchAddress("depositedEnergySolenoid",&depositedEnergySolenoid);

  t1->SetBranchAddress("depositedElecEnergyTotal",&depositedElecEnergyTotal);
  t1->SetBranchAddress("depositedElecEnergyWorld",&depositedElecEnergyWorld);
  t1->SetBranchAddress("depositedElecEnergyTiming_f",&depositedElecEnergyTiming_f);
  t1->SetBranchAddress("depositedElecEnergyTiming_r",&depositedElecEnergyTiming_r);
  t1->SetBranchAddress("depositedElecEnergyECAL_f",&depositedElecEnergyECAL_f);
  t1->SetBranchAddress("depositedElecEnergyECAL_r",&depositedElecEnergyECAL_r);
  t1->SetBranchAddress("depositedElecEnergyHCALAct",&depositedElecEnergyHCALAct);
  t1->SetBranchAddress("depositedElecEnergyHCALPas",&depositedElecEnergyHCALPas);
  t1->SetBranchAddress("depositedElecEnergyServices",&depositedElecEnergyServices);
  t1->SetBranchAddress("depositedElecEnergyTimingGap",&depositedElecEnergyTimingGap);
  t1->SetBranchAddress("depositedElecEnergyEcalGap",&depositedElecEnergyEcalGap);
  t1->SetBranchAddress("depositedElecEnergyEcalDet",&depositedElecEnergyEcalDet);
  t1->SetBranchAddress("depositedElecEnergySolenoid",&depositedElecEnergySolenoid);


  t1->SetBranchAddress("depositedIonEnergyTotal",&depositedIonEnergyTotal);
  t1->SetBranchAddress("depositedIonEnergyWorld",&depositedIonEnergyWorld);
  t1->SetBranchAddress("depositedIonEnergyTiming_f",&depositedIonEnergyTiming_f);
  t1->SetBranchAddress("depositedIonEnergyTiming_r",&depositedIonEnergyTiming_r);
  t1->SetBranchAddress("depositedIonEnergyECAL_f",&depositedIonEnergyECAL_f);
  t1->SetBranchAddress("depositedIonEnergyECAL_r",&depositedIonEnergyECAL_r);
  t1->SetBranchAddress("depositedIonEnergyHCALAct",&depositedIonEnergyHCALAct);
  t1->SetBranchAddress("depositedIonEnergyHCALPas",&depositedIonEnergyHCALPas);
  t1->SetBranchAddress("depositedIonEnergyServices",&depositedIonEnergyServices);
  t1->SetBranchAddress("depositedIonEnergyTimingGap",&depositedIonEnergyTimingGap);
  t1->SetBranchAddress("depositedIonEnergyEcalGap",&depositedIonEnergyEcalGap);
  t1->SetBranchAddress("depositedIonEnergyEcalDet",&depositedIonEnergyEcalDet);
  t1->SetBranchAddress("depositedIonEnergySolenoid",&depositedIonEnergySolenoid);

  t1->SetBranchAddress("tot_phot_cer_Timing_f",&tot_phot_cer_Timing_f);
  t1->SetBranchAddress("tot_phot_cer_Timing_r",&tot_phot_cer_Timing_r);
  t1->SetBranchAddress("tot_phot_cer_ECAL_f",&tot_phot_cer_ECAL_f);
  t1->SetBranchAddress("tot_phot_cer_ECAL_r",&tot_phot_cer_ECAL_r);
  t1->SetBranchAddress("tot_phot_cer_HCAL",&tot_phot_cer_HCAL);

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
    std::cout<<"timing front energy deposited is "<<depositedEnergyTiming_f<<std::endl;
    std::cout<<"timing rear energy deposited is "<<depositedEnergyTiming_r<<std::endl;
    std::cout<<"ECAL front energy deposited is "<<depositedEnergyECAL_f<<std::endl;
    std::cout<<"ECAL rear energy deposited is "<<depositedEnergyECAL_r<<std::endl;
    std::cout<<"HCAL act  energy deposited is "<<depositedEnergyHCALAct<<std::endl;
    std::cout<<"HCAL Pas  energy deposited is "<<depositedEnergyHCALPas<<std::endl;
    std::cout<<"escape energy deposited is "<<depositedEnergyEscapeWorld<<std::endl;
    std::cout<<"Services energy deposited is "<<depositedEnergyServices<<std::endl;
    std::cout<<"Timing Gap energy deposited is "<<depositedEnergyTimingGap<<std::endl;
    std::cout<<"Ecal gap energy deposited is "<<depositedEnergyEcalGap<<std::endl;
    std::cout<<"Ecal det energy deposited is "<<depositedEnergyEcalDet<<std::endl;
    std::cout<<"Solenoid energy deposited is "<<depositedEnergySolenoid<<std::endl;
    */
    float eee=depositedEnergyTiming_f+depositedEnergyTiming_r+depositedEnergyECAL_f+depositedEnergyECAL_r+depositedEnergyHCALAct+depositedEnergyHCALPas+depositedEnergyWorld+depositedEnergyServices+depositedEnergyTimingGap+depositedEnergyEcalGap+depositedEnergyEcalDet+depositedEnergySolenoid;
    float fff=depositedEnergyTotal+depositedEnergyEscapeWorld;
    float ecaltotal=depositedEnergyTiming_f+depositedEnergyTiming_r+depositedEnergyECAL_f+depositedEnergyECAL_r;
    //std::cout<<" sum in detectors is "<<eee<<std::endl;
    //std::cout<<" deposited plus escaped is "<<fff<<std::endl;

    float estE= (depositedEnergyTiming_f-depositedIonEnergyTiming_f)+
      (depositedEnergyTiming_r-depositedIonEnergyTiming_r)+
      (depositedEnergyECAL_f-depositedIonEnergyECAL_f)+
      (depositedEnergyECAL_r-depositedIonEnergyECAL_r)+
      sampling_fraction*(depositedEnergyHCALAct-depositedIonEnergyHCALAct);
    std::cout<<"est E is "<<estE<<std::endl;


    hTotalE->Fill(depositedEnergyTotal/trueE);
    hestE->Fill(estE/trueE);
    hHcalPasE->Fill(depositedEnergyHCALPas/trueE);
    hHcalActE->Fill((depositedEnergyHCALAct-depositedIonEnergyHCALAct)/trueE);
    hEcalHcal->Fill(depositedEnergyHCALPas/trueE,ecaltotal/trueE);
    float yyy=depositedEnergyECAL_f+depositedEnergyECAL_r;
    float yy2=depositedElecEnergyECAL_f+depositedElecEnergyECAL_r;
    hEcalEcalem->Fill(yyy/trueE,yy2/trueE);
    if(yyy>0) hEcalfEcalem->Fill(yyy/trueE,yy2/yyy);
    if(depositedEnergyHCALAct>0) hHCALActPas->Fill(depositedEnergyHCALPas/depositedEnergyHCALAct);

    hcerTimingf->Fill(tot_phot_cer_Timing_f);
    hcerTimingr->Fill(tot_phot_cer_Timing_r);
    hcerECALf->Fill(tot_phot_cer_ECAL_f);
    hcerECALr->Fill(tot_phot_cer_ECAL_r);
    hcerHCAL->Fill(tot_phot_cer_HCAL);

		    
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  hestE->Write();
  hHcalPasE->Write();
  hHcalActE->Write();
  hEcalHcal->Write();
  hEcalEcalem->Write();
  hEcalfEcalem->Write();
  hcerTimingf->Write();
  hcerTimingr->Write();
  hcerECALf->Write();
  hcerECALr->Write();
  hcerHCAL->Write();
  hHCALActPas->Write();
  out->Close();

}


void CC_Analyzer(bool debug) {
  if(debug) {  
    std::cout<<"running on temp.root"<<std::endl;
    CC_Ana("temp.root","temp_hists.root");
  }
  else {
    std::cout<<"running on normal files"<<std::endl;
    //  CC_Ana("pions_1GeV.root","pions_1GeV_hists.root");
    //CC_Ana("pions_5GeV.root","pions_5GeV_hists.root");
  CC_Ana("pions_10GeV.root","pions_10GeV_hists.root");
  CC_Ana("pions_20GeV.root","pions_20GeV_hists.root");
  CC_Ana("pions_50GeV.root","pions_50GeV_hists.root");
  //CC_Ana("pions_100GeV.root","pions_100GeV_hists.root");
  }
}
