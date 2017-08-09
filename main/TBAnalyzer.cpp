#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TApplication.h"


unsigned int CH0 = 0;
unsigned int CH1 = 1;
unsigned int CH2 = 2;
unsigned int CH3 = 3;
unsigned int CH8 = 4;

//=== focusing on ch1 for now (bottom left) ===
unsigned int CH = CH2;
unsigned int CHref = CH2;

struct TreeVars
{  
  float t_beamE;
  unsigned int t_nCh;
  
  float* t_beamX;
  float* t_beamY;

  float* t_chargeSig;
  float* t_chargeTot;
  
  float* t_ampMax;
  float* t_timeMax;
  float* t_chi2Max;
  float* t_timeCFD;

  int t_nSampl;
  int* t_wfCh;
  float* t_wfTime;
  float* t_wfVal;

 };

void InitTreeVars(TTree* chain1, TreeVars& treeVars)
{
  chain1 -> SetBranchStatus("Energy",1); chain1 -> SetBranchAddress("Energy",&treeVars.t_beamE);
  chain1 -> SetBranchStatus("n_channels",1); chain1 -> SetBranchAddress("n_channels",&treeVars.t_nCh);
  
  treeVars.t_beamX = new float[2];
  treeVars.t_beamY = new float[2];
  chain1 -> SetBranchStatus("X",1); chain1 -> SetBranchAddress("X",treeVars.t_beamX);
  chain1 -> SetBranchStatus("Y",1); chain1 -> SetBranchAddress("Y",treeVars.t_beamY);
  
  treeVars.t_chargeSig = new float[treeVars.t_nCh];
  treeVars.t_chargeTot = new float[treeVars.t_nCh];
  chain1 -> SetBranchStatus("charge_sig",1); chain1 -> SetBranchAddress("charge_sig",treeVars.t_chargeSig);
  chain1 -> SetBranchStatus("charge_tot",1); chain1 -> SetBranchAddress("charge_tot",treeVars.t_chargeTot);

  treeVars.t_ampMax = new float[treeVars.t_nCh];
  treeVars.t_timeMax = new float[treeVars.t_nCh];
  treeVars.t_chi2Max = new float[treeVars.t_nCh];
  treeVars.t_timeCFD = new float[treeVars.t_nCh];
  chain1 -> SetBranchStatus("amp_max",1); chain1 -> SetBranchAddress("amp_max",treeVars.t_ampMax);
  chain1 -> SetBranchStatus("time_max",1); chain1 -> SetBranchAddress("time_max",treeVars.t_timeMax);
  chain1 -> SetBranchStatus("chi2_max",1); chain1 -> SetBranchAddress("chi2_max",treeVars.t_chi2Max);
  chain1 -> SetBranchStatus("time",1); chain1 -> SetBranchAddress("time",treeVars.t_timeCFD);

  chain1 -> SetBranchStatus("WF_samples",1); chain1 -> SetBranchAddress("WF_samples",&treeVars.t_nSampl);

  // treeVars.t_wfTime = new float[treeVars.t_nSampl];
  // treeVars.t_wfVal = new float[treeVars.t_nSampl];
  // treeVars.t_wfCh = new int[treeVars.t_nSampl];
  treeVars.t_wfTime = new float[5120];
  treeVars.t_wfVal = new float[5120];
  treeVars.t_wfCh = new int[5120];
  chain1 -> SetBranchStatus("WF_time",1); chain1 -> SetBranchAddress("WF_time",treeVars.t_wfTime);
  chain1 -> SetBranchStatus("WF_val",1); chain1 -> SetBranchAddress("WF_val",treeVars.t_wfVal);
  chain1 -> SetBranchStatus("WF_ch",1); chain1 -> SetBranchAddress("WF_ch",treeVars.t_wfCh);

}


//=== cut on the hodo position ===
bool SelectEvent_hodo(TreeVars treeVars)
{
  //if( !(treeVars.t_beamX[0]>-13 && treeVars.t_beamX[0]<6 && treeVars.t_beamY[0]>-11 && treeVars.t_beamY[0]<7) ) return false;
  if( !(treeVars.t_beamX[0]>-6 && treeVars.t_beamX[0]<-2 && treeVars.t_beamY[0]>-4 && treeVars.t_beamY[0]<0) ) return false;
  return true;
}
bool SelectEnPoints(TreeVars treeVars)
{
  if( treeVars.t_beamE == 200 ) return false;
  return true;
}
bool SelectEvent_scint(TreeVars treeVars)
{
  if( treeVars.t_ampMax[CH] < 40. ) return false;
  if( treeVars.t_ampMax[CHref] < 200. ) return false;
  if( (treeVars.t_timeMax[CH]/0.4-treeVars.t_timeMax[CH8]/0.4) < -240 ) return false;
  return true;
}
bool SelectEvent_kov(TreeVars treeVars)
{
  if( (treeVars.t_timeMax[CH]/0.4-treeVars.t_timeMax[CH8]/0.4) > -280 ) return false;
  //if( treeVars.t_ampMax[CH] < 2200. ) return false;
  return true;
}




int main(int argc, char** argv)
{
  setTDRStyle();
  TApplication* theApp;
  theApp = new TApplication("App", &argc, argv);


  std::string fileName = "/eos/cms/store/group/dpg_hcal/comm_hcal/deguio/CeF3_TB_Jul2017/EnScan/ntuples_v2/ntuples_v2_hadd.root";
  //std::string fileName = "/eos/cms/store/group/dpg_hcal/comm_hcal/deguio/CeF3_TB_Jul2017/EnScan/ntuples_v3/ntuples_v3_hadd.root";
  //std::string fileName = "./ntuples/H4July2017_test_8412.root";

  TFile* outFile = new TFile("TBAnalysis.root", "RECREATE");
  outFile -> cd();


  
  //=== setup chains ===
  TChain* chain1 = new TChain("info","info");
  TChain* chain2 = new TChain("hodo","hodo");
  TChain* chain3 = new TChain("digi","digi");
  TChain* chain4 = new TChain("wf","wf");
  TChain* chain5 = new TChain("h4","h4");

  chain1->Add(fileName.c_str());
  chain2->Add(fileName.c_str());
  chain3->Add(fileName.c_str());
  chain4->Add(fileName.c_str());
  chain5->Add(fileName.c_str());

  chain2 -> BuildIndex("index");
  chain1 -> AddFriend("hodo");
  chain3 -> BuildIndex("index");
  chain1 -> AddFriend("digi");
  chain4 -> BuildIndex("index");
  chain1 -> AddFriend("wf");
  chain5 -> BuildIndex("index");
  chain1 -> AddFriend("h4");
  chain1 -> BuildIndex("index");

  std::cout << " Read " << chain1->GetEntries() << " total events in tree " << chain1->GetName() << std::endl;
  std::cout << " Read " << chain2->GetEntries() << " total events in tree " << chain2->GetName() << std::endl;
  std::cout << " Read " << chain3->GetEntries() << " total events in tree " << chain3->GetName() << std::endl;
  std::cout << " Read " << chain4->GetEntries() << " total events in tree " << chain4->GetName() << std::endl;
  std::cout << " Read " << chain5->GetEntries() << " total events in tree " << chain5->GetName() << std::endl;
  
  TreeVars treeVars;
  InitTreeVars(chain1,treeVars);

  //=== book histograms ===
  TProfile2D* h_chargeMap_ch0 = new TProfile2D("h_chargeMap_ch0","h_chargeMap_ch0",80,-20,20,80,-20,20);
  TProfile2D* h_chargeMap_ch1 = new TProfile2D("h_chargeMap_ch1","h_chargeMap_ch1",80,-20,20,80,-20,20);
  TProfile2D* h_chargeMap_ch2 = new TProfile2D("h_chargeMap_ch2","h_chargeMap_ch2",80,-20,20,80,-20,20);
  TProfile2D* h_chargeMap_ch3 = new TProfile2D("h_chargeMap_ch3","h_chargeMap_ch3",80,-20,20,80,-20,20);

  TH2F* h_timeVsAmpli = new TH2F("h_timeVsAmpli","h_timeVsAmpli",5000,0.,5000.,500,-450,50);
  TH1F* h_timeDiff = new TH1F("h_timeDiff","h_timeDiff",5000,-250.,250.);

  std::map<float,TH1F*> h_ampMax;
  std::map<float,TProfile*> h_wf;

  //=== loop over events ===
  int nEntries = chain1 -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  //for(int entry = 0; entry < 50000; ++entry)
  {
    if( entry%10000 == 0 ) std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);

    //=== fill histograms ===
    h_chargeMap_ch0->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[CH0]);
    h_chargeMap_ch1->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[CH1]);
    h_chargeMap_ch2->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[CH2]);
    h_chargeMap_ch3->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[CH3]);

    //=== make event selection ===
    if (!SelectEnPoints(treeVars)) continue;
    if (!SelectEvent_hodo(treeVars)) continue;
    if (!SelectEvent_scint(treeVars)) continue;
    //if (!SelectEvent_kov(treeVars)) continue;

    //timing plots
    float time = treeVars.t_timeMax[CH]/0.4-treeVars.t_timeMax[CH8]/0.4;
    h_timeVsAmpli->Fill(treeVars.t_ampMax[CH],time);
    
    float timeDiff = treeVars.t_timeMax[CH2]-treeVars.t_timeMax[CH3];
    h_timeDiff->Fill(timeDiff);
    
    char histoName1[100];
    sprintf(histoName1, "h_wf_%f", treeVars.t_beamE);
    if(h_wf.find(treeVars.t_beamE)==h_wf.end())
      h_wf[treeVars.t_beamE] = new TProfile(histoName1,histoName1,1100,-500.,600.);
    
    //=== loop over WF samples ===
    for(unsigned int sampl=0; sampl<treeVars.t_nSampl; ++sampl)
      if(treeVars.t_wfCh[sampl] == CH)
     	h_wf[treeVars.t_beamE]->Fill(treeVars.t_wfTime[sampl]/0.4-treeVars.t_timeCFD[CH8]/0.4,treeVars.t_wfVal[sampl]);
    
    
    char histoName2[100];	
    sprintf(histoName2, "h_ampMax_%f", treeVars.t_beamE);
    if(h_ampMax.find(treeVars.t_beamE)==h_ampMax.end())
      h_ampMax[treeVars.t_beamE] = new TH1F(histoName2,histoName2,500,0.,4000.);
    
    h_ampMax[treeVars.t_beamE]->Fill(treeVars.t_ampMax[CH]);
    
    
  }//end loop
  std::cout << std::endl;


  TF1* myGauss = new TF1("myGauss","gaus");


  //=== draw plots ===
  TCanvas* c_chargeMap = new TCanvas("c_chargeMap","c_chargeMap");
  c_chargeMap->Divide(2,2);
  c_chargeMap->cd(2);
  h_chargeMap_ch0->Draw("colz");
  c_chargeMap->cd(1);
  h_chargeMap_ch1->Draw("colz");
  c_chargeMap->cd(3);
  h_chargeMap_ch3->Draw("colz");
  c_chargeMap->cd(4);
  h_chargeMap_ch2->Draw("colz");

  TCanvas* c_timeVsAmpli = new TCanvas("c_timeVsAmpli","c_timeVsAmpli");
  c_timeVsAmpli->cd();
  h_timeVsAmpli->Draw("COLZ");

  TCanvas* c_timeDiff = new TCanvas("c_timeDiff","c_timeDiff");
  c_timeDiff->cd();
  h_timeDiff->Draw();
  h_timeDiff->Fit(myGauss);

  
  TCanvas* c_ampMax = new TCanvas("c_ampMax","c_ampMax");
  c_ampMax->Divide(2,3);
  TGraphErrors* resGraph = new TGraphErrors();
  unsigned int count = 0;
  for(std::map<float,TH1F*>::const_iterator iter=h_ampMax.begin(); iter!=h_ampMax.end(); ++iter)
    {
      c_ampMax->cd(count+1);
      iter->second->Draw();
      iter->second->Fit(myGauss);
      
      //subtract noise. sigma(noise) = 1.85
      float mySigma = sqrt(myGauss->GetParameter(2)*myGauss->GetParameter(2)-1.85*1.85);
      resGraph->SetPoint(count,iter->first, mySigma/myGauss->GetParameter(1));
      ++count;
    }

  
  TF1* myRes = new TF1("myRes","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/x/x)");
  TCanvas* c_resGraph = new TCanvas("c_resGraph","c_resGraph");
  c_resGraph->cd();
  resGraph->Draw("AP");
  resGraph->Fit(myRes);

  TCanvas* c_wf = new TCanvas("c_wf","c_wf");
  c_wf->cd();
  count = 0;
  for(std::map<float,TProfile*>::const_iterator iter=h_wf.begin(); iter!=h_wf.end(); ++iter)
    {
      if(count==0)
	iter->second->Draw();
      else
	iter->second->Draw("sames");
      ++count;
    }

  theApp -> Run();  
  return 0;
}



