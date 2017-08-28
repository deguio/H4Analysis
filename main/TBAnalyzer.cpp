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


//=== focusing on ch1 for now (bottom left) ===
struct TreeVars
{  
  float t_beamE;
  float t_HV;
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
  chain1 -> SetBranchStatus("HV",1); chain1 -> SetBranchAddress("HV",&treeVars.t_HV);
  chain1 -> SetBranchStatus("n_channels",1); chain1 -> SetBranchAddress("n_channels",&treeVars.t_nCh);
  
  treeVars.t_beamX = new float[2];
  treeVars.t_beamY = new float[2];
  chain1 -> SetBranchStatus("X",1); chain1 -> SetBranchAddress("X",treeVars.t_beamX);
  chain1 -> SetBranchStatus("Y",1); chain1 -> SetBranchAddress("Y",treeVars.t_beamY);
  
  treeVars.t_chargeSig = new float[5];
  treeVars.t_chargeTot = new float[5];
  chain1 -> SetBranchStatus("charge_sig",1); chain1 -> SetBranchAddress("charge_sig",treeVars.t_chargeSig);
  chain1 -> SetBranchStatus("charge_tot",1); chain1 -> SetBranchAddress("charge_tot",treeVars.t_chargeTot);

  treeVars.t_ampMax = new float[5];
  treeVars.t_timeMax = new float[5];
  treeVars.t_chi2Max = new float[5];
  treeVars.t_timeCFD = new float[5];
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
bool SelectEvent_hodo(TreeVars treeVars, std::vector<float> positionCuts)
{
  if( !(treeVars.t_beamX[0]>positionCuts[0] && treeVars.t_beamX[0]<positionCuts[1] && treeVars.t_beamY[0]>positionCuts[2] && treeVars.t_beamY[0]<positionCuts[3]) ) return false;
  return true;
}
bool SelectEnPoints(TreeVars treeVars, std::vector<float> enpoints)
{
  for(int en=0; en<enpoints.size(); ++en)
    if(treeVars.t_beamE == enpoints[en])
      return true;

  return false;
}
bool SelectHVPoints(TreeVars treeVars, std::vector<float> hvpoints)
{
  for(int hv=0; hv<hvpoints.size(); ++hv)
    if(treeVars.t_HV == hvpoints[hv])
      return true;

  return false;
}
bool SelectEvent_ref(TreeVars treeVars, int mainch, float mainthr, std::vector<float> refch, std::vector<float> refthr)
{
  if( treeVars.t_ampMax[mainch] < mainthr ) return false;

  for(int ch=0; ch<refch.size(); ++ch)
    if( treeVars.t_ampMax[int(refch[ch])] < refthr[ch] ) return false;

  return true;
}




int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  setTDRStyle();

  if( argc < 2 )
    {
      std::cerr << ">>>>> TBAnalyzer.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
      return 1;
    }

  //=== Take Params from Config ===
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string fileName = opts.GetOpt<std::string>("Input.inputFile");
  std::string outName = opts.GetOpt<std::string>("Input.outputFile");
  std::vector<float> positionCuts = opts.GetOpt<std::vector<float>>("Cuts.position");
  std::vector<float> enpoints     = opts.GetOpt<std::vector<float>>("Cuts.enpoints");
  std::vector<float> hvpoints     = opts.GetOpt<std::vector<float>>("Cuts.hvpoints");

  int mainch    = opts.GetOpt<int>("Input.mainch");
  float mainthr = opts.GetOpt<int>("Input.mainthr");
  std::vector<float> refch   = opts.GetOpt<std::vector<float>>("Input.refch");
  std::vector<float> refthr  = opts.GetOpt<std::vector<float>>("Input.refthr");
  int trgch  = opts.GetOpt<int>("Input.trgch");

  //=== outFile ===
  TFile* outFile = new TFile(outName.c_str(), "RECREATE");
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

  TApplication* theApp;
  theApp = new TApplication("App", &argc, argv);


  //=== loop over events ===
  int nEntries = chain1 -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  //for(int entry = 0; entry < 50000; ++entry)
  {
    if( entry%10000 == 0 ) std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);

    //=== fill histograms ===
    h_chargeMap_ch0->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[0]);
    h_chargeMap_ch1->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[1]);
    h_chargeMap_ch2->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[2]);
    h_chargeMap_ch3->Fill(treeVars.t_beamX[0],treeVars.t_beamY[0],treeVars.t_chargeSig[3]);


    //=== make event selection ===
    if (!SelectEnPoints(treeVars,enpoints)) continue;
    if (!SelectHVPoints(treeVars,hvpoints)) continue;
    if (!SelectEvent_hodo(treeVars,positionCuts)) continue;
    if (!SelectEvent_ref(treeVars,mainch,mainthr,refch,refthr)) continue;


    //timing plots
    float time = treeVars.t_timeMax[mainch]/0.4-treeVars.t_timeCFD[trgch]/0.4;
    h_timeVsAmpli->Fill(treeVars.t_ampMax[mainch],time);
    
    float timeDiff = treeVars.t_timeMax[2]-treeVars.t_timeMax[1];
    h_timeDiff->Fill(timeDiff);
    
    char histoName1[100];
    sprintf(histoName1, "h_wf_%f", treeVars.t_beamE);
    if(h_wf.find(treeVars.t_beamE)==h_wf.end())
      h_wf[treeVars.t_beamE] = new TProfile(histoName1,histoName1,1300,-500.,800.);
    
    //=== loop over WF samples ===
    for(unsigned int sampl=0; sampl<treeVars.t_nSampl; ++sampl)
      {
	//std::cout << sampl << " " << treeVars.t_wfCh[sampl] << " " << treeVars.t_wfTime[sampl]/0.4-treeVars.t_timeCFD[trgch]/0.4 << " " << treeVars.t_timeCFD[trgch]  << " " << treeVars.t_wfVal[sampl] << std::endl;


	if(treeVars.t_wfCh[sampl] == mainch)
	  h_wf[treeVars.t_beamE]->Fill(treeVars.t_wfTime[sampl]/0.4-treeVars.t_timeCFD[trgch]/0.4,treeVars.t_wfVal[sampl]);
      }

    char histoName2[100];	
    sprintf(histoName2, "h_ampMax_%f", treeVars.t_beamE);
    if(h_ampMax.find(treeVars.t_beamE)==h_ampMax.end())
      h_ampMax[treeVars.t_beamE] = new TH1F(histoName2,histoName2,80,0.,4000.);
    
    h_ampMax[treeVars.t_beamE]->Fill(treeVars.t_ampMax[mainch]);

  }//end loop
  std::cout << std::endl;


  TF1* myGauss = new TF1("myGauss","gaus");


  //=== draw plots ===
  TCanvas* c_chargeMap = new TCanvas("c_chargeMap","c_chargeMap");
  c_chargeMap->Divide(2,2);
  c_chargeMap->cd(1);
  h_chargeMap_ch0->Draw("colz");
  c_chargeMap->cd(2);
  h_chargeMap_ch1->Draw("colz");
  c_chargeMap->cd(4);
  h_chargeMap_ch2->Draw("colz");
  c_chargeMap->cd(3);
  h_chargeMap_ch3->Draw("colz");

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
      //float mySigma = sqrt(myGauss->GetParameter(2)*myGauss->GetParameter(2)-1.85*1.85);
      float mySigma = myGauss->GetParameter(2);
      float mySigmaErr = myGauss->GetParError(2);
      float myMean = myGauss->GetParameter(1);
      float myMeanErr = myGauss->GetParError(1);
      resGraph->SetPoint(count,iter->first, mySigma/myMean);
      float err2 = pow(mySigmaErr/myMean,2) + pow(myMeanErr*mySigma/myMean/myMean,2);
      resGraph->SetPointError(count, 0., sqrt(err2));

      iter->second->Write();

      ++count;
    }

  
  TF1* myRes = new TF1("myRes","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/x/x)");
  TCanvas* c_resGraph = new TCanvas("c_resGraph","c_resGraph");
  c_resGraph->cd();
  resGraph->Draw("AP");
  resGraph->Fit(myRes);

  TCanvas* c_wf = new TCanvas("c_wf","c_wf");
  c_wf->Divide(2,3);
  count = 0;
  for(std::map<float,TProfile*>::const_iterator iter=h_wf.begin(); iter!=h_wf.end(); ++iter)
    {
      c_wf->cd(count+1);
      iter->second->Draw();

      ++count;
    }

  theApp -> Run();  
  outFile -> Close();
  

  return 0;
}



