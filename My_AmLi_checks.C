// Standard C++ headers                                                                                                                                                             
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT Headers                                                                                                                                                                     
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TTimeStamp.h"



void My_AmLi_checks(char *infiles) {

  TH1F *h_spes = new TH1F("h_spes"," ; pulse area (phd) ; ", 200, 0., 20.);
  TH2F *h_spes_width = new TH2F("h_spe_width", " ; pulse area (phd) ; pulse width (ns) ", 200,0.,20.,200,0.,1200.);
  TH1F *h_ses = new TH1F("h_ses"," ; pulse area (phd) ; ", 200, 0., 500.);
  TH2F *h_logS1S2 = new TH2F("h_logS1S2", " ; log(S1 area) (phd) ; log(S2 area) (phd) ",500,0.,1.E5,500,0.,1.E7);
  TH1F *h_skinarea = new TH1F("h_skinarea", "; skin area (phd) ; ", 200, 0., 2000.);
  TH1F *h_odarea = new TH1F("h_odarea", "; od area (phd) ; ", 200, 0., 1000.);
  TH1F *h_skinsumarea = new TH1F("h_skinsumarea", "; skin area (phd) ; ", 200, 0., 1000.);
  TH1F *h_odsumarea = new TH1F("h_odsumarea", "; od area (phd) ; ", 201, -5., 2005.);
  
  TH2F* hrs2s1s1 = new TH2F("hrs2s1s1", "Single Scatter: S2/S1 Vs S1; S1; log10(S2/S1)", 6000, 0., 6000., 70, 0., 7.);
  TH2F* hS2dTss = new TH2F("hS2dTss", " ;Drift Time/#{mu}s;log10(S2 [phd])", 1000, 0., 1000., 100, 0., 10.);
  TH2F* hrs2s1s1_m = new TH2F("hrs2s1s1_m", "Multiple Scatter: S2/S1 Vs S1; S1; log10(S2/S1)", 51000, -1000., 50000., 700, 0., 7.);
  TH2F* hS2dTms = new TH2F("hS2dTms", " ;Drift Time/#{mu}s;log10(S2 [phd])", 1000, 0., 1000., 1000, 0., 10.);


  TChain *events = new TChain("Events");
  TChain *scatters = new TChain("Scatters");

  TString txtFileList(infiles);
  if (txtFileList.Contains(".txt")) {
    cout << "Loading file names from "<<txtFileList << " Int_to "<< events->GetName()<<endl;
    ifstream fileList(txtFileList);
    string file;
    if (fileList.is_open()) {
      while ( getline(fileList, file) ) {
	events->AddFile(file.c_str());
	scatters->AddFile(file.c_str());
      }
      fileList.close();
      }else{
      cout<<"The file "<< txtFileList <<" doesn't exist. Exiting !!"<<endl;
      exit(-1);
    }
  }
  else if (txtFileList.Contains(".root")) {
    events->Add(infiles);
    scatters->Add(infiles);
  }


    unsigned long runID {};
    unsigned long eventID {};

  Int_t nPulsesOD;
  std::vector<Float_t> odPulseArea_phd;
  std::vector<Int_t> odAft5;
  std::vector<Int_t> odCoincidence;
  std::vector<Int_t> odPulseStart;

  Int_t nPulsesTPC;
  std::vector<Float_t> tpcPulseArea_phd;
  std::vector<Int_t> tpcPulseStart, tpcPulseEnd;
  std::vector<Float_t> tpcPulseTBA;
  std::vector<Int_t>   tpcAft5;
  std::vector<Int_t>   tpcAft50;
  std::vector<Int_t>   tpcAft95;
  std::vector<Int_t> tpcCoincidence;
  std::vector<Float_t> speProb;
  std::vector<Float_t> s1Prob;
  std::vector<Float_t> s2Prob;
  std::vector<Float_t> seProb;
  std::vector<Float_t> s2Xpos_cm;
  std::vector<Float_t> s2Ypos_cm;

  Int_t nPulsesSkin;
  std::vector<Float_t> skinPulseArea_phd;
  std::vector<Int_t> skinPulseStart;
  std::vector<Int_t> skinAft5;
  std::vector<Int_t> skinCoincidence;


  events->SetMakeClass(1);
  events->SetBranchStatus("*",0);
  // header                                                                                                                                                                                 
  events->SetBranchStatus("eventHeader.eventID",                    1);
  events->SetBranchAddress("eventHeader.eventID",                   &eventID);
  events->SetBranchStatus("eventHeader.runID",                      1);
  events->SetBranchAddress("eventHeader.runID",                     &runID);
  //events->SetBranchStatus("eventHeader.rawFileName",                1);
  //events->SetBranchAddress("eventHeader.rawFileName",               &rawFileName);
  //events->SetBranchStatus("eventHeader.triggerTimeStamp_ns",        1);
  //events->SetBranchAddress("eventHeader.triggerTimeStamp_ns",       &triggerTimestamp_ns);
  //events->SetBranchStatus("eventHeader.triggerTimeStamp_s",         1);
  //events->SetBranchAddress("eventHeader.triggerTimeStamp_s",        &triggerTimestamp_s);
    
  //OD                                                                                                                                                                                     
  events->SetBranchStatus("pulsesODHG.nPulses",                     1);
  events->SetBranchAddress("pulsesODHG.nPulses",                    &nPulsesOD);
  events->SetBranchStatus("pulsesODHG.pulseArea_phd",               1);
  events->SetBranchAddress("pulsesODHG.pulseArea_phd",              &odPulseArea_phd);
  events->SetBranchStatus("pulsesODHG.pulseStartTime_ns",           1);
  events->SetBranchAddress("pulsesODHG.pulseStartTime_ns",          &odPulseStart);
  events->SetBranchStatus("pulsesODHG.areaFractionTime5_ns",        1);
  events->SetBranchAddress("pulsesODHG.areaFractionTime5_ns",       &odAft5);
  events->SetBranchStatus("pulsesODHG.coincidence",                 1);
  events->SetBranchAddress("pulsesODHG.coincidence",                &odCoincidence);
  //TPC                                                                                                                                                                                     
  events->SetBranchStatus("pulsesTPCHG.nPulses",                    1);
  events->SetBranchAddress("pulsesTPCHG.nPulses",                   &nPulsesTPC);
  events->SetBranchStatus("pulsesTPCHG.pulseArea_phd",              1);
  events->SetBranchAddress("pulsesTPCHG.pulseArea_phd",             &tpcPulseArea_phd);
  events->SetBranchStatus("pulsesTPCHG.areaFractionTime5_ns",       1);
  events->SetBranchAddress("pulsesTPCHG.areaFractionTime5_ns",      &tpcAft5);
  events->SetBranchStatus("pulsesTPCHG.areaFractionTime50_ns",      1);
  events->SetBranchAddress("pulsesTPCHG.areaFractionTime50_ns",     &tpcAft50);
  events->SetBranchStatus("pulsesTPCHG.areaFractionTime95_ns",      1);
  events->SetBranchAddress("pulsesTPCHG.areaFractionTime95_ns",     &tpcAft95);
  events->SetBranchStatus("pulsesTPCHG.pulseStartTime_ns",          1);
  events->SetBranchAddress("pulsesTPCHG.pulseStartTime_ns",         &tpcPulseStart);
  events->SetBranchStatus("pulsesTPCHG.pulseEndTime_ns",            1);
  events->SetBranchAddress("pulsesTPCHG.pulseEndTime_ns",           &tpcPulseEnd);
  events->SetBranchStatus("pulsesTPCHG.topBottomAsymmetry",         1);
  events->SetBranchAddress("pulsesTPCHG.topBottomAsymmetry",        &tpcPulseTBA);
  events->SetBranchStatus("pulsesTPCHG.singlePEprobability",        1);
  events->SetBranchAddress("pulsesTPCHG.singlePEprobability",       &speProb);
  events->SetBranchStatus("pulsesTPCHG.s1Probability",              1);
  events->SetBranchAddress("pulsesTPCHG.s1Probability",             &s1Prob);
  events->SetBranchStatus("pulsesTPCHG.s2Probability",              1);
  events->SetBranchAddress("pulsesTPCHG.s2Probability",             &s2Prob);
  events->SetBranchStatus("pulsesTPCHG.singleElectronProbability",  1);
  events->SetBranchAddress("pulsesTPCHG.singleElectronProbability", &seProb);
  events->SetBranchStatus("pulsesTPCHG.s2Xposition_cm",             1);
  events->SetBranchAddress("pulsesTPCHG.s2Xposition_cm",            &s2Xpos_cm);
  events->SetBranchStatus("pulsesTPCHG.s2Yposition_cm",             1);
  events->SetBranchAddress("pulsesTPCHG.s2Yposition_cm",            &s2Ypos_cm);
  events->SetBranchStatus("pulsesTPCHG.coincidence",                1);
  events->SetBranchAddress("pulsesTPCHG.coincidence",               &tpcCoincidence);

  // Skin                                                                                                                                                                                   
  events->SetBranchStatus("pulsesSkin.nPulses",                     1);
  events->SetBranchAddress("pulsesSkin.nPulses",                    &nPulsesSkin);
  events->SetBranchStatus("pulsesSkin.pulseArea_phd",               1);
  events->SetBranchAddress("pulsesSkin.pulseArea_phd",              &skinPulseArea_phd);
  events->SetBranchStatus("pulsesSkin.pulseStartTime_ns",           1);
  events->SetBranchAddress("pulsesSkin.pulseStartTime_ns",          &skinPulseStart);
  events->SetBranchStatus("pulsesSkin.areaFractionTime5_ns",        1);
  events->SetBranchAddress("pulsesSkin.areaFractionTime5_ns",       &skinAft5);
  events->SetBranchStatus("pulsesSkin.coincidence",                 1);
  events->SetBranchAddress("pulsesSkin.coincidence",                &skinCoincidence);

  //Scatters
  Int_t nSingleScatters;
  Int_t nMultipleScatters;
  Int_t nKr83mScatters;
  Int_t nPileUpScatters;
  Int_t nOtherScatters;
  scatters->SetMakeClass(1);
  scatters->SetBranchAddress("ss.nSingleScatters", &nSingleScatters);
  scatters->SetBranchAddress("ms.nMultipleScatters", &nMultipleScatters);
  scatters->SetBranchAddress("kr83m.nKr83mScatters", &nKr83mScatters);
  scatters->SetBranchAddress("pileUp.nPileUpScatters", &nPileUpScatters);
  scatters->SetBranchAddress("other.nOtherScatters", &nOtherScatters);
  Float_t ss_s1Area_phd;
  Float_t ms_s1Area_phd;
  Float_t ss_s2Area_phd;
  vector<Float_t> ms_s2Area_phd;
  Float_t ss_x_cm;
  Float_t ss_y_cm;
  Float_t ss_drift;
  vector<Float_t> ms_x_cm;
  vector<Float_t> ms_y_cm;
  Float_t ms_drift;
  Int_t nS2s;
  Int_t ss_s1ID;
  Int_t ms_s1ID;
  Int_t ss_s2ID;
  vector<Int_t> ms_s2ID;
  scatters->SetBranchAddress("ss.s1PulseID", &ss_s1ID);
  scatters->SetBranchAddress("ms.s1PulseID", &ms_s1ID);
  scatters->SetBranchAddress("ss.s2PulseID", &ss_s2ID);
  scatters->SetBranchAddress("ms.s2PulseIDs", &ms_s2ID);
  scatters->SetBranchAddress("ss.s1Area_phd", &ss_s1Area_phd);
  scatters->SetBranchAddress("ms.s1Area_phd", &ms_s1Area_phd);
  scatters->SetBranchAddress("ss.s2Area_phd", &ss_s2Area_phd);
  scatters->SetBranchAddress("ms.s2Area_phd", &ms_s2Area_phd);
  scatters->SetBranchAddress("ss.x_cm", &ss_x_cm);
  scatters->SetBranchAddress("ss.y_cm", &ss_y_cm);
  scatters->SetBranchAddress("ss.driftTime_ns", &ss_drift);
  scatters->SetBranchAddress("ms.x_cm", &ms_x_cm);
  scatters->SetBranchAddress("ms.y_cm", &ms_y_cm);
  scatters->SetBranchAddress("ms.weightedDriftTime_ns", &ms_drift);
  scatters->SetBranchAddress("ms.nS2s", &nS2s);

  // variable to define ratio
  Float_t rs2cs1c, rs2s1, rs2s1ss,rs2cs1css;
  Float_t rs2cs1c_m, rs2s1_m, rs2s1ss_m,rs2cs1css_m;
  //Float_t s1_s = ss_s1Area_phd;
  //Float_t s2_s =log10(ss_s2Area_phd);
  //Float_t dt_s = ss_drift/1000.;
  //Float_t s1_m = ms_s1Area_phd;
  //Float_t s2_m =log10((ms_s2Area_phd);
  //Float_t dt_m = ms_drift/1000.;

  //if (s1_s!=0) rs2s1 = (s2_s / s1_s); else rs2s1=0;
  //rs2s1ss = log10(rs2s1);


  Int_t nEvents = events->GetEntries();
  Int_t ss = 0;
  Int_t ms = 0;
  Int_t kr83m = 0;
  Int_t pileup = 0;
  Int_t  other = 0;
  Float_t skinsum;
  Float_t skin_sum;
  Float_t odsum;
  Float_t od_sum;
  Float_t odP;
  Int_t odC;
  Float_t od_PST;
  Float_t TPC_PST;
  Int_t skin = 0;
  Int_t od = 0;
  Int_t tpc = 0;
  Int_t skinandod = 0;
  Int_t odandtpc = 0;
  Int_t skinandtpc = 0;
  Int_t allthree = 0;
  Int_t skinonly = 0;
  Int_t odonly = 0;
  Int_t tpconly = 0;
  for (Int_t evt=0; evt < nEvents; evt++) {
   
    skinsum = 0;
    odsum = 0;
    skin_sum = 0;
    od_sum = 0;
    rs2s1 = 0;
    rs2s1_m = 0;
    
    //if (evt%10==0) cout << "Processing event number " << evt << " of " << nEvents << endl;
    events->GetEntry(evt);
    scatters->GetEntry(evt);

    if (nSingleScatters>0 && ss_s1Area_phd < 50)
     {
            Float_t s1_s = ss_s1Area_phd;
            Float_t s2_s =log10(ss_s2Area_phd);
            Float_t dt_s = ss_drift/1000.;
            ++ss;
            od_sum += odPulseArea_phd[nPulsesOD];
            skin_sum += skinPulseArea_phd[nPulsesSkin];
            if (ss_s1Area_phd!=0) rs2s1 = (ss_s2Area_phd / ss_s1Area_phd); else rs2s1=0;
              //rs2s1ss = log10(rs2s1);

            hrs2s1s1->Fill(s1_s,log10(rs2s1));
            hS2dTss->Fill(dt_s, s2_s);

            //cout << "S2ID = " << s2pIDSS  <<endl;
            //cout << "SS "<< ss << ": --> runID = " << runID << ",   " << "eventID = "<<eventID  <<  endl;
            cout << "SS S1 Area [phd] = " << s1_s << " S2 Area [phd] = " << s2_s<< "   Ratio = " << log10(rs2s1) << " Drift Time = "<< dt_s<< endl;
            //cout << "SS: SkinPulseArea [phd] = " <<skin_sum << "  ODPulseArea = "<<  od_sum << endl;
      }

    //if (nMultipleScatters>0)// && ms_s1Area_phd < 50)
    for(Int_t i = 0; i < nMultipleScatters; i++)
     {  
      if(ms_s1Area_phd < 50){
        ++ms;
        Float_t s1_m = ms_s1Area_phd;
        Float_t s2_m =log10(ms_s2Area_phd[nMultipleScatters]);
        Float_t dt_m = ms_drift/1000.;
        //od_sum += odPulseArea_phd[nPulsesOD];
        skin_sum += skinPulseArea_phd[nPulsesSkin];
        if (ms_s1Area_phd!=0) rs2s1_m = (ms_s2Area_phd[nMultipleScatters] / ms_s1Area_phd); else rs2s1_m=0;
           //rs2s1ss_m = log10(rs2s1_m);

        hrs2s1s1_m->Fill(s1_m,log10(rs2s1_m));
        hS2dTms->Fill(dt_m, s2_m);        

            //cout << "S2ID = " << s2pIDSS  <<endl;
            //cout << "MS "<< ms << ": --> runID = " << runID << ",   " << "eventID = "<<eventID  <<  endl;
            cout << "MS S1 Area [phd] = " << s1_m << " S2 Area [phd] = " << s2_m << "   Ratio = " << log10(rs2s1_m) << " Drift Time = "<< dt_m << endl;
            //cout << "MS: SkinPulseArea [phd] = " <<skin_sum << "  ODPulseArea = "<<  od_sum << endl;
        }
      }



    //if (nSingleScatters>0) ss++;
    //if (nMultipleScatters>0) ms++;
    if (nKr83mScatters>0) kr83m++;
    if (nPileUpScatters>0) pileup++;
    if (nOtherScatters>0) other++;

    if (nPulsesSkin>0) skin++;
    if (nPulsesOD>0) od++;
    if (nPulsesTPC>0) tpc++;
    if (nPulsesSkin > 0 && nPulsesOD > 0 && nPulsesTPC == 0) skinandod++;
    if (nPulsesSkin > 0 && nPulsesTPC > 0 && nPulsesOD == 0) skinandtpc++;
    if (nPulsesOD > 0 && nPulsesTPC > 0 && nPulsesSkin == 0) odandtpc++;
    if (nPulsesOD > 0 && nPulsesTPC > 0 && nPulsesSkin > 0) allthree++;
    if (nPulsesTPC>0 && nPulsesSkin == 0 && nPulsesOD == 0) tpconly++;
    if (nPulsesTPC== 0 && nPulsesSkin ==0 && nPulsesOD > 0) odonly++;
    if (nPulsesTPC == 0 && nPulsesSkin > 0 && nPulsesOD ==0) skinonly++;


//Skin
    for (Int_t p=0; p<nPulsesSkin; p++) {
      h_skinarea->Fill(skinPulseArea_phd[p]);  
      if (skinCoincidence[p] > 1) {
	skinsum += skinPulseArea_phd[p];
      }
    }

//OD  
    for(Int_t p=0; p<nPulsesOD; p++) {
      h_odarea->Fill(odPulseArea_phd[p]);
      if (odCoincidence[p] > 4)// To remove dark counts
       {
  odC = odCoincidence[p];      
  odP = odPulseArea_phd[p];     
	odsum += odPulseArea_phd[p];
  od_PST = odPulseStart[p];

  for (Int_t i=0; i<nPulsesTPC; i++) { 
  Float_t time_constant = tpcPulseStart[i] - odPulseStart[p]; 
  //cout << time_constant << endl;
      if (nSingleScatters>0 && ss_s1Area_phd < 50 && (odCoincidence[p] > 4) && odPulseArea_phd[p] >100 && (tpcPulseStart[i]/1000. < 30. + odPulseStart[p]/1000.)) // && od_sum >1000000)
           {cout << "SS: ODPulseArea = "<<  odPulseArea_phd[p] << "Time Constant = " << (tpcPulseStart[i] - odPulseStart[p])/1000. <<endl;
            //cout << time_constant/1000. << endl;
           }
      if (nMultipleScatters>0 && ms_s1Area_phd < 50 && (odCoincidence[p] > 4) && odPulseArea_phd[p]>100 && (tpcPulseStart[i]/1000. < 30. + odPulseStart[p]/1000.)) // && od_sum >30000)
           {cout << "MS: ODPulseArea = "<<  odPulseArea_phd[p]<< "Time Constant = " << (tpcPulseStart[i] - odPulseStart[p])/1000. << endl;
            //cout << time_constant/1000. << endl;
           
         }
          }
      }
    }  
/*
      if (nSingleScatters>0 && ss_s1Area_phd < 50) // && od_sum >1000000)
           {cout << "SS: ODPulseArea = "<<  odsum << "Time Constant = " << time_constant << endl;}
      if (nMultipleScatters>0 && ms_s1Area_phd < 50) // && od_sum >30000)
           {cout << "MS: ODPulseArea = "<<  odsum << "Time Constant = " << time_constant << endl;}*/        
    
	
    if (skinsum>0)h_skinsumarea->Fill(skinsum);
    if (odsum>0)  h_odsumarea->Fill(odsum);

//TPC   
        for (Int_t p=0; p<nPulsesTPC; p++) {
      //if (speProb[p] == 1) { 
      if (tpcCoincidence[p] ==1) {
	h_spes->Fill(tpcPulseArea_phd[p]);
	h_spes_width->Fill(tpcPulseArea_phd[p], tpcAft95[p]-tpcAft5[p]);
      }
      if (seProb[p]==1) {
        h_ses->Fill(tpcPulseArea_phd[p]);
      }
  TPC_PST = tpcPulseStart[p];
    }

    if (ss_s1Area_phd>0) {
      h_logS1S2->Fill(log10(ss_s1Area_phd),log10(ss_s2Area_phd));
    }

 }


  cout << "Done!" << endl;
  cout << " Of " << nEvents << " events: " << endl;
 
  cout << ss << " are single scatters" << endl;
  cout << ms << " are multiple scatters" << endl;
  cout << kr83m << " are Kr83m scatters" << endl;
  cout << pileup << " are Pile Up scatters" << endl;
  cout << other << " are other scatters" << endl;

  cout << tpc << " have TPC pulses, and there are " << tpconly << " only in the TPC" <<  endl;
  cout << od << " have OD pulses, and there are " << odonly << " only in the OD" << endl;
  cout << skin << " have Skin pulses, and there are " << skinonly << " only in the Skin" << endl;
  cout << skinandtpc << " are only in the skin and TPC" << endl;
  cout << odandtpc <<" are only in the OD and TPC"<< endl;
  cout << skinandod << " are only in the vetoes" << endl;
  cout << allthree << " are in all three." << endl;


  TFile *f = new TFile("My_mdc3_checks_try.root","RECREATE");

  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);

  h_spes->Write();
  h_spes_width->Write();
  h_skinsumarea->Write();
  h_odsumarea->Write();
  hrs2s1s1->Write();
  hS2dTss->Write();
  hrs2s1s1_m->Write();
  hS2dTms->Write();  
  f->Close();

/*
  TCanvas *c = new TCanvas();
  c->cd();
  h_spes->Draw();
  TCanvas *c2 = new TCanvas();
  c2->cd();
  h_spes_width->Draw();
  TCanvas *c3 = new TCanvas();
  //c3->Divide(1,2);
  //c3->cd(1);
  //h_skinarea->Draw();
  c3->cd();
  //gPad->SetLogy();
  h_skinsumarea->Draw();
  TCanvas *c4 = new TCanvas();
  //c4->Divide(1,2);
  //c4->cd(1);
  //h_odarea->Draw();
  //c4->cd(2);
  
  c4->cd();
  //gPad->SetLogy();
  h_odsumarea->Draw();*/


}

