//To find the rate from a (very very large) min bias sample,
//Runs off of the track object ntuplizer

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;
TH1D* GetCumulative(TH1D* plot);
void printThresh(TH1D* histo, float thresh);
string noZeroes(float value);

void trkMET_rate() {
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add("CheckingMET_MB_CMSSW1117.root");

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }


  int numEventsToRun = 300000; //anything greater than 0 will only run that number of events
  bool rerun_trueTkMET = true;
  bool rerun_recoTkMET = true;
  float TP_minPt = 2.0;
  float TP_maxEta = 2.4;

  //Delta z cuts
  const int numEtaReg = 7;
  float z0Thresholds[numEtaReg] = {0.37, 0.5, 0.6, 0.75, 1.0, 1.6};
  float etaRegions[numEtaReg] = {0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4};


  // define leafs & branches
  float trueMET = 0;
  float trueTkMET = 0;
  float trkMET = 0;
  float trkMETEmu = 0;
  vector<float>* pv_L1reco;
  vector<float>* pv_L1reco_emu;
  //gen particles
  vector<float>* gen_pt;
  vector<float>* gen_phi;
  vector<int>* gen_pdgid;

  // tracking particles
  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_dxy;
  vector<float>* tp_z0;
  vector<float>* tp_d0;
  vector<int>* tp_pdgid;
  vector<int>* tp_nmatch;
  vector<int>* tp_nstub;
  vector<int>* tp_eventid;
  vector<int>* tp_charge;

  // Matched tracks
  vector<float>* matchtrk_pt;
  vector<float>* matchtrk_eta;
  vector<float>* matchtrk_phi;
  vector<float>* matchtrk_d0;
  vector<float>* matchtrk_z0;
  vector<float>* matchtrk_chi2;
  vector<float>* matchtrk_chi2dof;
  vector<float>* matchtrk_chi2rz;
  vector<float>* matchtrk_chi2rphi;
  vector<float>* matchtrk_MVA1;
  vector<float>* matchtrk_bendchi2;
  vector<int>* matchtrk_nstub;
  vector<int>* matchtrk_seed;

  // all tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_z0;
  vector<float>* trk_chi2;
  vector<float>* trk_chi2dof;
  vector<float>* trk_chi2rz;
  vector<float>* trk_chi2rphi;
  vector<float>* trk_MVA1;
  vector<float>* trk_bendchi2;
  vector<int>* trk_nstub;
  vector<int>* trk_seed;
  vector<int>* trk_fake;

  TBranch* b_gen_pt;
  TBranch* b_gen_phi;
  TBranch* b_gen_pdgid;

  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_dxy;
  TBranch* b_tp_z0;
  TBranch* b_tp_d0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_nmatch;
  TBranch* b_tp_nstub;
  TBranch* b_tp_eventid;
  TBranch* b_tp_charge;

  TBranch* b_matchtrk_pt;
  TBranch* b_matchtrk_eta;
  TBranch* b_matchtrk_phi;
  TBranch* b_matchtrk_d0;
  TBranch* b_matchtrk_z0;
  TBranch* b_matchtrk_chi2;
  TBranch* b_matchtrk_chi2dof;
  TBranch* b_matchtrk_chi2rz;
  TBranch* b_matchtrk_chi2rphi;
  TBranch* b_matchtrk_MVA1;
  TBranch* b_matchtrk_bendchi2;
  TBranch* b_matchtrk_nstub;
  TBranch* b_matchtrk_seed;

  TBranch* b_trk_pt;
  TBranch* b_trk_eta;
  TBranch* b_trk_phi;
  TBranch* b_trk_z0;
  TBranch* b_trk_chi2;
  TBranch* b_trk_chi2dof;
  TBranch* b_trk_chi2rz;
  TBranch* b_trk_chi2rphi;
  TBranch* b_trk_MVA1;
  TBranch* b_trk_bendchi2;
  TBranch* b_trk_nstub;
  TBranch* b_trk_seed;
  TBranch* b_trk_fake;

  TBranch* b_pv_L1reco;
  TBranch* b_pv_L1reco_emu;
  TBranch* b_trueMET;
  TBranch* b_trueTkMET;
  TBranch* b_trkMET;
  TBranch* b_trkMETEmu;

  gen_pt = 0;
  gen_phi = 0;
  gen_pdgid = 0;

  tp_pt = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_z0 = 0;
  tp_d0 = 0;
  tp_pdgid = 0;
  tp_nmatch = 0;
  tp_nstub = 0;
  tp_eventid = 0;
  tp_charge = 0;

  matchtrk_pt = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_d0 = 0;
  matchtrk_z0 = 0;
  matchtrk_chi2 = 0;
  matchtrk_chi2dof = 0;
  matchtrk_chi2rz = 0;
  matchtrk_chi2rphi = 0;
  matchtrk_MVA1 = 0;
  matchtrk_bendchi2 = 0;
  matchtrk_nstub = 0;
  matchtrk_seed = 0;

  trk_pt = 0;
  trk_eta = 0;
  trk_phi = 0;
  trk_z0 = 0;
  trk_chi2 = 0;
  trk_chi2dof = 0;
  trk_chi2rz = 0;
  trk_chi2rphi = 0;
  trk_MVA1 = 0;
  trk_bendchi2 = 0;
  trk_nstub = 0;
  trk_seed = 0;
  trk_fake = 0;

  pv_L1reco = 0;
  pv_L1reco_emu = 0;

  tree->SetBranchAddress("pv_L1reco", &pv_L1reco, &b_pv_L1reco);
  tree->SetBranchAddress("pv_L1reco_emu", &pv_L1reco_emu, &b_pv_L1reco_emu);
  tree->SetBranchAddress("trueMET", &trueMET, &b_trueMET);
  tree->SetBranchAddress("trueTkMET", &trueTkMET, &b_trueTkMET);
  tree->SetBranchAddress("trkMET", &trkMET, &b_trkMET);
  tree->SetBranchAddress("trkMETEmu", &trkMETEmu, &b_trkMETEmu);
  tree->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
  tree->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
  tree->SetBranchAddress("gen_pdgid", &gen_pdgid, &b_gen_pdgid);

  tree->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
  tree->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
  tree->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
  tree->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
  tree->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
  tree->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
  tree->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
  tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
  tree->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
  tree->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
  tree->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);

  tree->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
  tree->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
  tree->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
  tree->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
  tree->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
  tree->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
  tree->SetBranchAddress("matchtrk_chi2dof", &matchtrk_chi2dof, &b_matchtrk_chi2dof);
  tree->SetBranchAddress("matchtrk_chi2rz", &matchtrk_chi2rz, &b_matchtrk_chi2rz);
  tree->SetBranchAddress("matchtrk_chi2rphi", &matchtrk_chi2rphi, &b_matchtrk_chi2rphi);
  tree->SetBranchAddress("matchtrk_MVA1", &matchtrk_MVA1, &b_matchtrk_MVA1);
  tree->SetBranchAddress("matchtrk_bendchi2", &matchtrk_bendchi2, &b_matchtrk_bendchi2);
  tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);

  tree->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
  tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
  tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
  tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_chi2dof", &trk_chi2dof, &b_trk_chi2dof);
  tree->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
  tree->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
  tree->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);
  tree->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
  tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
  tree->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);

  float numBins = 100.0;
  TH1D* trueTkMET_thresh = new TH1D("trueTkMET_thresh",";TrackMET threshold [GeV];Events",numBins,0,numBins); //using tracking particles
  TH1D* recoTkMET_thresh = new TH1D("recoTkMET_thresh", ";TrackMET threshold [GeV];Events", numBins, 0, numBins); //using tracks
  TH1D* recoTkMET_thresh2 = new TH1D("recoTkMET_thresh2", ";TrackMET threshold [GeV];Events", numBins, 0, numBins);
  TH1D* recoTkMET_thresh3 = new TH1D("recoTkMET_thresh3", ";TrackMET threshold [GeV];Events", numBins, 0, numBins);
  TH1D* recoTkMET_thresh4 = new TH1D("recoTkMET_thresh4", ";TrackMET threshold [GeV];Events", numBins, 0, numBins);
  TH1D* recoTkMET_thresh5 = new TH1D("recoTkMET_thresh5", ";TrackMET threshold [GeV];Events", numBins, 0, numBins);
  TH1D* recoTkMET_thresh6 = new TH1D("recoTkMET_thresh6", ";TrackMET threshold [GeV];Events", numBins, 0, numBins);


  /////////////////////////////
  /////  Start of events  /////
  /////////////////////////////
  int nevt = tree->GetEntries();
  if (numEventsToRun>0) std::cout<<"Only running "<<numEventsToRun<<"/"<<nevt<<" events..."<<std::endl;
  else std::cout<<"Running "<<nevt<<" events..."<<std::endl;
  for (int i = 0; i < nevt; i++) {
    tree->GetEntry(i);
    if (i%10000==0) {
      if (numEventsToRun>0) std::cout << i << "/" << numEventsToRun << std::endl;
      else std::cout << i << "/" << nevt << std::endl;
    }
    if (i>numEventsToRun) break;

    float trueTkMET_ntuple = 0;  //grab these from ntuple
    float recoTkMET_ntuple = 0;  //grab these from ntuple
    float recoTkMETEmu_ntuple = 0;  //grab these from ntuple

    if (!rerun_trueTkMET) trueTkMET_ntuple = trueTkMET;
    if (!rerun_recoTkMET) recoTkMET_ntuple = trkMET;
    recoTkMETEmu_ntuple = trkMETEmu;

    if (rerun_trueTkMET) {
      float trueTkMET_calc = 0; float trueTkMETx = 0; float trueTkMETy = 0;
      for (int it = 0; it < (int)tp_pt->size(); it++) {
        float this_tp_pt = tp_pt->at(it);
        float this_tp_eta = tp_eta->at(it);
        float this_tp_phi = tp_phi->at(it);
        int this_tp_signal = tp_eventid->at(it);
        // float deltaZ = fabs(tp_z0->at(it) - pv_L1reco->at(0));
        float deltaZ = fabs(tp_z0->at(it) - pv_L1reco_emu->at(0));

        // kinematic cuts
        if (tp_pt->at(it)<TP_minPt || fabs(this_tp_eta)>TP_maxEta || tp_nstub->at(it)<4 || fabs(tp_z0->at(it))>15 || tp_charge->at(it)==0) continue;
        // if (tp_dxy->at(it) > 1) continue;

        float deltaZ_cut = 3.0;  // cuts out PU
        for (unsigned int reg = 0; reg < numEtaReg; reg++) {
          if (fabs(this_tp_eta) >= etaRegions[reg] && fabs(this_tp_eta) < etaRegions[reg+1]) {
            deltaZ_cut = z0Thresholds[reg];
            break;
          }
        }

        if (deltaZ > deltaZ_cut) continue;
        trueTkMETx += this_tp_pt * cos(this_tp_phi);
        trueTkMETy += this_tp_pt * sin(this_tp_phi);
      }  // end of tracking particle loop

      trueTkMET_calc = sqrt(trueTkMETx * trueTkMETx + trueTkMETy * trueTkMETy);
      trueTkMET_ntuple = trueTkMET_calc;
    }  //re-run trueTkMET


    float recoTkMET_calc = 0; float recoTkMET_calc2 = 0; float recoTkMET_calc3 = 0; float recoTkMET_calc4 = 0; float recoTkMET_calc5 = 0; float recoTkMET_calc6 = 0;
    if (rerun_recoTkMET) {
      float recoTkMETx = 0; float recoTkMETy = 0;
      float recoTkMETx2 = 0; float recoTkMETy2 = 0;
      float recoTkMETx3 = 0; float recoTkMETy3 = 0;
      float recoTkMETx4 = 0; float recoTkMETy4 = 0;
      float recoTkMETx5 = 0; float recoTkMETy5 = 0;
      float recoTkMETx6 = 0; float recoTkMETy6 = 0;

      for (int it = 0; it < (int)trk_pt->size(); it++) {
        float thisTrk_pt = trk_pt->at(it);
        float thisTrk_eta = trk_eta->at(it);
        int thisTrk_nstub = trk_nstub->at(it);
        float thisTrk_chi2dof = trk_chi2dof->at(it);
        float thisTrk_chi2rz = trk_chi2rz->at(it);
        float thisTrk_chi2rphi = trk_chi2rphi->at(it);
        float thisTrk_bendchi2 = trk_bendchi2->at(it);
        float thisTrk_phi = trk_phi->at(it);
        float thisTrk_MVA = trk_MVA1->at(it);
        int thisTrk_fake = trk_fake->at(it);
        float thisTrk_z0 = trk_z0->at(it);
        // float deltaZ = thisTrk_z0 - pv_L1reco->at(0);
        float deltaZ = thisTrk_z0 - pv_L1reco_emu->at(0);
        if (thisTrk_pt < TP_minPt || fabs(thisTrk_eta) > TP_maxEta || thisTrk_nstub < 4 || fabs(thisTrk_z0) > 15) continue;

        float deltaZ_cut = 3.0;  // cuts out PU
        for (unsigned int reg = 0; reg < numEtaReg; reg++) {
          if (fabs(thisTrk_eta) >= etaRegions[reg] && fabs(thisTrk_eta) < etaRegions[reg+1]) {
            deltaZ_cut = z0Thresholds[reg];
            break;
          }
        }

        recoTkMETx += thisTrk_pt * cos(thisTrk_phi);
        recoTkMETy += thisTrk_pt * sin(thisTrk_phi);

        if (fabs(deltaZ)<3.0) {
          recoTkMETx2 += thisTrk_pt * cos(thisTrk_phi);
          recoTkMETy2 += thisTrk_pt * sin(thisTrk_phi);
        }

        if (fabs(deltaZ)<1.0) {
          recoTkMETx3 += thisTrk_pt * cos(thisTrk_phi);
          recoTkMETy3 += thisTrk_pt * sin(thisTrk_phi);
        }

        if (fabs(deltaZ)<deltaZ_cut && thisTrk_chi2rz<5.0 && thisTrk_chi2rphi<20.0 && thisTrk_bendchi2 < 2.25) {
          recoTkMETx4 += thisTrk_pt * cos(thisTrk_phi);
          recoTkMETy4 += thisTrk_pt * sin(thisTrk_phi);
        }

        if (fabs(deltaZ)<3.0 && thisTrk_chi2rz<5.0 && thisTrk_chi2rphi<20.0 && thisTrk_bendchi2 < 2.25) {
          recoTkMETx5 += thisTrk_pt * cos(thisTrk_phi);
          recoTkMETy5 += thisTrk_pt * sin(thisTrk_phi);
        }

        if (fabs(deltaZ)<3.0 && thisTrk_chi2rz<20.0 && thisTrk_chi2rphi<40.0 && thisTrk_bendchi2 < 4.0) {
          recoTkMETx6 += thisTrk_pt * cos(thisTrk_phi);
          recoTkMETy6 += thisTrk_pt * sin(thisTrk_phi);
        }
      }  //end loop over tracks

      recoTkMET_calc = sqrt(recoTkMETx * recoTkMETx + recoTkMETy * recoTkMETy);
      recoTkMET_calc2 = sqrt(recoTkMETx2 * recoTkMETx2 + recoTkMETy2 * recoTkMETy2);
      recoTkMET_calc3 = sqrt(recoTkMETx3 * recoTkMETx3 + recoTkMETy3 * recoTkMETy3);
      recoTkMET_calc4 = sqrt(recoTkMETx4 * recoTkMETx4 + recoTkMETy4 * recoTkMETy4);
      recoTkMET_calc5 = sqrt(recoTkMETx5 * recoTkMETx5 + recoTkMETy5 * recoTkMETy5);
      recoTkMET_calc6 = sqrt(recoTkMETx6 * recoTkMETx6 + recoTkMETy6 * recoTkMETy6);
      recoTkMET_ntuple = recoTkMET_calc;
    }  //re-run reco tkMET

    trueTkMET_thresh->Fill(trueTkMET_ntuple);
    recoTkMET_thresh->Fill(recoTkMET_ntuple);
    recoTkMET_thresh2->Fill(recoTkMET_calc2);
    recoTkMET_thresh3->Fill(recoTkMET_calc3);
    recoTkMET_thresh4->Fill(recoTkMET_calc4);
    recoTkMET_thresh5->Fill(recoTkMET_calc5);
    // recoTkMET_thresh6->Fill(recoTkMET_calc6);
    recoTkMET_thresh6->Fill(recoTkMETEmu_ntuple);
  }  // end event loop

  // -------------------------------------------------------------------------------------------
  trueTkMET_thresh->SetLineColor(kBlack);
  recoTkMET_thresh->SetLineColor(kRed);
  recoTkMET_thresh2->SetLineColor(kBlue);
  recoTkMET_thresh3->SetLineColor(kGreen+1);
  recoTkMET_thresh4->SetLineColor(kCyan);
  recoTkMET_thresh5->SetLineColor(kMagenta);
  recoTkMET_thresh6->SetLineColor(kGray);

  TH1D* fCumulative_true = GetCumulative(trueTkMET_thresh);
  TH1D* fCumulative_reco = GetCumulative(recoTkMET_thresh);
  TH1D* fCumulative_reco2 = GetCumulative(recoTkMET_thresh2);
  TH1D* fCumulative_reco3 = GetCumulative(recoTkMET_thresh3);
  TH1D* fCumulative_reco4 = GetCumulative(recoTkMET_thresh4);
  TH1D* fCumulative_reco5 = GetCumulative(recoTkMET_thresh5);
  TH1D* fCumulative_reco6 = GetCumulative(recoTkMET_thresh6);

  TCanvas* can = new TCanvas("can", "can");
  can->SetGridx(); can->SetGridy(); can->SetLogy();
  // fCumulative_true->Draw();
  // fCumulative_reco->Draw("same");
  // fCumulative_reco2->Draw("same");
  // fCumulative_reco3->Draw("same");
  fCumulative_reco4->Draw();
  // fCumulative_reco5->Draw("same");
  fCumulative_reco6->Draw("same");

  TLegend* leg = new TLegend(0.48, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextSize(0.04); leg->SetTextFont(42);
  // leg->AddEntry(fCumulative_true, "trueTkMET", "lp");
  leg->AddEntry(fCumulative_reco4, "recoTkMET, sim", "lp");
  leg->AddEntry(fCumulative_reco6, "recoTkMET, emu", "lp");
  leg->Draw("same");


  TLine* line_rate = new TLine(0, 35, numBins, 35);
  line_rate->SetLineColor(kBlack); line_rate->SetLineWidth(2);
  line_rate->Draw("SAME");
  can->SaveAs("trkMET_minBias_RatePlot.root");

  // printThresh(fCumulative_true, 35);
  // printThresh(fCumulative_reco, 35);
  // printThresh(fCumulative_reco2, 35);
  // printThresh(fCumulative_reco3, 35);
  printThresh(fCumulative_reco4, 35);
  // printThresh(fCumulative_reco5, 35);
  printThresh(fCumulative_reco6, 35);
}

TH1D* GetCumulative(TH1D* plot) {
  std::string newName = Form("cumulative_%s", plot->GetName());
  TH1D* temp = (TH1D*)plot->Clone(newName.c_str());
  temp->SetDirectory(0);
  for (int i = 0; i < plot->GetNbinsX()+1; i++) {
    double content(0.0), error2(0.0);
    for (int j = i; j < plot->GetNbinsX()+1; j++) {
      content+=plot->GetBinContent(j);
      error2+=plot->GetBinError(j)*plot->GetBinError(j);
    }
    temp->SetBinContent(i, content);
    temp->SetBinError(i, TMath::Sqrt(error2));
  }
  temp->SetMarkerSize(0.7); temp->SetMarkerStyle(20);
  temp->SetMarkerColor(temp->GetLineColor());
  temp->Scale(3.0e4/temp->GetBinContent(1)); //Or should this be 3*10^4? Or 4*10^3?
  temp->GetYaxis()->SetRangeUser(1E-01, 1E05);
  temp->SetTitle("; Track MET [GeV]; L1 Rate [kHz]");
  temp->SetStats(0);
  return temp;
}

void printThresh(TH1D* histo, float thresh) {
  std::cout << "Track MET threshold at "<<noZeroes(thresh)<<" kHz: " << histo->GetBinLowEdge(histo->FindLastBinAbove(thresh)) <<" GeV"<< std::endl;
}

string noZeroes(float value) {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2) << value;
	std::string str = ss.str();
	if(str.find('.') != std::string::npos) {
		str = str.substr(0, str.find_last_not_of('0')+1);
		if(str.find('.') == str.size()-1) {
			str = str.substr(0, str.size()-1);
		}
	}
	return str;
}
