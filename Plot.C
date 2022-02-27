//Runs track efficiency and fake rate
//Runs trackMET resolution, residual/bias, and trigger efficiency
//TrackMET rate run in different script

//Input sample needs to be full GTT ntuplizer, or at least contain the tracks, tracking particles, genMET (or gen particles), and PV
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

string noZeroes(float value);
void makePrettyHisto(TH1F* &h, Color_t color);
std::vector< std::vector<float> > resSigAndErr(std::vector<TH1F*> histos);
// ----------------------------------------------------------------------------------------------------------------
// Main script
void Plot(TString type="Stop200PU") {
	gROOT->SetBatch();
	gErrorIgnoreLevel = kWarning;

	TChain* tree = new TChain("L1TrackNtuple/eventTree");
	if (type == "TTBar200PU") tree->Add("../ForEmail/ttbarPU200_D49_extended.root");
	else if (type == "DispMu200PU") tree->Add("../ForEmail/DispMu_PU200_D49_extended_trunc.root");
	else if (type == "TTBar0PU") tree->Add("../ForEmail/TTbar_PU0_D49.root");
	else if (type == "FlatMu0PU") tree->Add("../ForEmail/SingleMuFlatPt2To100_PU0_D49.root");
	else if (type=="Stop200PU") tree->Add("../ForEmail/CheckingMET_stopPU200_CMSSW1117.root");

	if (tree->GetEntries() == 0) {
		cout << "File doesn't exist or is empty, returning..." << endl;
		return;
	}

	//bools for which plots to run
	bool doEfficiency = true;
	bool doFakeRate = true;
	bool doTurnOn = true;
	bool doTkMETRes = true;
	bool doDeltaZRes = true;
	bool fullGTTFile = true;


	// Quality cuts
	float TP_minPt = 2.0;
	float TP_maxEta = 2.4;
	//float chi2dof_cut = 10.0;
	float chi2rz_cut = 5.0;
	float chi2rphi_cut = 20.0;
	float bendchi2_cut = 2.25;
	int numEventsToRun = 8000; //anything greater than 0 will only run that number of events

	//Delta z cuts
	const int numEtaReg = 7;
	float z0Thresholds[numEtaReg] = {0.37, 0.5, 0.6, 0.75, 1.0, 1.6};
	float etaRegions[numEtaReg] = {0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4};

	//Bins for trkMET resolution
	float bins[]={0,20,30,40,55,75,120,150,250};
	int binnum=8;


	//gen particles
	vector<float>* gen_pt;
	vector<float>* gen_phi;
	vector<float>* gen_pdgid;
	vector<float>* gen_z0;

	// tracking particles
	vector<float>* tp_pt;
	vector<float>* tp_eta;
	vector<float>* tp_phi;
	vector<float>* tp_dxy;
	vector<float>* tp_z0;
	vector<float>* tp_d0;
	vector<int>*   tp_pdgid;
	vector<int>*   tp_nmatch;
	vector<int>*   tp_nstub;
	vector<int>*   tp_eventid;
	vector<int>*   tp_charge;

	// Matched tracks
	vector<float>* matchtrk_pt;
	vector<float>* matchtrk_eta;
	vector<float>* matchtrk_phi;
	vector<float>* matchtrk_d0;
	vector<float>* matchtrk_z0;
	vector<float>* matchtrk_chi2;
	vector<float>* matchtrk_chi2rphi;
	vector<float>* matchtrk_chi2rz;
	vector<float>* matchtrk_chi2dof;
	vector<float>* matchtrk_bendchi2;
	vector<float>* matchtrk_MVA1;
	vector<int>*   matchtrk_nstub;
	vector<int>*   matchtrk_seed;

	// all tracks
	vector<float>* trk_pt;
	vector<float>* trk_eta;
	vector<float>* trk_phi;
	vector<float>* trk_z0;
	vector<float>* trk_chi2;
	vector<float>* trk_chi2rphi;
	vector<float>* trk_chi2rz;
	vector<float>* trk_chi2dof;
	vector<float>* trk_bendchi2;
	vector<float>* trk_MVA1;
	vector<int>*   trk_nstub;
	vector<int>*   trk_seed;
	vector<int>*   trk_fake;
	vector<int>*   trk_genuine;

	float trueMET; // from gen neutrinos and any SUSY particles (for SUSY stop sample)
	float trkMET; // from ntuple
	float trkMETEmu; // from ntuple
	vector<float>* pv_L1reco; // From L1 tracks, simulated PV
	vector<float>* pv_L1reco_emu; // From L1 tracks, emulated PV
	vector<float>* pv_MC; // from gen particles

	TBranch* b_gen_pt;
	TBranch* b_gen_phi;
	TBranch* b_gen_pdgid;
	TBranch* b_gen_z0;

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
	TBranch* b_matchtrk_chi2rz;
	TBranch* b_matchtrk_chi2rphi;
	TBranch* b_matchtrk_chi2dof;
	TBranch* b_matchtrk_bendchi2;
	TBranch* b_matchtrk_MVA1;
	TBranch* b_matchtrk_nstub;

	TBranch* b_trk_pt;
	TBranch* b_trk_eta;
	TBranch* b_trk_phi;
	TBranch* b_trk_z0;
	TBranch* b_trk_chi2;
	TBranch* b_trk_chi2rz;
	TBranch* b_trk_chi2rphi;
	TBranch* b_trk_chi2dof;
	TBranch* b_trk_bendchi2;
	TBranch* b_trk_MVA1;
	TBranch* b_trk_nstub;
	TBranch* b_trk_fake;

	TBranch* b_trueMET;
	TBranch* b_trkMET;
	TBranch* b_trkMETEmu;
	TBranch* b_pv_L1reco;
	TBranch* b_pv_L1reco_emu;
	TBranch* b_pv_MC;

	gen_pt = 0;
	gen_phi = 0;
	gen_pdgid = 0;
	gen_z0 = 0;

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
	matchtrk_chi2rz = 0;
	matchtrk_chi2rphi = 0;
	matchtrk_chi2dof = 0;
	matchtrk_bendchi2 = 0;
	matchtrk_MVA1 = 0;
	matchtrk_nstub = 0;

	trk_pt = 0;
	trk_eta = 0;
	trk_phi = 0;
	trk_z0 = 0;
	trk_chi2 = 0;
	trk_chi2rz = 0;
	trk_chi2rphi = 0;
	trk_chi2dof = 0;
	trk_bendchi2 = 0;
	trk_MVA1 = 0;
	trk_nstub = 0;
	trk_fake = 0;

	trueMET = 0;
	trkMET = 0;
	trkMETEmu = 0;
	pv_L1reco = 0;
	pv_L1reco_emu = 0;
	pv_MC = 0;

	if (fullGTTFile) {
		tree->SetBranchAddress("trueMET",   &trueMET,   &b_trueMET);
		tree->SetBranchAddress("trkMET",   &trkMET,   &b_trkMET);
		tree->SetBranchAddress("trkMETEmu",   &trkMETEmu,   &b_trkMETEmu);
		tree->SetBranchAddress("pv_MC",   &pv_MC,   &b_pv_MC);
		tree->SetBranchAddress("pv_L1reco",   &pv_L1reco,   &b_pv_L1reco);
		tree->SetBranchAddress("pv_L1reco_emu",   &pv_L1reco_emu,   &b_pv_L1reco_emu);
		tree->SetBranchAddress("gen_pt",     &gen_pt,     &b_gen_pt);
		tree->SetBranchAddress("gen_phi",     &gen_phi,     &b_gen_phi);
		tree->SetBranchAddress("gen_pdgid",     &gen_pdgid,     &b_gen_pdgid);
		tree->SetBranchAddress("gen_z0",     &gen_z0,     &b_gen_z0);
	}

	tree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
	tree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
	tree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
	tree->SetBranchAddress("tp_dxy",    &tp_dxy,    &b_tp_dxy);
	tree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
	tree->SetBranchAddress("tp_d0",     &tp_d0,     &b_tp_d0);
	tree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
	tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
	tree->SetBranchAddress("tp_nstub",      &tp_nstub,      &b_tp_nstub);
	tree->SetBranchAddress("tp_eventid",    &tp_eventid,    &b_tp_eventid);
	tree->SetBranchAddress("tp_charge",    &tp_charge,    &b_tp_charge);

	tree->SetBranchAddress("matchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
	tree->SetBranchAddress("matchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
	tree->SetBranchAddress("matchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
	tree->SetBranchAddress("matchtrk_d0",    &matchtrk_d0,    &b_matchtrk_d0);
	tree->SetBranchAddress("matchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
	tree->SetBranchAddress("matchtrk_chi2",  &matchtrk_chi2,  &b_matchtrk_chi2);
	tree->SetBranchAddress("matchtrk_chi2rz",  &matchtrk_chi2rz,  &b_matchtrk_chi2rz);
	tree->SetBranchAddress("matchtrk_chi2rphi",  &matchtrk_chi2rphi,  &b_matchtrk_chi2rphi);
	tree->SetBranchAddress("matchtrk_chi2dof",  &matchtrk_chi2dof,  &b_matchtrk_chi2dof);
	tree->SetBranchAddress("matchtrk_bendchi2",  &matchtrk_bendchi2,  &b_matchtrk_bendchi2);
	tree->SetBranchAddress("matchtrk_MVA1",  &matchtrk_MVA1,  &b_matchtrk_MVA1);
	tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);

	tree->SetBranchAddress("trk_pt",   &trk_pt,   &b_trk_pt);
	tree->SetBranchAddress("trk_eta",  &trk_eta,  &b_trk_eta);
	tree->SetBranchAddress("trk_phi",  &trk_phi,  &b_trk_phi);
	tree->SetBranchAddress("trk_z0",  &trk_z0,  &b_trk_z0);
	tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
	tree->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
	tree->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
	tree->SetBranchAddress("trk_chi2dof", &trk_chi2dof, &b_trk_chi2dof);
	tree->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
	tree->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);
	tree->SetBranchAddress("trk_nstub",   &trk_nstub,   &b_trk_nstub);
	tree->SetBranchAddress("trk_fake",    &trk_fake,    &b_trk_fake);

	// ----------------------------------------------------------------------------------------------------------------
	// for efficiencies
	TH1F* h_tp_pt   = new TH1F("tp_pt", ";Matched track p_{T} [GeV]; Matched tracks", 200,  0,   200.0);
	TH1F* h_tp_pt_zoom= new TH1F("tp_pt_zoom", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_tp_eta  = new TH1F("tp_eta", ";Matched track #eta; Matched tracks",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt   = new TH1F("match_tp_pt", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_eta  = new TH1F("match_tp_eta", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);

	//numerators
	TH1F* h_match_tp_pt_0   = new TH1F("match_tp_pt_0", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_eta_0  = new TH1F("match_tp_eta_0", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_zoom_0   = new TH1F("match_tp_pt_zoom_0", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_pt_1   = new TH1F("match_tp_pt_1", ";Matched track p_{T} [GeV]; Matched tracks", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_1  = new TH1F("match_tp_pt_zoom_1", ";Matched track p_{T} [GeV]; Matched tracks", 100, 0, 8.0);
	TH1F* h_match_tp_eta_1  = new TH1F("match_tp_eta_1", ";Matched track #eta; Matched tracks",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_2   = new TH1F("match_tp_pt_2", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_2  = new TH1F("match_tp_pt_zoom_2", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_eta_2  = new TH1F("match_tp_eta_2", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_3   = new TH1F("match_tp_pt_3", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_3  = new TH1F("match_tp_pt_zoom_3", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_eta_3  = new TH1F("match_tp_eta_3", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_4   = new TH1F("match_tp_pt_4", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_4  = new TH1F("match_tp_pt_zoom_4", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_eta_4  = new TH1F("match_tp_eta_4", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_5   = new TH1F("match_tp_pt_5", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_5  = new TH1F("match_tp_pt_zoom_5", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_eta_5  = new TH1F("match_tp_eta_5", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);
	TH1F* h_match_tp_pt_6   = new TH1F("match_tp_pt_6", ";Tracking particle p_{T} [GeV]; Tracking particles", 200,  0,   200.0);
	TH1F* h_match_tp_pt_zoom_6  = new TH1F("match_tp_pt_zoom_6", ";Tracking particle p_{T} [GeV]; Tracking particles", 100, 0, 8.0);
	TH1F* h_match_tp_eta_6  = new TH1F("match_tp_eta_6", ";Tracking particle #eta; Tracking particles",             50, -2.4,   2.4);

	// Fake rate, denoms
	TH1F* h_trk_pt0 = new TH1F("trk_pt0", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta0 = new TH1F("trk_eta0",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt1 = new TH1F("trk_pt1", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta1 = new TH1F("trk_eta1",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt2 = new TH1F("trk_pt2", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta2 = new TH1F("trk_eta2",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt3 = new TH1F("trk_pt3", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta3 = new TH1F("trk_eta3",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt4 = new TH1F("trk_pt4", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta4 = new TH1F("trk_eta4",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt5 = new TH1F("trk_pt5", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta5 = new TH1F("trk_eta5",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trk_pt6 = new TH1F("trk_pt6", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trk_eta6 = new TH1F("trk_eta6",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);

	//Fake rate, numerators
	TH1F* h_trkFake_pt_0 = new TH1F("trk_fake_pt_0", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_0 = new TH1F("trk_fake_eta_0",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_1 = new TH1F("trk_fake_pt_1", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_1 = new TH1F("trk_fake_eta_1",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_2 = new TH1F("trk_fake_pt_2", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_2 = new TH1F("trk_fake_eta_2",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_3 = new TH1F("trk_fake_pt_3", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_3 = new TH1F("trk_fake_eta_3",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_4 = new TH1F("trk_fake_pt_4", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_4 = new TH1F("trk_fake_eta_4",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_5 = new TH1F("trk_fake_pt_5", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_5 = new TH1F("trk_fake_eta_5",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);
	TH1F* h_trkFake_pt_6 = new TH1F("trk_fake_pt_6", ";L1 track p_{T} [GeV]; L1 tracks", 200,  0, 200.0);
	TH1F* h_trkFake_eta_6 = new TH1F("trk_fake_eta_6",";L1 track #eta; L1 tracks", 50, -2.4, 2.4);

	// MET
	TH1F* h_trueMET = new TH1F("trueMET", ";trueMET [GeV]; Events", 40, 0, 800.0); //used for trigger efficiency, so different binning
	TH1F* h_trueTkMET = new TH1F("trueTkMET", ";trueTkMET [GeV]; Events",  300, 0, 300.0);
	TH1F* h_trueTkMET_tightZ = new TH1F("trueTkMET_tightZ", ";trueTkMET [GeV]; Events",  300, 0, 300.0);
	TH1F* h_recoTkMET_1 = new TH1F("recoTkMET_1", ";recoTkMET [GeV]; Events", 50, 0, 300.0);
	TH1F* h_recoTkMET_2 = new TH1F("recoTkMET_2", ";recoTkMET [GeV]; Events", 50, 0, 300.0);
	TH1F* h_recoTkMET_3 = new TH1F("recoTkMET_3", ";recoTkMET [GeV]; Events", 50, 0, 300.0);
	TH1F* h_recoTkMET_4 = new TH1F("recoTkMET_4", ";recoTkMET [GeV]; Events", 50, 0, 300.0);
	TH1F* h_recoTkMET_5 = new TH1F("recoTkMET_5", ";recoTkMET [GeV]; Events", 50, 0, 300.0);
	TH1F* h_recoTkMET_6 = new TH1F("recoTkMET_6", ";recoTkMET [GeV]; Events", 50, 0, 300.0);

	// recoTkMET resolution, binned in trueTkkMET
	TH1F* h_recoTkMET_1_bin1 = new TH1F("recoTkMET_1_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin2 = new TH1F("recoTkMET_1_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin3 = new TH1F("recoTkMET_1_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin4 = new TH1F("recoTkMET_1_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin5 = new TH1F("recoTkMET_1_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin6 = new TH1F("recoTkMET_1_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin7 = new TH1F("recoTkMET_1_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_1_bin8 = new TH1F("recoTkMET_1_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	TH1F* h_recoTkMET_2_bin1 = new TH1F("recoTkMET_2_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin2 = new TH1F("recoTkMET_2_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin3 = new TH1F("recoTkMET_2_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin4 = new TH1F("recoTkMET_2_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin5 = new TH1F("recoTkMET_2_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin6 = new TH1F("recoTkMET_2_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin7 = new TH1F("recoTkMET_2_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_2_bin8 = new TH1F("recoTkMET_2_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	TH1F* h_recoTkMET_3_bin1 = new TH1F("recoTkMET_3_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin2 = new TH1F("recoTkMET_3_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin3 = new TH1F("recoTkMET_3_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin4 = new TH1F("recoTkMET_3_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin5 = new TH1F("recoTkMET_3_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin6 = new TH1F("recoTkMET_3_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin7 = new TH1F("recoTkMET_3_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_3_bin8 = new TH1F("recoTkMET_3_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	TH1F* h_recoTkMET_4_bin1 = new TH1F("recoTkMET_4_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin2 = new TH1F("recoTkMET_4_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin3 = new TH1F("recoTkMET_4_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin4 = new TH1F("recoTkMET_4_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin5 = new TH1F("recoTkMET_4_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin6 = new TH1F("recoTkMET_4_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin7 = new TH1F("recoTkMET_4_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_4_bin8 = new TH1F("recoTkMET_4_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	TH1F* h_recoTkMET_5_bin1 = new TH1F("recoTkMET_5_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin2 = new TH1F("recoTkMET_5_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin3 = new TH1F("recoTkMET_5_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin4 = new TH1F("recoTkMET_5_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin5 = new TH1F("recoTkMET_5_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin6 = new TH1F("recoTkMET_5_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin7 = new TH1F("recoTkMET_5_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_5_bin8 = new TH1F("recoTkMET_5_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	TH1F* h_recoTkMET_6_bin1 = new TH1F("recoTkMET_6_bin1", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin2 = new TH1F("recoTkMET_6_bin2", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin3 = new TH1F("recoTkMET_6_bin3", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin4 = new TH1F("recoTkMET_6_bin4", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin5 = new TH1F("recoTkMET_6_bin5", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin6 = new TH1F("recoTkMET_6_bin6", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin7 = new TH1F("recoTkMET_6_bin7", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);
	TH1F* h_recoTkMET_6_bin8 = new TH1F("recoTkMET_6_bin8", ";recoTkMET - trueTkMET [GeV]; Events", 50, -100, 100.0);

	// trigger efficiency numerators (all denominators are trueMET)
	TH1F* h_trueTkMET_turnon = new TH1F("trueTkMET_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_trueTkMET_tightZ_turnon = new TH1F("trueTkMET_tightZ_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_1_turnon = new TH1F("recoTkMET_1_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_3_turnon = new TH1F("recoTkMET_3_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_2_turnon = new TH1F("recoTkMET_2_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_4_turnon = new TH1F("recoTkMET_4_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_5_turnon = new TH1F("recoTkMET_5_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_6_turnon = new TH1F("recoTkMET_6_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_500_turnon = new TH1F("recoTkMET_500_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);
	TH1F* h_recoTkMET_200_turnon = new TH1F("recoTkMET_200_turnon", ";trueMET [GeV]; Events", 40, 0, 800.0);

	// delta Z
	TH1F *h_trk_deltaz = new TH1F("trk_deltaz",";pv_MC-trk_z0;All tracks", 50, -5,5);
	TH1F *h_tp_deltaz = new TH1F("tp_deltaz",";pv_MC-trk_z0;All tps", 50, -0.5,0.5);

	TH1F *h_trk_deltaz_eta0to0p7 = new TH1F("trk_deltaz_eta0to0p7",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);
	TH1F *h_trk_deltaz_eta0p7to1 = new TH1F("trk_deltaz_eta0p7to1",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);
	TH1F *h_trk_deltaz_eta1to1p2 = new TH1F("trk_deltaz_eta1to1p2",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);
	TH1F *h_trk_deltaz_eta1p2to1p6 = new TH1F("trk_deltaz_eta1p2to1p6",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);
	TH1F *h_trk_deltaz_eta1p6to2 = new TH1F("trk_deltaz_eta1p6to2",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);
	TH1F *h_trk_deltaz_eta2to2p4 = new TH1F("trk_deltaz_eta2to2p4",";pv_L1reco-matchtrk_z0;All tracks", 50, -2,2);

	TH1F *h_tp_deltaz_eta0to0p7 = new TH1F("tp_deltaz_eta0to0p7",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);
	TH1F *h_tp_deltaz_eta0p7to1 = new TH1F("tp_deltaz_eta0p7to1",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);
	TH1F *h_tp_deltaz_eta1to1p2 = new TH1F("tp_deltaz_eta1to1p2",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);
	TH1F *h_tp_deltaz_eta1p2to1p6 = new TH1F("tp_deltaz_eta1p2to1p6",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);
	TH1F *h_tp_deltaz_eta1p6to2 = new TH1F("tp_deltaz_eta1p6to2",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);
	TH1F *h_tp_deltaz_eta2to2p4 = new TH1F("tp_deltaz_eta2to2p4",";pv_MC-tp_z0;All tracks", 50, -0.3,0.3);

	std::vector<TH1F*> tp_dz_histos = {h_tp_deltaz_eta0to0p7,h_tp_deltaz_eta0p7to1,h_tp_deltaz_eta1to1p2,h_tp_deltaz_eta1p2to1p6,h_tp_deltaz_eta1p6to2,h_tp_deltaz_eta2to2p4};
	std::vector<TH1F*> trk_dz_histos = {h_trk_deltaz_eta0to0p7,h_trk_deltaz_eta0p7to1,h_trk_deltaz_eta1to1p2,h_trk_deltaz_eta1p2to1p6,h_trk_deltaz_eta1p6to2,h_trk_deltaz_eta2to2p4};
	// ----------------------------------------------------------------------------------------------------------------
	//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
	// ----------------------------------------------------------------------------------------------------------------
	int nevt = tree->GetEntries();
	cout << "Number of events = " << nevt << endl;

	// ----------------------------------------------------------------------------------------------------------------
	// event loop
	for (int i=0; i<nevt; i++) {
		if (numEventsToRun>0 && i>numEventsToRun) break;
		tree->GetEntry(i,0);

		//Gen-level MET
		float this_trueMET = 0;
		float trueMETx = 0.0; float trueMETy = 0.0; float trueMET_calc = 0.0;
    for (size_t i = 0; i < gen_pt->size(); ++i) {
			//only gen parts saved are status==1, so don't need that cut
      int id = gen_pdgid->at(i);
      float pt = gen_pt->at(i);
      float phi = gen_phi->at(i);
      bool isNeutrino = false;
      if ((fabs(id) == 12 || fabs(id) == 14 || fabs(id) == 16)) isNeutrino = true;
      if ((isNeutrino || id == 1000022)) {  //MET = neutrino or SUSY stop parent
        trueMETx += pt * cos(phi);
        trueMETy += pt * sin(phi);
      }
    }
    trueMET_calc = sqrt(trueMETx * trueMETx + trueMETy * trueMETy);
    this_trueMET = trueMET_calc;

		float trueTkMET = 0; float trueTkMETx = 0; float trueTkMETy = 0;
		float trueTkMET_tightZ = 0; float trueTkMETx_tightZ = 0; float trueTkMETy_tightZ = 0;
		float recoTkMET_1 = 0; float recoTkMETx_1 = 0; float recoTkMETy_1 = 0;
		float recoTkMET_2 = 0; float recoTkMETx_2 = 0; float recoTkMETy_2 = 0;
		float recoTkMET_3 = 0; float recoTkMETx_3 = 0; float recoTkMETy_3 = 0;
		float recoTkMET_4 = 0; float recoTkMETx_4 = 0; float recoTkMETy_4 = 0;
		float recoTkMET_5 = 0; float recoTkMETx_5 = 0; float recoTkMETy_5 = 0;
		float recoTkMET_6 = 0; float recoTkMETx_6 = 0; float recoTkMETy_6 = 0;

		// tracking particle loop
		for (int it=0; it<(int)tp_pt->size(); it++) {
			float this_tp_pt = tp_pt->at(it);
			float this_tp_eta = tp_eta->at(it);
			float this_tp_phi = tp_phi->at(it);
			int this_tp_eventId = tp_eventid->at(it);
			float deltaZ = tp_z0->at(it)-pv_L1reco_emu->at(0);

			// kinematic cuts
			if (this_tp_eventId!=0) continue; //Signal only
			// if (tp_dxy->at(it) > 1) continue;
			if (tp_charge->at(it)==0 || fabs(this_tp_eta)>TP_maxEta || tp_pt->at(it)<TP_minPt || tp_nstub->at(it)<4 || fabs(tp_z0->at(it))>15) continue;

			//DeltaZ cuts, if you want to apply to the tracking particles
			float deltaZ_cut = 3.0;
			for (unsigned int reg = 0; reg < numEtaReg; reg++) {
				if (fabs(this_tp_eta) >= etaRegions[reg] && fabs(this_tp_eta) < etaRegions[reg+1]) {
					deltaZ_cut = z0Thresholds[reg];
					tp_dz_histos[reg]->Fill(deltaZ);
					break;
				}
			}

			//Fill efficiency denominator
			h_tp_pt->Fill(tp_pt->at(it));
			h_tp_pt_zoom->Fill(tp_pt->at(it));
			h_tp_eta->Fill(tp_eta->at(it));

			trueTkMETx += this_tp_pt*cos(this_tp_phi);
			trueTkMETy += this_tp_pt*sin(this_tp_phi);

			if (fabs(deltaZ)<0.5) {
				trueTkMETx_tightZ += this_tp_pt*cos(this_tp_phi);
				trueTkMETy_tightZ += this_tp_pt*sin(this_tp_phi);
			}

			//Access matched tracks
			if (tp_nmatch->at(it) < 1) continue;
			int nStubs = matchtrk_nstub->at(it);
			float thisTrk_bendchi2 = matchtrk_bendchi2->at(it);
			float chi2 = matchtrk_chi2->at(it);
			float chi2rz = matchtrk_chi2rz->at(it)/(nStubs-2);
			float chi2rphi = matchtrk_chi2rphi->at(it)/(nStubs-2);
			float chi2dof = matchtrk_chi2dof->at(it);
			float deltaZ_match = matchtrk_z0->at(it)-pv_L1reco_emu->at(0);
			float matchtrkEta = matchtrk_eta->at(it);
			float matchtrkz0 = matchtrk_z0->at(it);
			float matchtrkpT = matchtrk_pt->at(it);
			float matchtrk_MVA = matchtrk_MVA1->at(it);
			if (matchtrkpT<0.0) continue;


			float deltaZmatch_cut = 3.0;
			for (unsigned int reg = 0; reg < numEtaReg; reg++) {
				if (fabs(matchtrkEta) >= etaRegions[reg] && fabs(matchtrkEta) < etaRegions[reg+1]) {
					deltaZmatch_cut = z0Thresholds[reg];
					trk_dz_histos[reg]->Fill(deltaZ_match);
					break;
				}
			}


			//Fill efficiency numerators
			h_match_tp_pt_0->Fill(tp_pt->at(it));
			h_match_tp_pt_zoom_0->Fill(tp_pt->at(it));
			h_match_tp_eta_0->Fill(this_tp_eta);

			h_match_tp_pt_1->Fill(this_tp_pt);
			h_match_tp_pt_zoom_1->Fill(this_tp_pt);
			h_match_tp_eta_1->Fill(this_tp_eta);


			// if (fabs(deltaZ_match)>deltaZmatch_cut) continue;
			if (fabs(deltaZ_match)<3.0) {
				h_match_tp_pt_2->Fill(this_tp_pt);
				h_match_tp_pt_zoom_2->Fill(this_tp_pt);
				h_match_tp_eta_2->Fill(this_tp_eta);
			}

			if (fabs(deltaZ_match)<deltaZmatch_cut) {
				h_match_tp_pt_3->Fill(this_tp_pt);
				h_match_tp_pt_zoom_3->Fill(this_tp_pt);
				h_match_tp_eta_3->Fill(this_tp_eta);
			}

			if (fabs(deltaZ_match)<deltaZmatch_cut && chi2rz<5.0 && chi2rphi<20.0 && thisTrk_bendchi2<2.25) {
				h_match_tp_pt_4->Fill(this_tp_pt);
				h_match_tp_pt_zoom_4->Fill(this_tp_pt);
				h_match_tp_eta_4->Fill(this_tp_eta);
			}

			if (fabs(deltaZ_match)<3.0 && chi2rz<5.0 && chi2rphi<20.0 && thisTrk_bendchi2<2.25) {
				h_match_tp_pt_5->Fill(this_tp_pt);
				h_match_tp_pt_zoom_5->Fill(this_tp_pt);
				h_match_tp_eta_5->Fill(this_tp_eta);
			}

			if (fabs(deltaZ_match)<3.0 && chi2rz<chi2rz_cut && chi2rphi<chi2rphi_cut && thisTrk_bendchi2<bendchi2_cut) {
				h_match_tp_pt_6->Fill(this_tp_pt);
				h_match_tp_pt_zoom_6->Fill(this_tp_pt);
				h_match_tp_eta_6->Fill(this_tp_eta);
			}

		} // end of tracking particle loop


		// loop over all tracks
		for (int it=0; it<(int)trk_pt->size(); it++) {
			int thisTrk_nstub = trk_nstub->at(it);
			float thisTrk_pt = trk_pt->at(it);
			float thisTrk_eta = trk_eta->at(it);
			float thisTrk_chi2rz = trk_chi2rz->at(it);
			float thisTrk_chi2rphi = trk_chi2rphi->at(it);
			float thisTrk_chi2dof = trk_chi2dof->at(it);
			float thisTrk_MVA = trk_MVA1->at(it);
			float thisTrk_bendchi2 = trk_bendchi2->at(it);
			float thisTrk_phi = trk_phi->at(it);
			int thisTrk_fake = trk_fake->at(it);
			float thisTrk_z0 = trk_z0->at(it);
			float deltaZ = thisTrk_z0-pv_L1reco_emu->at(0);

			if (fabs(thisTrk_eta)>TP_maxEta || thisTrk_pt<TP_minPt || thisTrk_nstub<4 || fabs(thisTrk_z0)>15) continue;
			if (thisTrk_chi2dof>=100) thisTrk_chi2dof=99.9;
			if (deltaZ>=5.0) deltaZ=4.99; if (deltaZ<=-5.0) deltaZ=-4.99;

			if (thisTrk_fake==0) h_trk_deltaz->Fill(deltaZ);
			float deltaZ_cut = 3.0;
			for (unsigned int reg = 0; reg < numEtaReg; reg++) {
				if (fabs(thisTrk_eta) >= etaRegions[reg] && fabs(thisTrk_eta) < etaRegions[reg+1]) {
					deltaZ_cut = z0Thresholds[reg];
					break;
				}
			}

			h_trk_pt1->Fill(thisTrk_pt); h_trk_eta1->Fill(thisTrk_eta);
			recoTkMETx_1 += thisTrk_pt*cos(thisTrk_phi); recoTkMETy_1 += thisTrk_pt*sin(thisTrk_phi);
			if (thisTrk_fake == 0) {
				h_trkFake_pt_1->Fill(thisTrk_pt); h_trkFake_eta_1->Fill(thisTrk_eta);
			}

			// if (abs(deltaZ) > deltaZ_cut) continue;
			if (abs(deltaZ)<3.0) {
				h_trk_pt2->Fill(thisTrk_pt); h_trk_eta2->Fill(thisTrk_eta);
				recoTkMETx_2 += thisTrk_pt*cos(thisTrk_phi); recoTkMETy_2 += thisTrk_pt*sin(thisTrk_phi);
				if (thisTrk_fake == 0) {
					h_trkFake_pt_2->Fill(thisTrk_pt); h_trkFake_eta_2->Fill(thisTrk_eta);
				}
			}

			if (abs(deltaZ)<deltaZ_cut) {
				h_trk_pt3->Fill(thisTrk_pt); h_trk_eta3->Fill(thisTrk_eta);
				recoTkMETx_3 += thisTrk_pt*cos(thisTrk_phi); recoTkMETy_3 += thisTrk_pt*sin(thisTrk_phi);
				if (thisTrk_fake == 0) {
					h_trkFake_pt_3->Fill(thisTrk_pt); h_trkFake_eta_3->Fill(thisTrk_eta);
				}
			}

			if (abs(deltaZ)<deltaZ_cut && thisTrk_chi2rz<5.0 && thisTrk_chi2rphi<20.0 && thisTrk_bendchi2<2.25) {
				h_trk_pt4->Fill(thisTrk_pt); h_trk_eta4->Fill(thisTrk_eta);
				recoTkMETx_4 += thisTrk_pt*cos(thisTrk_phi); recoTkMETy_4 += thisTrk_pt*sin(thisTrk_phi);
				if (thisTrk_fake == 0) {
					h_trkFake_pt_4->Fill(thisTrk_pt); h_trkFake_eta_4->Fill(thisTrk_eta);
				}
			}

			if (abs(deltaZ)<3.0  && thisTrk_chi2rz<5.0 && thisTrk_chi2rphi<20.0 && thisTrk_bendchi2<2.25) {
				h_trk_pt5->Fill(thisTrk_pt);
				h_trk_eta5->Fill(thisTrk_eta);
				recoTkMETx_5 += thisTrk_pt*cos(thisTrk_phi);
				recoTkMETy_5 += thisTrk_pt*sin(thisTrk_phi);
				if (thisTrk_fake == 0) {
					h_trkFake_pt_5->Fill(thisTrk_pt);
					h_trkFake_eta_5->Fill(thisTrk_eta);
				}
			}

			if (abs(deltaZ)<3.0  && thisTrk_chi2rz<20.0 && thisTrk_chi2rphi<40.0 && thisTrk_bendchi2<4.0) {
				h_trk_pt6->Fill(thisTrk_pt); h_trk_eta6->Fill(thisTrk_eta);
				recoTkMETx_6 += thisTrk_pt*cos(thisTrk_phi); recoTkMETy_6 += thisTrk_pt*sin(thisTrk_phi);
				if (thisTrk_fake == 0) {
					h_trkFake_pt_6->Fill(thisTrk_pt); h_trkFake_eta_6->Fill(thisTrk_eta);
				}
			}

		} //end loop over tracks



		//calculating MET
		if (this_trueMET !=-999 && doTkMETRes) {
			trueTkMET = sqrt(trueTkMETx*trueTkMETx + trueTkMETy*trueTkMETy);
			// trueTkMET_tightZ = sqrt(trueTkMETx_tightZ*trueTkMETx_tightZ + trueTkMETy_tightZ*trueTkMETy_tightZ);
			recoTkMET_1 = sqrt(recoTkMETx_1*recoTkMETx_1 + recoTkMETy_1*recoTkMETy_1);
			recoTkMET_2 = sqrt(recoTkMETx_2*recoTkMETx_2 + recoTkMETy_2*recoTkMETy_2);
			recoTkMET_3 = sqrt(recoTkMETx_3*recoTkMETx_3 + recoTkMETy_3*recoTkMETy_3);
			recoTkMET_4 = sqrt(recoTkMETx_4*recoTkMETx_4 + recoTkMETy_4*recoTkMETy_4);
			recoTkMET_5 = sqrt(recoTkMETx_5*recoTkMETx_5 + recoTkMETy_5*recoTkMETy_5);
			recoTkMET_6 = sqrt(recoTkMETx_6*recoTkMETx_6 + recoTkMETy_6*recoTkMETy_6);


			h_trueTkMET->Fill(trueTkMET);
			// h_trueTkMET_tightZ->Fill(trueTkMET_tightZ);
			h_trueMET->Fill(this_trueMET);

			//thresholds from rate plot, at 35 kHz
			//if reco MET passes, fill with trueMET
			if (trueTkMET>66)  h_trueTkMET_turnon->Fill(this_trueMET);
			// if (trueTkMET_tightZ>60)  h_trueTkMET_tightZ_turnon->Fill(this_trueMET);
			if (recoTkMET_1>98) h_recoTkMET_1_turnon->Fill(this_trueMET);
			if (recoTkMET_2>95) h_recoTkMET_2_turnon->Fill(this_trueMET);
			if (recoTkMET_3>94) h_recoTkMET_3_turnon->Fill(this_trueMET);
			if (recoTkMET_4>66) h_recoTkMET_4_turnon->Fill(this_trueMET);
			if (recoTkMET_5>80) h_recoTkMET_5_turnon->Fill(this_trueMET);
			if (trkMETEmu>66) h_recoTkMET_6_turnon->Fill(this_trueMET);
			// if (recoTkMET_6>84) h_recoTkMET_6_turnon->Fill(this_trueMET);


			float METdiff_1 = recoTkMET_1-trueTkMET;
			float METdiff_2 = recoTkMET_2-trueTkMET;
			float METdiff_3 = recoTkMET_3-trueTkMET;
			float METdiff_4 = recoTkMET_4-trueTkMET;
			float METdiff_5 = recoTkMET_5-trueTkMET;
			float METdiff_6 = trkMETEmu-trueTkMET;
			// float METdiff_6 = recoTkMET_6-trueTkMET;


			if (trueTkMET>=bins[0] && trueTkMET<bins[1]) {
				h_recoTkMET_1_bin1->Fill(METdiff_1);
				h_recoTkMET_2_bin1->Fill(METdiff_2);
				h_recoTkMET_3_bin1->Fill(METdiff_3);
				h_recoTkMET_4_bin1->Fill(METdiff_4);
				h_recoTkMET_5_bin1->Fill(METdiff_5);
				h_recoTkMET_6_bin1->Fill(METdiff_6);
			}
			else if (trueTkMET>=bins[1] && trueTkMET<bins[2]) {
				h_recoTkMET_1_bin2->Fill(METdiff_1);
				h_recoTkMET_2_bin2->Fill(METdiff_2);
				h_recoTkMET_3_bin2->Fill(METdiff_3);
				h_recoTkMET_4_bin2->Fill(METdiff_4);
				h_recoTkMET_5_bin2->Fill(METdiff_5);
				h_recoTkMET_6_bin2->Fill(METdiff_6);
			}
			else if (trueTkMET>=bins[2] && trueTkMET<bins[3]) {
				h_recoTkMET_1_bin3->Fill(METdiff_1);
				h_recoTkMET_2_bin3->Fill(METdiff_2);
				h_recoTkMET_3_bin3->Fill(METdiff_3);
				h_recoTkMET_4_bin3->Fill(METdiff_4);
				h_recoTkMET_5_bin3->Fill(METdiff_5);
				h_recoTkMET_6_bin3->Fill(METdiff_6);
			}
			else if (trueTkMET>=bins[3] && trueTkMET<bins[4]) {
				h_recoTkMET_1_bin4->Fill(METdiff_1);
				h_recoTkMET_2_bin4->Fill(METdiff_2);
				h_recoTkMET_3_bin4->Fill(METdiff_3);
				h_recoTkMET_4_bin4->Fill(METdiff_4);
				h_recoTkMET_5_bin4->Fill(METdiff_5);
				h_recoTkMET_6_bin4->Fill(METdiff_6);
			}
			else if (trueTkMET>=bins[4] && trueTkMET<bins[5]) {
				h_recoTkMET_1_bin5->Fill(METdiff_1);
				h_recoTkMET_2_bin5->Fill(METdiff_2);
				h_recoTkMET_3_bin5->Fill(METdiff_3);
				h_recoTkMET_4_bin5->Fill(METdiff_4);
				h_recoTkMET_5_bin5->Fill(METdiff_5);
				h_recoTkMET_6_bin5->Fill(METdiff_6);
			}

			if (type.Contains("TTBar")) {
				if (trueTkMET>=bins[5]) {
					h_recoTkMET_1_bin6->Fill(METdiff_1);
					h_recoTkMET_2_bin6->Fill(METdiff_2);
					h_recoTkMET_3_bin6->Fill(METdiff_3);
					h_recoTkMET_4_bin6->Fill(METdiff_4);
					h_recoTkMET_5_bin6->Fill(METdiff_5);
					h_recoTkMET_6_bin6->Fill(METdiff_6);
				}
			}

			else if (type.Contains("Stop")) {
				if (trueTkMET>=bins[5] && trueTkMET<bins[6]) {
					h_recoTkMET_1_bin6->Fill(METdiff_1);
					h_recoTkMET_2_bin6->Fill(METdiff_2);
					h_recoTkMET_3_bin6->Fill(METdiff_3);
					h_recoTkMET_4_bin6->Fill(METdiff_4);
					h_recoTkMET_5_bin6->Fill(METdiff_5);
					h_recoTkMET_6_bin6->Fill(METdiff_6);
				}
				else if (trueTkMET>=bins[6] && trueTkMET<bins[7]) {
					h_recoTkMET_1_bin7->Fill(METdiff_1);
					h_recoTkMET_2_bin7->Fill(METdiff_2);
					h_recoTkMET_3_bin7->Fill(METdiff_3);
					h_recoTkMET_4_bin7->Fill(METdiff_4);
					h_recoTkMET_5_bin7->Fill(METdiff_5);
					h_recoTkMET_6_bin7->Fill(METdiff_6);
				}
				else if (trueTkMET>=bins[7]) {
					h_recoTkMET_1_bin8->Fill(METdiff_1);
					h_recoTkMET_2_bin8->Fill(METdiff_2);
					h_recoTkMET_3_bin8->Fill(METdiff_3);
					h_recoTkMET_4_bin8->Fill(METdiff_4);
					h_recoTkMET_5_bin8->Fill(METdiff_5);
					h_recoTkMET_6_bin8->Fill(METdiff_6);
				}
			}

			//for overflow bin
			if (recoTkMET_1>=400) recoTkMET_1=399.99;
			if (recoTkMET_2>=400) recoTkMET_2=399.99;
			if (recoTkMET_3>=400) recoTkMET_3=399.99;
			if (recoTkMET_4>=400) recoTkMET_4=399.99;
			if (recoTkMET_5>=400) recoTkMET_5=399.99;
			if (recoTkMET_6>=400) recoTkMET_6=399.99;
			h_recoTkMET_1->Fill(recoTkMET_1);
			h_recoTkMET_2->Fill(recoTkMET_2);
			h_recoTkMET_3->Fill(recoTkMET_3);
			h_recoTkMET_4->Fill(recoTkMET_4);
			h_recoTkMET_5->Fill(recoTkMET_5);
			h_recoTkMET_6->Fill(trkMETEmu);
			// h_recoTkMET_6->Fill(recoTkMET_6);
		}

	} // end of event loop

	// -------------------------------------------------------------------------------------------
	// output file for histograms
	TString foutName = "test_"+type+".root";
	TFile* fout = new TFile(foutName,"recreate");
	std::cout<<"Output file: "<<foutName<<std::endl;

	TDirectory *cdeff = fout->mkdir("Efficiency");
	TDirectory *cdfake = fout->mkdir("FakeRate");
	TDirectory *cdres = fout->mkdir("TkMETResolution");
	TDirectory *cdturnon = fout->mkdir("TurnOn");
	TDirectory *cdresZ = fout->mkdir("ZResolution");


	// efficiency plots
	if (doEfficiency) {
		cdeff->cd();

		h_tp_pt->Rebin(4); h_tp_pt_zoom->Rebin(4);	h_tp_eta->Rebin(2);
		h_match_tp_pt_0->Rebin(4); h_match_tp_pt_zoom_0->Rebin(4); h_match_tp_eta_0->Rebin(2);
		h_match_tp_pt_1->Rebin(4); h_match_tp_pt_zoom_1->Rebin(4); h_match_tp_eta_1->Rebin(2);
		h_match_tp_pt_2->Rebin(4); h_match_tp_pt_zoom_2->Rebin(4); h_match_tp_eta_2->Rebin(2);
		h_match_tp_pt_3->Rebin(4); h_match_tp_pt_zoom_3->Rebin(4); h_match_tp_eta_3->Rebin(2);
		h_match_tp_pt_4->Rebin(4); h_match_tp_pt_zoom_4->Rebin(4); h_match_tp_eta_4->Rebin(2);
		h_match_tp_pt_5->Rebin(4); h_match_tp_pt_zoom_5->Rebin(4); h_match_tp_eta_5->Rebin(2);
		h_match_tp_pt_6->Rebin(4); h_match_tp_pt_zoom_6->Rebin(4); h_match_tp_eta_6->Rebin(2);


		//sumw2
		h_tp_pt_zoom->Sumw2(); h_tp_pt->Sumw2(); h_tp_eta->Sumw2();
		h_match_tp_pt_0->Sumw2(); h_match_tp_pt_zoom_0->Sumw2(); h_match_tp_eta_0->Sumw2();
		h_match_tp_pt_1->Sumw2(); h_match_tp_pt_zoom_1->Sumw2(); h_match_tp_eta_1->Sumw2();
		h_match_tp_pt_2->Sumw2(); h_match_tp_pt_zoom_2->Sumw2(); h_match_tp_eta_2->Sumw2();
		h_match_tp_pt_3->Sumw2(); h_match_tp_pt_zoom_3->Sumw2(); h_match_tp_eta_3->Sumw2();
		h_match_tp_pt_4->Sumw2(); h_match_tp_pt_zoom_4->Sumw2(); h_match_tp_eta_4->Sumw2();
		h_match_tp_pt_5->Sumw2(); h_match_tp_pt_zoom_5->Sumw2(); h_match_tp_eta_5->Sumw2();
		h_match_tp_pt_6->Sumw2(); h_match_tp_pt_zoom_6->Sumw2(); h_match_tp_eta_6->Sumw2();

		// calculate the effeciency
		TH1F* h_eff_pt_0 = (TH1F*) h_match_tp_pt_0->Clone("eff_pt_0");
		h_eff_pt_0->Divide(h_match_tp_pt_0, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_0, kRed);

		TH1F* h_eff_pt_zoom_0 = (TH1F*) h_match_tp_pt_zoom_0->Clone("eff_pt_zoom_0");
		h_eff_pt_zoom_0->Divide(h_match_tp_pt_zoom_0, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_0, kRed);

		TH1F* h_eff_eta_0 = (TH1F*) h_match_tp_eta_0->Clone("eff_eta_0");
		h_eff_eta_0->Divide(h_match_tp_eta_0, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_0, kRed);

		TH1F* h_eff_pt_1 = (TH1F*) h_match_tp_pt_1->Clone("eff_pt_1");
		h_eff_pt_1->Divide(h_match_tp_pt_1, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_1, kRed);

		TH1F* h_eff_pt_zoom_1 = (TH1F*) h_match_tp_pt_zoom_1->Clone("eff_pt_zoom_1");
		h_eff_pt_zoom_1->Divide(h_match_tp_pt_zoom_1, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_1, kRed);

		TH1F* h_eff_eta_1 = (TH1F*) h_match_tp_eta_1->Clone("eff_eta_1");
		h_eff_eta_1->Divide(h_match_tp_eta_1, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_1, kBlue);

		TH1F* h_eff_pt_2 = (TH1F*) h_match_tp_pt_2->Clone("eff_pt_2");
		h_eff_pt_2->Divide(h_match_tp_pt_2, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_2, kBlue);

		TH1F* h_eff_pt_zoom_2 = (TH1F*) h_match_tp_pt_zoom_2->Clone("eff_pt_zoom_2");
		h_eff_pt_zoom_2->Divide(h_match_tp_pt_zoom_2, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_2, kBlue);

		TH1F* h_eff_eta_2 = (TH1F*) h_match_tp_eta_2->Clone("eff_eta_2");
		h_eff_eta_2->Divide(h_match_tp_eta_2, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_2, kBlue);

		TH1F* h_eff_pt_3 = (TH1F*) h_match_tp_pt_3->Clone("eff_pt_3");
		h_eff_pt_3->Divide(h_match_tp_pt_3, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_3, kGreen+1);

		TH1F* h_eff_pt_zoom_3 = (TH1F*) h_match_tp_pt_zoom_3->Clone("eff_pt_zoom_3");
		h_eff_pt_zoom_3->Divide(h_match_tp_pt_zoom_3, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_3, kGreen+1);

		TH1F* h_eff_eta_3 = (TH1F*) h_match_tp_eta_3->Clone("eff_eta_3");
		h_eff_eta_3->Divide(h_match_tp_eta_3, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_3, kGreen+1);

		TH1F* h_eff_pt_4 = (TH1F*) h_match_tp_pt_4->Clone("eff_pt_4");
		h_eff_pt_4->Divide(h_match_tp_pt_4, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_4, kCyan);

		TH1F* h_eff_pt_zoom_4 = (TH1F*) h_match_tp_pt_zoom_4->Clone("eff_pt_zoom_4");
		h_eff_pt_zoom_4->Divide(h_match_tp_pt_zoom_4, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_4, kCyan);

		TH1F* h_eff_eta_4 = (TH1F*) h_match_tp_eta_4->Clone("eff_eta_4");
		h_eff_eta_4->Divide(h_match_tp_eta_4, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_4, kCyan);

		TH1F* h_eff_pt_5 = (TH1F*) h_match_tp_pt_5->Clone("eff_pt_5");
		h_eff_pt_5->Divide(h_match_tp_pt_5, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_5, kMagenta);

		TH1F* h_eff_pt_zoom_5 = (TH1F*) h_match_tp_pt_zoom_5->Clone("eff_pt_zoom_5");
		h_eff_pt_zoom_5->Divide(h_match_tp_pt_zoom_5, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_5, kMagenta);

		TH1F* h_eff_eta_5 = (TH1F*) h_match_tp_eta_5->Clone("eff_eta_5");
		h_eff_eta_5->Divide(h_match_tp_eta_5, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_5, kMagenta);

		TH1F* h_eff_pt_6 = (TH1F*) h_match_tp_pt_6->Clone("eff_pt_6");
		h_eff_pt_6->Divide(h_match_tp_pt_6, h_tp_pt, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_6, kGray);

		TH1F* h_eff_pt_zoom_6 = (TH1F*) h_match_tp_pt_zoom_6->Clone("eff_pt_zoom_6");
		h_eff_pt_zoom_6->Divide(h_match_tp_pt_zoom_6, h_tp_pt_zoom, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_pt_zoom_6, kGray);

		TH1F* h_eff_eta_6 = (TH1F*) h_match_tp_eta_6->Clone("eff_eta_6");
		h_eff_eta_6->Divide(h_match_tp_eta_6, h_tp_eta, 1.0, 1.0, "B");
		makePrettyHisto(h_eff_eta_6, kGray);

		// set the axis range
		h_eff_pt_0->SetAxisRange(0,1.1,"Y"); h_eff_eta_0->SetAxisRange(0,1.1,"Y");
		h_eff_pt_1->SetAxisRange(0,1.1,"Y"); h_eff_eta_1->SetAxisRange(0,1.1,"Y");
		h_eff_pt_2->SetAxisRange(0,1.1,"Y"); h_eff_eta_2->SetAxisRange(0,1.1,"Y");
		h_eff_pt_3->SetAxisRange(0,1.1,"Y"); h_eff_eta_3->SetAxisRange(0,1.1,"Y");
		h_eff_pt_4->SetAxisRange(0,1.1,"Y"); h_eff_eta_4->SetAxisRange(0,1.1,"Y");
		h_eff_pt_5->SetAxisRange(0,1.1,"Y"); h_eff_eta_5->SetAxisRange(0,1.1,"Y");
		h_eff_pt_6->SetAxisRange(0,1.1,"Y"); h_eff_eta_6->SetAxisRange(0,1.1,"Y");

		// Final canvas
		TLegend* leg_eff = new TLegend(0.60,0.7,0.93,0.91);
	  leg_eff->SetBorderSize(0);
		// leg_eff->AddEntry(h_eff_eta_0,"No cuts (matching efficiency)","lp");
		leg_eff->AddEntry(h_eff_eta_1, "1", "lp");
		leg_eff->AddEntry(h_eff_eta_2, "2", "lp");
		leg_eff->AddEntry(h_eff_eta_3, "3", "lp");
		leg_eff->AddEntry(h_eff_eta_4, "4", "lp");
		leg_eff->AddEntry(h_eff_eta_5, "5", "lp");
		leg_eff->AddEntry(h_eff_eta_6, "6", "lp");


		TCanvas *c_eff_pt = new TCanvas("c_eff_pt","c_eff_pt", 600,600);
		h_eff_pt_0->SetStats(0); h_eff_pt_0->SetMinimum(0.7); h_eff_pt_0->SetMaximum(1.0);
		h_eff_pt_0->SetTitle(";Tracking particle p_{T} [GeV]; Efficiency");
		h_eff_pt_1->SetStats(0); h_eff_pt_1->SetMinimum(0.7); h_eff_pt_1->SetMaximum(1.0);
		h_eff_pt_1->SetTitle(";Tracking particle p_{T} [GeV]; Efficiency");
		// h_eff_pt_0->Draw();
		h_eff_pt_1->Draw("");
		h_eff_pt_2->Draw("same");
		h_eff_pt_3->Draw("same");
		h_eff_pt_4->Draw("same");
		h_eff_pt_5->Draw("same");
		h_eff_pt_6->Draw("same");
		leg_eff->Draw("same");
		c_eff_pt->Write();

		TCanvas *c_eff_pt_zoom = new TCanvas("c_eff_pt_zoom","c_eff_pt_zoom", 600,600);
		h_eff_pt_zoom_0->SetStats(0); h_eff_pt_zoom_0->SetMinimum(0.7); h_eff_pt_zoom_0->SetMaximum(1.0);
		h_eff_pt_zoom_0->SetTitle(";Tracking particle p_{T} [GeV]; Efficiency");
		h_eff_pt_zoom_1->SetStats(0); h_eff_pt_zoom_1->SetMinimum(0.7); h_eff_pt_zoom_1->SetMaximum(1.0);
		h_eff_pt_zoom_1->SetTitle(";Tracking particle p_{T} [GeV]; Efficiency");
		// h_eff_pt_zoom_0->Draw();
		h_eff_pt_zoom_1->Draw("");
		h_eff_pt_zoom_2->Draw("same");
		h_eff_pt_zoom_3->Draw("same");
		h_eff_pt_zoom_4->Draw("same");
		h_eff_pt_zoom_5->Draw("same");
		h_eff_pt_zoom_6->Draw("same");
		leg_eff->Draw("same");
		c_eff_pt_zoom->Write();

		TCanvas *c_eff_eta = new TCanvas("c_eff_eta","c_eff_eta", 600,600);
		h_eff_eta_0->SetStats(0); h_eff_eta_0->SetMinimum(0.7); h_eff_eta_0->SetMaximum(1.0);
		h_eff_eta_0->SetTitle("; Tracking particle #eta;Efficiency");
		h_eff_eta_1->SetStats(0);h_eff_eta_1->SetMinimum(0.7);h_eff_eta_1->SetMaximum(1.0);
		h_eff_eta_1->SetTitle("; Tracking particle #eta;Efficiency");
		// h_eff_eta_0->Draw();
		h_eff_eta_1->Draw("");
		h_eff_eta_2->Draw("same");
		h_eff_eta_3->Draw("same");
		h_eff_eta_4->Draw("same");
		h_eff_eta_5->Draw("same");
		h_eff_eta_6->Draw("same");
		leg_eff->Draw("same");
		c_eff_eta->Write();

		// draw and write individual plots
		h_eff_pt_0->Draw(); h_eff_pt_0->Write();
		h_eff_pt_1->Draw(); h_eff_pt_1->Write();
		h_eff_pt_2->Draw(); h_eff_pt_2->Write();
		h_eff_pt_3->Draw(); h_eff_pt_3->Write();
		h_eff_pt_4->Draw(); h_eff_pt_4->Write();
		h_eff_pt_5->Draw(); h_eff_pt_5->Write();
		h_eff_pt_6->Draw(); h_eff_pt_6->Write();

		h_eff_eta_0->Draw(); h_eff_eta_0->Write();
		h_eff_eta_1->Draw(); h_eff_eta_1->Write();
		h_eff_eta_2->Draw(); h_eff_eta_2->Write();
		h_eff_eta_3->Draw(); h_eff_eta_3->Write();
		h_eff_eta_4->Draw(); h_eff_eta_4->Write();
		h_eff_eta_5->Draw(); h_eff_eta_5->Write();
		h_eff_eta_6->Draw(); h_eff_eta_6->Write();
	} //end if doEfficiency


	// fake rate plots
	if (doFakeRate) {
		cdfake->cd();

		h_trk_pt0->Rebin(4); h_trk_eta0->Rebin(2);
		h_trk_pt1->Rebin(4); h_trk_eta1->Rebin(2);
		h_trk_pt2->Rebin(4); h_trk_eta2->Rebin(2);
		h_trk_pt3->Rebin(4); h_trk_eta3->Rebin(2);
		h_trk_pt4->Rebin(4); h_trk_eta4->Rebin(2);
		h_trk_pt5->Rebin(4); h_trk_eta5->Rebin(2);
		h_trk_pt6->Rebin(4); h_trk_eta6->Rebin(2);

		h_trkFake_pt_0->Rebin(4); h_trkFake_eta_0->Rebin(2);
		h_trkFake_pt_1->Rebin(4); h_trkFake_eta_1->Rebin(2);
		h_trkFake_pt_2->Rebin(4); h_trkFake_eta_2->Rebin(2);
		h_trkFake_pt_3->Rebin(4); h_trkFake_eta_3->Rebin(2);
		h_trkFake_pt_4->Rebin(4); h_trkFake_eta_4->Rebin(2);
		h_trkFake_pt_5->Rebin(4); h_trkFake_eta_5->Rebin(2);
		h_trkFake_pt_6->Rebin(4); h_trkFake_eta_6->Rebin(2);

		// Fake rate
		h_trk_pt0->Sumw2(); h_trk_eta0->Sumw2();
		h_trk_pt1->Sumw2(); h_trk_eta1->Sumw2();
		h_trk_pt2->Sumw2(); h_trk_eta2->Sumw2();
		h_trk_pt3->Sumw2(); h_trk_eta3->Sumw2();
		h_trk_pt4->Sumw2(); h_trk_eta4->Sumw2();
		h_trk_pt5->Sumw2(); h_trk_eta5->Sumw2();
		h_trk_pt6->Sumw2(); h_trk_eta6->Sumw2();

		h_trkFake_pt_0->Sumw2(); h_trkFake_eta_0->Sumw2();
		h_trkFake_pt_1->Sumw2(); h_trkFake_eta_1->Sumw2();
		h_trkFake_pt_2->Sumw2(); h_trkFake_eta_2->Sumw2();
		h_trkFake_pt_3->Sumw2(); h_trkFake_eta_3->Sumw2();
		h_trkFake_pt_4->Sumw2(); h_trkFake_eta_4->Sumw2();
		h_trkFake_pt_5->Sumw2(); h_trkFake_eta_5->Sumw2();
		h_trkFake_pt_6->Sumw2(); h_trkFake_eta_6->Sumw2();


		// Plot the fake rate
		TH1F* h_fakeInt_pt_0 = (TH1F*) h_trkFake_pt_0->Clone("fake_pt_0");
		h_fakeInt_pt_0->Divide(h_trkFake_pt_0, h_trk_pt0, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_0, kRed);

		TH1F* h_fakeInt_eta_0 = (TH1F*) h_trkFake_eta_0->Clone("fake_eta_0");
		h_fakeInt_eta_0->Divide(h_trkFake_eta_0, h_trk_eta0, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_0, kRed);

		TH1F* h_fakeInt_pt_1 = (TH1F*) h_trkFake_pt_1->Clone("fake_pt_1");
		h_fakeInt_pt_1->Divide(h_trkFake_pt_1, h_trk_pt1, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_1, kRed);

		TH1F* h_fakeInt_eta_1 = (TH1F*) h_trkFake_eta_1->Clone("fake_eta_1");
		h_fakeInt_eta_1->Divide(h_trkFake_eta_1, h_trk_eta1, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_1, kRed);

		TH1F* h_fakeInt_pt_2 = (TH1F*) h_trkFake_pt_2->Clone("fake_pt_2");
		h_fakeInt_pt_2->Divide(h_trkFake_pt_2, h_trk_pt2, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_2, kBlue);

		TH1F* h_fakeInt_eta_2 = (TH1F*) h_trkFake_eta_2->Clone("fake_eta_2");
		h_fakeInt_eta_2->Divide(h_trkFake_eta_2, h_trk_eta2, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_2, kBlue);

		TH1F* h_fakeInt_pt_3 = (TH1F*) h_trkFake_pt_3->Clone("fake_pt_3");
		h_fakeInt_pt_3->Divide(h_trkFake_pt_3, h_trk_pt3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_3, kGreen+1);

		TH1F* h_fakeInt_eta_3 = (TH1F*) h_trkFake_eta_3->Clone("fake_eta_3");
		h_fakeInt_eta_3->Divide(h_trkFake_eta_3, h_trk_eta3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_3, kGreen+1);

		TH1F* h_fakeInt_pt_4 = (TH1F*) h_trkFake_pt_4->Clone("fake_pt_4");
		h_fakeInt_pt_4->Divide(h_trkFake_pt_4, h_trk_pt3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_4, kCyan);

		TH1F* h_fakeInt_eta_4 = (TH1F*) h_trkFake_eta_4->Clone("fake_eta_4");
		h_fakeInt_eta_4->Divide(h_trkFake_eta_4, h_trk_eta3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_4, kCyan);

		TH1F* h_fakeInt_pt_5 = (TH1F*) h_trkFake_pt_5->Clone("fake_pt_5");
		h_fakeInt_pt_5->Divide(h_trkFake_pt_5, h_trk_pt3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_5, kMagenta);

		TH1F* h_fakeInt_eta_5 = (TH1F*) h_trkFake_eta_5->Clone("fake_eta_5");
		h_fakeInt_eta_5->Divide(h_trkFake_eta_5, h_trk_eta3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_5, kMagenta);

		TH1F* h_fakeInt_pt_6 = (TH1F*) h_trkFake_pt_6->Clone("fake_pt_6");
		h_fakeInt_pt_6->Divide(h_trkFake_pt_6, h_trk_pt3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_pt_6, kGray);

		TH1F* h_fakeInt_eta_6 = (TH1F*) h_trkFake_eta_6->Clone("fake_eta_6");
		h_fakeInt_eta_6->Divide(h_trkFake_eta_6, h_trk_eta3, 1.0, 1.0, "B");
		makePrettyHisto(h_fakeInt_eta_6, kGray);

		// Final canvas, fake rate
		TLegend* leg_fake = new TLegend(0.60,0.7,0.93,0.91);
	  leg_fake->SetBorderSize(0);
		// leg_fake->AddEntry(h_fakeInt_pt_0,"No cuts","lp");
		leg_fake->AddEntry(h_fakeInt_pt_1, "1","lp");
		leg_fake->AddEntry(h_fakeInt_pt_2, "2","lp");
		leg_fake->AddEntry(h_fakeInt_pt_3, "3","lp");
		leg_fake->AddEntry(h_fakeInt_pt_4, "4","lp");
		leg_fake->AddEntry(h_fakeInt_pt_5, "5","lp");
		leg_fake->AddEntry(h_fakeInt_pt_6, "6","lp");

		TCanvas *c_Intfake_pt = new TCanvas("c_Intfake_pt","c_Intfake_pt", 600,600);
		h_fakeInt_pt_0->SetStats(0); h_fakeInt_pt_0->SetMinimum(0.0); h_fakeInt_pt_0->SetMaximum(0.7);
		h_fakeInt_pt_0->SetTitle("; Track p_{T} [GeV];Fake rate");
		h_fakeInt_pt_1->SetStats(0); h_fakeInt_pt_1->SetMinimum(0.0); h_fakeInt_pt_1->SetMaximum(0.7);
		h_fakeInt_pt_1->SetTitle("; Track p_{T} [GeV];Fake rate");
		// h_fakeInt_pt_0->Draw();
		h_fakeInt_pt_1->Draw("");
		h_fakeInt_pt_2->Draw("same");
		h_fakeInt_pt_3->Draw("same");
		h_fakeInt_pt_4->Draw("same");
		h_fakeInt_pt_5->Draw("same");
		h_fakeInt_pt_6->Draw("same");
		leg_fake->Draw("same");
		c_Intfake_pt->Write();

		TCanvas *c_Intfake_eta = new TCanvas("c_Intfake_eta","c_Intfake_eta", 600,600);
		h_fakeInt_eta_0->SetStats(0); h_fakeInt_eta_0->SetMinimum(0.0); h_fakeInt_eta_0->SetMaximum(0.05);
		h_fakeInt_eta_0->SetTitle("; Track #eta;Fake rate");
		h_fakeInt_eta_1->SetStats(0); h_fakeInt_eta_1->SetMinimum(0.0); h_fakeInt_eta_1->SetMaximum(0.05);
		h_fakeInt_eta_1->SetTitle("; Track #eta;Fake rate");
		// h_fakeInt_eta_0->Draw();
		h_fakeInt_eta_1->Draw("");
		h_fakeInt_eta_2->Draw("same");
		h_fakeInt_eta_3->Draw("same");
		h_fakeInt_eta_4->Draw("same");
		h_fakeInt_eta_5->Draw("same");
		h_fakeInt_eta_6->Draw("same");
		leg_fake->Draw("same");
		c_Intfake_eta->Write();

		// Draw and write individual plots
		h_fakeInt_pt_0->Draw(); h_fakeInt_pt_0->Write();
		h_fakeInt_pt_1->Draw(); h_fakeInt_pt_1->Write();
		h_fakeInt_pt_2->Draw(); h_fakeInt_pt_2->Write();
		h_fakeInt_pt_3->Draw(); h_fakeInt_pt_3->Write();
		h_fakeInt_eta_0->Draw(); h_fakeInt_eta_0->Write();
		h_fakeInt_eta_1->Draw(); h_fakeInt_eta_1->Write();
		h_fakeInt_eta_2->Draw(); h_fakeInt_eta_2->Write();
		h_fakeInt_eta_3->Draw(); h_fakeInt_eta_3->Write();
	} //end if doFakeRate


	// tkMET turn-on curve
	if (doTurnOn) {
		cdturnon->cd();

		// sumw2
		h_trueMET->Sumw2(); h_trueTkMET->Sumw2();
		h_trueTkMET_turnon->Sumw2();
		// h_trueTkMET_tightZ_turnon->Sumw2();
		h_recoTkMET_1_turnon->Sumw2();
		h_recoTkMET_2_turnon->Sumw2();
		h_recoTkMET_3_turnon->Sumw2();
		h_recoTkMET_4_turnon->Sumw2();
		h_recoTkMET_5_turnon->Sumw2();
		h_recoTkMET_6_turnon->Sumw2();

		//Make trigger efficiency plots
		TH1F* h_turnon = (TH1F*) h_trueTkMET_turnon->Clone("turnon");
		h_turnon->Divide(h_trueTkMET_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_turnon, kBlack);

		// TH1F* h_turnon_tightZ = (TH1F*) h_trueTkMET_tightZ_turnon->Clone("turnon_tightZ");
		// h_turnon_tightZ->Divide(h_trueTkMET_tightZ_turnon, h_trueMET, 1.0, 1.0, "B");
		// makePrettyHisto(h_turnon_tightZ, kBlack);

		TH1F* h_1_turnon = (TH1F*) h_recoTkMET_1_turnon->Clone();
		h_1_turnon->Divide(h_recoTkMET_1_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_1_turnon, kRed);

		TH1F* h_2_turnon = (TH1F*) h_recoTkMET_2_turnon->Clone("reco_turnon_2");
		h_2_turnon->Divide(h_recoTkMET_2_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_2_turnon, kBlue);

		TH1F* h_3_turnon = (TH1F*) h_recoTkMET_3_turnon->Clone("reco_turnon_3");
		h_3_turnon->Divide(h_recoTkMET_3_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_3_turnon, kGreen+1);

		TH1F* h_4_turnon = (TH1F*) h_recoTkMET_4_turnon->Clone("reco_turnon_4");
		h_4_turnon->Divide(h_recoTkMET_4_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_4_turnon, kCyan);

		TH1F* h_5_turnon = (TH1F*) h_recoTkMET_5_turnon->Clone("reco_turnon_5");
		h_5_turnon->Divide(h_recoTkMET_5_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_5_turnon, kMagenta);

		TH1F* h_6_turnon = (TH1F*) h_recoTkMET_6_turnon->Clone("reco_turnon_6");
		h_6_turnon->Divide(h_recoTkMET_6_turnon, h_trueMET, 1.0, 1.0, "B");
		makePrettyHisto(h_6_turnon, kGray);


		// plot final tkMET turn-on curve
		gStyle->SetOptStat(0);
		TCanvas *c_turnon = new TCanvas("c_turnon","c_turnon", 600,600);
		h_turnon->SetStats(0); h_turnon->SetMinimum(0); h_turnon->SetMaximum(1.4);

		TLegend* leg_turnon = new TLegend(0.60,0.7,0.93,0.91);
	  leg_turnon->SetBorderSize(0);
		// leg_turnon->AddEntry(h_turnon,"trueTkMET, signal only","lp");
		//leg_turnon->AddEntry(h_turnon_tightZ,"trueTkMET, #Deltaz<0.5 cm","lp");
		// leg_turnon->AddEntry(h_1_turnon,"recoTkMET, 1","lp");
		// leg_turnon->AddEntry(h_2_turnon,"recoTkMET, 2","lp");
		// leg_turnon->AddEntry(h_3_turnon,"recoTkMET, 3","lp");
		leg_turnon->AddEntry(h_4_turnon,"trackMET, sim","lp");
		// leg_turnon->AddEntry(h_5_turnon,"recoTkMET, 5","lp");
		leg_turnon->AddEntry(h_6_turnon,"trackMET, emu","lp");

		h_turnon->SetTitle(";Gen-level MET [GeV];Efficiency");
		// h_turnon->Draw();
		//h_turnon_tightZ->Draw();
		// h_1_turnon->Draw("same");
		// h_3_turnon->Draw("same");
		// h_2_turnon->Draw("same");
		h_4_turnon->Draw();
		// h_5_turnon->Draw("same");
		h_6_turnon->Draw("same");
		leg_turnon->Draw("same");
		c_turnon->Write();


		float trueTkMET_95eff = h_turnon->GetBinLowEdge(h_turnon->FindFirstBinAbove(0.95));
		// float recoTkMET_95eff_1 = h_1_turnon->GetBinLowEdge(h_1_turnon->FindFirstBinAbove(0.95));
		// float recoTkMET_95eff_2 = h_2_turnon->GetBinLowEdge(h_2_turnon->FindFirstBinAbove(0.95));
		// float recoTkMET_95eff_3 = h_3_turnon->GetBinLowEdge(h_3_turnon->FindFirstBinAbove(0.95));
		float recoTkMET_95eff_4 = h_4_turnon->GetBinLowEdge(h_4_turnon->FindFirstBinAbove(0.95));
		// float recoTkMET_95eff_5 = h_5_turnon->GetBinLowEdge(h_5_turnon->FindFirstBinAbove(0.95));
		float recoTkMET_95eff_6 = h_6_turnon->GetBinLowEdge(h_6_turnon->FindFirstBinAbove(0.95));
		std::cout << "TrueTkMET: 95 point at " << trueTkMET_95eff << " GeV (offline)"<< std::endl;
		// std::cout << "RecoTkMET1: 95 point at " << recoTkMET_95eff_1 << " GeV (offline)"<< std::endl;
		// std::cout << "RecoTkMET2: 95 point at " << recoTkMET_95eff_2 << " GeV (offline)"<< std::endl;
		// std::cout << "RecoTkMET3: 95 point at " << recoTkMET_95eff_3 << " GeV (offline)"<< std::endl;
		std::cout << "RecoTkMET4: 95 point at " << recoTkMET_95eff_4 << " GeV (offline)"<< std::endl;
		// std::cout << "RecoTkMET5: 95 point at " << recoTkMET_95eff_5 << " GeV (offline)"<< std::endl;
		std::cout << "RecoTkMET6: 95 point at " << recoTkMET_95eff_6 << " GeV (offline)"<< std::endl;

		//Write trkMET plots
		h_trueMET->Draw(); 	h_trueMET->Write();
		h_trueTkMET_turnon->Draw(); h_trueTkMET_turnon->Write();
		// h_trueTkMET_tightZ_turnon->Draw(); h_trueTkMET_tightZ_turnon->Write();
		h_recoTkMET_1_turnon->Draw(); h_recoTkMET_1_turnon->Write();
		h_recoTkMET_2_turnon->Draw(); h_recoTkMET_2_turnon->Write();
		h_recoTkMET_3_turnon->Draw(); h_recoTkMET_3_turnon->Write();
		h_recoTkMET_4_turnon->Draw(); h_recoTkMET_4_turnon->Write();
		h_recoTkMET_5_turnon->Draw(); h_recoTkMET_5_turnon->Write();
		h_recoTkMET_6_turnon->Draw(); h_recoTkMET_6_turnon->Write();

		//Draw individual plots
		h_turnon->Draw(); h_turnon->Write();
		// h_turnon_tightZ->Draw(); h_turnon_tightZ->Write();
		h_1_turnon->Draw(); h_1_turnon->Write();
		h_2_turnon->Draw(); h_2_turnon->Write();
		h_3_turnon->Draw(); h_3_turnon->Write();
		h_4_turnon->Draw(); h_4_turnon->Write();
		h_5_turnon->Draw(); h_5_turnon->Write();
		h_6_turnon->Draw(); h_6_turnon->Write();
	} //end if doTurnOn


	// TrackMET resolution
	if (doTkMETRes) {
		gPad->SetGridx(0); gPad->SetGridy(0);
		cdres->cd();

		std::vector<TH1F*> tkMET1_histos = {h_recoTkMET_1_bin1,h_recoTkMET_1_bin2,h_recoTkMET_1_bin3,h_recoTkMET_1_bin4,h_recoTkMET_1_bin5,h_recoTkMET_1_bin6,h_recoTkMET_1_bin7,h_recoTkMET_1_bin8};
		std::vector<TH1F*> tkMET2_histos = {h_recoTkMET_2_bin1,h_recoTkMET_2_bin2,h_recoTkMET_2_bin3,h_recoTkMET_2_bin4,h_recoTkMET_2_bin5,h_recoTkMET_2_bin6,h_recoTkMET_2_bin7,h_recoTkMET_2_bin8};
		std::vector<TH1F*> tkMET3_histos = {h_recoTkMET_3_bin1,h_recoTkMET_3_bin2,h_recoTkMET_3_bin3,h_recoTkMET_3_bin4,h_recoTkMET_3_bin5,h_recoTkMET_3_bin6,h_recoTkMET_3_bin7,h_recoTkMET_3_bin8};
		std::vector<TH1F*> tkMET4_histos = {h_recoTkMET_4_bin1,h_recoTkMET_4_bin2,h_recoTkMET_4_bin3,h_recoTkMET_4_bin4,h_recoTkMET_4_bin5,h_recoTkMET_4_bin6,h_recoTkMET_4_bin7,h_recoTkMET_4_bin8};
		std::vector<TH1F*> tkMET5_histos = {h_recoTkMET_5_bin1,h_recoTkMET_5_bin2,h_recoTkMET_5_bin3,h_recoTkMET_5_bin4,h_recoTkMET_5_bin5,h_recoTkMET_5_bin6,h_recoTkMET_5_bin7,h_recoTkMET_5_bin8};
		std::vector<TH1F*> tkMET6_histos = {h_recoTkMET_6_bin1,h_recoTkMET_6_bin2,h_recoTkMET_6_bin3,h_recoTkMET_6_bin4,h_recoTkMET_6_bin5,h_recoTkMET_6_bin6,h_recoTkMET_6_bin7,h_recoTkMET_6_bin8};

		// Fitting delta(tkMET)
		std::vector< std::vector<float> >  tkMET_1 = resSigAndErr(tkMET1_histos);
		std::vector< std::vector<float> >  tkMET_2 = resSigAndErr(tkMET2_histos);
		std::vector< std::vector<float> >  tkMET_3 = resSigAndErr(tkMET3_histos);
		std::vector< std::vector<float> >  tkMET_4 = resSigAndErr(tkMET4_histos);
		std::vector< std::vector<float> >  tkMET_5 = resSigAndErr(tkMET5_histos);
		std::vector< std::vector<float> >  tkMET_6 = resSigAndErr(tkMET6_histos);

		// Final resolution plots
		TH1F *trkMETResolution_1 = new TH1F("trkMETResolution_1", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);
		TH1F *trkMETResolution_2 = new TH1F("trkMETResolution_2", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);
		TH1F *trkMETResolution_3 = new TH1F("trkMETResolution_3", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);
		TH1F *trkMETResolution_4 = new TH1F("trkMETResolution_4", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);
		TH1F *trkMETResolution_5 = new TH1F("trkMETResolution_5", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);
		TH1F *trkMETResolution_6 = new TH1F("trkMETResolution_6", ";MET (GeV);#sigma_{MET} (GeV)", binnum, bins);

		std::vector<float> trkMET_1_sig = tkMET_1[0]; std::vector<float> trkMET_1_err = tkMET_1[1];
		std::vector<float> trkMET_2_sig = tkMET_2[0]; std::vector<float> trkMET_2_err = tkMET_2[1];
		std::vector<float> trkMET_3_sig = tkMET_3[0]; std::vector<float> trkMET_3_err = tkMET_3[1];
		std::vector<float> trkMET_4_sig = tkMET_4[0]; std::vector<float> trkMET_4_err = tkMET_4[1];
		std::vector<float> trkMET_5_sig = tkMET_5[0]; std::vector<float> trkMET_5_err = tkMET_5[1];
		std::vector<float> trkMET_6_sig = tkMET_6[0]; std::vector<float> trkMET_6_err = tkMET_6[1];

		for(int i=0;i<=binnum;i++) {
			float avgMET = (bins[i]+bins[i+1])/2;
			//to show a tkMET spread
			//float avgMET = 1.0;
			trkMETResolution_1->SetBinContent(i+1,trkMET_1_sig[i]/avgMET); trkMETResolution_1->SetBinError(i+1,trkMET_1_err[i]/avgMET);
			trkMETResolution_2->SetBinContent(i+1,trkMET_2_sig[i]/avgMET); trkMETResolution_2->SetBinError(i+1,trkMET_2_err[i]/avgMET);
			trkMETResolution_3->SetBinContent(i+1,trkMET_3_sig[i]/avgMET); trkMETResolution_3->SetBinError(i+1,trkMET_3_err[i]/avgMET);
			trkMETResolution_4->SetBinContent(i+1,trkMET_4_sig[i]/avgMET); trkMETResolution_4->SetBinError(i+1,trkMET_4_err[i]/avgMET);
			trkMETResolution_5->SetBinContent(i+1,trkMET_5_sig[i]/avgMET); trkMETResolution_5->SetBinError(i+1,trkMET_5_err[i]/avgMET);
			trkMETResolution_6->SetBinContent(i+1,trkMET_6_sig[i]/avgMET); trkMETResolution_6->SetBinError(i+1,trkMET_6_err[i]/avgMET);
		}

		// Resolution final canvas
		TCanvas *cFinalRes = new TCanvas("cFinalRes","cFinalRes", 600,600);
		trkMETResolution_1->SetTitle("recoTkMET resolution; trueTkMET (GeV); #sigma_{TkMET} (GeV)");
		trkMETResolution_1->SetMaximum(3.0); trkMETResolution_1->SetMinimum(0);
		trkMETResolution_1->SetStats(0); makePrettyHisto(trkMETResolution_1, kRed);
		trkMETResolution_2->SetStats(0); makePrettyHisto(trkMETResolution_2, kBlue);
		trkMETResolution_3->SetStats(0); makePrettyHisto(trkMETResolution_3, kGreen+1);
		trkMETResolution_4->SetStats(0); makePrettyHisto(trkMETResolution_4, kCyan);
		trkMETResolution_5->SetStats(0); makePrettyHisto(trkMETResolution_5, kMagenta);
		trkMETResolution_6->SetStats(0); makePrettyHisto(trkMETResolution_6, kGray);

		trkMETResolution_1->Draw();
		trkMETResolution_2->Draw("SAME");
		trkMETResolution_3->Draw("SAME");
		trkMETResolution_4->Draw("SAME");
		trkMETResolution_5->Draw("SAME");
		trkMETResolution_6->Draw("SAME");

		TLegend* leg_res = new TLegend(0.60,0.7,0.93,0.91);
	  leg_res->SetBorderSize(0);
		leg_res->AddEntry(trkMETResolution_1,"Tracks, 1","l");
		leg_res->AddEntry(trkMETResolution_2,"Tracks, 2","l");
		leg_res->AddEntry(trkMETResolution_3,"Tracks, 3","l");
		leg_res->AddEntry(trkMETResolution_4,"Tracks, 4","l");
		leg_res->AddEntry(trkMETResolution_5,"Tracks, 5","l");
		leg_res->AddEntry(trkMETResolution_6,"Tracks, 6","l");

		leg_res->Draw("SAME");
		cFinalRes->Write();

		// Draw and write plots
		h_trueTkMET->Draw(); h_trueTkMET->Write();
		h_recoTkMET_1->Draw(); h_recoTkMET_1->Write();
		h_recoTkMET_2->Draw(); h_recoTkMET_2->Write();
		h_recoTkMET_3->Draw(); h_recoTkMET_3->Write();
		h_recoTkMET_4->Draw(); h_recoTkMET_4->Write();
		h_recoTkMET_5->Draw(); h_recoTkMET_5->Write();
		h_recoTkMET_6->Draw(); h_recoTkMET_6->Write();

		trkMETResolution_1->Write();
		trkMETResolution_2->Write();
		trkMETResolution_3->Write();
		trkMETResolution_4->Write();
		trkMETResolution_5->Write();
		trkMETResolution_6->Write();


		// Divided canvas with delta(tkMET) distributions
		TCanvas *c_div_recoTkMET_1RES = new TCanvas("c_div_recoTkMET_1RES","c_div_recoTkMET_1RES", 600,600);
		c_div_recoTkMET_1RES->Divide(4,2);
		for (int i=0; i<tkMET1_histos.size(); i++) {
			c_div_recoTkMET_1RES->cd(i+1); tkMET1_histos[i]->Draw();
		}
		c_div_recoTkMET_1RES->Write();

		TCanvas *c_div_recoTkMET_2RES = new TCanvas("c_div_recoTkMET_2RES","c_div_recoTkMET_2RES", 600,600);
		c_div_recoTkMET_2RES->Divide(4,2);
		for (int i=0; i<tkMET2_histos.size(); i++) {
			c_div_recoTkMET_2RES->cd(i+1); tkMET2_histos[i]->Draw();
		}
		// c_div_recoTkMET_2RES->cd(1); h_recoTkMET_2_bin1->Draw();
		// c_div_recoTkMET_2RES->cd(2); h_recoTkMET_2_bin2->Draw();
		// c_div_recoTkMET_2RES->cd(3); h_recoTkMET_2_bin3->Draw();
		// c_div_recoTkMET_2RES->cd(4); h_recoTkMET_2_bin4->Draw();
		// c_div_recoTkMET_2RES->cd(5); h_recoTkMET_2_bin5->Draw();
		// c_div_recoTkMET_2RES->cd(6); h_recoTkMET_2_bin6->Draw();
		// c_div_recoTkMET_2RES->cd(7); h_recoTkMET_2_bin7->Draw();
		// c_div_recoTkMET_2RES->cd(8); h_recoTkMET_2_bin8->Draw();
		c_div_recoTkMET_2RES->Write();

		TCanvas *c_div_recoTkMET_3RES = new TCanvas("c_div_recoTkMET_3RES","c_div_recoTkMET_3RES", 600,600);
		c_div_recoTkMET_3RES->Divide(4,2);
		for (int i=0; i<tkMET3_histos.size(); i++) {
			c_div_recoTkMET_3RES->cd(i+1); tkMET3_histos[i]->Draw();
		}
		c_div_recoTkMET_3RES->Write();

		TCanvas *c_div_recoTkMET_4RES = new TCanvas("c_div_recoTkMET_4RES","c_div_recoTkMET_4RES", 600,600);
		c_div_recoTkMET_4RES->Divide(4,2);
		for (int i=0; i<tkMET4_histos.size(); i++) {
			c_div_recoTkMET_4RES->cd(i+1); tkMET4_histos[i]->Draw();
		}
		c_div_recoTkMET_4RES->Write();

		TCanvas *c_div_recoTkMET_5RES = new TCanvas("c_div_recoTkMET_5RES","c_div_recoTkMET_5RES", 600,600);
		c_div_recoTkMET_5RES->Divide(4,2);
		for (int i=0; i<tkMET5_histos.size(); i++) {
			c_div_recoTkMET_5RES->cd(i+1); tkMET5_histos[i]->Draw();
		}
		c_div_recoTkMET_5RES->Write();

		TCanvas *c_div_recoTkMET_6RES = new TCanvas("c_div_recoTkMET_6RES","c_div_recoTkMET_6RES", 600,600);
		c_div_recoTkMET_6RES->Divide(4,2);
		for (int i=0; i<tkMET6_histos.size(); i++) {
			c_div_recoTkMET_6RES->cd(i+1); tkMET6_histos[i]->Draw();
		}
		c_div_recoTkMET_6RES->Write();
	} // end ifDoTkMETRes


	// ----------------------------------------------------------------------------------------------------------------
	// delta Z plots
	if (doDeltaZRes) {
		cdresZ->cd();


		//deltaZ resolution
		// std::vector<TH1F*> tp_deltaz_histos = {h_tp_deltaz_eta0to0p7, h_tp_deltaz_eta0p7to1, h_tp_deltaz_eta1to1p2, h_tp_deltaz_eta1p2to1p6, h_tp_deltaz_eta1p6to2, h_tp_deltaz_eta2to2p4};
		// std::vector<TH1F*> trk_deltaz_histos = {h_trk_deltaz_eta0to0p7, h_trk_deltaz_eta0p7to1, h_trk_deltaz_eta1to1p2, h_trk_deltaz_eta1p2to1p6, h_trk_deltaz_eta1p6to2, h_trk_deltaz_eta2to2p4};
		std::vector< std::vector<float> >  tpdz_fit = resSigAndErr(tp_dz_histos);
		std::vector< std::vector<float> >  trkdz_fit = resSigAndErr(trk_dz_histos);

		float z_bins[]={0,0.7,1.0,1.2,1.6,2.0,2.4};
		int z_binnum = 6;
		TH1F *z_tpMETResolution = new TH1F("z_tpMETResolution", ";track #eta;#sigma_{z} (cm)", z_binnum, z_bins);
		TH1F *z_trkMETResolution = new TH1F("z_trkMETResolution", ";track #eta;#sigma_{z} (cm)", z_binnum, z_bins);

		std::vector<float> tp_dz_sig = tpdz_fit[0]; std::vector<float> tp_dz_err = tpdz_fit[1];
		std::vector<float> trk_dz_sig = trkdz_fit[0]; std::vector<float> trk_dz_err = trkdz_fit[1];

		for(int i=0;i<binnum;i++) {
			z_tpMETResolution->SetBinContent(i+1,tp_dz_sig[i]); z_tpMETResolution->SetBinError(i+1,tp_dz_err[i]);
			z_trkMETResolution->SetBinContent(i+1,trk_dz_sig[i]); z_trkMETResolution->SetBinError(i+1,trk_dz_err[i]);
		}

		TCanvas *z_cFinalRes = new TCanvas("z_cFinalRes","z_cFinalRes", 600,600);
		z_tpMETResolution->SetTitle("; Track #eta; #sigma_{z} [cm]");
		z_tpMETResolution->SetMaximum(1.0); z_tpMETResolution->SetMinimum(0);
		z_tpMETResolution->SetStats(0); makePrettyHisto(z_tpMETResolution, kRed);
		z_trkMETResolution->SetStats(0); makePrettyHisto(z_trkMETResolution, kBlack);
		z_tpMETResolution->Draw();
		z_trkMETResolution->Draw("SAME");

		TLegend* leg_res_z = new TLegend(0.60,0.7,0.93,0.91);
	  leg_res_z->SetBorderSize(0);
		leg_res_z->AddEntry(z_tpMETResolution,"Tracking particles","lp");
		leg_res_z->AddEntry(z_trkMETResolution,"Matched tracks","lp");
		leg_res_z->Draw("SAME");
		z_cFinalRes->Write();

		//Draw and write individual plots
		z_tpMETResolution->Write();
		z_trkMETResolution->Write();

		h_trk_deltaz->Draw(); h_trk_deltaz->Write();
		h_tp_deltaz->Draw(); 	h_tp_deltaz->Write();

		TLegend* leg_deltaz = new TLegend(0.60,0.7,0.93,0.91);
		leg_deltaz->SetBorderSize(0);
		leg_deltaz->AddEntry(h_trk_deltaz,"All tracks","l");
		leg_deltaz->AddEntry(h_tp_deltaz,"signal only tracking particles","l");

		TCanvas *c_deltaz = new TCanvas("c_deltaz","c_deltaz", 600,600);
		c_deltaz->SetLogy();
		h_trk_deltaz->SetLineColor(kBlue);
		h_tp_deltaz->SetLineColor(kRed);
		h_trk_deltaz->Draw();
		h_tp_deltaz->Draw("same");
		leg_deltaz->Draw("same");
		c_deltaz->Write();

		TCanvas *c_div_TrkdeltazRes = new TCanvas("c_div_TrkdeltazRes","c_div_TrkdeltazRes", 600,600);
		c_div_TrkdeltazRes->Divide(3,2);
		for (int i=0; i<trk_dz_histos.size(); i++) {
			c_div_TrkdeltazRes->cd(i+1); trk_dz_histos[i]->Draw();
		}
		// c_div_TrkdeltazRes->cd(1); h_trk_deltaz_eta0to0p7->Draw();
		// c_div_TrkdeltazRes->cd(2); h_trk_deltaz_eta0p7to1->Draw();
		// c_div_TrkdeltazRes->cd(3); h_trk_deltaz_eta1to1p2->Draw();
		// c_div_TrkdeltazRes->cd(4); h_trk_deltaz_eta1p2to1p6->Draw();
		// c_div_TrkdeltazRes->cd(5); h_trk_deltaz_eta1p6to2->Draw();
		// c_div_TrkdeltazRes->cd(6); 	h_trk_deltaz_eta2to2p4->Draw();
		c_div_TrkdeltazRes->Write();

		TCanvas *c_div_TPdeltazRes = new TCanvas("c_div_TPdeltazRes","c_div_TPdeltazRes", 600,600);
		c_div_TPdeltazRes->Divide(3,2);
		for (int i=0; i<trk_dz_histos.size(); i++) {
			c_div_TPdeltazRes->cd(i+1); tp_dz_histos[i]->Draw();
		}
		c_div_TPdeltazRes->Write();
	} // if doDeltaZRes


	fout->Close();
}


//Change a float to a string, removing any trailing zeroes
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

void makePrettyHisto(TH1F* &h, Color_t color) {
	h->SetMarkerStyle(20); h->SetMarkerSize(0.7);
	h->SetMarkerColor(color); 	h->SetLineColor(color);
}

std::vector< std::vector<float> > resSigAndErr(std::vector<TH1F*> histos) {
	std::vector<float> sigmas(histos.size(), 0.0); std::vector<float> errors(histos.size(), 0.0);
	for (int i=0; i<histos.size(); i++) {
		if (histos[i]->GetEntries()!=0) {
			histos[i]->Fit("gaus","Q");
			TF1 *this_fit = histos[i]->GetFunction("gaus");
			sigmas[i] = this_fit->GetParameter(2);
			errors[i] = this_fit->GetParError(2);
		}
	}
	return {sigmas,errors};
}
