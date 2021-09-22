#include <cmath>
#include <fstream>
#include <iostream>

#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TRandom.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLegend.h"
#include "TF1.h"

#include "TVector3.h"
#include "TText.h"

#include "TGraph.h"
#include "TPad.h"
#include "TStyle.h"

#include "THStack.h"
#include "TMinuit.h"

#include <math.h>
#include <fstream>
bool DEBUG = false;
using namespace std;

// This script is for the pre-processing of the MCTruth Hits
// Usage:
// 1. Change the Input MCNTuple File Directory
// 2. Change the "surface" to define coverage
// 3. if blurIt == true, blur the hits with effective resolutions
// Change the line below:
// "CHANGE HERE TO ADD MORE SILICON"
// to add more silicon layers.

bool IsInSiliconArea(float x, float y, Int_t surface)
{
	// Definition of the coverage of the silicon inner tracker
	if (surface == 3)
	{ // Upgrade 1, Inner Tracker Silicon area, 3 3
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 3 and y > -100 * 3)))
			{ // Upgrade 1, Inner Tracker Silicon area, 3 3
				return true;
			}
		}
		return false;
	}
	if (surface == 4)
	{ // Upgrade 1, Inner Tracker Silicon area, 4 4
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 4 and y > -100 * 4)))
			{ // Upgrade 1, Inner Tracker Silicon area, 4 4
				return true;
			}
		}
		return false;
	}
	if (surface == 5)
	{ // Upgrade 1, Inner Tracker Silicon area, 5 5
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 5 and y > -100 * 5)))
			{ // Upgrade 1, Inner Tracker Silicon area, 5 5
				return true;
			}
		}
		return false;
	}
	if (surface == 6)
	{ // Upgrade 1, Inner Tracker Silicon area, 2 3 3 2
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 3 and y > -100 * 3)) ||
				((x < 540 * 2 and x > -540 * 2) and (y < 100 * 2 and y > -100 * 2)))
			{ // Upgrade 1, Inner Tracker Silicon area, 2 3 3 2
				return true;
			}
		}
		return false;
	}

	{
		cerr << "Error: Surface type not defined." << endl;
		return false;
	}
}

bool IsInBox(float x, float y)
{
	// Define the big box
	if ((x < -100 or x > 100) or (y < -100 or y > 100))
	{
		if (
			((x < 1500 and x > -1500) and (y < 600 and y > -600)))
		{ // Upgrade 1, Inner Tracker Silicon area, 3 3
			return true;
		}
		else
			return false;
	}
	else
	{
		return false;
	}
}

class hit
{
public:
	hit *previous_layer_hit;
	hit *next_layer_hit;
	float x, y, z;
	TVector3 vec;
	hit(float _x, float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
		vec = TVector3(x, y, z);
		previous_layer_hit = NULL;
		next_layer_hit = NULL;
	}

	void set_previous_layer(hit *_hit)
	{
		previous_layer_hit = _hit;
	}

	void set_next_layer(hit *_hit)
	{
		next_layer_hit = _hit;
	}

	int which_layer()
	{

		if (z > 7816 && z < 7836)
			return 1;
		if (z > 8024 && z < 8047)
			return 2;
		if (z > 8497 && z < 8518)
			return 3;
		if (z > 8706 && z < 8728)
			return 4;
		if (z > 9183 && z < 9203)
			return 5;
		if (z > 9391 && z < 9414)
			return 6;
		else
			return -1;
	}

private:
};

std::shared_ptr<hit> recoverhit(float posz, std::vector<std::shared_ptr<hit>> inputs)
{

	float besthit1y = 0;
	float besthit1x = 0;

	if (inputs.size() == 3)
	{
		float _z[] = {inputs[0]->z, inputs[1]->z, inputs[2]->z};
		float _y[] = {inputs[0]->y, inputs[1]->y, inputs[2]->y};
		float _x[] = {inputs[0]->x, inputs[1]->x, inputs[2]->x};

		TGraph gy(3, _z, _y);
		TGraph gx(3, _z, _x);
		TFitResultPtr fpx(0);
		TFitResultPtr fpy(0);
		fpx = gx.Fit("pol1", "Q");
		fpy = gy.Fit("pol1", "Q");

		TF1 *fx = gx.GetFunction("pol1");
		TF1 *fy = gy.GetFunction("pol1");
		besthit1y = fy->Eval(posz);
		besthit1x = fx->Eval(posz);
		posz = posz + besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1x = besthit1x + fx->GetParameter(1) * besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1y = besthit1y * (1 / (1 - 0.003601 * fy->GetParameter(1)));
	}

	if (inputs.size() == 2)
	{
		float _z[] = {inputs[0]->z, inputs[1]->z};
		float _y[] = {inputs[0]->y, inputs[1]->y};
		float _x[] = {inputs[0]->x, inputs[1]->x};

		TGraph gy(2, _z, _y);
		TGraph gx(2, _z, _x);
		TFitResultPtr fpx(0);
		TFitResultPtr fpy(0);
		fpx = gx.Fit("pol1", "Q");
		fpy = gy.Fit("pol1", "Q");

		TF1 *fx = gx.GetFunction("pol1");
		TF1 *fy = gy.GetFunction("pol1");
		besthit1y = fy->Eval(posz);
		besthit1x = fx->Eval(posz);
		posz = posz + besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1x = besthit1x + fx->GetParameter(1) * besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1y = besthit1y * (1 / (1 - 0.003601 * fy->GetParameter(1)));
	}

	if (inputs.size() == 1)
	{
		float _z[] = {0, inputs[0]->z};
		float _y[] = {0, inputs[0]->y};
		float _x[] = {0, inputs[0]->x};

		TGraph gy(2, _z, _y);
		TGraph gx(2, _z, _x);
		TFitResultPtr fpx(0);
		TFitResultPtr fpy(0);
		fpx = gx.Fit("pol1", "Q");
		fpy = gy.Fit("pol1", "Q");

		TF1 *fx = gx.GetFunction("pol1");
		TF1 *fy = gy.GetFunction("pol1");
		besthit1y = fy->Eval(posz);
		besthit1x = fx->Eval(posz);

		posz = posz + besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1x = besthit1x + fx->GetParameter(1) * besthit1y / (1 / 0.003601 - fy->GetParameter(1));
		besthit1y = besthit1y * (1 / (1 - 0.003601 * fy->GetParameter(1)));
	}

	if (inputs.size() > 3 or inputs.size() < 1)
		cout << "Recovering error!!! Inputs size > 3 or < 1" << endl;
	if (besthit1x == 0 and besthit1y == 0)
		cout << "Recovering error!!! ERROR !!! it is going to explode... run for your life" << endl;

	std::shared_ptr<hit> _hit(new hit(besthit1x, besthit1y, posz));
	return _hit;
}

void blurHits(Int_t surface = 3, TString inputFile = "../../standaloneanddeadareas/DataSets/original/MCtracks_U1b_small.root", TString outputFile = "outputFile.root", bool blurIt = true)
{
	gStyle->SetOptStat(1110);

	TFile f(inputFile);
	TDirectoryFile *dir = (TDirectoryFile *)gDirectory->Get("MCParticleNTuple");
	TTree *tracks = (TTree *)dir->Get("Tracks");

	int nentries;
	ULong64_t eventNumber;
	float hitzpos[30], hitypos[30], hitxpos[30];
	double p, pz;
	double qop;
	Double_t OriginVertexTime = -1;
	Double_t PrimaryVertexTime = -1;
	Int_t particleID = 0;

	tracks->SetBranchAddress("eventNumber", &eventNumber);
	tracks->SetBranchAddress("HitZpos", hitzpos);
	tracks->SetBranchAddress("HitXpos", hitxpos);
	tracks->SetBranchAddress("HitYpos", hitypos);
	tracks->SetBranchAddress("p", &p);
	tracks->SetBranchAddress("pz", &pz);
	tracks->SetBranchAddress("qop", &qop);
	tracks->SetBranchAddress("OriginVertexTime", &OriginVertexTime);
	tracks->SetBranchAddress("PrimaryVertexTime", &PrimaryVertexTime);
	TBranch *pid_branch_exist_test = (TBranch *)tracks->GetListOfBranches()->FindObject("particleID");
	if (pid_branch_exist_test)
		tracks->SetBranchAddress("particleID", &particleID);

	nentries = tracks->GetEntries();

	cout << "Total number of entries in Tree: " << nentries << endl;

	int previous_eventnumber = -1;

	int eventCounter = 0;

	TFile *hfile = TFile::Open(outputFile, "RECREATE");

	float hitPos[6][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	bool bHitReco[6] = {false, false, false, false, false, false};
	bool bPrompt = false;
	bool bSciFiHit[6] = {false, false, false, false, false, false};

	TTree *output_tree2 = new TTree("RTracks", "Tracks with recovered hits");
	output_tree2->Branch("eventNumber", &eventNumber, "eventNumber/l");
	output_tree2->Branch("hit1x", &hitPos[0][0], "hit1x/F");
	output_tree2->Branch("hit1y", &hitPos[0][1], "hit1y/F");
	output_tree2->Branch("hit1z", &hitPos[0][2], "hit1z/F");
	output_tree2->Branch("hit2x", &hitPos[1][0], "hit2x/F");
	output_tree2->Branch("hit2y", &hitPos[1][1], "hit2y/F");
	output_tree2->Branch("hit2z", &hitPos[1][2], "hit2z/F");
	output_tree2->Branch("hit3x", &hitPos[2][0], "hit3x/F");
	output_tree2->Branch("hit3y", &hitPos[2][1], "hit3y/F");
	output_tree2->Branch("hit3z", &hitPos[2][2], "hit3z/F");
	output_tree2->Branch("hit4x", &hitPos[3][0], "hit4x/F");
	output_tree2->Branch("hit4y", &hitPos[3][1], "hit4y/F");
	output_tree2->Branch("hit4z", &hitPos[3][2], "hit4z/F");
	output_tree2->Branch("hit5x", &hitPos[4][0], "hit5x/F");
	output_tree2->Branch("hit5y", &hitPos[4][1], "hit5y/F");
	output_tree2->Branch("hit5z", &hitPos[4][2], "hit5z/F");
	output_tree2->Branch("hit6x", &hitPos[5][0], "hit6x/F");
	output_tree2->Branch("hit6y", &hitPos[5][1], "hit6y/F");
	output_tree2->Branch("hit6z", &hitPos[5][2], "hit6z/F");
	output_tree2->Branch("p", &p, "p/D");
	output_tree2->Branch("pz", &pz, "pz/D");
	output_tree2->Branch("bhit1reco", &bHitReco[0], "bhit1reco/O");
	output_tree2->Branch("bhit2reco", &bHitReco[1], "bhit2reco/O");
	output_tree2->Branch("bhit3reco", &bHitReco[2], "bhit3reco/O");
	output_tree2->Branch("bhit4reco", &bHitReco[3], "bhit4reco/O");
	output_tree2->Branch("bhit5reco", &bHitReco[4], "bhit5reco/O");
	output_tree2->Branch("bhit6reco", &bHitReco[5], "bhit6reco/O");
	output_tree2->Branch("bprompt", &bPrompt, "bprompt/O");
	output_tree2->Branch("bscifihit1", &bSciFiHit[0], "bscifihit1/O");
	output_tree2->Branch("bscifihit2", &bSciFiHit[1], "bscifihit2/O");
	output_tree2->Branch("bscifihit3", &bSciFiHit[2], "bscifihit3/O");
	output_tree2->Branch("bscifihit4", &bSciFiHit[3], "bscifihit4/O");
	output_tree2->Branch("bscifihit5", &bSciFiHit[4], "bscifihit5/O");
	output_tree2->Branch("bscifihit6", &bSciFiHit[5], "bscifihit6/O");

	if (pid_branch_exist_test)
		output_tree2->Branch("particleID", &particleID, "particleID/I");

	TRandom rndgen = TRandom(time(0));

	std::vector<TH1F *> v_h1_hitx_recovered;
	std::vector<TH1F *> v_h1_hitx_recovered_and_mc;
	std::vector<TH1F *> v_h1_hitx_mc;
	std::vector<TH2F *> v_h2_hit_mc;
	std::vector<TH2F *> v_h2_hit_recovered;
	std::vector<TH2F *> v_h2_hit_recovered_prompt;
	std::vector<TH2F *> v_h2_hit_recovered_notprompt;
	std::vector<TH2F *> v_h2_hit_recovered_and_mc;

	std::vector<TH1F *> v_h1_hitx_predictionerror_p5e1;
	std::vector<TH1F *> v_h1_hity_predictionerror_p5e1;

	std::vector<TH1F *> v_h1_hitx_predictionerror_p20e1;
	std::vector<TH1F *> v_h1_hity_predictionerror_p20e1;

	std::vector<TH1F *> v_h1_hitx_predictionerror;
	std::vector<TH1F *> v_h1_hity_predictionerror;

	std::vector<TH1F *> v_h1_hitx_predictionerror_qm;
	std::vector<TH1F *> v_h1_hity_predictionerror_qm;

	std::vector<TH1F *> v_h1_hitx_predictionerror_qp;
	std::vector<TH1F *> v_h1_hity_predictionerror_qp;

	std::vector<TH2F *> v_h2_hitx_p_predictionerror;
	std::vector<TH2F *> v_h2_hity_p_predictionerror;

	std::vector<TH2F *> v_h2_hitx_p_predictionerror_qp;
	std::vector<TH2F *> v_h2_hity_p_predictionerror_qp;

	std::vector<TH2F *> v_h2_hitx_p_predictionerror_qm;
	std::vector<TH2F *> v_h2_hity_p_predictionerror_qm;

	std::vector<TH2F *> v_h2_hit_blured;

	for (unsigned int i = 1; i < 7; i++)
	{

		v_h1_hitx_predictionerror_p5e1.push_back(new TH1F(Form("h1_hit%dx_predictionerror_p5e1", i), Form("hit%d x prediction error p>5GeV in acc", i), 100, -1, 1));
		v_h1_hity_predictionerror_p5e1.push_back(new TH1F(Form("h1_hit%dy_predictionerror_p5e1", i), Form("hit%d y prediction error p>5GeV in acc", i), 100, -1, 1));

		v_h1_hitx_predictionerror_p20e1.push_back(new TH1F(Form("h1_hit%dx_predictionerror_p20e1", i), Form("hit%d x prediction error p>20GeV in acc", i), 100, -1, 1));
		v_h1_hity_predictionerror_p20e1.push_back(new TH1F(Form("h1_hit%dy_predictionerror_p20e1", i), Form("hit%d y prediction error p>20GeV in acc", i), 100, -1, 1));

		v_h1_hitx_predictionerror.push_back(new TH1F(Form("h1_hit%dx_predictionerror", i), Form("hit%d x prediction error", i), 100, -1, 1));
		v_h1_hity_predictionerror.push_back(new TH1F(Form("h1_hit%dy_predictionerror", i), Form("hit%d y prediction error", i), 100, -1, 1));

		v_h1_hitx_predictionerror_qm.push_back(new TH1F(Form("h1_hit%dx_predictionerror_qm", i), Form("hit%d x prediction error (q<0)", i), 100, -1, 1));
		v_h1_hity_predictionerror_qm.push_back(new TH1F(Form("h1_hit%dy_predictionerror_qm", i), Form("hit%d y prediction error (q<0)", i), 100, -0.6, 0.6));
		v_h1_hitx_predictionerror_qp.push_back(new TH1F(Form("h1_hit%dx_predictionerror_qp", i), Form("hit%d x prediction error (q>0)", i), 100, -1, 1));
		v_h1_hity_predictionerror_qp.push_back(new TH1F(Form("h1_hit%dy_predictionerror_qp", i), Form("hit%d y prediction error (q>0)", i), 100, -0.6, 0.6));

		v_h2_hitx_p_predictionerror.push_back(new TH2F(Form("h2_hit%dx_p_predictionerror", i), Form("hit%d x prediction error by p", i), 20, 0, 20000, 100, -5, 5));
		v_h2_hity_p_predictionerror.push_back(new TH2F(Form("h2_hit%dy_p_predictionerror", i), Form("hit%d y prediction error by p", i), 20, 0, 20000, 120, -0.6, 0.6));
		v_h1_hitx_recovered.push_back(new TH1F(Form("h1_hit%dx_recovered", i), Form("hit%d x recovered", i), 4 * 540, -4 * 540, 4 * 540));
		v_h1_hitx_recovered_and_mc.push_back(new TH1F(Form("h1_hit%dx_recovered_and_mc", i), Form("hit%d x recovered and mc", i), 4 * 1080, -4 * 540, 4 * 540));
		v_h1_hitx_mc.push_back(new TH1F(Form("h1_hit%dx_mc", i), Form("hit%d x mc", i), 4 * 1080, -4 * 540, 4 * 540));
		v_h2_hit_mc.push_back(new TH2F(Form("h2_hit%d_mc", i), Form("hit%d mc", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));
		v_h2_hit_recovered.push_back(new TH2F(Form("h2_hit%d_recovered", i), Form("hit%d recovered", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));
		v_h2_hit_recovered_prompt.push_back(new TH2F(Form("h2_hit%d_recovered_prompt", i), Form("hit%d recovered prompt", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));
		v_h2_hit_recovered_notprompt.push_back(new TH2F(Form("h2_hit%d_recovered_notprompt", i), Form("hit%d recovered not prompt", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));
		v_h2_hit_recovered_and_mc.push_back(new TH2F(Form("h2_hit%d_recovered_and_mc", i), Form("hit%d recovered and mc", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));

		v_h2_hit_blured.push_back(new TH2F(Form("h2_hit %d_blured", i), Form("h2_hit %d_blured", i), 4 * 1080, -4 * 540, 4 * 540, 1400, -700, 700));

		v_h1_hitx_predictionerror[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h1_hity_predictionerror[i - 1]->GetXaxis()->SetTitle("y [mm]");
		v_h2_hitx_p_predictionerror[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hitx_p_predictionerror[i - 1]->GetYaxis()->SetTitle("x [mm]");
		v_h2_hity_p_predictionerror[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hity_p_predictionerror[i - 1]->GetYaxis()->SetTitle("y [mm]");
		v_h1_hitx_recovered[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h1_hitx_recovered_and_mc[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h1_hitx_mc[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_mc[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_mc[i - 1]->GetYaxis()->SetTitle("y [mm]");
		v_h2_hit_recovered[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_recovered[i - 1]->GetYaxis()->SetTitle("y [mm]");
		v_h2_hit_recovered_prompt[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_recovered_prompt[i - 1]->GetYaxis()->SetTitle("y [mm]");
		v_h2_hit_recovered_notprompt[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_recovered_notprompt[i - 1]->GetYaxis()->SetTitle("y [mm]");
		v_h2_hit_recovered_and_mc[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_recovered_and_mc[i - 1]->GetYaxis()->SetTitle("y [mm]");

		v_h2_hit_blured[i - 1]->GetXaxis()->SetTitle("x [mm]");
		v_h2_hit_blured[i - 1]->GetYaxis()->SetTitle("y [mm]");

		v_h2_hitx_p_predictionerror_qp.push_back(new TH2F(Form("h2_hit%dx_p_predictionerror_qp", i), Form("hit%d x prediction error by p (q>0)", i), 20, 0, 20000, 100, -5, 5));
		v_h2_hity_p_predictionerror_qp.push_back(new TH2F(Form("h2_hit%dy_p_predictionerror_qp", i), Form("hit%d y prediction error by p (q>0)", i), 20, 0, 20000, 120, -0.6, 0.6));
		v_h2_hitx_p_predictionerror_qm.push_back(new TH2F(Form("h2_hit%dx_p_predictionerror_qm", i), Form("hit%d x prediction error by p (q<0)", i), 20, 0, 20000, 100, -5, 5));
		v_h2_hity_p_predictionerror_qm.push_back(new TH2F(Form("h2_hit%dy_p_predictionerror_qm", i), Form("hit%d y prediction error by p (q<0)", i), 20, 0, 20000, 120, -0.6, 0.6));

		v_h2_hitx_p_predictionerror_qp[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hitx_p_predictionerror_qp[i - 1]->GetYaxis()->SetTitle("#Delta x [mm]");
		v_h2_hity_p_predictionerror_qp[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hity_p_predictionerror_qp[i - 1]->GetYaxis()->SetTitle("#Delta y [mm]");

		v_h2_hitx_p_predictionerror_qm[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hitx_p_predictionerror_qm[i - 1]->GetYaxis()->SetTitle("#Delta x [mm]");
		v_h2_hity_p_predictionerror_qm[i - 1]->GetXaxis()->SetTitle("p [MeV/c]");
		v_h2_hity_p_predictionerror_qm[i - 1]->GetYaxis()->SetTitle("#Delta y [mm]");
	}

	for (int i = 0; i < nentries; i++)
	{
		for (int iLayer = 0; iLayer < 6; iLayer++)
		{
			hitPos[iLayer][0] = 0;
			hitPos[iLayer][1] = 0;
			hitPos[iLayer][2] = 0;
			bSciFiHit[iLayer] = 0;
			bHitReco[iLayer] = 0;
		}

		if (previous_eventnumber == -1)
			previous_eventnumber = eventNumber;

		if (i % 100000 == 0)
			cout << "Progress: " << i * 100. / nentries << endl;

		if (DEBUG)
			cout << "Entry number: " << i << endl;
		tracks->GetEntry(i);

		bPrompt = (OriginVertexTime == PrimaryVertexTime);
		float eta = 1. / 2 * log((p + pz) / (p - pz)); // log is in e base in C++

		std::shared_ptr<hit> hit1 = NULL;
		std::shared_ptr<hit> hit12u = NULL;
		std::shared_ptr<hit> hit12v = NULL;
		std::shared_ptr<hit> hit2 = NULL;

		std::shared_ptr<hit> hit3 = NULL;
		std::shared_ptr<hit> hit34u = NULL;
		std::shared_ptr<hit> hit34v = NULL;
		std::shared_ptr<hit> hit4 = NULL;

		std::shared_ptr<hit> hit5 = NULL;
		std::shared_ptr<hit> hit56u = NULL;
		std::shared_ptr<hit> hit56v = NULL;
		std::shared_ptr<hit> hit6 = NULL;

		std::vector<std::shared_ptr<hit>> v_hits;

		for (unsigned int j = 0; j < 12; j++)
			v_hits.push_back(NULL);

		if (DEBUG)
			cout << "v_hists filled with NULL" << endl;

		for (int j = 29; j >= 0; j--)
		{

			if (DEBUG)
				cout << "Creating temporary hit: " << hitxpos[j] << " " << hitypos[j] << " " << hitzpos[j] << endl;
			std::shared_ptr<hit> _hit(new hit(hitxpos[j], hitypos[j], hitzpos[j]));

			if (DEBUG)
				cout << "Checking if hit is in the big box" << endl;

			if (not IsInBox(_hit->x, _hit->y))
				continue;

			if (_hit->z > 7816 && _hit->z < 7836)
				bSciFiHit[0] = true;
			if (_hit->z > 8024 && _hit->z < 8047)
				bSciFiHit[1] = true;
			if (_hit->z > 8497 && _hit->z < 8518)
				bSciFiHit[2] = true;
			if (_hit->z > 8706 && _hit->z < 8728)
				bSciFiHit[3] = true;
			if (_hit->z > 9183 && _hit->z < 9203)
				bSciFiHit[4] = true;
			if (_hit->z > 9391 && _hit->z < 9414)
				bSciFiHit[5] = true;

			if (DEBUG)
				cout << "Hit in the box area" << endl;

			if (_hit->z > 7816 && _hit->z < 7836)
			{
				hit1 = _hit;
				v_hits[0] = _hit;
				hitPos[0][0] = hit1->x;
				hitPos[0][1] = hit1->y;
				hitPos[0][2] = hit1->z;
				if (DEBUG)
					cout << "Found hit1" << endl;
			}
			if (_hit->z > 7880 && _hit->z < 7910)
			{
				hit12u = _hit;

				v_hits[1] = _hit;
				if (DEBUG)
					cout << "Found hit12u" << endl;
			}
			if (_hit->z > 7950 && _hit->z < 7980)
			{
				hit12v = _hit;

				v_hits[2] = _hit;
				if (DEBUG)
					cout << "Found hit12v" << endl;
			}
			if (_hit->z > 8024 && _hit->z < 8047)
			{
				hit2 = _hit;

				v_hits[3] = _hit;
				hitPos[1][0] = hit2->x;
				hitPos[1][1] = hit2->y;
				hitPos[1][2] = hit2->z;
				if (DEBUG)
					cout << "Found hit2" << endl;
			}
			if (_hit->z > 8497 && _hit->z < 8518)
			{
				hit3 = _hit;

				v_hits[4] = _hit;
				hitPos[2][0] = hit3->x;
				hitPos[2][1] = hit3->y;
				hitPos[2][2] = hit3->z;
				if (DEBUG)
					cout << "Found hit3" << endl;
			}
			if (_hit->z > 8560 && _hit->z < 8600)
			{
				hit34u = _hit;

				v_hits[5] = _hit;
				if (DEBUG)
					cout << "Found hit34u" << endl;
			}
			if (_hit->z > 8630 && _hit->z < 8670)
			{
				hit34v = _hit;

				v_hits[6] = _hit;
				if (DEBUG)
					cout << "Found hit34v" << endl;
			}
			if (_hit->z > 8706 && _hit->z < 8728)
			{
				hit4 = _hit;

				v_hits[7] = _hit;
				hitPos[3][0] = hit4->x;
				hitPos[3][1] = hit4->y;
				hitPos[3][2] = hit4->z;
				if (DEBUG)
					cout << "Found hit4" << endl;
			}
			if (_hit->z > 9183 && _hit->z < 9203)
			{
				hit5 = _hit;

				v_hits[8] = _hit;
				hitPos[4][0] = hit5->x;
				hitPos[4][1] = hit5->y;
				hitPos[4][2] = hit5->z;
				if (DEBUG)
					cout << "Found hit5" << endl;
			}
			if (_hit->z > 9240 && _hit->z < 9280)
			{
				hit56u = _hit;

				v_hits[9] = _hit;
				if (DEBUG)
					cout << "Found hit56u" << endl;
			}
			if (_hit->z > 9310 && _hit->z < 9350)
			{
				hit56v = _hit;

				v_hits[10] = _hit;
				if (DEBUG)
					cout << "Found hit56v" << endl;
			}
			if (_hit->z > 9391 && _hit->z < 9414)
			{
				hit6 = _hit;

				v_hits[11] = _hit;
				hitPos[5][0] = hit6->x;
				hitPos[5][1] = hit6->y;
				hitPos[5][2] = hit6->z;
				if (DEBUG)
					cout << "Found hit6" << endl;
			}
		}

		if (DEBUG)
			cout << "Starting Hits recover algorithm" << endl;

		for (unsigned int k = 1; k < 7; k++)
		{
			// loop over the k MT layers 1-6

			std::vector<std::shared_ptr<hit>> inputs;
			int z = ((k) / 2) * 4;
			if (k % 2 == 0)
				z -= 1;
			//cout << "z: " << z << endl;
			if (k % 2 != 0)
			{
				if (v_hits[z + 1] and v_hits[z + 2] and v_hits[z + 3])
				{
					inputs.push_back(v_hits[z + 1]);
					inputs.push_back(v_hits[z + 2]);
					inputs.push_back(v_hits[z + 3]);
				}
				else if (v_hits[z + 1] and v_hits[z + 2])
				{
					inputs.push_back(v_hits[z + 1]);
					inputs.push_back(v_hits[z + 2]);
				}
				else if (v_hits[z + 1] and v_hits[z + 3])
				{
					inputs.push_back(v_hits[z + 1]);
					inputs.push_back(v_hits[z + 3]);
				}
				else if (v_hits[z + 2] and v_hits[z + 3])
				{
					inputs.push_back(v_hits[z + 2]);
					inputs.push_back(v_hits[z + 3]);
				}

				//cout << "z+3: " << z+3 << endl;
			}
			else
			{
				if (v_hits[z - 1] and v_hits[z - 2] and v_hits[z - 3])
				{
					inputs.push_back(v_hits[z - 1]);
					inputs.push_back(v_hits[z - 2]);
					inputs.push_back(v_hits[z - 3]);
				}
				else if (v_hits[z - 1] and v_hits[z - 2])
				{
					inputs.push_back(v_hits[z - 1]);
					inputs.push_back(v_hits[z - 2]);
				}
				else if (v_hits[z - 1] and v_hits[z - 3])
				{
					inputs.push_back(v_hits[z - 1]);
					inputs.push_back(v_hits[z - 3]);
				}
				else if (v_hits[z - 2] and v_hits[z - 3])
				{
					inputs.push_back(v_hits[z - 2]);
					inputs.push_back(v_hits[z - 3]);
				}

				//cout << "z-3: " << z-3 << endl;
			}

			if (inputs.size() > 0)
			{
				std::shared_ptr<hit> _hit(NULL);

				if (k == 1)
					_hit = recoverhit(7825.9, inputs);
				if (k == 2)
					_hit = recoverhit(8035.8, inputs);
				if (k == 3)
					_hit = recoverhit(8507.8, inputs);
				if (k == 4)
					_hit = recoverhit(8717.9, inputs);
				if (k == 5)
					_hit = recoverhit(9192.9, inputs);
				if (k == 6)
					_hit = recoverhit(9402.8, inputs);

				std::shared_ptr<hit> _mhit = NULL;

				if (v_hits[z])
					_mhit = v_hits[z];

				//cout << "mhit assigned" << endl;

				if (_mhit)
				{
					if (IsInSiliconArea(_hit->x, _hit->y, surface) && IsInSiliconArea(_mhit->x, _mhit->y, surface))
					{
						if (DEBUG)
							cout << "mhit in silicon area" << endl;
						//cout << "Processing info for " << k << " layer" << endl;
						v_h1_hitx_predictionerror[k - 1]->Fill((_mhit->x - _hit->x));
						v_h1_hity_predictionerror[k - 1]->Fill((_mhit->y - _hit->y));

						if (p > 5000. and eta > 2 and eta < 5)
							v_h1_hitx_predictionerror_p5e1[k - 1]->Fill((_mhit->x - _hit->x));
						if (p > 5000. and eta > 2 and eta < 5)
							v_h1_hity_predictionerror_p5e1[k - 1]->Fill((_mhit->y - _hit->y));

						if (p > 20000. and eta > 2 and eta < 5)
							v_h1_hitx_predictionerror_p20e1[k - 1]->Fill((_mhit->x - _hit->x));
						if (p > 20000. and eta > 2 and eta < 5)
							v_h1_hity_predictionerror_p20e1[k - 1]->Fill((_mhit->y - _hit->y));

						v_h2_hitx_p_predictionerror[k - 1]->Fill(p, (_mhit->x - _hit->x));
						v_h2_hity_p_predictionerror[k - 1]->Fill(p, (_mhit->y - _hit->y));

						v_h2_hit_mc[k - 1]->Fill(_mhit->x, _mhit->y);
						v_h2_hit_recovered_and_mc[k - 1]->Fill(_mhit->x, _mhit->y);
						v_h1_hitx_recovered_and_mc[k - 1]->Fill(_mhit->x);
						v_h1_hitx_mc[k - 1]->Fill(_hit->x);
						//cout << "Filled mc-recovered hit info" << endl;

						if (qop > 0)
						{
							v_h2_hitx_p_predictionerror_qp[k - 1]->Fill(p, (_mhit->x - _hit->x));
							v_h2_hity_p_predictionerror_qp[k - 1]->Fill(p, (_mhit->y - _hit->y));
							v_h1_hity_predictionerror_qp[k - 1]->Fill(_mhit->y - _hit->y);
							v_h1_hitx_predictionerror_qp[k - 1]->Fill(_mhit->x - _hit->x);
						}
						if (qop < 0)
						{
							v_h2_hitx_p_predictionerror_qm[k - 1]->Fill(p, (_mhit->x - _hit->x));
							v_h2_hity_p_predictionerror_qm[k - 1]->Fill(p, (_mhit->y - _hit->y));
							v_h1_hity_predictionerror_qm[k - 1]->Fill(_mhit->y - _hit->y);
							v_h1_hitx_predictionerror_qm[k - 1]->Fill(_mhit->x - _hit->x);
						}
					}
					else if (IsInBox(_mhit->x, _mhit->y))
					{
						if (DEBUG)
							cout << "mhit not in silicon area, in box" << endl;

						v_h2_hit_mc[k - 1]->Fill(_mhit->x, _mhit->y);
						v_h2_hit_recovered_and_mc[k - 1]->Fill(_mhit->x, _mhit->y);
						v_h1_hitx_recovered_and_mc[k - 1]->Fill(_mhit->x);
						v_h1_hitx_mc[k - 1]->Fill(_hit->x);
					}
					else if (not IsInBox(_mhit->x, _mhit->y))
					{
						cout << "WTF" << endl;
					}
				}
				else
				{
					// Only recover the hits in silicon area and proper layer
					// CHANGE HERE TO ADD MORE SILICON
					if (IsInSiliconArea(_hit->x, _hit->y, surface) and (k == 3 or k == 4))
					{
						if (DEBUG)
							cout << "recovered hit in silicon area" << endl;
						v_h2_hit_recovered_and_mc[k - 1]->Fill(_hit->x, _hit->y);
						v_h2_hit_recovered[k - 1]->Fill(_hit->x, _hit->y);

						if (OriginVertexTime == PrimaryVertexTime)
							v_h2_hit_recovered_prompt[k - 1]->Fill(_hit->x, _hit->y);
						else
							v_h2_hit_recovered_notprompt[k - 1]->Fill(_hit->x, _hit->y);

						v_h1_hitx_recovered_and_mc[k - 1]->Fill(_hit->x);
						v_h1_hitx_recovered[k - 1]->Fill(_hit->x);

						//flag hit on scifi after the recovery algorithm as well

						if (_hit->z > 7816 && _hit->z < 7836)
							bSciFiHit[0] = true;
						if (_hit->z > 8024 && _hit->z < 8047)
							bSciFiHit[1] = true;
						if (_hit->z > 8497 && _hit->z < 8518)
							bSciFiHit[2] = true;
						if (_hit->z > 8706 && _hit->z < 8728)
							bSciFiHit[3] = true;
						if (_hit->z > 9183 && _hit->z < 9203)
							bSciFiHit[4] = true;
						if (_hit->z > 9391 && _hit->z < 9414)
							bSciFiHit[5] = true;
						//hit1 = _hit;
						hitPos[k - 1][0] = _hit->x;
						hitPos[k - 1][1] = _hit->y;
						hitPos[k - 1][2] = _hit->z;
						bHitReco[k - 1] = true;
					}
				}
			}
		}

		for (int iLayer = 0; iLayer < 6; iLayer++)
		{
			if (bSciFiHit[iLayer])
			{
				if (blurIt)
				{
					// Layer 2 and 3 are Hybrid Layers
					// CHANGE HERE TO ADD MORE SILICON
					if ((iLayer == 2 or iLayer == 3) and IsInSiliconArea(hitPos[iLayer][0], hitPos[iLayer][1], surface))
					{
						// Silicon Pixel Resolution: x50um*y150um /sqrt(12). According to FTDR
						hitPos[iLayer][0] = hitPos[iLayer][0] + rndgen.Gaus(0.0, 0.050 / 3.4641016);
						hitPos[iLayer][1] = hitPos[iLayer][1] + rndgen.Gaus(0.0, 0.150 / 3.4641016);
					}
					else
					{
						// SciFi Resolution: 70um in x. Hence in y is 1134um. According to SciFi Resolution in LHCb-PUB-2015-008.
						hitPos[iLayer][0] = hitPos[iLayer][0] + rndgen.Gaus(0.0, 0.070);
						hitPos[iLayer][1] = hitPos[iLayer][1] + rndgen.Gaus(0.0, 1.134);
					}
					// Bigger resolution Test
					// hitPos[iLayer][0] = hitPos[iLayer][0] + rndgen.Gaus(0.0, 0.25);
					// hitPos[iLayer][1] = hitPos[iLayer][1] + rndgen.Gaus(0.0, 4);
				}
				v_h2_hit_blured[iLayer]->Fill(hitPos[iLayer][0], hitPos[iLayer][1]);
			}
		}
		output_tree2->Fill();
		v_hits.clear();
	}

	TCanvas *c2 = new TCanvas("c2", "c2", 1800, 600);
	TCanvas *c1 = new TCanvas();
	c1->cd();

	for (unsigned int i = 1; i < 7; i++)
	{
		c1->cd();
		gStyle->SetOptStat(1110);

		v_h2_hitx_p_predictionerror_qm[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dx_p_predictionerror_qm.png", i));
		v_h2_hity_p_predictionerror_qm[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dy_p_predictionerror_qm.png", i));

		v_h2_hitx_p_predictionerror_qp[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dx_p_predictionerror_qp.png", i));
		v_h2_hity_p_predictionerror_qp[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dy_p_predictionerror_qp.png", i));

		v_h1_hitx_predictionerror[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_predictionerror.png", i));
		v_h1_hity_predictionerror[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dy_predictionerror.png", i));

		v_h1_hity_predictionerror_qm[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dy_predictionerror_qm.png", i));
		v_h1_hity_predictionerror_qp[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dy_predictionerror_qp.png", i));
		v_h1_hitx_predictionerror_qm[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_predictionerror_qm.png", i));
		v_h1_hitx_predictionerror_qp[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_predictionerror_qp.png", i));

		v_h1_hitx_predictionerror_p5e1[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_predictionerror_p5e1.png", i));
		v_h1_hity_predictionerror_p5e1[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dy_predictionerror_p5e1.png", i));

		v_h1_hitx_predictionerror_p20e1[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_predictionerror_p20e1.png", i));
		v_h1_hity_predictionerror_p20e1[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dy_predictionerror_p20e1.png", i));

		v_h2_hitx_p_predictionerror[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dx_p_predictionerror.png", i));
		v_h2_hity_p_predictionerror[i - 1]->Draw("COLZ");
		c1->SaveAs(Form("h2_hit%dy_p_predictionerror.png", i));
		v_h1_hitx_recovered[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_recovered.png", i));
		v_h1_hitx_recovered_and_mc[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_recovered_and_mc.png", i));
		v_h1_hitx_mc[i - 1]->Draw();
		c1->SaveAs(Form("h1_hit%dx_mc.png", i));

		gStyle->SetOptStat(0);
		c2->cd();
		v_h2_hit_mc[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("h2_hit%d_mc.png", i));
		v_h2_hit_recovered[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("h2_hit%d_recovered.png", i));
		v_h2_hit_recovered_and_mc[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("h2_hit%d_recovered_and_mc.png", i));

		v_h2_hit_blured[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("v_h2_hit%d_blured.png", i));

		v_h2_hit_recovered_prompt[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("h2_hit%d_recovered_prompt.png", i));
		v_h2_hit_recovered_notprompt[i - 1]->Draw("COLZ");
		c2->SaveAs(Form("h2_hit%d_recovered_notprompt.png", i));
	}

	hfile->Write();
}
