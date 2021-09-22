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

// #define debug 1
using namespace std;

bool IsInSiliconArea(float x, float y, Int_t surface)
{
	if (surface == 0)
	{ // Upgrade 2, Full Mighty Tracker Silicon area, 2 3 4 5 5 4 3 2
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 5 and y > -100 * 5)) ||
				((x < 540 * 2 and x > -540 * 2) and (y < 100 * 4 and y > -100 * 4)) ||
				((x < 540 * 3 and x > -540 * 3) and (y < 100 * 3 and y > -100 * 3)) ||
				((x < 540 * 4 and x > -540 * 4) and (y < 100 * 2 and y > -100 * 2)))
			{ // Upgrade 2, Full Mighty Tracker Silicon area, 2 3 4 5 5 4 3 2
				return true;
			}
		}
		return false;
	}
	if (surface == 1)
	{ // Upgrade 2, Modest Mighty Tracker Silicon area, 1 2 3 4 4 3 2 1
		if ((x < -100 or x > 100) or (y < -100 or y > 100))
		{
			if (
				((x < 540 * 1 and x > -540 * 1) and (y < 100 * 4 and y > -100 * 4)) ||
				((x < 540 * 2 and x > -540 * 2) and (y < 100 * 3 and y > -100 * 3)) ||
				((x < 540 * 3 and x > -540 * 3) and (y < 100 * 2 and y > -100 * 2)) ||
				((x < 540 * 4 and x > -540 * 4) and (y < 100 * 1 and y > -100 * 1)))
			{ // Upgrade 2, Modest Mighty Tracker Silicon area, 1 2 3 4 4 3 2 1
				// cout << "True" << endl;
				return true;
			}
		}
		return false;
	}
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
	{ // Upgrade 1, Inner Tracker Silicon area, 1 3 3 1
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
		// if (
		// 	((x < 1500 and x > -1500) and (y < 600 and y > -600)))
		// {
		// All the hits not in the beam region are allowed.
		return true;
		// }
		// else
		// 	return false;
	}
	else
	{
		return false;
	}
}

bool uUpOrDown(float x, float y)
{
	// Decide Upper/Lower Plane of the u layer
	if ((y + x * TMath::Tan(5.0 / 180)) > 0)
	{
		// Up
		return true;
	}
	else
	{
		// Down
		return false;
	}
}

bool vUpOrDown(float x, float y)
{
	// Decide Upper/Lower Plane of the v layer
	if ((y - x * TMath::Tan(5.0 / 180)) > 0)
	{
		// Up
		return true;
	}
	else
	{
		// Down
		return false;
	}
}

class hit
{
public:
	hit *previous_layer_hit = NULL;
	hit *next_layer_hit = NULL;
	float x, y, z;
	hit(float _x, float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
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

		// Tilt Modification
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

		// Tilt Modification
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

		// Tilt Modification
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
// void samplePrep(Int_t layers = 2, Int_t surface = 3, TString inputFile = "../DataSets/original/MCtracks_U1b_small.root", TString outputFile = "2LayerSmall.root", bool blurIt = true)

void samplePrep(Int_t layers = 2, Int_t surface = 3, TString inputFile = "../DataSets/original/MCtracks_U1b-5000ev.root", TString outputFile = "3LayerBig.root", bool blurIt = true)
{
	gStyle->SetOptStat(1110);

	TFile f(inputFile);
	TDirectoryFile *dir = (TDirectoryFile *)gDirectory->Get("MCParticleNTuple");
	TTree *tracks = (TTree *)dir->Get("Tracks");

	// Input Variables
	int nentries;
	ULong64_t eventNumber = 0;
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

	TFile *hfile = TFile::Open(outputFile, "RECREATE");

	// Output Variables
	// T1 x-u-v-x hits position
	float hitT1up[4] = {0, 0, 0, 0};
	float hitT1down[4] = {0, 0, 0, 0};
	// T2 Silicon Hits position
	float hitT2[2][3] = {{0, 0, 0}, {0, 0, 0}};
	// T3 x-u-v-x hits position
	float hitT3up[4] = {0, 0, 0, 0};
	float hitT3down[4] = {0, 0, 0, 0};

	double eta = 0;

	bool bT1HitUp[4] = {false, false, false, false};
	bool bT1HitDown[4] = {false, false, false, false};
	bool bT2Hit[2] = {false, false};
	bool bT3HitUp[4] = {false, false, false, false};
	bool bT3HitDown[4] = {false, false, false, false};

	TTree *output_tree2 = new TTree("RTracks", "Tracks with recovered hits");
	output_tree2->Branch("eventNumber", &eventNumber, "eventNumber/l");

	output_tree2->Branch("t1x1u", &hitT1up[0], "t1x1u/F");
	output_tree2->Branch("t1x1d", &hitT1down[0], "t1x1d/F");
	output_tree2->Branch("t1uu", &hitT1up[1], "t1uu/F");
	output_tree2->Branch("t1ud", &hitT1down[1], "t1ud/F");
	output_tree2->Branch("t1vu", &hitT1up[2], "t1vu/F");
	output_tree2->Branch("t1vd", &hitT1down[2], "t1vd/F");
	output_tree2->Branch("t1x2u", &hitT1up[3], "t1x2u/F");
	output_tree2->Branch("t1x2d", &hitT1down[3], "t1x2d/F");

	output_tree2->Branch("t2x1", &hitT2[0][0], "t2x1/F");
	output_tree2->Branch("t2y1", &hitT2[0][1], "t2y1/F");
	output_tree2->Branch("t2z1", &hitT2[0][2], "t2z1/F");
	output_tree2->Branch("t2x2", &hitT2[1][0], "t2x2/F");
	output_tree2->Branch("t2y2", &hitT2[1][1], "t2y2/F");
	output_tree2->Branch("t2z2", &hitT2[1][2], "t2z2/F");

	output_tree2->Branch("t3x1u", &hitT3up[0], "t3x1u/F");
	output_tree2->Branch("t3x1d", &hitT3down[0], "t3x1d/F");
	output_tree2->Branch("t3uu", &hitT3up[1], "t3uu/F");
	output_tree2->Branch("t3ud", &hitT3down[1], "t3ud/F");
	output_tree2->Branch("t3vu", &hitT3up[2], "t3vu/F");
	output_tree2->Branch("t3vd", &hitT3down[2], "t3vd/F");
	output_tree2->Branch("t3x2u", &hitT3up[3], "t3x2u/F");
	output_tree2->Branch("t3x2d", &hitT3down[3], "t3x2d/F");

	output_tree2->Branch("p", &p, "p/D");
	output_tree2->Branch("eta", &eta, "eta/D");

	if (pid_branch_exist_test)
		output_tree2->Branch("particleID", &particleID, "particleID/I");

	TRandom rndgen = TRandom(time(0));

	for (int i = 0; i < nentries; i++)
	{
		hitT1up[0] = 0.;
		hitT1up[1] = 0.;
		hitT1up[2] = 0.;
		hitT1up[3] = 0.;
		hitT1down[0] = 0.;
		hitT1down[1] = 0.;
		hitT1down[2] = 0.;
		hitT1down[3] = 0.;

		hitT2[0][0] = 0.;
		hitT2[0][1] = 0.;
		hitT2[0][2] = 0.;
		hitT2[1][0] = 0.;
		hitT2[1][1] = 0.;
		hitT2[1][2] = 0.;

		hitT3up[0] = 0.;
		hitT3up[1] = 0.;
		hitT3up[2] = 0.;
		hitT3up[3] = 0.;
		hitT3down[0] = 0.;
		hitT3down[1] = 0.;
		hitT3down[2] = 0.;
		hitT3down[3] = 0.;

		eta = 0.;

		bT1HitUp[0] = false;
		bT1HitUp[1] = false;
		bT1HitUp[2] = false;
		bT1HitUp[3] = false;

		bT1HitDown[0] = false;
		bT1HitDown[1] = false;
		bT1HitDown[2] = false;
		bT1HitDown[3] = false;

		bT2Hit[0] = false;
		bT2Hit[1] = false;

		bT3HitUp[0] = false;
		bT3HitUp[1] = false;
		bT3HitUp[2] = false;
		bT3HitUp[3] = false;

		bT3HitDown[0] = false;
		bT3HitDown[1] = false;
		bT3HitDown[2] = false;
		bT3HitDown[3] = false;

		if (i % 100000 == 0)
			cout << "Progress: " << i * 100. / nentries << endl;

		tracks->GetEntry(i);

		eta = 1. / 2 * log((p + pz) * 1.0 / (p - pz)); // log is in e base in C++

		std::vector<std::shared_ptr<hit>> v_hits;

		for (unsigned int j = 0; j < 12; j++)
			v_hits.push_back(NULL);

#ifdef debug
		cout << "v_hists filled with NULL" << endl;
#endif

		for (int j = 29; j >= 0; j--)
		{
			// cout << "Checking if hit is in the big box" << endl;
			if (not IsInBox(hitxpos[j], hitypos[j]))
				continue;

#ifdef debug
			cout << "Hit in the box area" << endl;
			cout << "Creating temporary hit: " << hitxpos[j] << " " << hitypos[j] << " " << hitzpos[j] << endl;
#endif
			std::shared_ptr<hit> _hit(new hit(hitxpos[j], hitypos[j], hitzpos[j]));

			if (_hit->z > 7816 && _hit->z < 7836)
			{
				// hit1 = _hit;
				v_hits[0] = _hit;
				if (_hit->y > 0)
				{
					hitT1up[0] = _hit->x;
					bT1HitUp[0] = true;
				}
				else
				{
					hitT1down[0] = _hit->x;
					bT1HitDown[0] = true;
				}
#ifdef debug
				cout << "Found T1 x" << endl;
#endif
			}
			else if (_hit->z > 7880 && _hit->z < 7910)
			{
				// hit12u = _hit;
				v_hits[1] = _hit;
				if (uUpOrDown(_hit->x, _hit->y))
				{
					hitT1up[1] = _hit->x * cos(5.0 / 180) + _hit->y * sin(5.0 / 180);
					bT1HitUp[1] = true;
				}
				else
				{
					hitT1down[1] = _hit->x * cos(5.0 / 180) + _hit->y * sin(5.0 / 180);
					bT1HitDown[1] = true;
				}

#ifdef debug
				cout << "Found hit12u" << endl;
#endif
			}
			else if (_hit->z > 7950 && _hit->z < 7980)
			{
				// hit12v = _hit;
				v_hits[2] = _hit;
				if (vUpOrDown(_hit->x, _hit->y))
				{
					hitT1up[2] = _hit->x * cos(5.0 / 180) - _hit->y * sin(5.0 / 180);
					bT1HitUp[2] = true;
				}
				else
				{
					hitT1down[2] = _hit->x * cos(5.0 / 180) - _hit->y * sin(5.0 / 180);
					bT1HitDown[2] = true;
				}

#ifdef debug
				cout << "Found hit12v" << endl;
#endif
			}
			else if (_hit->z > 8024 && _hit->z < 8047)
			{
				// hit2 = _hit;
				v_hits[3] = _hit;
				if (_hit->y > 0)
				{
					hitT1up[3] = _hit->x;
					bT1HitUp[3] = true;
				}
				else
				{
					hitT1down[3] = _hit->x;
					bT1HitDown[3] = true;
				}
#ifdef debug
				cout << "Found hit2" << endl;
#endif
			}
			else if (_hit->z > 8497 && _hit->z < 8518)
			{
				// hit3 = _hit;
				v_hits[4] = _hit;
				hitT2[0][0] = _hit->x;
				hitT2[0][1] = _hit->y;
				hitT2[0][2] = _hit->z;
				bT2Hit[0] = true;
#ifdef debug
				cout << "Found hit3" << endl;
#endif
			}
			if (_hit->z > 8560 && _hit->z < 8600)
			{
				// hit34u = _hit;
				v_hits[5] = _hit;
#ifdef debug
				cout << "Found hit34u" << endl;
#endif
			}
			if (_hit->z > 8630 && _hit->z < 8670)
			{
				// hit34v = _hit;
				v_hits[6] = _hit;
#ifdef debug
				cout << "Found hit34v" << endl;
#endif
			}
			if (_hit->z > 8706 && _hit->z < 8728)
			{
				// hit4 = _hit;
				v_hits[7] = _hit;
				hitT2[1][0] = _hit->x;
				hitT2[1][1] = _hit->y;
				hitT2[1][2] = _hit->z;
				bT2Hit[1] = true;
#ifdef debug
				cout << "Found hit4" << endl;
#endif
			}
			if (_hit->z > 9183 && _hit->z < 9203)
			{
				// hit5 = _hit;
				v_hits[8] = _hit;
				if (_hit->y > 0)
				{
					hitT3up[0] = _hit->x;
					bT3HitUp[0] = true;
				}
				else
				{
					hitT3down[0] = _hit->x;
					bT3HitDown[0] = true;
				}

#ifdef debug
				cout << "Found hit5" << endl;
#endif
			}
			if (_hit->z > 9240 && _hit->z < 9280)
			{
				// hit56u = _hit;
				v_hits[9] = _hit;
				if (uUpOrDown(_hit->x, _hit->y))
				{
					hitT3up[1] = _hit->x * cos(5.0 / 180) + _hit->y * sin(5.0 / 180);
					bT3HitUp[1] = true;
				}
				else
				{
					hitT3down[1] = _hit->x * cos(5.0 / 180) + _hit->y * sin(5.0 / 180);
					bT3HitDown[1] = true;
				}
#ifdef debug
				cout << "Found hit56u" << endl;
#endif
			}
			if (_hit->z > 9310 && _hit->z < 9350)
			{
				// hit56v = _hit;
				v_hits[10] = _hit;
				if (vUpOrDown(_hit->x, _hit->y))
				{
					hitT3up[2] = _hit->x * cos(5.0 / 180) - _hit->y * sin(5.0 / 180);
					bT3HitUp[2] = true;
				}
				else
				{
					hitT3down[2] = _hit->x * cos(5.0 / 180) - _hit->y * sin(5.0 / 180);
					bT3HitDown[2] = true;
				}

#ifdef debug
				cout << "Found hit56v" << endl;
#endif
			}
			if (_hit->z > 9391 && _hit->z < 9414)
			{
				// hit6 = _hit;
				v_hits[11] = _hit;
				if (_hit->y > 0)
				{
					hitT3up[3] = _hit->x;
					bT3HitUp[3] = true;
				}
				else
				{
					hitT3down[3] = _hit->x;
					bT3HitDown[3] = true;
				}

#ifdef debug
				cout << "Found hit6" << endl;
#endif
			}
		}
#ifdef debug
		cout << "Starting Hits recover algorithm" << endl;
#endif
		for (unsigned int k = 1; k < 7; k++)
		{
			// Only recover the silicon layer hits
			if (not(k == 3 or k == 4))
				continue;
			// loop over the k MT layers 1-6
			std::vector<std::shared_ptr<hit>> inputs;
			int z = ((k) / 2) * 4;
			if (k % 2 == 0)
				z -= 1;
			// z point to the index of x
			//cout << "z: " << z << endl;
			if (k % 2 != 0)
			{
				if (v_hits[z + 1] and v_hits[z + 2] and v_hits[z + 3])
				{
					inputs.push_back(v_hits[z + 1]);
					inputs.push_back(v_hits[z + 2]);
					inputs.push_back(v_hits[z + 3]);
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
			}
			else
			{
				if (v_hits[z - 1] and v_hits[z - 2] and v_hits[z - 3])
				{
					inputs.push_back(v_hits[z - 1]);
					inputs.push_back(v_hits[z - 2]);
					inputs.push_back(v_hits[z - 3]);
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
				{
					_mhit = v_hits[z];
					//cout << "mhit assigned" << endl;
				}
				if (_mhit)
				{
					// Draw Figure
				}
				else
				{
					// Original Hit don't exist, recover it
					// Only recover the hits in silicon area and proper layer and not
					if (IsInSiliconArea(_hit->x, _hit->y, surface))
					{
						hitT2[k - 3][0] = _hit->x;
						hitT2[k - 3][0] = _hit->y;
						hitT2[k - 3][0] = _hit->z;
						bT2Hit[k - 3] = true;
					}
				}
			}
		} // Recover hit of silicon

		if (blurIt)
		{
			for (int j = 0; j < 4; j++)
			{
				if (bT1HitUp[j])
				{
					hitT1up[j] += rndgen.Gaus(0.0, 0.070);
				}
				if (bT3HitUp[j])
				{
					hitT3up[j] += rndgen.Gaus(0.0, 0.070);
				}
				if (bT1HitDown[j])
				{
					hitT1down[j] += rndgen.Gaus(0.0, 0.070);
				}
				if (bT3HitDown[j])
				{
					hitT3down[j] += rndgen.Gaus(0.0, 0.070);
				}
				// cout << bT1HitUp[j] + bT3HitUp[j] + bT1HitDown[j] + bT3HitDown[j] << endl;
			}
			if (bT2Hit[0])
			{
				if (IsInSiliconArea(hitT2[0][0], hitT2[0][1], surface))
				{
					hitT2[0][0] += rndgen.Gaus(0.0, 0.050 / sqrt(12));
					hitT2[0][1] += rndgen.Gaus(0.0, 0.150 / sqrt(12));
				}
				else
				{
					hitT2[0][0] += rndgen.Gaus(0.0, 0.070);
					hitT2[0][1] += rndgen.Gaus(0.0, 1.134);
				}
			}
			if (bT2Hit[1])
			{
				if (IsInSiliconArea(hitT2[1][0], hitT2[1][1], surface))
				{
					hitT2[1][0] += rndgen.Gaus(0.0, 0.050 / sqrt(12));
					hitT2[1][1] += rndgen.Gaus(0.0, 0.150 / sqrt(12));
				}
				else
				{
					hitT2[1][0] += rndgen.Gaus(0.0, 0.070);
					hitT2[1][1] += rndgen.Gaus(0.0, 1.134);
				}
			}
		}
		output_tree2->Fill();
		v_hits.clear();
	} // Loop on entries
	hfile->Write();
}
