#include <cmath>
#include <fstream>
#include <iostream>

#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLegend.h"
#include "TF1.h"
#include "TText.h"

#include "TVector3.h"

#include "TGraph.h"
#include "TPad.h"
#include "TFitResult.h"
#include "TStyle.h"
#include <chrono>
#include <math.h>

using namespace std;

// #define modfication
#define layerStudy

class track;

class hit3
{
public:
    float x, y, z;
    float eventID = -1;
    double hP = -1;
    double hEta = -1;
    int particleID = 0;
    int used = 0;
    bool flag = false;

    std::vector<track *> relatedTracks;
    track *motherTrack = NULL;

    bool specialFlag = false;
    hit3(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    hit3(float _eventID, float _hP, float _hEta, float _x, float _y, float _z)
    {
        hEta = _hEta;
        hP = _hP;
        eventID = _eventID;
        x = _x;
        y = _y;
        z = _z;
    }
    hit3(float _eventID, float _hP, float _hEta, int _particleID, float _x, float _y, float _z)
    {
        particleID = _particleID;
        hEta = _hEta;
        hP = _hP;
        eventID = _eventID;
        x = _x;
        y = _y;
        z = _z;
    }

    ~hit3()
    {
    }

    int operator==(const hit3 *anotherHit)
    {
        if ((anotherHit->x == this->x) and (anotherHit->y == this->y) and (anotherHit->z == this->z) and (anotherHit->hP == this->hP) and (anotherHit->hEta == this->hEta))
            return 1;
        else
            return 0;
    }

private:
};

class hit1
{
public:
    // Position along the x/u/v axis
    float s;

    float eventID = -1;
    double hP = -1;
    double hEta = -1;
    int particleID = 0;
    int used = 0;
    bool flag = false;

    std::vector<track *> relatedTracks;
    track *motherTrack = NULL;

    bool specialFlag = false;
    hit1(float _s)
    {
        s = _s;
    }
    hit1(float _eventID, float _hP, float _hEta, float _s)
    {
        hEta = _hEta;
        hP = _hP;
        eventID = _eventID;
        s = _s;
    }
    hit1(float _eventID, float _hP, float _hEta, int _particleID, float _s)
    {
        particleID = _particleID;
        _hEta = _hEta;
        hP = _hP;
        eventID = _eventID;
        s = _s;
    }

    ~hit1()
    {
    }

    int operator==(const hit1 *anotherHit)
    {
        if ((anotherHit->s == this->s) and (anotherHit->hP == this->hP) and (anotherHit->hEta == this->hEta))
            return 1;
        else
            return 0;
    }

private:
};

class SiliconLayer
{
public:
    std::vector<hit3 *> hits;
};

class SciFiLayer
{
public:
    std::vector<hit1 *> hits;
};

class track
{
public:
    hit1 *hitT1Set[4] = {NULL, NULL, NULL, NULL};
    hit3 *hitT2Set[2] = {NULL, NULL};
    hit1 *hitT3Set[4] = {NULL, NULL, NULL, NULL};

    float xMeasure[6];
    float yMeasure[6];

    float pChi2X = 0;
    float pChi2Y = 0;

    float fChi2X = 0;
    float fChi2Y = 0;

    float dFitX[6] = {0, 0, 0, 0, 0, 0};
    float dFitY[6] = {0, 0, 0, 0, 0, 0};

    bool match[9];
    bool trackMatch;

    bool trackAccept = true;
    int nTrueHits = 0;
    int theTrueHit = 0;
    float purity = 0.0;
    track *motherTrack = NULL;

    // Only for motherTracks
    bool particleMatched = false;
    bool bFriendlyFired = false;

    track(hit1 **_hitT1Set, hit3 **_hitT2Set, hit1 **_hitT3Set, bool _match)
    {
        for (Int_t j = 0; j < 4; j++)
        {

            hitT1Set[j] = _hitT1Set[j];
            _hitT1Set[j]->used++;
            _hitT1Set[j]->relatedTracks.push_back(this);
            if (j == 0)
            {
                xMeasure[j] = hitT1Set[j]->s;
            }
            else if (j == 3)
            {
                xMeasure[j - 2] = hitT1Set[j]->s;
            }

            hitT3Set[j] = _hitT3Set[j];
            _hitT3Set[j]->used++;
            _hitT3Set[j]->relatedTracks.push_back(this);
            if (j == 0)
            {
                xMeasure[j + 4] = hitT3Set[j]->s;
            }
            else if (j == 3)
            {
                xMeasure[j + 2] = hitT3Set[j]->s;
            }

            if (j < 2)
            {
                hitT2Set[j] = _hitT2Set[j];
                _hitT2Set[j]->used++;
                _hitT2Set[j]->relatedTracks.push_back(this);
                xMeasure[j + 2] = hitT2Set[j]->x;
                yMeasure[j + 2] = hitT2Set[j]->y;
            }
        }

        float zFit1[6] = {7825.9, 8035.8, 8507.8, 8717.9, 9192.9, 9402.8};

        TGraph gx = TGraph(6, zFit1, xMeasure);
        TFitResultPtr fpx = gx.Fit("pol2", "SQ");
        fChi2X = fpx->Chi2();
        TF1 *fx = gx.GetFunction("pol2");

        yMeasure[0] = (hitT1Set[1]->s - fx->Eval(7895) * cos(5.0 / 180)) / sin(5.0 / 180);
        yMeasure[1] = (hitT1Set[2]->s - fx->Eval(7965) * cos(5.0 / 180)) / sin(-5.0 / 180);
        yMeasure[4] = (hitT3Set[1]->s - fx->Eval(9260) * cos(5.0 / 180)) / sin(5.0 / 180);
        yMeasure[5] = (hitT3Set[2]->s - fx->Eval(9330) * cos(5.0 / 180)) / sin(-5.0 / 180);

        float zFit2[6] = {7895, 7965, 8507.8, 8717.9, 9260, 9330};
        TGraph gy = TGraph(6, zFit2, yMeasure);
        TFitResultPtr fpy = gy.Fit("pol1", "SQ");
        fChi2Y = fpy->Chi2();

        trackMatch = _match;
        if (not trackMatch)
        {
            for (Int_t j = 0; j < 4; j++)
            {
                if (nTrueHits < 5)
                {
                    if ((hitT1Set[j]->motherTrack != NULL) && (*(hitT1Set[j]->motherTrack) == *this) > nTrueHits)
                    {
                        nTrueHits = (*(hitT1Set[j]->motherTrack) == *this);
                        motherTrack = hitT1Set[j]->motherTrack;
                        theTrueHit = j;
                    }
                    if (j < 2)
                    {
                        if ((hitT2Set[j]->motherTrack != NULL) && (*(hitT2Set[j]->motherTrack) == *this) > nTrueHits)
                        {
                            nTrueHits = (*(hitT2Set[j]->motherTrack) == *this);
                            motherTrack = hitT2Set[j]->motherTrack;
                            theTrueHit = j + 4;
                        }
                    }
                    if ((hitT3Set[j]->motherTrack != NULL) && (*(hitT3Set[j]->motherTrack) == *this) > nTrueHits)
                    {
                        nTrueHits = (*(hitT3Set[j]->motherTrack) == *this);
                        motherTrack = hitT3Set[j]->motherTrack;
                        theTrueHit = j + 6;
                    }
                }
                else
                    break;
            }
        }
        else
        {
            nTrueHits = 10;
            theTrueHit = 0;
            motherTrack = hitT2Set[0]->motherTrack;
        }
        if (motherTrack)
        {
            purity = nTrueHits * 1.0 / motherTrack->nTrueHits;
            if (purity >= 0.7)
                trackMatch = true;
        }
    }

    track(hit1 **_hitT1Set, hit3 **_hitT2Set, hit1 **_hitT3Set, int _nTrueHits)
    {
        for (Int_t j = 0; j < 4; j++)
        {
            if (_hitT1Set[j] != NULL)
            {
                hitT1Set[j] = _hitT1Set[j];
            }
            if (_hitT3Set[j] != NULL)
            {
                hitT3Set[j] = _hitT3Set[j];
            }
            if (j < 2)
            {
                if (_hitT2Set[j] != NULL)
                {
                    hitT2Set[j] = _hitT2Set[j];
                }
            }
        }
        nTrueHits = _nTrueHits;
    }

    int GetUsedSum()
    {
        int usedSum = 0;
        for (Int_t j = 0; j < 4; j++)
        {
            if (hitT1Set[j] != NULL)
            {
                usedSum += hitT1Set[j]->used;
            }
            if (hitT3Set[j] != NULL)
            {
                usedSum += hitT3Set[j]->used;
            }
            if (j < 2)
            {
                if (hitT2Set[j] != NULL)
                {
                    usedSum += hitT2Set[j]->used;
                }
            }
        }
        return usedSum;
    }

    int operator==(const track anotherTrack)
    {
        int sum = 0;
        for (Int_t j = 0; j < 4; j++)
        {
            if ((this->hitT1Set[j] != NULL) && (anotherTrack.hitT1Set[j] != NULL))
            {
                sum += (this->hitT1Set[j] == anotherTrack.hitT1Set[j]);
            }
            if ((this->hitT3Set[j] != NULL) && (anotherTrack.hitT3Set[j] != NULL))
            {
                sum += (this->hitT3Set[j] == anotherTrack.hitT3Set[j]);
            }
            if (j < 2)
            {
                if ((this->hitT2Set[j] != NULL) && (anotherTrack.hitT2Set[j] != NULL))
                {
                    sum += (this->hitT2Set[j] == anotherTrack.hitT2Set[j]);
                }
            }
        }
        // Return the number of hits matched
        return sum;
    }
};

bool IsInDeadArea(float x, float y)
{
    if (x < 100 and x > -100)
        if (y < 100 and y > -100)
            return false; // beam area

    // if (fabs(remainder((x / 20), 1)) < 0.1 / 20)
    // {
    //     if (fabs(remainder(((x - 20) / 20), 2)) < 2 * 0.1 / 20)
    //         return false;
    //     return true;
    // }

    // if (fabs(remainder((y / 20), 1)) < 0.1 / 20)
    // {
    //     return true;
    // }
    return false;
}

bool uUpOrDown(float x, float y)
{
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

int linearFit(float *linFitX, float *linFitZ, float *parabolaParams, float *fitResult)
{
    // Use hit 3, 4, 5 and 6, predict hit1x and hit2x.
    TGraph gx = TGraph(4, linFitZ, linFitX);

    TFitResultPtr fpx(0);
    fpx = gx.Fit("pol1", "QS");
    TF1 *fx = gx.GetFunction("pol1");

    parabolaParams[0] = fx->GetParameter(0);
    parabolaParams[1] = fx->GetParameter(1);

    fitResult[0] = fx->Eval(7825.9);
    // fitResult[1] = fx->Eval(7895);
    // fitResult[2] = fx->Eval(7965);
    fitResult[3] = fx->Eval(8035.8);
    return true;
}

int parabolaFit(float *linFitX, float *linFitZ, float *parabolaParams, float *fitResult)
{
    // Use hit 3, 4, 5 and 6, predict hit1x and hit2x.
    TGraph gx = TGraph(4, linFitZ, linFitX);

    TFitResultPtr fpx(0);
    fpx = gx.Fit("pol2", "QS");
    TF1 *fx = gx.GetFunction("pol2");

    parabolaParams[0] = fx->GetParameter(0);
    parabolaParams[1] = fx->GetParameter(1);
    parabolaParams[2] = fx->GetParameter(2);

    fitResult[0] = fx->Eval(7825.9);
    // fitResult[1] = fx->Eval(7895);
    // fitResult[2] = fx->Eval(7965);
    fitResult[3] = fx->Eval(8035.8);
    return true;
}

int cubicFit(float *linFitX, float *linFitZ, float *cubicParams, float *fitResult)
{
    // Use hit 3, 4, 5 and 6, predict hit1x and hit2x.
    TGraph gx = TGraph(4, linFitZ, linFitX);

    TFitResultPtr fpx(0);
    fpx = gx.Fit("pol3", "QS");
    TF1 *fx = gx.GetFunction("pol3");

    cubicParams[0] = fx->GetParameter(0);
    cubicParams[1] = fx->GetParameter(1);
    cubicParams[2] = fx->GetParameter(2);
    cubicParams[3] = fx->GetParameter(3);

    fitResult[0] = fx->Eval(7825.9);
    // fitResult[1] = fx->Eval(7895);
    // fitResult[2] = fx->Eval(7965);
    fitResult[3] = fx->Eval(8035.8);
    return true;
}

void make_inefficiency_window_plot(TCanvas *c1, TH1F *h1_good, TH1F *h1_bad, TString outputplotsfolder)
{
    h1_good->Scale(1. / h1_good->Integral());
    h1_bad->Scale(1. / h1_bad->Integral());
    h1_good->SetLineWidth(3);
    h1_good->SetLineColor(kGreen - 3);
    h1_bad->SetLineWidth(3);
    h1_bad->SetLineColor(kRed);
    c1->SetGrid();
    c1->SetLogy(1);
    h1_good->GetCumulative(false)->Draw();
    float threshold999 = 0;
    float threshold99 = 0;
    float threshold98 = 0;
    float threshold95 = 0;
    for (unsigned int i = 0; i < h1_good->GetXaxis()->GetNbins(); i++)
    {
        if (h1_good->GetCumulative(false)->GetBinContent(i) > 1e-4)
        {
            threshold999 = h1_good->GetBinLowEdge(i);
        }
        if (h1_good->GetCumulative(false)->GetBinContent(i) > 1e-3)
        {
            threshold99 = h1_good->GetBinLowEdge(i);
        }
        if (h1_good->GetCumulative(false)->GetBinContent(i) > 2e-3)
        {
            threshold98 = h1_good->GetBinLowEdge(i);
        }
        if (h1_good->GetCumulative(false)->GetBinContent(i) > 5e-3)
        {
            threshold95 = h1_good->GetBinLowEdge(i);
        }
    }

    h1_bad->GetCumulative(false)->Draw("same");
    TText *xlabel = new TText();
    xlabel->SetNDC();
    xlabel->SetTextFont(1);
    xlabel->SetTextColor(1);
    xlabel->SetTextSize(0.03);
    //xlabel->SetTextAlign(0.8);
    //xlabel->SetTextAngle(0);
    xlabel->DrawText(0.15, 0.3, Form("Threshold for Eff.(99.99\%): %.3f[mm]", threshold999));
    xlabel->DrawText(0.15, 0.25, Form("Threshold for Eff.(99.9\%): %.3f[mm]", threshold99));
    xlabel->DrawText(0.15, 0.20, Form("Threshold for Eff.(99.8\%): %.3f[mm]", threshold98));
    xlabel->DrawText(0.15, 0.15, Form("Threshold for Eff.(99.5\%): %.3f[mm]", threshold95));
    c1->SaveAs(outputplotsfolder + +"/ineff/" + h1_good->GetName() + "_inefficiency.png");

    h1_good->GetYaxis()->SetTitle("%");
    h1_bad->GetYaxis()->SetTitle("%");
}

void stereoSeeding(TString settings = "1Casesp20", TString inputfilestring = "2LayerBig.root", TString outputplotsfolder = "2LayerU1Big/10HitDef/")
{
    // system("rm -rf " + outputplotsfolder);
    outputplotsfolder = outputplotsfolder + settings + '/';
    system("mkdir -p " + outputplotsfolder);
    system("mkdir -p " + outputplotsfolder + "h1");
    system("mkdir -p " + outputplotsfolder + "h2");
    system("mkdir -p " + outputplotsfolder + "ineff");

    cout << "Settings: " << settings << endl;
    ofstream logOut;
    logOut.open(outputplotsfolder + "run.log", ios_base::app);
    logOut << "Settings: " << settings << endl;
    logOut << "inputfilestring: " << inputfilestring << endl;

    bool conf_eta = true;
    float conf_min_p = 5000.;
    float conf_max_p = 7.1e6;
    float conf_max_sx = -1.;
    float conf_max_sy = -1.;
    float conf_deadarea = 0.;
    float conf_correctionT2 = 0;
    float conf_alone = 0;

    // No p0p5 tracks is reconstructed in the first case
    // // Finished windows for MCTruth Hits
    // Float_t MAXDX[5][5] = {{20, 0.5, 0.5, 0.5, 0.5}, {40, 3, 0.8, 1.5, 2}, {80, 10, 3, 5, 5}, {120, 30, 5, 7, 7}, {150, 70, 6, 10, 10}};
    // Float_t MAXDY[5][5] = {{3, 0.5, 0.5, 0.5, 0.5}, {10, 0.8, 0.8, 0.8, 0.8}, {20, 2, 2, 2, 2}, {30, 3, 3, 3, 3}, {50, 15, 6, 5, 5}};

    // Finished windows for Blured Hits
    // // Few p0p5 tracks is reconstructed in the first case
    // Float_t MAXDX[5][5] = {{20, 4, 0.8, 2.4, 1.6}, {50, 4.8, 0.8, 1.6, 2}, {80, 10, 5, 5, 8}, {120, 30, 7, 14, 11}, {160, 120, 7, 40, 30}};
    // Float_t MAXDY[5][5] = {{4.9, 12.5, 6.5, 11.1, 13}, {14, 12, 5.5, 11, 13.5}, {40, 18, 7, 14, 15}, {46, 20, 8, 15, 18}, {50, 26, 16, 20, 20}};
    // Float_t MAXDX[5][5] = {{20, 4, 0.8, 2.4, 1.6}, {50, 4.8, 0.8, 1.6, 2}, {80, 10, 5, 5, 8}, {120, 30, 7, 14, 11}, {160, 80, 7, 14, 11}};
    // Float_t MAXDY[5][5] = {{4.9, 12.5, 6.5, 11.1, 13}, {14, 12, 5.5, 11, 13.5}, {40, 18, 7, 14, 15}, {46, 20, 8, 15, 18}, {50, 20, 8, 15, 18}};

    // T2 Windows
    float DXT2[5] = {50, 50, 110, 160, 220};
    float DYT2[5] = {3, 20, 30, 50, 100};
    // T3 Windows
    float DX1T3[5] = {2, 10, 20, 50, 80};
    float DX2T3[5] = {0.5, 3, 10, 18, 25};
    float DUT3[5] = {0.5, 1.25, 1.75, 2, 2.5};
    float DVT3[5] = {0.5, 1.25, 1.75, 2, 2.5};
    // T1 Windows
    float DX2T1[5] = {2, 10, 30, 45, 60};
    float DVT1[5] = {0.5, 0.5, 1.75, 2, 2};
    float DUT1[5] = {0.5, 0.5, 1.75, 2, 2};
    float DX1T1[5] = {2, 15, 35, 50, 60};

    int correctionT2[5] = {5, 10, 35, 40, 45};
    int nShareHits[5] = {9, 8, 7, 6, 5};
    int nGlobalShareHits = 4;

    if (settings.Contains("dead1"))
        conf_deadarea = 1.;
    if (settings.Contains("eta0"))
        conf_eta = 0.;

    if (settings.Contains("p0"))
        conf_min_p = 0.;
    if (settings.Contains("p1"))
        conf_min_p = 1000.;
    if (settings.Contains("p5"))
        conf_min_p = 5000.;
    if (settings.Contains("p20"))
        conf_min_p = 20000.;

    int iCases = 0;
    int nCases = 5;

    if (settings.Contains("Case5"))
    {
        iCases = 4;
        nCases = 5;
    }
    else if (settings.Contains("Case4"))
    {
        iCases = 3;
        nCases = 4;
    }
    else if (settings.Contains("Case3"))
    {
        iCases = 2;
        nCases = 5;
    }
    else if (settings.Contains("Case2"))
    {
        iCases = 1;
        nCases = 2;
    }
    else if (settings.Contains("Case1"))
    {
        iCases = 0;
        nCases = 1;
    }

    if (settings.Contains("5Cases"))
    {
        iCases = 0;
        nCases = 5;
    }
    else if (settings.Contains("4Cases"))
    {
        iCases = 0;
        nCases = 4;
    }
    else if (settings.Contains("3Cases"))
    {
        iCases = 0;
        nCases = 3;
    }
    else if (settings.Contains("2Cases"))
    {
        iCases = 0;
        nCases = 2;
    }
    else if (settings.Contains("1Cases"))
    {
        iCases = 0;
        nCases = 1;
    }

    if (settings.Contains("p0") && settings.Contains("p5"))
    {
        conf_max_p = 5000;
        conf_min_p = 0.;
    }
    if (settings.Contains("p5") && settings.Contains("p20"))
    {
        conf_max_p = 20000;
        conf_min_p = 5000.;
    }

    if (settings.Contains("sx10"))
        conf_max_sx = 10.;
    if (settings.Contains("sx15"))
        conf_max_sx = 15.;
    if (settings.Contains("sx20"))
        conf_max_sx = 20.;

    if (settings.Contains("sy2"))
        conf_max_sy = 2.;
    if (settings.Contains("sy3"))
        conf_max_sy = 3.;
    if (settings.Contains("sy4"))
        conf_max_sy = 4.;

    // By default remove the fake tracks
    if (settings.Contains("alone"))
        conf_alone = 1;

    TCanvas *c1 = new TCanvas();
    c1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0101);

    // Track reconstructed counter
    Int_t nReconstructible = 0;

    Int_t nTrueReconstructed = 0;

    Int_t nClone = 0;
    Int_t notFound = 0;
    Int_t nFakeReconstructed = 0;
    Int_t nPerfectFake = 0;
    Int_t nTotalReconstructed = 0;

    Int_t nFriendlyFire = 0;
    Int_t nNotFoundAsFired = 0;
    Int_t nRejectedTracks = 0;

    // Input file directory
    TFile f(inputfilestring);

    // Pointers to original and recovered hits track trees respectively
    TTree *rTracks = NULL;

    // Check if the recovered tracks available
    rTracks = (TTree *)gDirectory->Get("RTracks");
    if (not rTracks)
        cout << "rTracks not found" << endl;

    // Hit position

    // T1 x-u-v-x hits
    float hitT1up[4] = {0, 0, 0, 0};
    float hitT1down[4] = {0, 0, 0, 0};
    // Silicon Hits of T2
    float hitT2[2][3] = {{0, 0, 0}, {0, 0, 0}};
    // T3 x-u-v-x hits
    float hitT3up[4] = {0, 0, 0, 0};
    float hitT3down[4] = {0, 0, 0, 0};

    // Kinetics Variables
    double eta = 0;
    double p = 0;

    // Number of entries in the tree
    int nEntries;

    // Event number of the tracks
    ULong64_t eventNumber = -1;
    int eventCounter = 0;

    // PID branch
    TBranch *pidBranchExistTest = NULL;
    Int_t particleID = 0;

    if (rTracks)
    {
        // Set the branches for recovered hits
        cout << "rTracks found" << endl;

        rTracks->SetBranchAddress("eventNumber", &eventNumber);
        rTracks->SetBranchAddress("t1x1u", &hitT1up[0]);
        rTracks->SetBranchAddress("t1uu", &hitT1up[1]);
        rTracks->SetBranchAddress("t1vu", &hitT1up[2]);
        rTracks->SetBranchAddress("t1x2u", &hitT1up[3]);
        rTracks->SetBranchAddress("t1x1d", &hitT1down[0]);
        rTracks->SetBranchAddress("t1ud", &hitT1down[1]);
        rTracks->SetBranchAddress("t1vd", &hitT1down[2]);
        rTracks->SetBranchAddress("t1x2d", &hitT1down[3]);

        rTracks->SetBranchAddress("t2x1", &hitT2[0][0]);
        rTracks->SetBranchAddress("t2y1", &hitT2[0][1]);
        rTracks->SetBranchAddress("t2z1", &hitT2[0][2]);
        rTracks->SetBranchAddress("t2x2", &hitT2[1][0]);
        rTracks->SetBranchAddress("t2y2", &hitT2[1][1]);
        rTracks->SetBranchAddress("t2z2", &hitT2[1][2]);

        rTracks->SetBranchAddress("t3x1u", &hitT3up[0]);
        rTracks->SetBranchAddress("t3uu", &hitT3up[1]);
        rTracks->SetBranchAddress("t3vu", &hitT3up[2]);
        rTracks->SetBranchAddress("t3x2u", &hitT3up[3]);
        rTracks->SetBranchAddress("t3x1d", &hitT3down[0]);
        rTracks->SetBranchAddress("t3ud", &hitT3down[1]);
        rTracks->SetBranchAddress("t3vd", &hitT3down[2]);
        rTracks->SetBranchAddress("t3x2d", &hitT3down[3]);

        rTracks->SetBranchAddress("p", &p);
        rTracks->SetBranchAddress("eta", &eta);

        pidBranchExistTest = (TBranch *)rTracks->GetListOfBranches()->FindObject("particleID");
        if (pidBranchExistTest)
            rTracks->SetBranchAddress("particleID", &particleID);
        nEntries = rTracks->GetEntries();
    }
    else
    {
        cout << "Error: Recovered Tracks NOT found." << endl;
    }
    logOut << "Total number of entries in Tree: " << nEntries << endl;

    TH1F *h1_hitPurity = new TH1F("h1_hitPurity", "h1_hitPurity", 100, 0, 1.1);
    TH1F *h1_whereFake = new TH1F("h1_whereFake", "h1_whereFake", 100, 0, 11);
    TH1F *h1_notFound = new TH1F("h1_notFound", "h1_notFound", 100, 0, 11);
    TH1F *h1_whereNotFound = new TH1F("h1_whereNotFound", "h1_whereNotFound", 100, 0, 10);
    TH1F *h1_rejectPurity = new TH1F("h1_rejectPurity", "h1_rejectPurity", 100, 0, 1.1);

    // T1 Figures
    TH1F *h1ErrT1x1Good = new TH1F("h1ErrT1x1Good", "h1ErrT1x1Good", 160, 0, 160);
    h1ErrT1x1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT1x1Bad = new TH1F("h1ErrT1x1Bad", "h1ErrT1x1Bad", 160, 0, 160);
    h1ErrT1x1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT1x2Good = new TH1F("h1ErrT1x2Good", "h1ErrT1x2Good", 160, 0, 160);
    h1ErrT1x2Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT1x2Bad = new TH1F("h1ErrT1x2Bad", "h1ErrT1x2Bad", 160, 0, 160);
    h1ErrT1x2Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT1u1Good = new TH1F("h1ErrT1u1Good", "h1ErrT1u1Good", 160, 0, 80);
    h1ErrT1u1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT1u1Bad = new TH1F("h1ErrT1u1Bad", "h1ErrT1u1Bad", 160, 0, 80);
    h1ErrT1u1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT1v1Good = new TH1F("h1ErrT1v1Good", "h1ErrT1v1Good", 160, 0, 80);
    h1ErrT1v1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT1v1Bad = new TH1F("h1ErrT1v1Bad", "h1ErrT1v1Bad", 160, 0, 80);
    h1ErrT1v1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    // T2 Figures
    TH1F *h1ErrT2x2xGood = new TH1F("h2ErrT2x2xGood", "h2ErrT2x2xGood", 500, 0, 500);
    h1ErrT2x2xGood->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT2x2xBad = new TH1F("h1ErrT2x2xBad", "h1ErrT2x2xBad", 500, 0, 500);
    h1ErrT2x2xBad->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT2x2yGood = new TH1F("h1ErrT2x2yGood", "h1ErrT2x2yGood", 500, 0, 500);
    h1ErrT2x2yGood->GetXaxis()->SetTitle("Y_{pred}-Y [mm]");
    TH1F *h1ErrT2x2yBad = new TH1F("h2ErrT2x2yBad", "h2ErrT2x2yBad", 500, 0, 500);
    h1ErrT2x2yBad->GetYaxis()->SetTitle("Y_{pred}-Y [mm]");

    TH2F *h2ErrT2x2Good = new TH2F("h2ErrT2x2Good", "h2ErrT2x2Good", 400, -200, 200, 80, -40, 40);
    h2ErrT2x2Good->GetYaxis()->SetTitle("Y_{pred}-Y [mm]");
    h2ErrT2x2Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH2F *h2ErrT2x2Bad = new TH2F("h2ErrT2x2Bad", "h2ErrT2x2Bad", 400, -200, 200, 80, -40, 40);
    h2ErrT2x2Bad->GetYaxis()->SetTitle("Y_{pred}-Y [mm]");
    h2ErrT2x2Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    // T3 Figures
    TH1F *h1ErrT3x1Good = new TH1F("h1ErrT3x1Good", "h1ErrT3x1Good", 200, 0, 200);
    h1ErrT3x1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT3x1Bad = new TH1F("h1ErrT3x1Bad", "h1ErrT3x1Bad", 200, 0, 200);
    h1ErrT3x1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT3x2Good = new TH1F("h1ErrT3x2Good", "h1ErrT3x2Good", 200, 0, 200);
    h1ErrT3x2Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT3x2Bad = new TH1F("h1ErrT3x2Bad", "h1ErrT3x2Bad", 200, 0, 200);
    h1ErrT3x2Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT3v1Good = new TH1F("h1ErrT3v1Good", "h1ErrT3v1Good", 200, 0, 80);
    h1ErrT3v1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT3v1Bad = new TH1F("h1ErrT3v1Bad", "h1ErrT3v1Bad", 160, 0, 80);
    h1ErrT3v1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    TH1F *h1ErrT3u1Good = new TH1F("h1ErrT3u1Good", "h1ErrT3u1Good", 160, 0, 80);
    h1ErrT3u1Good->GetXaxis()->SetTitle("X_{pred}-X [mm]");
    TH1F *h1ErrT3u1Bad = new TH1F("h1ErrT3u1Bad", "h1ErrT3u1Bad", 160, 0, 80);
    h1ErrT3u1Bad->GetXaxis()->SetTitle("X_{pred}-X [mm]");

    // Container of hits
    SciFiLayer *T1Layers[4];
    std::vector<SciFiLayer> T1LayersUp;
    std::vector<SciFiLayer> T1LayersDown;
    std::vector<SiliconLayer> T2Layers;
    SciFiLayer *T3Layers[4];
    std::vector<SciFiLayer> T3LayersDown;
    std::vector<SciFiLayer> T3LayersUp;
    for (Int_t j = 0; j < 4; j++)
    {
        T1LayersUp.push_back(SciFiLayer());
        T1LayersDown.push_back(SciFiLayer());
        T3LayersUp.push_back(SciFiLayer());
        T3LayersDown.push_back(SciFiLayer());
        if (j < 2)
        {
            T2Layers.push_back(SiliconLayer());
        }
    }

    std::vector<track *> motherTracks;
    std::vector<track *> caseTracks;
    std::vector<track *> allTracks;

    // counters for tracks
    unsigned int nReconstructableElectron = 0, nTrueReconstructedElectron = 0;
    unsigned int nReconstructableP0P5Eta1 = 0, nTrueReconstructedP0P5Eta1 = 0;
    unsigned int nReconstructableP5P20Eta1 = 0, nTrueReconstructedP5P20Eta1 = 0;
    unsigned int nReconstructableP20Eta1 = 0, nTrueReconstructedP20Eta1 = 0;
    unsigned int nReconstructable10Hit = 0;
    unsigned int nReconstructableSilicon = 0;

    // Just for read the file
    int previousEventnumber = -1;
    int cubicCount = 0;
    int parabolaCount = 0;
    int linearCount = 0;

    // Output stream for debug
    ofstream goodOut;
    ofstream fakeOut;
    goodOut.open(outputplotsfolder + "goodTracks.txt", ios_base::app);
    fakeOut.open(outputplotsfolder + "fakeTracks.txt", ios_base::app);

    // hit1 *blankhit1=new hit1(0);
    // hit3 *blankhit3=new hit3(0);

    chrono::time_point<chrono::high_resolution_clock> timingStart = chrono::high_resolution_clock::now();

    for (Int_t i = 0; i < nEntries; i++)
    {
        previousEventnumber = eventNumber;

        if (i % 1000000 == 0)
            cout << "Progress: " << i * 100. / nEntries << endl;

        // Read a entry to memory
        rTracks->GetEntry(i);
        // if(eventCounter>10)
        // break;
        if (previousEventnumber != -1 and (previousEventnumber != eventNumber or i == nEntries - 1))
        {
            // Run the reconstruction in the end of an event or in the end of the readout
            //             // print the hits on every layer
            cout << "Event summary by layers: T1 " << T1LayersUp[0].hits.size() << " " << T1LayersUp[1].hits.size() << " " << T1LayersUp[2].hits.size() << " " << T1LayersUp[3].hits.size() << " T2: " << T2Layers[0].hits.size() << " " << T2Layers[1].hits.size() << " T3: " << T3LayersDown[0].hits.size() << " " << T3LayersDown[1].hits.size() << " " << T3LayersDown[2].hits.size() << " " << T3LayersDown[3].hits.size() << " " << endl;

            //             // Perform Reconstruction
            for (int iCase = iCases; iCase < nCases; iCase++)
            {
                for (vector<hit3 *>::iterator iT2x1 = T2Layers[0].hits.begin(); iT2x1 != T2Layers[0].hits.end(); iT2x1++)
                {
                    // New track variables
                    // Prediction
                    float xPreT1[4] = {0.0, 0.0, 0.0, 0.0};
                    float yPreT1[2] = {0.0, 0.0};
                    float uPreT1 = 0.0, vPreT1 = 0.0;

                    float xPreT2[2] = {0.0, 0.0};
                    float yPreT2[2] = {0.0, 0.0};

                    float xPreT3[4] = {0.0, 0.0, 0.0, 0.0};
                    float yPreT3[2] = {0.0, 0.0};
                    float uPreT3 = 0.0, vPreT3 = 0.0;

                    float sErrT1[4] = {0.0, 0.0, 0.0, 0.0};

                    float xErrT2[2] = {0.0, 0.0};
                    float yErrT2[2] = {0.0, 0.0};

                    float sErrT3[4] = {0.0, 0.0, 0.0, 0.0};

                    bool match[9] = {false, false, false, false, false, false, false, false, false};

                    bool bFound[10];

                    // Hits of this track
                    hit1 *trackT1[4] = {NULL, NULL, NULL, NULL};
                    hit3 *trackT2[2] = {NULL, NULL};
                    hit1 *trackT3[4] = {NULL, NULL, NULL, NULL};

                    // Start the seed from the T2x1 layer
                    trackT2[0] = *iT2x1;
                    if (trackT2[0]->flag)
                        continue;
                    // Variable Fill
                    xPreT2[0] = trackT2[0]->x;
                    yPreT2[0] = trackT2[0]->y;

// Tilt modification
#ifdef modfication
                    xPreT2[1] = xPreT2[1] + yPreT2[1] / ((trackT2[0]->z / trackT2[0]->x) * (1.0 / 0.003601 - trackT2[0]->y / trackT2[0]->z));
                    yPreT2[1] = yPreT2[1] / (1.0 - 0.003601 * trackT2[0]->y / trackT2[0]->z);
#endif
                    for (vector<hit3 *>::iterator iT2x2 = T2Layers[1].hits.begin(); iT2x2 != T2Layers[1].hits.end(); iT2x2++)
                    {
                        // Select a hit from the T2x2 layer
                        trackT2[1] = *iT2x2;
                        if (trackT2[1]->flag)
                            continue;
                        // Rough Estimation with Respect to the Origin
                        xPreT2[1] = trackT2[0]->x * (trackT2[1]->z / trackT2[0]->z);
                        if (conf_correctionT2)
                        {
                            if (xPreT2[1] > 0)
                                xPreT2[1] += correctionT2[iCase];
                            else
                                xPreT2[1] -= correctionT2[iCase];
                        }
                        yPreT2[1] = trackT2[0]->y * (trackT2[1]->z / trackT2[0]->z);
                        // Prediction Error
                        xErrT2[1] = xPreT2[1] - trackT2[1]->x;
                        yErrT2[1] = yPreT2[1] - trackT2[1]->y;

                        if ((trackT2[0]->motherTrack != NULL) and ((trackT2[0]->motherTrack->hitT2Set[1] != NULL) and (trackT2[0]->motherTrack->hitT2Set[1] == trackT2[1])))
                            match[0] = true;
                        else
                            match[0] = false;

                        if (match[0])
                        {
                            h2ErrT2x2Good->Fill(xErrT2[1], yErrT2[1]);
                            h1ErrT2x2xGood->Fill(TMath::Abs(xErrT2[1]));
                            h1ErrT2x2yGood->Fill(TMath::Abs(yErrT2[1]));
                        }
                        else
                        {
                            h2ErrT2x2Bad->Fill(xErrT2[1], yErrT2[1]);
                            h1ErrT2x2xBad->Fill(TMath::Abs(xErrT2[1]));
                            h1ErrT2x2yBad->Fill(TMath::Abs(yErrT2[1]));
                        }

                        // First Selection use window 0
                        if (TMath::Abs(yErrT2[1]) > DYT2[iCase])
                            continue;
                        else if (TMath::Abs(xErrT2[1]) > DXT2[iCase])
                            continue;

                        // Use hit T2x1 and hit T2x2, Linear Fit to predict hit T3x1
                        // yPreT3[0] = trackT2[1]->y * 685.1 / 210.1 - trackT2[0]->y * 475 / 210.1;
                        xPreT3[0] = trackT2[1]->x * 685.1 / 210.1 - trackT2[0]->x * 475 / 210.1;
                        yPreT3[0] = trackT2[1]->y * 752.2 / 210.1 - trackT2[0]->y * 542.1 / 210.1;
                        yPreT3[1] = trackT2[1]->y * 822.2 / 210.1 - trackT2[0]->y * 612.1 / 210.1;
#ifdef modfication
                        xPreT3[0] = xPreT3[0] + yPreT3[0] * ((trackT2[1]->x - trackT2[0]->x) / 210.1) / ((1.0 / 0.003601) - (trackT2[1]->y - trackT2[0]->y) / 210.1);
                        yPreT3[0] = yPreT3[0] / (1 - 0.003601 * (trackT2[1]->y - trackT2[0]->y) / 210.1);
#endif

                        if (yPreT3[0] > 0)
                        {
                            T3Layers[0] = &T3LayersUp[0];
                            T3Layers[3] = &T3LayersUp[3];
                        }
                        else
                        {
                            T3Layers[0] = &T3LayersDown[0];
                            T3Layers[3] = &T3LayersDown[3];
                        }
                        for (vector<hit1 *>::iterator iT3x1 = T3Layers[0]->hits.begin(); iT3x1 != T3Layers[0]->hits.end(); iT3x1++)
                        {
                            // Select a hit from the T3x1 layer
                            trackT3[0] = *iT3x1;
                            if (trackT3[0]->flag)
                                continue;

                            // Prediction Error
                            sErrT3[0] = TMath::Abs(xPreT3[0] - trackT3[0]->s);

                            if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT3Set[0] != NULL) and (trackT2[0]->motherTrack->hitT3Set[0] == trackT3[0]))
                                match[1] = true;
                            else
                                match[1] = false;

                            if (match[0] and match[1])
                            {
                                h1ErrT3x1Good->Fill(sErrT3[0]);
                            }
                            else if (match[0] and !match[1])
                            {
                                h1ErrT3x1Bad->Fill(sErrT3[0]);
                            }

                            if (sErrT3[0] > DX1T3[iCase])
                                continue;

                            xPreT3[3] = trackT3[0]->s * 895.0 / 685.1 - trackT2[0]->x * 209.9 / 685.1;

                            for (Int_t iT3x2 = 0; iT3x2 < T3Layers[3]->hits.size(); iT3x2++)
                            {
                                // Select a hit from the T3x2 layer
                                trackT3[3] = T3Layers[3]->hits[iT3x2];
                                if (trackT3[3]->flag)
                                    continue;
                                // Prediction Error
                                sErrT3[3] = TMath::Abs(xPreT3[3] - trackT3[3]->s);

                                if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT3Set[3] != NULL) and (trackT2[0]->motherTrack->hitT3Set[3] == trackT3[3]))
                                    match[2] = true;
                                else
                                    match[2] = false;

                                if (match[0] and match[1] and match[2])
                                {
                                    h1ErrT3x2Good->Fill(sErrT3[3]);
                                }
                                else if (match[0] and match[1] and !match[2])
                                {
                                    h1ErrT3x2Bad->Fill(sErrT3[3]);
                                }

                                if (sErrT3[3] > DX2T3[iCase])
                                    continue;

                                // Prediction of y on u/v of T3

                                // Prediction of x on u/v of T3
                                xPreT3[1] = trackT3[3]->s * 67.1 / 209.9 + trackT3[0]->s * 142.8 / 209.9;
                                xPreT3[2] = trackT3[3]->s * 137.1 / 209.9 + trackT3[0]->s * 72.8 / 209.9;

                                // Convert Prediction in y to u and v
                                uPreT3 = xPreT3[1] * cos(5.0 / 180) + yPreT3[0] * sin(5.0 / 180);
                                vPreT3 = xPreT3[2] * cos(5.0 / 180) - yPreT3[1] * sin(5.0 / 180);

                                if (uUpOrDown(xPreT3[1], yPreT3[1]))
                                {
                                    T3Layers[1] = &T3LayersUp[1];
                                }
                                else
                                {
                                    T3Layers[1] = &T3LayersDown[1];
                                }
                                for (Int_t iT3u1 = 0; iT3u1 < T3Layers[1]->hits.size(); iT3u1++)
                                {
                                    // Select a hit from the T3u1 layer
                                    trackT3[1] = T3Layers[1]->hits[iT3u1];
                                    if (trackT3[1]->flag)
                                        continue;
                                    // Prediction Error
                                    sErrT3[1] = uPreT3 - trackT3[1]->s;

                                    if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT3Set[1] != NULL) and (trackT2[0]->motherTrack->hitT3Set[1] == trackT3[1]))
                                        match[3] = true;
                                    else
                                        match[3] = false;

                                    if (match[0] and match[1] and match[2] and match[3])
                                    {
                                        h1ErrT3u1Good->Fill(TMath::Abs(sErrT3[1]));
                                    }
                                    else if (match[0] and match[1] and match[2] and !match[3])
                                    {
                                        h1ErrT3u1Bad->Fill(TMath::Abs(sErrT3[1]));
                                    }

                                    if (TMath::Abs(sErrT3[1]) > DUT3[iCase])
                                        continue;

                                    if (vUpOrDown(xPreT3[2], yPreT3[1]))
                                    {
                                        T3Layers[2] = &T3LayersUp[2];
                                    }
                                    else
                                    {
                                        T3Layers[2] = &T3LayersDown[2];
                                    }
                                    for (Int_t iT3v1 = 0; iT3v1 < T3Layers[2]->hits.size(); iT3v1++)
                                    {
                                        // Select a hit from the T3u1 layer
                                        trackT3[2] = T3Layers[2]->hits[iT3v1];
                                        if (trackT3[2]->flag)
                                            continue;
                                        // Prediction Error
                                        sErrT3[2] = TMath::Abs(vPreT3 - trackT3[2]->s);

                                        if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT3Set[2] != NULL) and (trackT2[0]->motherTrack->hitT3Set[2] == trackT3[2]))
                                            match[4] = true;
                                        else
                                            match[4] = false;

                                        if (match[0] and match[1] and match[2] and match[3] and match[4])
                                        {
                                            h1ErrT3v1Good->Fill(sErrT3[2]);
                                        }
                                        else if (match[0] and match[1] and match[2] and match[3] and !match[4])
                                        {
                                            h1ErrT3v1Bad->Fill(sErrT3[2]);
                                        }

                                        if (sErrT3[2] > DVT3[iCase])
                                            continue;

                                        // T1, T3 DONE!
                                        // Reconstruct the T3 coordinates in x-z and y-z to fit on T1

                                        float parabolaParams[4] = {0, 0, 0, 0};
                                        float fitResult[4] = {0.0, 0.0, 0.0, 0.0};

                                        float linFitX[4] = {trackT2[0]->x, trackT2[1]->x, xPreT3[0], xPreT3[3]};
                                        // float linFitY[4] = {trackT2[0]->y, trackT2[1]->y, yPreT3[0], yPreT3[1]};

                                        float linFitZ[4] = {trackT2[0]->z, trackT2[1]->z, 9192.9, 9402.8};
                                        // if (iCase == 0)
                                        // {
                                        //     linearFit(linFitX, linFitZ, parabolaParams, fitResult);
                                        //     linearCount++;
                                        // }
                                        // else
                                        {
                                            parabolaFit(linFitX, linFitZ, parabolaParams, fitResult);
                                            parabolaCount++;
                                        }
                                        // else
                                        // {
                                        //     cubicFit(linFitX, linFitZ, parabolaParams, fitResult);
                                        //     cubicCount++;
                                        // }

                                        // Predict x at x of T1
                                        xPreT1[0] = fitResult[0];
                                        // xPreT1[1] = fitResult[1];
                                        // xPreT1[2] = fitResult[2];
                                        xPreT1[3] = fitResult[3];
                                        // Predict y at U/V of T1
                                        yPreT1[0] = trackT2[0]->y * 822.9 / 210.1 - trackT2[1]->y * 612.8 / 210.1;
                                        yPreT1[1] = trackT2[0]->y * 752.9 / 210.1 - trackT2[1]->y * 542.8 / 210.1;

                                        if (yPreT1[1] > 0)
                                        {
                                            T1Layers[0] = &T1LayersUp[0];
                                            T1Layers[3] = &T1LayersUp[3];
                                        }
                                        else
                                        {
                                            T1Layers[0] = &T1LayersDown[0];
                                            T1Layers[3] = &T1LayersDown[3];
                                        }
                                        for (Int_t iT1x2 = 0; iT1x2 < T1Layers[3]->hits.size(); iT1x2++)
                                        {
                                            trackT1[3] = T1Layers[3]->hits[iT1x2];
                                            if (trackT1[3]->flag)
                                                continue;
                                            sErrT1[3] = TMath::Abs(xPreT1[3] - trackT1[3]->s);

                                            if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT1Set[3] != NULL) and (trackT2[0]->motherTrack->hitT1Set[3] == trackT1[3]))
                                                match[5] = true;
                                            else
                                                match[5] = false;

                                            if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5])
                                            {
                                                h1ErrT1x2Good->Fill(sErrT1[3]);
                                            }
                                            else if (match[0] and match[1] and match[2] and match[3] and match[4] and !match[5])
                                            {
                                                h1ErrT1x2Bad->Fill(sErrT1[3]);
                                            }

                                            if (sErrT1[3] > DX2T1[iCase])
                                                continue;

                                            for (Int_t iT1x1 = 0; iT1x1 < T1Layers[0]->hits.size(); iT1x1++)
                                            {
                                                trackT1[0] = T1Layers[0]->hits[iT1x1];
                                                if (trackT1[0]->flag)
                                                    continue;
                                                sErrT1[0] = TMath::Abs(xPreT1[0] - trackT1[0]->s);

                                                if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT1Set[0] != NULL) and (trackT2[0]->motherTrack->hitT1Set[0] == trackT1[0]))
                                                    match[6] = true;
                                                else
                                                    match[6] = false;

                                                if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6])
                                                {
                                                    h1ErrT1x1Good->Fill(sErrT1[0]);
                                                }
                                                else if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and !match[6])
                                                {
                                                    h1ErrT1x1Bad->Fill(sErrT1[0]);
                                                }

                                                if (sErrT1[0] > DX1T1[iCase])
                                                    continue;

                                                // Prediction of x on u/v of T1
                                                xPreT1[1] = trackT1[3]->s * 69.1 / 209.9 + trackT1[0]->s * 140.8 / 209.9;
                                                xPreT1[2] = trackT1[3]->s * 139.1 / 209.9 + trackT1[0]->s * 70.8 / 209.9;

                                                // Convert Prediction in y to u and v
                                                uPreT1 = xPreT1[1] * cos(5.0 / 180) + yPreT1[0] * sin(5.0 / 180);
                                                vPreT1 = xPreT1[2] * cos(5.0 / 180) - yPreT1[1] * sin(5.0 / 180);

                                                if (uUpOrDown(xPreT1[1], yPreT1[0]))
                                                {
                                                    T1Layers[1] = &T1LayersUp[1];
                                                }
                                                else
                                                {
                                                    T1Layers[1] = &T1LayersDown[1];
                                                }
                                                for (Int_t iT1u1 = 0; iT1u1 < T1Layers[1]->hits.size(); iT1u1++)
                                                {
                                                    trackT1[1] = T1Layers[1]->hits[iT1u1];
                                                    if (trackT1[1]->flag)
                                                        continue;
                                                    sErrT1[1] = TMath::Abs(uPreT1 - trackT1[1]->s);

                                                    if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT1Set[1] != NULL) and (trackT2[0]->motherTrack->hitT1Set[1] == trackT1[1]))
                                                        match[7] = true;
                                                    else
                                                        match[7] = false;

                                                    if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6] and match[7])
                                                    {
                                                        h1ErrT1u1Good->Fill(sErrT1[1]);
                                                    }
                                                    else if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6] and !match[7])
                                                    {
                                                        h1ErrT1u1Bad->Fill(sErrT1[1]);
                                                    }

                                                    if (sErrT1[1] > DUT1[iCase])
                                                        continue;

                                                    if (vUpOrDown(xPreT1[2], yPreT1[1]))
                                                    {
                                                        T1Layers[2] = &T1LayersUp[2];
                                                    }
                                                    else
                                                    {
                                                        T1Layers[2] = &T1LayersDown[2];
                                                    }
                                                    for (Int_t iT1v1 = 0; iT1v1 < T1Layers[2]->hits.size(); iT1v1++)
                                                    {
                                                        trackT1[2] = T1Layers[2]->hits[iT1v1];
                                                        if (trackT1[2]->flag)
                                                            continue;
                                                        sErrT1[2] = TMath::Abs(vPreT1 - trackT1[2]->s);

                                                        if ((trackT2[0]->motherTrack != NULL) and (trackT2[0]->motherTrack->hitT1Set[2] != NULL) and (trackT2[0]->motherTrack->hitT1Set[2] == trackT1[2]))
                                                            match[8] = true;
                                                        else
                                                            match[8] = false;

                                                        if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6] and match[7] and match[8])
                                                        {
                                                            h1ErrT1v1Good->Fill(sErrT1[2]);
                                                        }
                                                        else if (match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6] and match[7] and !match[8])
                                                        {
                                                            h1ErrT1v1Bad->Fill(sErrT1[2]);
                                                        }

                                                        if (sErrT1[2] > DVT1[iCase])
                                                            continue;

                                                        // Form the Track Candidate
                                                        // cout << "Form the Track Candidate" << endl;
                                                        track *trackCandidate = new track(trackT1, trackT2, trackT3, match[0] and match[1] and match[2] and match[3] and match[4] and match[5] and match[6] and match[7] and match[8]);
                                                        caseTracks.push_back(trackCandidate);
                                                        allTracks.push_back(trackCandidate);
                                                    } // Loop on T1v1 Layer
                                                }     // Loop on T1u1 Layer
                                            }         // Loop on T1x1 Layer
                                        }             // Loop on T1x2 Layer
                                    }                 // Loop on T3v1 Layer
                                }                     // Loop on T3u1 Layer
                            }                         // Loop on T3x2 Layer
                        }                             // Loop on T3x1 Layer
                    }                                 // Loop on T2x2 Layer
                }                                     // Loop on T2x1 Layer

                // Local Purge
                int usedSum;
                vector<track *>::iterator _t;
                for (_t = caseTracks.begin(); _t != caseTracks.end(); _t++)
                {
                    if (!conf_alone)
                    {
                        usedSum = (*_t)->GetUsedSum();
                        // Purge the tracks
                        if (usedSum > 10)
                        {
                            // Check every hit, if it is used more that once, purge it
                            for (Int_t j = 0; j < 4; j++)
                            {
                                if ((*_t)->trackAccept == false)
                                {
                                    break;
                                }

                                if ((*_t)->hitT1Set[j]->used != 1)
                                {
                                    // Loop on every other related tracks.
                                    for (auto p : (*_t)->hitT1Set[j]->relatedTracks)
                                    {
                                        // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                        if (((*p) == (**_t)) < nShareHits[iCase] or p->trackAccept == false)
                                        {
                                            continue;
                                        }

                                        // If the chi2 of this track is larger than other clone track, reject this track and break.
                                        if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                        {
                                            (*_t)->trackAccept = false;
                                            // See how many friendly fire
                                            h1_rejectPurity->Fill((*_t)->purity);
                                            nRejectedTracks++;
                                            if ((*_t)->trackMatch)
                                            {
                                                // Record the friendly fire rate
                                                (*_t)->motherTrack->bFriendlyFired = true;
                                                nFriendlyFire++;
                                            }
                                            break;
                                        }
                                    }
                                }

                                if (j < 2 and (*_t)->hitT2Set[j]->used != 1)
                                {
                                    if ((*_t)->trackAccept == false)
                                    {
                                        break;
                                    }
                                    // Loop on every other related tracks.
                                    for (auto p : (*_t)->hitT2Set[j]->relatedTracks)
                                    {
                                        // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                        if (((*p) == (**_t)) < nShareHits[iCase] or p->trackAccept == false)
                                        {
                                            continue;
                                        }

                                        // If the chi2 of this track is larger than other clone track, reject this track and break.
                                        if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                        {
                                            (*_t)->trackAccept = false;
                                            // See how many friendly fire
                                            h1_rejectPurity->Fill((*_t)->purity);
                                            nRejectedTracks++;
                                            if ((*_t)->trackMatch)
                                            {
                                                // Record the friendly fire rate
                                                (*_t)->motherTrack->bFriendlyFired = true;
                                                nFriendlyFire++;
                                            }
                                            break;
                                        }
                                    }
                                }

                                if ((*_t)->hitT3Set[j]->used != 1)
                                {
                                    if ((*_t)->trackAccept == false)
                                    {
                                        break;
                                    }
                                    // Loop on every other related tracks.
                                    for (auto p : (*_t)->hitT3Set[j]->relatedTracks)
                                    {
                                        // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                        if (((*p) == (**_t)) < nShareHits[iCase] or p->trackAccept == false)
                                        {
                                            continue;
                                        }

                                        // If the chi2 of this track is larger than other clone track, reject this track and break.
                                        if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                        {
                                            (*_t)->trackAccept = false;
                                            // See how many friendly fire
                                            h1_rejectPurity->Fill((*_t)->purity);
                                            nRejectedTracks++;
                                            if ((*_t)->trackMatch)
                                            {
                                                // Record the friendly fire rate
                                                (*_t)->motherTrack->bFriendlyFired = true;
                                                nFriendlyFire++;
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            // TODO How to reject usedSum = 6 hits fake tracks ? (*_t)->trackAccept = true;
                        }
                        // Take care of the hits and counters of the rejected tracks
                        if ((*_t)->trackAccept == false)
                        {
                            for (Int_t j = 0; j < 4; j++)
                            {
                                (*_t)->hitT1Set[j]->used--;
                                (*_t)->hitT3Set[j]->used--;
                                if (j < 2)
                                {
                                    (*_t)->hitT2Set[j]->used--;
                                }
                            }
                        }
                        else
                        {

                            for (Int_t j = 0; j < 4; j++)
                            {
                                (*_t)->hitT1Set[j]->flag = true;
                                (*_t)->hitT3Set[j]->flag = true;
                                if (j < 2)
                                {
                                    (*_t)->hitT2Set[j]->flag = true;
                                }
                            }
                        }
                    }
                }
                caseTracks.clear();
            }
            int usedSum;
            vector<track *>::iterator _t;
            for (_t = allTracks.begin(); _t != allTracks.end(); _t++)
            {
                usedSum = (*_t)->GetUsedSum();
                if (!conf_alone)
                {
                    // If this track is already rejected, just continue
                    if ((*_t)->trackAccept == false)
                    {
                        continue;
                    }
                    // Purge the tracks
                    if (usedSum > 10)
                    {
                        // Check every hit, if it is used more that once, purge it
                        for (Int_t j = 0; j < 4; j++)
                        {
                            if ((*_t)->trackAccept == false)
                            {
                                break;
                            }

                            if ((*_t)->hitT1Set[j]->used != 1)
                            {
                                // Loop on every other related tracks.
                                for (auto p : (*_t)->hitT1Set[j]->relatedTracks)
                                {
                                    // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                    if (((*p) == (**_t)) < nGlobalShareHits or p->trackAccept == false)
                                    {
                                        continue;
                                    }

                                    // If the chi2 of this track is larger than other clone track, reject this track and break.
                                    if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                    {
                                        (*_t)->trackAccept = false;
                                        // See how many friendly fire
                                        h1_rejectPurity->Fill((*_t)->purity);
                                        nRejectedTracks++;
                                        if ((*_t)->trackMatch)
                                        {
                                            // Record the friendly fire rate
                                            (*_t)->motherTrack->bFriendlyFired = true;
                                            nFriendlyFire++;
                                        }
                                        break;
                                    }
                                }
                            }

                            if (j < 2 and (*_t)->hitT2Set[j]->used != 1)
                            {
                                if ((*_t)->trackAccept == false)
                                {
                                    break;
                                }
                                // Loop on every other related tracks.
                                for (auto p : (*_t)->hitT2Set[j]->relatedTracks)
                                {
                                    // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                    if (((*p) == (**_t)) < nGlobalShareHits or p->trackAccept == false)
                                    {
                                        continue;
                                    }

                                    // If the chi2 of this track is larger than other clone track, reject this track and break.
                                    if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                    {
                                        (*_t)->trackAccept = false;
                                        // See how many friendly fire
                                        h1_rejectPurity->Fill((*_t)->purity);
                                        nRejectedTracks++;
                                        if ((*_t)->trackMatch)
                                        {
                                            // Record the friendly fire rate
                                            (*_t)->motherTrack->bFriendlyFired = true;
                                            nFriendlyFire++;
                                        }
                                        break;
                                    }
                                }
                            }

                            if ((*_t)->hitT3Set[j]->used != 1)
                            {
                                if ((*_t)->trackAccept == false)
                                {
                                    break;
                                }
                                // Loop on every other related tracks.
                                for (auto p : (*_t)->hitT3Set[j]->relatedTracks)
                                {
                                    // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                    if (((*p) == (**_t)) < nGlobalShareHits or p->trackAccept == false)
                                    {
                                        continue;
                                    }

                                    // If the chi2 of this track is larger than other clone track, reject this track and break.
                                    if ((*_t)->fChi2X > p->fChi2X or (*_t)->fChi2Y > p->fChi2Y)
                                    {
                                        (*_t)->trackAccept = false;
                                        // See how many friendly fire
                                        h1_rejectPurity->Fill((*_t)->purity);
                                        nRejectedTracks++;
                                        if ((*_t)->trackMatch)
                                        {
                                            // Record the friendly fire rate
                                            (*_t)->motherTrack->bFriendlyFired = true;
                                            nFriendlyFire++;
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // TODO How to reject usedSum = 6 hits fake tracks ? (*_t)->trackAccept = true;
                    }
                    // Take care of the hits and counters of the rejected tracks
                    if ((*_t)->trackAccept == false)
                    {
                        for (Int_t j = 0; j < 4; j++)
                        {
                            (*_t)->hitT1Set[j]->used--;
                            (*_t)->hitT3Set[j]->used--;
                            if (j < 2)
                            {
                                (*_t)->hitT2Set[j]->used--;
                            }
                        }
                    }
                }
                // After purging, a track is either rejected or accepted (default)
                // No purging, all tracks are accepted.
            } // Loop on all the tracks, purge and tally the result

            // Tally The Result

            // int usedSum;

            for (_t = allTracks.begin(); _t != allTracks.end(); _t++)
            {
                // Tally the result
                if ((*_t)->trackAccept)
                {
                    nTotalReconstructed++;
                    if ((*_t)->trackMatch)
                    {
                        p = (*_t)->hitT2Set[0]->hP;
                        eta = (*_t)->hitT2Set[0]->hEta; // log is in e base in C++
                        if ((*_t)->theTrueHit < 4)
                        {
                            if ((*_t)->hitT1Set[(*_t)->theTrueHit]->motherTrack->particleMatched)
                            {
                                nClone++;
                                delete *_t;
                                continue;
                            }
                            else
                            {
                                (*_t)->hitT1Set[(*_t)->theTrueHit]->motherTrack->particleMatched = true;
                            }
                        }
                        else if ((*_t)->theTrueHit < 6)
                        {
                            if ((*_t)->hitT2Set[(*_t)->theTrueHit - 4]->motherTrack->particleMatched)
                            {
                                nClone++;
                                delete *_t;
                                continue;
                            }
                            else
                            {
                                (*_t)->hitT2Set[(*_t)->theTrueHit - 4]->motherTrack->particleMatched = true;
                            }
                        }
                        else
                        {
                            if ((*_t)->hitT3Set[(*_t)->theTrueHit - 6]->motherTrack->particleMatched)
                            {
                                nClone++;
                                delete *_t;
                                continue;
                            }
                            else
                                (*_t)->hitT3Set[(*_t)->theTrueHit - 6]->motherTrack->particleMatched = true;
                        }

                        if (eta > 2 and eta < 5)
                        {
                            if (p < 5000)
                            {
                                nTrueReconstructedP0P5Eta1++;
                            }
                            else if (p >= 5000 and p < 20000)
                            {
                                nTrueReconstructedP5P20Eta1++;
                            }
                            else if (p >= 20000)
                            {
                                nTrueReconstructedP20Eta1++;
                            }
                            nTrueReconstructed++;
                        }
                    }
                    else
                    {
                        nFakeReconstructed++;
                        if ((*_t)->GetUsedSum() == 10)
                            nPerfectFake++;

                        h1_hitPurity->Fill(((*_t)->purity));
                        for (int j = 0; j < 9; j++)
                        {
                            if (!(*_t)->match[j])
                                h1_whereFake->Fill(j);
                        }
                    }
                }
                delete *_t;
            }

            allTracks.clear();

            for (auto p : motherTracks)
            {
                if (!p->particleMatched)
                {
                    h1_notFound->Fill(p->nTrueHits);
                    for (Int_t j = 0; j < 4; j++)
                    {
                        if (!p->hitT1Set[j])
                            h1_whereNotFound->Fill(j);
                        if (j < 2 and !p->hitT2Set[j])
                            h1_whereNotFound->Fill(j + 4);
                        if (!p->hitT3Set[j])
                            h1_whereNotFound->Fill(j + 6);
                    }
                    notFound++;
                    if (p->bFriendlyFired)
                        nNotFoundAsFired++;
                }
                delete p;
            }
            motherTracks.clear();

            // Clear the layers after finish the task
            for (Int_t j = 0; j < 4; j++)
            {
                // delete every hit
                for (auto p : T1LayersUp[j].hits)
                {
                    delete p;
                }
                T1LayersUp[j].hits.clear();
                for (auto p : T1LayersDown[j].hits)
                {
                    delete p;
                }
                T1LayersDown[j].hits.clear();
                for (auto p : T3LayersUp[j].hits)
                {
                    delete p;
                }
                T3LayersUp[j].hits.clear();
                for (auto p : T3LayersUp[j].hits)
                {
                    delete p;
                }
                T3LayersDown[j].hits.clear();

                if (j < 2)
                {
                    for (auto p : T2Layers[j].hits)
                    {
                        delete p;
                    }
                    T2Layers[j].hits.clear();
                }
                // delete the layer
            }
            eventCounter++;
        } // Run the reconstruction in the end of an event or in the end of the readout
        else
        {
            // Already got an entry before else
            int bT1XSum = 0;
            int bT1UVSum = 0;
            int bT2Sum = 0;
            int bT3XSum = 0;
            int bT3UVSum = 0;

            hit1 *_hitT1Set[4] = {NULL, NULL, NULL, NULL};
            hit3 *_hitT2Set[2] = {NULL, NULL};
            hit1 *_hitT3Set[4] = {NULL, NULL, NULL, NULL};

            // Don't read this track if not satisfies the cut
            if (p > conf_max_p or p < conf_min_p or ((eta < 2 or eta > 5) and (conf_eta)))
                continue;

            if (rTracks)
            {
                // Create and Count the Hits

                for (Int_t j = 0; j < 10; j++)
                {
                    if (j < 4)
                    {
                        if (not(hitT1up[j] == 0))
                        {
                            _hitT1Set[j] = new hit1(hitT1up[j]);
                            _hitT1Set[j]->hP = p;
                            _hitT1Set[j]->hEta = eta;
                            _hitT1Set[j]->eventID = eventNumber;
                            T1LayersUp[j].hits.push_back(_hitT1Set[j]);
                            // cout << "pushbackT1" << endl;
                            if (j == 0 or j == 3)
                            {
                                bT1XSum++;
                            }
                            else
                            {
                                bT1UVSum++;
                            }
                        }
                        if (not(hitT1down[j] == 0))
                        {
                            _hitT1Set[j] = new hit1(hitT1down[j]);
                            _hitT1Set[j]->hP = p;
                            _hitT1Set[j]->hEta = eta;
                            _hitT1Set[j]->eventID = eventNumber;
                            T1LayersDown[j].hits.push_back(_hitT1Set[j]);
                            // cout << "pushbackT1" << endl;
                            if (j == 0 or j == 3)
                            {
                                bT1XSum++;
                            }
                            else
                            {
                                bT1UVSum++;
                            }
                        }
                    }
                    else if (j == 4 or j == 5)
                    {
                        if (not(hitT2[j - 4][0] == 0 and hitT2[j - 4][1] == 0 and hitT2[j - 4][2] == 0))
                        {
                            _hitT2Set[j - 4] = new hit3(hitT2[j - 4][0], hitT2[j - 4][1], hitT2[j - 4][2]);
                            _hitT2Set[j - 4]->hP = p;
                            _hitT2Set[j - 4]->hEta = eta;
                            _hitT2Set[j - 4]->eventID = eventNumber;
                            T2Layers[j - 4].hits.push_back(_hitT2Set[j - 4]);
                            // cout << "pushbackT2" << endl;
                            bT2Sum++;
                        }
                    }
                    else if (j > 5)
                    {
                        if (not(hitT3up[j - 6] == 0))
                        {
                            _hitT3Set[j - 6] = new hit1(hitT3up[j - 6]);
                            _hitT3Set[j - 6]->hP = p;
                            _hitT3Set[j - 6]->hEta = eta;
                            _hitT3Set[j - 6]->eventID = eventNumber;
                            T3LayersUp[j - 6].hits.push_back(_hitT3Set[j - 6]);
                            // cout << "pushbackT1" << endl;
                            if (j == 6 or j == 9)
                            {
                                bT3XSum++;
                            }
                            else
                            {
                                bT3UVSum++;
                            }
                        }
                        if (not(hitT3down[j - 6] == 0))
                        {
                            _hitT3Set[j - 6] = new hit1(hitT3down[j - 6]);
                            _hitT3Set[j - 6]->hP = p;
                            _hitT3Set[j - 6]->hEta = eta;
                            _hitT3Set[j - 6]->eventID = eventNumber;
                            T3LayersDown[j - 6].hits.push_back(_hitT3Set[j - 6]);
                            // cout << "pushbackT1" << endl;
                            if (j == 6 or j == 9)
                            {
                                bT3XSum++;
                            }
                            else
                            {
                                bT3UVSum++;
                            }
                        }
                    }
                    // If no hit, just continue
                }
                if ((bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum) == 0)
                {
                    continue;
                }
            }
            else
            {
                cout << "rTracks not found" << endl;
            }
            // Define Reconstrutible as At least One Complete hit in each Station.

            // if (not ((bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum )== 10))
            if (not(bT1XSum and bT1UVSum and bT2Sum and bT3XSum and bT3UVSum and ((bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum) == 10)))
            // if (not(bT1XSum and bT1UVSum and bT2Sum and bT3XSum and bT3UVSum and ((bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum) <= 10)))
            {
                continue;
            }

            // if( (bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum ) >10)
            // {
            //     cout << i << endl;
            //     cout << (bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum ) << endl;
            //     cout<< bT1XSum <<" " <<bT1UVSum  <<" " << bT2Sum  <<" " << bT3XSum  <<" " << bT3UVSum  << endl;
            //     cout << hitT3down[0] << hitT3down[1] << hitT3down[2] << hitT3down[3] << endl;
            //     cout << hitT3up[0] << hitT3up[1] << hitT3up[2] << hitT3up[3] << endl;
            // }

            // cout << "bT1XSum:" << bT1XSum << " bT1UVSum:" << bT1UVSum << " bT2Sum:" << bT2Sum << " bT3Sum:" << bT3XSum << " bT3UVSum:" << bT3UVSum << endl;

            if (eta > 2 and eta < 5)
            {
                if (p < 5000)
                {
                    nReconstructableP0P5Eta1++;
                }
                else if (p > 5000 and p < 20000)
                {
                    nReconstructableP5P20Eta1++; // Always count tracks with momentum > 5000 and in eta acceptance
                }
                else if (p > 20000)
                {
                    nReconstructableP20Eta1++; // Always count tracks with momentum > 5000 and in eta acceptance
                }
            }
            if (p > conf_min_p and p < conf_max_p and ((eta > 2 and eta < 5) or (conf_eta == 0)))
            {
                nReconstructible++;
                if ((bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum) == 10)
                    nReconstructable10Hit++;

                if (particleID == -11 or particleID == 11)
                    nReconstructableElectron++;
            }

            // Set the motherTrack of hits
            track *_motherTrack = new track(_hitT1Set, _hitT2Set, _hitT3Set, bT1XSum + bT1UVSum + bT2Sum + bT3XSum + bT3UVSum);
            motherTracks.push_back(_motherTrack);
            for (Int_t j = 0; j < 4; j++)
            {
                if (_hitT1Set[j] != NULL)
                {
                    _hitT1Set[j]->motherTrack = _motherTrack;
                }
                if (_hitT3Set[j] != NULL)
                {
                    _hitT3Set[j]->motherTrack = _motherTrack;
                }
                if (j < 2)
                {
                    if (_hitT2Set[j] != NULL)
                    {
                        _hitT2Set[j]->motherTrack = _motherTrack;
                    }
                }
            }
        } //continue reading events and filling the information for the layers
    }

    fakeOut.close();
    goodOut.close();
    logOut << "---------- Total ----------" << endl;
    logOut << "# of Events:" << eventCounter << endl;
    logOut << "Time/Event (ms):" << 0.001 * chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - timingStart).count() / eventCounter << endl;
    logOut << "Electron fraction (reconstructable): " << (nReconstructableElectron * 1. / nReconstructible) << " ( " << nReconstructableElectron << " / " << nReconstructible << " ) " << endl;
    if (nReconstructableElectron > 0)
        logOut << "Electron reconstructed rate: " << (nTrueReconstructedElectron * 1. / nReconstructableElectron) << " ( " << nTrueReconstructedElectron << " / " << nReconstructableElectron << " ) " << endl;

    if (nReconstructableP0P5Eta1 > 0)
        logOut << "Eff. (p0p5eta1): " << (nTrueReconstructedP0P5Eta1 * 1. / nReconstructableP0P5Eta1) << " ( " << nTrueReconstructedP0P5Eta1 << " / " << nReconstructableP0P5Eta1 << " ) " << endl;
    if (nReconstructableP5P20Eta1 > 0)
        logOut << "Eff. (p5p20eta1): " << (nTrueReconstructedP5P20Eta1 * 1. / nReconstructableP5P20Eta1) << " ( " << nTrueReconstructedP5P20Eta1 << " / " << nReconstructableP5P20Eta1 << " ) " << endl;
    if (nReconstructableP20Eta1 > 0)
        logOut << "Eff. (p20eta1): " << (nTrueReconstructedP20Eta1 * 1. / nReconstructableP20Eta1) << " ( " << nTrueReconstructedP20Eta1 << " / " << nReconstructableP20Eta1 << " ) " << endl;

    logOut << "Eff.: " << (nTrueReconstructed * 1.0 / nReconstructible) << " ( " << nTrueReconstructed << " / " << nReconstructible << " ) " << endl;
    logOut << "nReconstructable10Hit: " << nReconstructable10Hit << endl;
    logOut << "Ghost rate(70\% hits true): " << (nTotalReconstructed - nTrueReconstructed) * 1.0 / nTotalReconstructed << " ( " << (nTotalReconstructed - nTrueReconstructed) << " ) " << endl;
    logOut << "Ghost rate (70\% hits true alternative. counting): " << (nFakeReconstructed)*1.0 / nTotalReconstructed << " ( " << nFakeReconstructed << " ) " << endl;
    logOut << "nPerfectFake: " << nPerfectFake << endl;

    logOut << "Friendly Fire Rate:" << nFriendlyFire * 1.0 / nRejectedTracks << " ( " << nFriendlyFire << "/" << nRejectedTracks << " ) " << endl;
    logOut << "nClones:" << nClone << endl;
    logOut << "notFound:" << notFound << endl;
    logOut << "nNotFoundAsFired" << nNotFoundAsFired << endl;
    logOut << "Parabola Count:" << parabolaCount << endl;
    logOut << "Linear Count:" << linearCount << endl;
    logOut << "Cubic Count:" << cubicCount << endl;
    logOut << "---------- SWs ----------" << endl;
    logOut << "DX1T1 " << DX1T1[0] << " " << DX1T1[1] << " " << DX1T1[2] << " " << DX1T1[3] << " " << DX1T1[4] << endl;
    logOut << "DX2T1 " << DX2T1[0] << " " << DX2T1[1] << " " << DX2T1[2] << " " << DX2T1[3] << " " << DX2T1[4] << endl;
    logOut << "DUT1 " << DUT1[0] << " " << DUT1[1] << " " << DUT1[2] << " " << DUT1[3] << " " << DUT1[4] << endl;
    logOut << "DVT1 " << DVT1[0] << " " << DVT1[1] << " " << DVT1[2] << " " << DVT1[3] << " " << DVT1[4] << endl;
    logOut << "DYT2 " << DYT2[0] << " " << DYT2[1] << " " << DYT2[2] << " " << DYT2[3] << " " << DYT2[4] << endl;
    logOut << "DXT2 " << DXT2[0] << " " << DXT2[1] << " " << DXT2[2] << " " << DXT2[3] << " " << DXT2[4] << endl;
    logOut << "DX1T3 " << DX1T3[0] << " " << DX1T3[1] << " " << DX1T3[2] << " " << DX1T3[3] << " " << DX1T3[4] << endl;
    logOut << "DX2T3 " << DX2T3[0] << " " << DX2T3[1] << " " << DX2T3[2] << " " << DX2T3[3] << " " << DX2T3[4] << endl;
    logOut << "DUT3 " << DUT3[0] << " " << DUT3[1] << " " << DUT3[2] << " " << DUT3[3] << " " << DUT3[4] << endl;
    logOut << "DVT3 " << DVT3[0] << " " << DVT3[1] << " " << DUT3[2] << " " << DUT3[3] << " " << DVT3[4] << endl;

    c1->SetLogy(0);

    h1_hitPurity->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1_hitPurity.png");

    h1_rejectPurity->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1_rejectPurity.png");

    h1_whereFake->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1_whereFake.png");

    h1_notFound->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1_notFound.png");

    h1_whereNotFound->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1_whereNotFound.png");

    h1ErrT1x1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1x1Good.png");
    h1ErrT1x1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1x1Bad.png");

    h1ErrT1u1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1u1Good.png");
    h1ErrT1u1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1u1Bad.png");

    h1ErrT1v1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1v1Good.png");
    h1ErrT1v1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1v1Bad.png");

    h1ErrT1x2Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1x2Good.png");
    h1ErrT1x2Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT1x2Bad.png");
    // T2 Figures
    h1ErrT2x2xGood->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT2x2xGood.png");
    h1ErrT2x2xBad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT2x2xBad.png");
    h1ErrT2x2yGood->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT2x2yGood.png");
    h1ErrT2x2yBad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT2x2yBad.png");

    h2ErrT2x2Good->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "/h2/" + "h2ErrT2x2Good.png");
    h2ErrT2x2Bad->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "/h2/" + "h2ErrT2x2Bad.png");

    // T3 Figures
    h1ErrT3x1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3x1Good.png");
    h1ErrT3x1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3x1bad.png");

    h1ErrT3u1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3u1Good.png");
    h1ErrT3u1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3u1Bad.png");

    h1ErrT3v1Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3v1Good.png");
    h1ErrT3v1Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3v1Bad.png");

    h1ErrT3x2Good->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3x2Good.png");
    h1ErrT3x2Bad->Draw();
    c1->SaveAs(outputplotsfolder + "/h1/" + "h1ErrT3x2Bad.png");

    c1->SetLogy(1);

    make_inefficiency_window_plot(c1, h1ErrT1x1Good, h1ErrT1x1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT1u1Good, h1ErrT1u1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT1v1Good, h1ErrT1v1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT1x2Good, h1ErrT1x2Bad, outputplotsfolder);

    make_inefficiency_window_plot(c1, h1ErrT2x2xGood, h1ErrT2x2xBad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT2x2yGood, h1ErrT2x2yBad, outputplotsfolder);

    make_inefficiency_window_plot(c1, h1ErrT3x1Good, h1ErrT3x1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT3u1Good, h1ErrT3u1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT3v1Good, h1ErrT3v1Bad, outputplotsfolder);
    make_inefficiency_window_plot(c1, h1ErrT3x2Good, h1ErrT3x2Bad, outputplotsfolder);
}
