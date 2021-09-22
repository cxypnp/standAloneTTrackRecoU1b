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

#define modfication
#define layerStudy
class track;
class hit
{
public:
    float x, y, z;
    float eventID = -1;
    double hP = -1;
    double hPz = -1;
    int particleID = 0;
    int used = 0;
    bool flag = false;

    hit *preLayerHit = NULL;
    hit *nextLayerHit = NULL;
    std::vector<track *> relatedTracks;
    std::shared_ptr<track> motherTrack = NULL;

    bool specialFlag = false;
    hit(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
        preLayerHit = NULL;
        nextLayerHit = NULL;
    }
    hit(float _eventID, float _hP, float _hPz, float _x, float _y, float _z)
    {
        hPz = _hPz;
        hP = _hP;
        eventID = _eventID;
        x = _x;
        y = _y;
        z = _z;
        preLayerHit = NULL;
        nextLayerHit = NULL;
    }
    hit(float _eventID, float _hP, float _hPz, int _particleID, float _x, float _y, float _z)
    {
        particleID = _particleID;
        hPz = _hPz;
        hP = _hP;
        eventID = _eventID;
        x = _x;
        y = _y;
        z = _z;
        preLayerHit = NULL;
        nextLayerHit = NULL;
    }

    ~hit()
    {
        motherTrack.reset();
    }

    void setPreLayer(hit *_hit)
    {
        preLayerHit = _hit;
    }

    void setNextLayer(hit *_hit)
    {
        nextLayerHit = _hit;
    }

    int whichLayer()
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

    bool operator==(const hit *anotherHit)
    {
        if ((anotherHit->x == this->x) and (anotherHit->y == this->y) and (anotherHit->z == this->z) and (anotherHit->hP == this->hP) and (anotherHit->hPz == this->hPz))
            return true;
        else
            return false;
    }

private:
};

class TLayer
{
public:
    std::vector<hit *> hits;
};

class track
{
public:
    std::vector<hit *> hits;

    float xPre[6];
    float yPre[6];

    float dRoughX = 0;
    float dRoughY = 0;

    float pChi2X = 0;
    float pChi2Y = 0;

    float fChi2X = 0;
    float fChi2Y = 0;

    float dFitX[6] = {0, 0, 0, 0, 0, 0};
    float dFitY[6] = {0, 0, 0, 0, 0, 0};

    int matchSum = 0;
    bool trackMatch;

    bool trackAccept = true;
    int nTrueHits = 0;
    int theTrueHit = 0;
    track *motherTrack = NULL;

    // Only for motherTracks
    bool particleMatched = false;

    track(hit **trackHit,
          float _xErr[6], float _yErr[6], bool _match[5])
    {

        for (Int_t j = 0; j < 6; j++)
        {
            if (_match[j] && j != 5)
            {
                matchSum++;
            }
            hits.push_back(trackHit[j]);
            // Get chi2X value, ignore the seed layer error
            if (j != 3)
            {
                pChi2X += _xErr[j] * _xErr[j];
                pChi2Y += _yErr[j] * _yErr[j];
            }
            trackHit[j]->used++;
            trackHit[j]->relatedTracks.push_back(this);
        }

        trackMatch = _match[0] and _match[1] and _match[2] and _match[3] and _match[4];

        if (not trackMatch)
        {
            for (Int_t j = 0; j < 6; j++)
            {
                if (nTrueHits < 3)
                {
                    if (((hits[j]->motherTrack != NULL) && *(hits[j]->motherTrack.get()) == *this) > nTrueHits)
                    {
                        nTrueHits = (*(hits[j]->motherTrack.get()) == *this);
                        motherTrack = hits[j]->motherTrack.get();
                        theTrueHit = j;
                    }
                }
                else
                    break;
            }
        }
        else
        {
            nTrueHits = 6;
            theTrueHit = 0;
        }

        float xFit[6], yFit[6], zFit[6];
        for (Int_t j = 0; j < 6; j++)
        {
            xFit[j] = hits[j]->x;
            yFit[j] = hits[j]->y;
            zFit[j] = hits[j]->z;
        }
        TGraph gy = TGraph(6, zFit, yFit);
        TFitResultPtr fpy = gy.Fit("pol1", "SQ");
        fChi2Y = fpy->Chi2();
        TGraph gx = TGraph(6, zFit, xFit);
        TFitResultPtr fpx = gx.Fit("pol3", "SQ");
        fChi2X = fpx->Chi2();
    }

    track(hit **trackHit)
    {
        for (Int_t j = 0; j < 6; j++)
        {
            hits.push_back(trackHit[j]);
        }
    }

    int GetUsedSum()
    {
        int usedSum = 0;
        for (Int_t j = 0; j < 6; j++)
        {
            usedSum = usedSum + hits[j]->used;
        }
        return usedSum;
    }

    int operator==(const track anotherTrack)
    {
        // Return the number of hits matched
        return int((anotherTrack.hits[0] == this->hits[0]) + (anotherTrack.hits[1] == this->hits[1]) + (anotherTrack.hits[2] == this->hits[2]) + (anotherTrack.hits[3] == this->hits[3]) + (anotherTrack.hits[4] == this->hits[4]) + (anotherTrack.hits[5] == this->hits[5]));
    }
};

bool IsInDeadArea(float x, float y)
{
    if (x < 100 and x > -100)
        if (y < 100 and y > -100)
            return false; // beam area

    if (fabs(remainder((x / 20), 1)) < 0.1 / 20)
    {
        if (fabs(remainder(((x - 20) / 20), 2)) < 2 * 0.1 / 20)
            return false;
        return true;
    }

    if (fabs(remainder((y / 20), 1)) < 0.1 / 20)
    {
        return true;
    }
    return false;
}

int parabolaFit(float *linFitX, float *linFitY, float *linFitZ, float *parabolaParams, float *fitResult)
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
    fitResult[1] = fx->Eval(8035.8);
    return true;
}

int cubicFit(float *linFitX, float *linFitY, float *linFitZ, float *cubicParams, float *fitResult)
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
    fitResult[1] = fx->Eval(8035.8);
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
    float threshold99 = 0;
    float threshold98 = 0;
    float threshold95 = 0;
    for (unsigned int i = 0; i < h1_good->GetXaxis()->GetNbins(); i++)
    {
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
    xlabel->DrawText(0.15, 0.25, Form("Threshold for Eff.(99.9\%): %.3f[mm]", threshold99));
    xlabel->DrawText(0.15, 0.20, Form("Threshold for Eff.(99.8\%): %.3f[mm]", threshold98));
    xlabel->DrawText(0.15, 0.15, Form("Threshold for Eff.(99.5\%): %.3f[mm]", threshold95));
    c1->SaveAs(outputplotsfolder + h1_good->GetName() + "_inefficiency.pdf");

    h1_good->GetYaxis()->SetTitle("%");
    h1_bad->GetYaxis()->SetTitle("%");
}

void patternReco(TString settings = "settings", TString inputfilestring = "inputFile.root", TString outputplotsfolder = "outputFolder/")
// void patternReco(TString settings = "Case1p5eff99", TString inputfilestring = "../DataSets/bigTestNoBlur.root", TString outputplotsfolder = "bigTestNoBlur/")
{
    // system("rm -rf " + outputplotsfolder);
    outputplotsfolder = outputplotsfolder + settings + '/';
    system("mkdir -p " + outputplotsfolder);
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
    float conf_alone = 0;

    // No p0p5 tracks is reconstructed in the first case
    // // Finished windows for MCTruth Hits
    // Float_t MAXDX[5][5] = {{20, 0.5, 0.5, 0.5, 0.5}, {40, 3, 0.8, 1.5, 2}, {80, 10, 3, 5, 5}, {120, 30, 5, 7, 7}, {150, 70, 6, 10, 10}};
    // Float_t MAXDY[5][5] = {{3, 0.5, 0.5, 0.5, 0.5}, {10, 0.8, 0.8, 0.8, 0.8}, {20, 2, 2, 2, 2}, {30, 3, 3, 3, 3}, {50, 15, 6, 5, 5}};

    // Finished windows for Blured Hits
    // Few p0p5 tracks is reconstructed in the first case
    Float_t MAXDX[5][5] = {{20, 4, 0.8, 2.4, 1.6}, {50, 4.8, 0.8, 1.6, 2}, {80, 10, 5, 5, 8}, {120, 30, 7, 14, 11}, {160, 80, 7, 14, 11}};
    Float_t MAXDY[5][5] = {{4.9, 12.5, 6.5, 11.1, 13}, {14, 12, 5.5, 11, 13.5}, {40, 18, 7, 14, 15}, {46, 20, 8, 15, 18}, {50, 20, 8, 15, 18}};

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

    TH1F *h1_dy_good_5 = new TH1F("h1_dy_good_5", "Prediction Inefficiency- Layers 5", 300, 0, 30);
    h1_dy_good_5->GetYaxis()->SetTitle("Inefficiency (%)");
    h1_dy_good_5->GetXaxis()->SetTitle("Y_{5Pred}-Y_{5} [mm]");
    TH1F *h1_dy_good_4 = new TH1F("h1_dy_good_4", "Prediction Inefficiency- Layers 4", 500, 0, 50);
    h1_dy_good_4->GetYaxis()->SetTitle("Inefficiency (%)");
    h1_dy_good_4->GetXaxis()->SetTitle("Y_{4Pred}-Y_{4} [mm]");
    TH1F *h1_dy_good_6 = new TH1F("h1_dy_good_6", "Prediction Inefficiency- Layers 6", 500, 0, 50);
    h1_dy_good_6->GetYaxis()->SetTitle("Inefficiency (%)");
    h1_dy_good_4->GetXaxis()->SetTitle("Y_{6Pred}-Y_{6} [mm]");
    TH1F *h1_dy_good_2 = new TH1F("h1_dy_good_2", "Prediction Inefficiency- Layers 2", 500, 0, 50);
    h1_dy_good_2->GetYaxis()->SetTitle("Inefficiency (%)");
    h1_dy_good_4->GetXaxis()->SetTitle("Y_{2Pred}-Y_{2} [mm]");
    TH1F *h1_dy_good_1 = new TH1F("h1_dy_good_1", "Prediction Inefficiency- Layers 1", 500, 0, 50);
    h1_dy_good_1->GetYaxis()->SetTitle("Inefficiency (%)");
    h1_dy_good_4->GetXaxis()->SetTitle("Y_{1Pred}-Y_{1} [mm]");

    TH1F *h1_dx_good_5 = new TH1F("h1_dx_good_5", "x Prediction Error -Layer 5", 500, 0, 200);
    TH1F *h1_dx_good_4 = new TH1F("h1_dx_good_4", "h1_dx_good_4", 500, 0, 200);
    TH1F *h1_dx_good_6 = new TH1F("h1_dx_good_6", "h1_dx_good_6", 500, 0, 200);
    TH1F *h1_dx_good_2 = new TH1F("h1_dx_good_2", "h1_dx_good_2", 500, 0, 200);
    TH1F *h1_dx_good_1 = new TH1F("h1_dx_good_1", "h1_dx_good_1", 500, 0, 200);

    TH1F *h1_dy_bad_5 = new TH1F("h1_dy_bad_5", "h1_dy_bad_5", 500, 0, 200);
    TH1F *h1_dy_bad_4 = new TH1F("h1_dy_bad_4", "h1_dy_bad_4", 500, 0, 200);
    TH1F *h1_dy_bad_6 = new TH1F("h1_dy_bad_6", "h1_dy_bad_6", 500, 0, 200);
    TH1F *h1_dy_bad_2 = new TH1F("h1_dy_bad_2", "h1_dy_bad_2", 500, 0, 200);
    TH1F *h1_dy_bad_1 = new TH1F("h1_dy_bad_1", "h1_dy_bad_1", 500, 0, 200);

    TH1F *h1_dx_bad_5 = new TH1F("h1_dx_bad_5", "h1_dx_bad_5", 500, 0, 200);
    TH1F *h1_dx_bad_4 = new TH1F("h1_dx_bad_4", "h1_dx_bad_4", 500, 0, 200);
    TH1F *h1_dx_bad_6 = new TH1F("h1_dx_bad_6", "h1_dx_bad_6", 500, 0, 200);
    TH1F *h1_dx_bad_2 = new TH1F("h1_dx_bad_2", "h1_dx_bad_2", 500, 0, 200);
    TH1F *h1_dx_bad_1 = new TH1F("h1_dx_bad_1", "h1_dx_bad_1", 500, 0, 200);

    TH1F *h1_good_match = new TH1F("h1_good_match", "h1_good_match", 10, 0, 5);
    TH1F *h1_bad_match = new TH1F("h1_bad_match", "h1_bad_match", 10, 0, 5);

    TH1F *h1_matchSum = new TH1F("h1_matchSum", "h1_matchSum", 12, 0, 6);

    TH1F *h1_fake_used = new TH1F("h1_fake_used", "h1_fake_used", 300, 0, 150);
    TH1F *h1_true_used = new TH1F("h1_true_used", "h1_true_used", 300, 0, 150);

    TH1F *h1_fake_chi2x = new TH1F("h1_fake_chi2x", "h1_fake_chi2x", 300, 0, 150);
    TH1F *h1_true_chi2x = new TH1F("h1_true_chi2x", "h1_true_chi2x", 300, 0, 150);
    TH1F *h1_fake_chi2y = new TH1F("h1_fake_chi2y", "h1_fake_chi2y", 300, 0, 150);
    TH1F *h1_true_chi2y = new TH1F("h1_true_chi2y", "h1_true_chi2y", 300, 0, 150);

    TH1F *h1_in_a = new TH1F("h1_in_a", "h1_in_a", 200, 0, 1);
    TH1F *h1_in_b = new TH1F("h1_in_b", "h1_in_b", 200, -200, 200);
    TH1F *h1_in_c = new TH1F("h1_in_c", "h1_in_c", 200, 0, 10000);

    TH1F *h1_out_a = new TH1F("h1_out_a", "h1_out_a", 200, 0, 1);
    TH1F *h1_out_b = new TH1F("h1_out_b", "h1_out_b", 200, -200, 200);
    TH1F *h1_out_c = new TH1F("h1_out_c", "h1_out_c", 200, 0, 10000);

    TH2F *h2_dx_dy_good_5 = new TH2F("h2_dx_dy_good_5", "(X_{pre}-X5, Y_{pre}-Y5) True Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_good_5->GetYaxis()->SetTitle("Y5_{pred}-Y5 [mm]");
    h2_dx_dy_good_5->GetXaxis()->SetTitle("X5_{pred}-X5 [mm]");
    TH2F *h2_dx_dy_bad_5 = new TH2F("h2_dx_dy_bad_5", "(X_{pre}-X5, Y_{pre}-Y5) Fake Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_bad_5->GetYaxis()->SetTitle("Y5_{pred}-Y5 [mm]");
    h2_dx_dy_bad_5->GetXaxis()->SetTitle("X5_{pred}-X5 [mm]");

    TH2F *h2_dx_dy_good_4 = new TH2F("h2_dx_dy_good_4", "(X_{pre}-X4, Y_{pre}-Y4) True Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_good_4->GetYaxis()->SetTitle("Y4_{pred}-Y4 [mm]");
    h2_dx_dy_good_4->GetXaxis()->SetTitle("X4_{pred}-X4 [mm]");
    TH2F *h2_dx_dy_bad_4 = new TH2F("h2_dx_dy_bad_4", "(X_{pre}-X4, Y_{pre}-Y4) Fake Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_bad_4->GetYaxis()->SetTitle("Y4_{pred}-Y4 [mm]");
    h2_dx_dy_bad_4->GetXaxis()->SetTitle("X4_{pred}-X4 [mm]");

    TH2F *h2_dx_dy_good_6 = new TH2F("h2_dx_dy_good_6", "(X_pre-X6 Y_pre-Y6) True Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_good_6->GetYaxis()->SetTitle("Y6_{pred}-Y3 [mm]");
    h2_dx_dy_good_6->GetXaxis()->SetTitle("X6_{pred}-X3 [mm]");
    TH2F *h2_dx_dy_bad_6 = new TH2F("h2_dx_dy_bad_6", "(X_pre-X6 Y_pre-Y6) Fake Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_bad_6->GetYaxis()->SetTitle("Y6_{pred}-Y6 [mm]");
    h2_dx_dy_bad_6->GetXaxis()->SetTitle("X6_{pred}-X6 [mm]");

    TH2F *h2_dx_dy_good_2 = new TH2F("h2_dx_dy_good_2", "(X_pre-X2 Y_pre-Y2) True Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_good_2->GetYaxis()->SetTitle("Y2_{pred}-Y2 [mm]");
    h2_dx_dy_good_2->GetXaxis()->SetTitle("X2_{pred}-X2 [mm]");
    TH2F *h2_dx_dy_bad_2 = new TH2F("h2_dx_dy_bad_2", "(X_pre-X2 Y_pre-Y2) Fake Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_bad_2->GetYaxis()->SetTitle("Y2_{pred}-Y2 [mm]");
    h2_dx_dy_bad_2->GetXaxis()->SetTitle("X2_{pred}-X2 [mm]");

    TH2F *h2_dx_dy_good_1 = new TH2F("h2_dx_dy_good_1", "(X_pre-X1 Y_pre-Y1) True Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_good_1->GetYaxis()->SetTitle("Y1_{pred}-Y1 [mm]");
    h2_dx_dy_good_1->GetXaxis()->SetTitle("X1_{pred}-X1 [mm]");
    TH2F *h2_dx_dy_bad_1 = new TH2F("h2_dx_dy_bad_1", "(X_pre-X1 Y_pre-Y1) Fake Tracks", 160, -80, 80, 80, -40, 40);
    h2_dx_dy_bad_1->GetYaxis()->SetTitle("Y1_{pred}-Y1 [mm]");
    h2_dx_dy_bad_1->GetXaxis()->SetTitle("X1_{pred}-X1 [mm]");

    TH2F *h2_chi2_good = new TH2F("h2_chi2_good", "h2_chi2_good", 50, 0, 5, 50, 0, 5);
    h2_chi2_good->GetYaxis()->SetTitle("Chi2 Y");
    h2_chi2_good->GetXaxis()->SetTitle("Chi2 X");
    TH2F *h2_chi2_bad = new TH2F("h2_chi2_bad", "h2_chi2_bad", 50, 0, 5, 20, 0, 5);
    h2_chi2_bad->GetYaxis()->SetTitle("Chi2 Y");
    h2_chi2_bad->GetXaxis()->SetTitle("Chi2 X");

    TH2F *h2_bad_origin = new TH2F("h2_bad_origin", "h2_bad_origin", 300, -1000, 1000, 300, -1000, 1000);
    TH2F *h2_good_origin = new TH2F("h2_good_origin", "h2_good_origin", 300, -1000, 1000, 300, -1000, 1000);
    h2_bad_origin->GetYaxis()->SetTitle("Y_{origin} [mm]");
    h2_good_origin->GetXaxis()->SetTitle("X_{origin} [mm]");

    // Track reconstructed counter
    Int_t nReconstructable6HitTrack = 0;
    Int_t nClones = 0;
    Int_t nTrueReconstructed = 0;
    Int_t n5HitReconstructed = 0;
    Int_t n4HitReconstructed = 0;

    Int_t nFakeReconstructed = 0;
    Int_t nTotalReconstructed = 0;

    Int_t nFriendlyFire = 0;
    Int_t nRejectedTracks = 0;

    // Input file directory
    TFile f(inputfilestring);

    // Pointers to original and recovered hits track trees respectively
    TTree *oTracks = NULL;
    TTree *rTracks = NULL;

    // Check if the recovered tracks available
    rTracks = (TTree *)gDirectory->Get("RTracks");
    if (not rTracks)
        cout << "rTracks not found" << endl;

    // Hit position
    float hitPos[6][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    bool bRecover[6] = {false, false, false, false, false, false};

    // Kinetics Variables
    double eta, p, pz;

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
        rTracks->SetBranchAddress("hit1x", &hitPos[0][0]);
        rTracks->SetBranchAddress("hit1y", &hitPos[0][1]);
        rTracks->SetBranchAddress("hit1z", &hitPos[0][2]);
        rTracks->SetBranchAddress("hit2x", &hitPos[1][0]);
        rTracks->SetBranchAddress("hit2y", &hitPos[1][1]);
        rTracks->SetBranchAddress("hit2z", &hitPos[1][2]);
        rTracks->SetBranchAddress("hit3x", &hitPos[2][0]);
        rTracks->SetBranchAddress("hit3y", &hitPos[2][1]);
        rTracks->SetBranchAddress("hit3z", &hitPos[2][2]);
        rTracks->SetBranchAddress("hit4x", &hitPos[3][0]);
        rTracks->SetBranchAddress("hit4y", &hitPos[3][1]);
        rTracks->SetBranchAddress("hit4z", &hitPos[3][2]);
        rTracks->SetBranchAddress("hit5x", &hitPos[4][0]);
        rTracks->SetBranchAddress("hit5y", &hitPos[4][1]);
        rTracks->SetBranchAddress("hit5z", &hitPos[4][2]);
        rTracks->SetBranchAddress("hit6x", &hitPos[5][0]);
        rTracks->SetBranchAddress("hit6y", &hitPos[5][1]);
        rTracks->SetBranchAddress("hit6z", &hitPos[5][2]);
        rTracks->SetBranchAddress("bhit1reco", &bRecover[0]);
        rTracks->SetBranchAddress("bhit2reco", &bRecover[1]);
        rTracks->SetBranchAddress("bhit3reco", &bRecover[2]);
        rTracks->SetBranchAddress("bhit4reco", &bRecover[3]);
        rTracks->SetBranchAddress("bhit5reco", &bRecover[4]);
        rTracks->SetBranchAddress("bhit6reco", &bRecover[5]);
        rTracks->SetBranchAddress("p", &p);
        rTracks->SetBranchAddress("pz", &pz);
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

    // Container of hits
    std::vector<TLayer> MTLayers;
    for (Int_t j = 0; j < 6; j++)
    {
        MTLayers.push_back(TLayer());
    }

    std::vector<track *> caseTracks;
    std::vector<track *> allTracks;

    // counters for tracks
    unsigned int nReconstructableElectron = 0, nTrueReconstructedElectron = 0;
    unsigned int nReconstructableP0P5Eta1 = 0, nTrueReconstructedP0P5Eta1 = 0;
    unsigned int nReconstructableP5P20Eta1 = 0, nTrueReconstructedP5P20Eta1 = 0;
    unsigned int nReconstructableP20Eta1 = 0, nTrueReconstructedP20Eta1 = 0;

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

    chrono::time_point<chrono::high_resolution_clock> timingStart = chrono::high_resolution_clock::now();
    for (Int_t i = 0; i < nEntries; i++)
    {
        previousEventnumber = eventNumber;

        if (i % 1000000 == 0)
            cout << "Progress: " << i * 100. / nEntries << endl;

        // Read a entry to memory
        if (not rTracks)
            oTracks->GetEntry(i);
        else
            rTracks->GetEntry(i);

        if (previousEventnumber != -1 and (previousEventnumber != eventNumber or i == nEntries - 1))
        { // Run the reconstruction in the end of an event or in the end of the readout
            // print the hits on every layer
            // cout << "Event summary by layers 6-5-4-3-2-1: " << MTLayers[5].hits.size() << " " << MTLayers[4].hits.size() << " " << MTLayers[3].hits.size() << " " << MTLayers[2].hits.size() << " " << MTLayers[1].hits.size() << " " << MTLayers[0].hits.size() << endl;

            // Perform Reconstruction
            for (int iCase = iCases; iCase < nCases; iCase++)
            {
                for (Int_t i3 = 0; i3 < MTLayers[2].hits.size(); i3++)
                {
                    // New track variables
                    // Prediction
                    float xPre[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                    float yPre[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                    // Prediction Error
                    float xErr[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                    float yErr[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                    // Match Result
                    bool match[5] = {false, false, false, false, false};

                    // Hits of this track
                    hit *trackHit[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

                    // Start the seed from the 3rd layer
                    trackHit[2] = MTLayers[2].hits[i3];
                    if (trackHit[2]->flag == true)
                    {
                        continue;
                    }

                    yPre[2] = trackHit[2]->y;
                    xPre[2] = trackHit[2]->x;

                    // Try every hit on the 4th layer
                    for (Int_t i4 = 0; i4 < MTLayers[3].hits.size(); i4++)
                    {
                        // Choose a hit from the 4th layer
                        trackHit[3] = MTLayers[3].hits[i4];
                        if (trackHit[3]->flag == true)
                        {
                            continue;
                        }
                        // A rough estimation for hit 4
                        yPre[3] = trackHit[2]->y * (trackHit[3]->z / trackHit[2]->z);
                        xPre[3] = trackHit[2]->x * (trackHit[3]->z / trackHit[2]->z);
#ifdef modfication
                        xPre[3] = xPre[3] + yPre[3] / ((trackHit[2]->z / trackHit[2]->x) * (1.0 / 0.003601 - trackHit[2]->y / trackHit[2]->z));
                        yPre[3] = yPre[3] / (1.0 - 0.003601 * trackHit[2]->y / trackHit[2]->z);
#endif
                        xErr[3] = xPre[3] - trackHit[3]->x;
                        yErr[3] = yPre[3] - trackHit[3]->y;

                        // Forward Tracking
                        if ((trackHit[3]->preLayerHit != NULL) && (trackHit[3]->preLayerHit == trackHit[2]))
                            match[2] = true;
                        else
                            match[2] = false;

                        if (match[2])
                        {
                            h2_dx_dy_good_4->Fill(xErr[3], yErr[3]);
                            h1_dy_good_4->Fill(yErr[3]);
                            h1_dx_good_4->Fill(xErr[3]);
                        }
                        else if (!match[2])
                        {
                            h2_dx_dy_bad_4->Fill(xErr[3], yErr[3]);
                            h1_dy_bad_4->Fill(yErr[3]);
                            h1_dx_bad_4->Fill(xErr[3]);
                        }

                        // First Selection use window 0
                        if (TMath::Abs(yErr[3]) > MAXDY[iCase][0])
                            continue;
                        else if (TMath::Abs(xErr[3]) > MAXDX[iCase][0])
                            continue;

                        // Use hit 3 and hit 4, Linear Fit to predict hit5
                        yPre[4] = trackHit[3]->y * 685.1 / 210.1 - trackHit[2]->y * 475 / 210.1;
                        xPre[4] = trackHit[3]->x * 685.1 / 210.1 - trackHit[2]->x * 475 / 210.1;
#ifdef modfication
                        xPre[4] = xPre[4] + yPre[4] * ((trackHit[3]->x - trackHit[2]->x) / 210.1) / ((1.0 / 0.003601) - (trackHit[3]->y - trackHit[2]->y) / 210.1);
                        yPre[4] = yPre[4] / (1 - 0.003601 * (trackHit[3]->y - trackHit[2]->y) / 210.1);
#endif

                        for (Int_t i5 = 0; i5 < MTLayers[4].hits.size(); i5++)
                        {
                            // Choose a hit from the 5th layer
                            trackHit[4] = MTLayers[4].hits[i5];
                            if (trackHit[4]->flag == true)
                            {
                                continue;
                            }
                            xErr[4] = xPre[4] - trackHit[4]->x;
                            yErr[4] = yPre[4] - trackHit[4]->y;

                            // Forward Tracking
                            if ((trackHit[4]->preLayerHit != NULL) && (trackHit[4]->preLayerHit == trackHit[3]))
                                match[3] = true;
                            else
                                match[3] = false;

                            if (match[2] && match[3])
                            {
                                h2_dx_dy_good_5->Fill(xErr[4], yErr[4]);
                                h1_dy_good_5->Fill(yErr[4]);
                                h1_dx_good_5->Fill(xErr[4]);
                            }
                            else if (match[2] && !match[3])
                            {
                                h2_dx_dy_bad_5->Fill(xErr[4], yErr[4]);
                                h1_dy_bad_5->Fill(yErr[4]);
                                h1_dx_bad_5->Fill(xErr[4]);
                            }

                            // Second Selection use window 1
                            if (TMath::Abs(yErr[4]) > MAXDY[iCase][1])
                                continue;
                            else if (TMath::Abs(xErr[4]) > MAXDX[iCase][1])
                                continue;

                            yPre[5] = trackHit[4]->y * 684.9 / 475.0 - trackHit[3]->y * 209.9 / 475.0;
                            xPre[5] = trackHit[4]->x * 684.9 / 475.0 - trackHit[3]->x * 209.9 / 475.0;
#ifdef modfication
                            xPre[5] = xPre[5] + yPre[5] * (trackHit[4]->x - trackHit[3]->x) / 475.0 / (1.0 / 0.003601 - (trackHit[4]->y - trackHit[3]->y) / 475.0);
                            yPre[5] = yPre[5] / (1 - 0.003601 * (trackHit[4]->y - trackHit[3]->y) / 475.0);
#endif

                            for (Int_t i6 = 0; i6 < MTLayers[5].hits.size(); i6++)
                            {
                                // Choose a hit from the 5th layer
                                trackHit[5] = MTLayers[5].hits[i6];
                                if (trackHit[5]->flag == true)
                                {
                                    continue;
                                }
                                // Linear Fit to predict hit5, 3 dots
                                // Use hit 3, hit 4 and hit 5, Linear Fit to predict hit6

                                xErr[5] = xPre[5] - trackHit[5]->x;
                                yErr[5] = yPre[5] - trackHit[5]->y;

                                if ((trackHit[5]->preLayerHit != NULL) && (trackHit[5]->preLayerHit == trackHit[4]))
                                    match[4] = true;
                                else
                                    match[4] = false;

                                if (match[2] && match[3] && match[4])
                                {
                                    h2_dx_dy_good_6->Fill(xErr[5], yErr[5]);
                                    h1_dy_good_6->Fill(yErr[5]);
                                    h1_dx_good_6->Fill(xErr[5]);
                                }
                                else if (match[2] && match[3] && !match[4])
                                {
                                    h2_dx_dy_bad_6->Fill(xErr[5], yErr[5]);
                                    h1_dy_bad_6->Fill(yErr[5]);
                                    h1_dx_bad_6->Fill(xErr[5]);
                                }

                                // Third Selection use window 3
                                if (TMath::Abs(yErr[5]) > MAXDY[iCase][2])
                                    continue;
                                else if (TMath::Abs(xErr[5]) > MAXDX[iCase][2])
                                    continue;

                                // Use hit 3, 4, 5, 6, Linear Fit to predict hit2 and hit1 y.
                                float parabolaParams[4] = {0, 0, 0, 0};
                                float paraResult[2] = {0.0, 0.0};

                                float linFitX[4] = {trackHit[2]->x, trackHit[3]->x, trackHit[4]->x, trackHit[5]->x};
                                float linFitY[4] = {trackHit[2]->y, trackHit[3]->y, trackHit[4]->y, trackHit[5]->y};

                                float linFitZ[4] = {trackHit[2]->z, trackHit[3]->z, trackHit[4]->z, trackHit[5]->z};

                                float fitResult[2] = {0, 0};

                                bool bPara = false;
                                // Use the parabola fit only when t
                                // if (iCase == 0)
                                // {
                                //     yPre[0] = trackHit[2]->y * 892.0 / 210.1 - trackHit[3]->y * 681.9 / 210.1;
                                //     yPre[1] = trackHit[2]->y * 682.1 / 210.1 - trackHit[3]->y * 472.0 / 210.1;

                                //     xPre[0] = trackHit[2]->x * 892.0 / 210.1 - trackHit[3]->x * 681.9 / 210.1;
                                //     xPre[1] = trackHit[2]->x * 682.1 / 210.1 - trackHit[3]->x * 472.0 / 210.1;

                                //     #ifdef modfication
                                //     xPre[0] = xPre[0] + yPre[0]*((trackHit[3]->x-trackHit[2]->x)/210.1)/((1.0/0.003601)-(trackHit[3]->y-trackHit[2]->y)/210.1);
                                //     yPre[0] = yPre[0] / (1 - 0.003601 * (trackHit[3]->y-trackHit[2]->y)/210.1);
                                //     xPre[1] = xPre[1] + yPre[1]*((trackHit[3]->x-trackHit[2]->x)/210.1)/((1.0/0.003601)-(trackHit[3]->y-trackHit[2]->y)/210.1);
                                //     yPre[1] = yPre[1] / (1 - 0.003601 * (trackHit[3]->y-trackHit[2]->y)/210.1);
                                //     #endif
                                //     linearCount++;
                                // }
                                // else
                                if (iCase < 2)
                                {
                                    bPara = parabolaFit(linFitX, linFitY, linFitZ, parabolaParams, paraResult);
                                    xPre[0] = paraResult[0];
                                    xPre[1] = paraResult[1];

                                    yPre[0] = trackHit[2]->y * 892.0 / 210.1 - trackHit[3]->y * 681.9 / 210.1;
                                    yPre[1] = trackHit[2]->y * 682.1 / 210.1 - trackHit[3]->y * 472.0 / 210.1;

#ifdef modfication
                                    xPre[0] = xPre[0] + yPre[0] * (2 * parabolaParams[2] * 7825.9 + parabolaParams[1]) / ((1.0 / 0.003601) - (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    xPre[1] = xPre[1] + yPre[1] * (2 * parabolaParams[2] * 8035.8 + parabolaParams[1]) / ((1.0 / 0.003601) - (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    yPre[0] = yPre[0] / (1 - 0.003601 * (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    yPre[1] = yPre[1] / (1 - 0.003601 * (trackHit[3]->y - trackHit[2]->y) / 210.1);
#endif

                                    parabolaCount++;
                                }
                                else
                                {
                                    cubicFit(linFitX, linFitY, linFitZ, parabolaParams, paraResult);
                                    xPre[0] = paraResult[0];
                                    xPre[1] = paraResult[1];
                                    yPre[0] = trackHit[2]->y * 892.0 / 210.1 - trackHit[3]->y * 681.9 / 210.1;
                                    yPre[1] = trackHit[2]->y * 682.1 / 210.1 - trackHit[3]->y * 472.0 / 210.1;

#ifdef modfication
                                    xPre[0] = xPre[0] + yPre[0] * (3 * parabolaParams[3] * 7825.9 * 7825.9 + 2 * parabolaParams[2] * 7825.9 + parabolaParams[1]) / ((1.0 / 0.003601) - (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    xPre[1] = xPre[1] + yPre[1] * (3 * parabolaParams[3] * 8035.8 * 8035.8 + 2 * parabolaParams[2] * 8035.8 + parabolaParams[1]) / ((1.0 / 0.003601) - (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    yPre[0] = yPre[0] / (1 - 0.003601 * (trackHit[3]->y - trackHit[2]->y) / 210.1);
                                    yPre[1] = yPre[1] / (1 - 0.003601 * (trackHit[3]->y - trackHit[2]->y) / 210.1);
#endif
                                    cubicCount++;
                                }

                                for (Int_t i2 = 0; i2 < MTLayers[1].hits.size(); i2++)
                                {
                                    // Choose a hit from the 5th layer
                                    trackHit[1] = MTLayers[1].hits[i2];
                                    if (trackHit[1]->flag == true)
                                    {
                                        continue;
                                    }
                                    xErr[1] = xPre[1] - trackHit[1]->x;
                                    yErr[1] = yPre[1] - trackHit[1]->y;
                                    // Backward Tracking
                                    if ((trackHit[2]->preLayerHit != NULL) && (trackHit[2]->preLayerHit == trackHit[1]))
                                        match[1] = true;
                                    else
                                        match[1] = false;

                                    if (match[1] && match[2] && match[3] && match[4])
                                    {
                                        h2_dx_dy_good_2->Fill(xErr[1], yErr[1]);
                                        h1_dy_good_2->Fill(yErr[1]);
                                        h1_dx_good_2->Fill(xErr[1]);
                                    }
                                    else if (!match[1] && match[2] && match[3] && match[4])
                                    {
                                        h2_dx_dy_bad_2->Fill(xErr[1], yErr[1]);
                                        h1_dy_bad_2->Fill(yErr[1]);
                                        h1_dx_bad_2->Fill(xErr[1]);
                                    }

                                    // Fourth Selection use window 4
                                    if (TMath::Abs(yErr[1]) > MAXDY[iCase][3])
                                        continue;
                                    else if (TMath::Abs(xErr[1]) > MAXDX[iCase][3])
                                    {
                                        if (bPara)
                                        {
                                            h1_out_a->Fill(parabolaParams[2]);
                                            h1_out_b->Fill(parabolaParams[1]);
                                            h1_out_c->Fill(parabolaParams[0]);
                                        }
                                        continue;
                                    }
                                    else if (bPara)
                                    {
                                        h1_in_a->Fill(parabolaParams[2]);
                                        h1_in_b->Fill(parabolaParams[1]);
                                        h1_in_c->Fill(parabolaParams[0]);
                                    }

                                    for (Int_t i1 = 0; i1 < MTLayers[0].hits.size(); i1++)
                                    {
                                        // Choose a hit from the 5th layer
                                        trackHit[0] = MTLayers[0].hits[i1];
                                        if (trackHit[0]->flag == true)
                                        {
                                            continue;
                                        }
                                        xErr[0] = xPre[0] - trackHit[0]->x;
                                        yErr[0] = yPre[0] - trackHit[0]->y;

                                        if ((trackHit[1]->preLayerHit != NULL) && (trackHit[1]->preLayerHit == trackHit[0]))
                                            match[0] = true;
                                        else
                                            match[0] = false;

                                        if (match[0] && match[1] && match[2] && match[3] && match[4])
                                        {
                                            h2_dx_dy_good_1->Fill(xErr[0], yErr[0]);
                                            h1_dy_good_1->Fill(yErr[0]);
                                            h1_dx_good_1->Fill(xErr[0]);
                                        }
                                        else if (!match[0] && match[1] && match[2] && match[3] && match[4])
                                        {
                                            h2_dx_dy_bad_1->Fill(xErr[0], yErr[0]);
                                            h1_dy_bad_1->Fill(yErr[0]);
                                            h1_dx_bad_1->Fill(xErr[0]);
                                        }

                                        // Fifth Selection use window 5
                                        if (TMath::Abs(yErr[0]) > MAXDY[iCase][4])
                                            continue;
                                        else if (TMath::Abs(xErr[0]) > MAXDX[iCase][4])
                                            continue;
                                        track *trackCandidate = new track(trackHit, xErr, yErr, match);
                                        caseTracks.push_back(trackCandidate); // Store all the tracks reconstructed
                                    }                                         // Loop on the 1st layer
                                }                                             // Loop on the 2nd layer
                            }                                                 // Loop on the 6th layer
                        }                                                     // Loop on the 5th layer
                    }                                                         // Loop on the 4th layer
                }                                                             // Loop on the 3th layer

                int usedSum;
                vector<track *>::iterator _t;
                for (_t = caseTracks.begin(); _t != caseTracks.end(); _t++)
                {
                    usedSum = (*_t)->GetUsedSum();
                    if (!conf_alone)
                    {
                        // Purge the tracks
                        if (usedSum > 6)
                        {
                            // Check every hit, if it is used more that once, purge it
                            for (Int_t j = 0; j < 6; j++)
                            {
                                if ((*_t)->trackMatch)
                                {
                                    h2_chi2_good->Fill((*_t)->fChi2X, (*_t)->fChi2Y);
                                }
                                else
                                {
                                    h2_chi2_bad->Fill((*_t)->fChi2X, (*_t)->fChi2Y);
                                }
                                // If this track is already rejected, just break the loop
                                if ((*_t)->trackAccept == false)
                                {
                                    break;
                                }
                                // if true, then it has related tracks.
                                if ((*_t)->hits[j]->used != 1)
                                {
                                    // Loop on every other related tracks.
                                    for (auto p : (*_t)->hits[j]->relatedTracks)
                                    {
                                        // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                        if (((*p) == (**_t)) < 1 or p->trackAccept == false)
                                        {
                                            continue;
                                        }

                                        // If the chi2 of this track is larger than other clone track, reject this track and break.
                                        if ((*_t)->fChi2X > p->fChi2X and (*_t)->fChi2Y > p->fChi2Y)
                                        {
                                            (*_t)->trackAccept = false;
                                            // See how many friendly fire
                                            h1_matchSum->Fill((*_t)->matchSum);
                                            nRejectedTracks++;
                                            if ((*_t)->matchSum == 5)
                                            {
                                                // Record the friendly fire rate
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
                            // if ((*_t)->trackMatch)
                            // {
                            //     h2_good_origin->Fill(((*_t)->hits[0]->x * (*_t)->hits[5]->z-(*_t)->hits[5]->x*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z),((*_t)->hits[0]->y * (*_t)->hits[5]->z-(*_t)->hits[5]->y*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z));
                            // }
                            // else
                            // {
                            //     h2_bad_origin->Fill(((*_t)->hits[0]->x * (*_t)->hits[5]->z-(*_t)->hits[5]->x*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z),((*_t)->hits[0]->y * (*_t)->hits[5]->z-(*_t)->hits[5]->y*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z));
                            // }
                        }
                        // Take care of the hits and counters of the rejected tracks
                        if ((*_t)->trackAccept == false)
                        {
                            for (Int_t j = 0; j < 6; j++)
                            {
                                (*_t)->hits[j]->used--;
                            }
                        }
                    }
                    // After purging, a track is either rejected or accepted (default)
                    // No purging, all tracks are accepted.
                    if ((*_t)->trackAccept)
                    {
                        // For all the accepted tracks, flag their hits used in this case, no longer use them.
                        for (Int_t j = 0; j < 6; j++)
                        {
                            (*_t)->hits[j]->flag = true;
                        }
                    }
                } // Loop on case tracks, purge.
                allTracks.insert(allTracks.end(), caseTracks.begin(), caseTracks.end());
                caseTracks.clear();

                // cout << "---------- Case " << iCase << " ----------" << endl;
                // if (nReconstructableP5Eta1 > 0)
                //     cout << "Eff. (p5eta1): " << (nTrueReconstructedP5Eta1 * 1. / nReconstructableP5Eta1) << " ( " << nTrueReconstructedP5Eta1 << " / " << nReconstructableP5Eta1 << " ) " << endl;
                // if (nReconstructableP20Eta1 > 0)
                //     cout << "Eff. (p20eta1): " << (nTrueReconstructedP20Eta1 * 1. / nReconstructableP20Eta1) << " ( " << nTrueReconstructedP20Eta1 << " / " << nReconstructableP20Eta1 << " ) " << endl;
                // cout << "Eff.: " << (nTrueReconstructed * 1.0 / nReconstructable6HitTrack) << " ( " << nTrueReconstructed << " / " << nReconstructable6HitTrack << " ) " << endl;

                // cout << "Fake rate: " << (nTotalReconstructed - nTrueReconstructed) * 1.0 / nTotalReconstructed << " ( " << (nTotalReconstructed - nTrueReconstructed) << " ) " << endl;
                // cout << "Fake rate (alternative. counting): " << (nFakeReconstructed)*1.0 / nTotalReconstructed << " ( " << nFakeReconstructed << " ) " << endl;

                // cout << "Parabola Count:" << parabolaCount << endl;
                // cout << "Linear Count:" << linearCount << endl;
                // cout << "Cubic Count:" << cubicCount << endl;
            } // Loop for 3 Cases

            // Global Purge
            int usedSum;
            vector<track *>::iterator _t;

            for (_t = allTracks.begin(); _t != allTracks.end(); _t++)
            {
                usedSum = (*_t)->GetUsedSum();
                if (!conf_alone)
                {
                    // Purge the tracks
                    if (usedSum > 6)
                    {
                        // Check every hit, if it is used more that once, purge it
                        for (Int_t j = 0; j < 6; j++)
                        {
                            if ((*_t)->trackMatch)
                            {
                                h2_chi2_good->Fill((*_t)->fChi2X, (*_t)->fChi2Y);
                            }
                            else
                            {
                                h2_chi2_bad->Fill((*_t)->fChi2X, (*_t)->fChi2Y);
                            }
                            // If this track is already rejected, just break the loop
                            if ((*_t)->trackAccept == false)
                            {
                                break;
                            }
                            // if true, then it has related tracks.
                            if ((*_t)->hits[j]->used != 1)
                            {
                                // Loop on every other related tracks.
                                for (auto p : (*_t)->hits[j]->relatedTracks)
                                {
                                    // At least share 70% of their hits in the T-stations seeding region, is a clone tracks
                                    if (((*p) == (**_t)) < 1 or p->trackAccept == false)
                                    {
                                        continue;
                                    }

                                    // If the chi2 of this track is larger than other clone track, reject this track and break.
                                    if ((*_t)->fChi2X > p->fChi2X and (*_t)->fChi2Y > p->fChi2Y)
                                    {
                                        (*_t)->trackAccept = false;
                                        // See how many friendly fire
                                        h1_matchSum->Fill((*_t)->matchSum);
                                        nRejectedTracks++;
                                        if ((*_t)->matchSum == 5)
                                        {
                                            // Record the friendly fire rate
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
                        // if ((*_t)->trackMatch)
                        // {
                        //     h2_good_origin->Fill(((*_t)->hits[0]->x * (*_t)->hits[5]->z-(*_t)->hits[5]->x*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z),((*_t)->hits[0]->y * (*_t)->hits[5]->z-(*_t)->hits[5]->y*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z));
                        // }
                        // else
                        // {
                        //     h2_bad_origin->Fill(((*_t)->hits[0]->x * (*_t)->hits[5]->z-(*_t)->hits[5]->x*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z),((*_t)->hits[0]->y * (*_t)->hits[5]->z-(*_t)->hits[5]->y*(*_t)->hits[0]->z)/( (*_t)->hits[5]->z- (*_t)->hits[0]->z));
                        // }
                    }
                    // Take care of the hits and counters of the rejected tracks
                    if ((*_t)->trackAccept == false)
                    {
                        for (Int_t j = 0; j < 6; j++)
                        {
                            (*_t)->hits[j]->used--;
                        }
                    }
                }
                // After purging, a track is either rejected or accepted (default)
                // No purging, all tracks are accepted.

                // Mark the motherTrack
                if ((*_t)->trackAccept && (*_t)->trackMatch)
                    (*_t)->hits[0]->motherTrack->particleMatched = true;
            } // Loop on all the tracks, purge and tally the result

            for (_t = allTracks.begin(); _t != allTracks.end(); _t++)
            {
                // Tally the result
                if ((*_t)->trackAccept)
                {
#ifdef layerStudy
                    // DONT COUNT RECOVERED HITS ON 1ST LAYER
                    if ((*_t)->hits[0]->specialFlag && (*_t)->trackMatch)
                    {
                        continue;
                    }
#endif
                    nTotalReconstructed++;
                    if ((*_t)->trackMatch)
                    {
                        p = (*_t)->hits[5]->hP;
                        pz = (*_t)->hits[5]->hPz;
                        eta = 1.0 / 2.0 * log((p + pz) / (p - pz)); // log is in e base in C++

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
                        if ((*_t)->nTrueHits == 5)
                        {
                            if (!((*_t)->hits[(*_t)->theTrueHit]->motherTrack->particleMatched))
                            {
                                n5HitReconstructed++;
                                (*_t)->hits[(*_t)->theTrueHit]->motherTrack->particleMatched = true;
                            }
                            else
                            {
                                nClones++;
                            }
                        }
                        else if ((*_t)->nTrueHits == 4)
                        {
                            if (!((*_t)->hits[(*_t)->theTrueHit]->motherTrack->particleMatched))
                            {
                                n4HitReconstructed++;
                                (*_t)->hits[(*_t)->theTrueHit]->motherTrack->particleMatched = true;
                            }
                            else
                            {
                                nClones++;
                            }
                        }
                    }
                }
                delete *_t;
            }

            allTracks.clear();

            // Clear the layers after finish the task
            for (Int_t j = 0; j < 6; j++)
            {
                // delete every hit
                for (auto p : MTLayers[j].hits)
                {
                    delete p;
                }
                // delete the layer
                MTLayers[j].hits.clear();
            }
            eventCounter++;
        } // Run the reconstruction in the end of an event or in the end of the readout
        else
        {
            // Already got an entry before else

            bool bCount[6] = {false, false, false, false, false, false};

            hit *hitSet[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

            eta = 1.0 / 2.0 * log((p + pz) / (p - pz)); // log is in e base in C++

            // Don't read this track if not satisfies the cut
            if (p > conf_max_p or p < conf_min_p or ((eta < 2 or eta > 5) and (conf_eta)))
                continue;

            if (rTracks)
            {
                for (Int_t j = 0; j < 6; j++)
                {
                    hitSet[j] = new hit(hitPos[j][0], hitPos[j][1], hitPos[j][2]);
                    hitSet[j]->hP = p;
                    hitSet[j]->hPz = pz;
                    hitSet[j]->eventID = eventNumber;
                    bCount[j] = (hitPos[j][0] != 0 and hitPos[j][1] != 0 and hitPos[j][2] != 0);
                    hitSet[j]->used = 0;
                }

                // At least one hit in six layer (does not take into account dead area)
                if (!(bCount[0] or bCount[1] or bCount[2] or bCount[3] or bCount[4] or bCount[5]))
                {
                    for (Int_t j = 0; j < 6; j++)
                    {
                        delete hitSet[j];
                    }
                    continue;
                }

                for (Int_t j = 0; j < 6; j++)
                {
                    if (((conf_deadarea == 0.) or !IsInDeadArea(hitSet[j]->x, hitSet[j]->y)) and bCount[j])
                        MTLayers[j].hits.push_back(hitSet[j]);
                    else
                        delete hitSet[j];
                }
            }
            else
            {
                cout << "rTracks not found" << endl;
            }
            // All the Reconstructible Tracks refer to tracks with one hit at every layer.
            if (!(bCount[0] and bCount[1] and bCount[2] and bCount[3] and bCount[4] and bCount[5]))
            {
                continue;
            }

#ifdef layerStudy
            if (!bRecover[0])
            {
#endif

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
                    nReconstructable6HitTrack++;

                    if (particleID == -11 or particleID == 11)
                        nReconstructableElectron++;
                }
#ifdef layerStudy
            }
            else
            {
                hitSet[0]->specialFlag = true;
            }
#endif
            // COMMENT HERE
            // Set the motherTrack of hits
            hitSet[0]->motherTrack.reset(new track(hitSet));
            for (int j = 0; j < 5; j++)
            {
                hitSet[j + 1]->motherTrack = hitSet[j]->motherTrack;
            }
            if (p > conf_min_p and p < conf_max_p and ((eta > 2 and eta < 5) or (conf_eta == 0)) and
                ((bCount[0] and bCount[1] and bCount[2] and bCount[3] and bCount[4] and bCount[5]))) //At least one hit considering dead area in each layer
            {
                for (Int_t j = 0; j < 5; j++)
                {
                    hitSet[j + 1]->setPreLayer(hitSet[j]);
                    hitSet[j]->setNextLayer(hitSet[j + 1]);
                }
            }
        } //continue reading events and filling the information for the layers
    }

    fakeOut.close();
    goodOut.close();
    logOut << "---------- Total ----------" << endl;
    logOut << "# of Events:" << eventCounter << endl;
    logOut << "Time/Event (ms):" << 0.001 * chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - timingStart).count() / eventCounter << endl;
    logOut << "Electron fraction (reconstructable): " << (nReconstructableElectron * 1. / nReconstructable6HitTrack) << " ( " << nReconstructableElectron << " / " << nReconstructable6HitTrack << " ) " << endl;
    if (nReconstructableElectron > 0)
        logOut << "Electron reconstructed rate: " << (nTrueReconstructedElectron * 1. / nReconstructableElectron) << " ( " << nTrueReconstructedElectron << " / " << nReconstructableElectron << " ) " << endl;

    if (nReconstructableP0P5Eta1 > 0)
        logOut << "Eff. (p0p5eta1): " << (nTrueReconstructedP0P5Eta1 * 1. / nReconstructableP0P5Eta1) << " ( " << nTrueReconstructedP0P5Eta1 << " / " << nReconstructableP0P5Eta1 << " ) " << endl;
    if (nReconstructableP5P20Eta1 > 0)
        logOut << "Eff. (p5p20eta1): " << (nTrueReconstructedP5P20Eta1 * 1. / nReconstructableP5P20Eta1) << " ( " << nTrueReconstructedP5P20Eta1 << " / " << nReconstructableP5P20Eta1 << " ) " << endl;
    if (nReconstructableP20Eta1 > 0)
        logOut << "Eff. (p20eta1): " << (nTrueReconstructedP20Eta1 * 1. / nReconstructableP20Eta1) << " ( " << nTrueReconstructedP20Eta1 << " / " << nReconstructableP20Eta1 << " ) " << endl;

    logOut << "Eff.: " << (nTrueReconstructed * 1.0 / nReconstructable6HitTrack) << " ( " << nTrueReconstructed << " / " << nReconstructable6HitTrack << " ) " << endl;
    logOut << "Eff(70\% hits true): " << (nTrueReconstructed + n4HitReconstructed + n5HitReconstructed) * 1.0 / nReconstructable6HitTrack << " ( " << (nTrueReconstructed + n4HitReconstructed + n5HitReconstructed) << " ) " << endl;

    logOut << "Ghost rate(6 hits true): " << (nTotalReconstructed - nTrueReconstructed) * 1.0 / nTotalReconstructed << " ( " << (nTotalReconstructed - nTrueReconstructed) << " ) " << endl;
    logOut << "Ghost rate (6 hits true alternative. counting): " << (nFakeReconstructed)*1.0 / nTotalReconstructed << " ( " << nFakeReconstructed << " ) " << endl;

    logOut << "Ghost rate(70\% hits true): " << (nTotalReconstructed - nTrueReconstructed - n5HitReconstructed - n4HitReconstructed) * 1.0 / nTotalReconstructed << " ( " << (nTotalReconstructed - nTrueReconstructed - n5HitReconstructed - n4HitReconstructed) << " ) " << endl;
    logOut << "Ghost rate (70\% hits true alternative. counting): " << (nFakeReconstructed - n5HitReconstructed - n4HitReconstructed) * 1.0 / nTotalReconstructed << " ( " << nFakeReconstructed - n5HitReconstructed - n4HitReconstructed << " ) " << endl;

    logOut << "Friendly Fire Rate:" << nFriendlyFire * 1.0 / nRejectedTracks << " ( " << nFriendlyFire << "/" << nRejectedTracks << " ) " << endl;
    logOut << "Clone Rate:" << nClones*1.0/(nTrueReconstructed + n4HitReconstructed + n5HitReconstructed+nClones) << endl;
    logOut << "Parabola Count:" << parabolaCount << endl;
    logOut << "Linear Count:" << linearCount << endl;
    logOut << "Cubic Count:" << cubicCount << endl;

    //Plotting run parameters
    for (int iCase = iCases; iCase < nCases; iCase++)
    {
        logOut << "---------- Case " << iCase << " ----------" << endl;
        logOut << "Window 5 params X " << MAXDX[iCase][4] << " Y " << MAXDY[iCase][4] << endl;
        logOut << "Window 4 params X " << MAXDX[iCase][3] << " Y " << MAXDY[iCase][3] << endl;
        logOut << "Window 3 params X " << MAXDX[iCase][2] << " Y " << MAXDY[iCase][2] << endl;
        logOut << "Window 2 params X " << MAXDX[iCase][1] << " Y " << MAXDY[iCase][1] << endl;
        logOut << "Window 1 params X " << MAXDX[iCase][0] << " Y " << MAXDY[iCase][0] << endl;
    }



    c1->SetLogy(0);
    h2_dx_dy_good_5->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_good_5.pdf");
    h2_dx_dy_bad_5->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_bad_5.pdf");

    h2_dx_dy_good_4->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_good_4.pdf");
    h2_dx_dy_bad_4->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_bad_4.pdf");

    h2_dx_dy_good_6->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_good_6.pdf");
    h2_dx_dy_bad_6->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_bad_6.pdf");

    h2_dx_dy_good_2->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_good_2.pdf");
    h2_dx_dy_bad_2->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_bad_2.pdf");

    h2_dx_dy_good_1->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_good_1.pdf");
    h2_dx_dy_bad_1->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_dx_dy_bad_1.pdf");

    h2_good_origin->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_good_origin.pdf");
    h2_bad_origin->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_bad_origin.pdf");

    h2_chi2_bad->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_chi2_bad.pdf");
    h2_chi2_good->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h2_chi2_good.pdf");

    h1_bad_match->Scale(1. / h1_bad_match->Integral());
    h1_bad_match->Draw();
    c1->SaveAs(outputplotsfolder + "h1_bad_match.pdf");
    h1_good_match->Scale(1. / h1_good_match->Integral());
    h1_good_match->Draw();
    c1->SaveAs(outputplotsfolder + "h1_good_match.pdf");

    h1_true_used->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_true_used.pdf");
    h1_fake_used->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_fake_used.pdf");
    h1_matchSum->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_matchSum.pdf");

    h1_true_chi2x->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_true_chi2x.pdf");
    h1_true_chi2y->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_true_chi2y.pdf");
    h1_fake_chi2x->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_fake_chi2x.pdf");
    h1_fake_chi2y->Draw("");
    c1->SaveAs(outputplotsfolder + "h1_fake_chi2y.pdf");

    h1_in_a->Draw("");
    h1_in_a->Scale(1.0 / h1_in_a->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_in_a.pdf");
    h1_in_b->Draw("");
    h1_in_b->Scale(1.0 / h1_in_b->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_in_b.pdf");
    h1_in_c->Draw("");
    h1_in_c->Scale(1.0 / h1_in_c->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_in_c.pdf");

    h1_out_a->Draw("");
    h1_out_a->Scale(1.0 / h1_out_a->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_out_a.pdf");
    h1_out_b->Draw("");
    h1_out_b->Scale(1.0 / h1_out_b->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_out_b.pdf");
    h1_out_c->Draw("");
    h1_out_c->Scale(1.0 / h1_out_c->GetEntries());
    c1->SaveAs(outputplotsfolder + "h1_out_c.pdf");

        c1->SetLogy(1);

    
    h1_dy_good_5->Draw("");
    h1_dy_bad_5->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dy_5.pdf");
    make_inefficiency_window_plot(c1, h1_dy_good_5, h1_dy_bad_5, outputplotsfolder);
    
    h1_dy_good_4->Draw("");
    h1_dy_bad_4->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dy_4.pdf");
    make_inefficiency_window_plot(c1, h1_dy_good_4, h1_dy_bad_4, outputplotsfolder);
    
    h1_dy_good_6->Draw("");
    h1_dy_bad_6->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dy_6.pdf");
    make_inefficiency_window_plot(c1, h1_dy_good_6, h1_dy_bad_6, outputplotsfolder);
    
    h1_dy_good_2->Draw("");
    h1_dy_bad_2->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dy_2.pdf");
    make_inefficiency_window_plot(c1, h1_dy_good_2, h1_dy_bad_2, outputplotsfolder);
    
    h1_dy_good_1->Draw("");
    h1_dy_bad_1->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dy_1.pdf");
    make_inefficiency_window_plot(c1, h1_dy_good_1, h1_dy_bad_1, outputplotsfolder);

    
    h1_dx_good_5->Draw("");
    h1_dx_bad_5->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dx_5.pdf");
    make_inefficiency_window_plot(c1, h1_dx_good_5, h1_dx_bad_5, outputplotsfolder);
    
    h1_dx_good_4->Draw("");
    h1_dx_bad_4->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dx_4.pdf");
    make_inefficiency_window_plot(c1, h1_dx_good_4, h1_dx_bad_4, outputplotsfolder);
    
    h1_dx_good_6->Draw("");
    h1_dx_bad_6->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dx_6.pdf");
    make_inefficiency_window_plot(c1, h1_dx_good_6, h1_dx_bad_6, outputplotsfolder);
    
    h1_dx_good_2->Draw("");
    h1_dx_bad_2->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dx_2.pdf");
    make_inefficiency_window_plot(c1, h1_dx_good_2, h1_dx_bad_2, outputplotsfolder);
    
    h1_dx_good_1->Draw("");
    h1_dx_bad_1->Draw("same");
    c1->SaveAs(outputplotsfolder + "h1_dx_1.pdf");
    make_inefficiency_window_plot(c1, h1_dx_good_1, h1_dx_bad_1, outputplotsfolder);
}
