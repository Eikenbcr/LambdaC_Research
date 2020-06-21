//Testing Different Fit Functions
#define FitFunctionAnalysis_cxx
#include "FitFunctionAnalysis.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"
#include "DGOneMuTwoTotals.C"
#include "DGOneMuOneTotal.C"
#include "DGTwoMuTwoTotals.C"
#include "DGTwoMuOneTotal.C"

using namespace std;

#include <vector>
#include <algorithm>
#include <functional>
#include <chrono>
#include <thread>
#include <array>

std::array<std::array<std::array<std::array<double, 50>, 50>, 50>, 50> fomKPi;

TH1D * TestHistogram = nullptr;

TFile * File = nullptr;
TCanvas * c1 = nullptr;

void FitFunctionAnalysis::Begin(TTree * /*tree*/)
{

   TString option = GetOption();
}

void FitFunctionAnalysis::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

   TestHistogram = new TH1D("Mass [MeV]", "Mass of LambdaC", 100, 2185, 2385);
   TestHistogram->GetXaxis()->SetTitle("MeV");
   TestHistogram->GetYaxis()->SetTitle("Events Per 2 MeV");

   File = new TFile("FitFunctionTest.root", "RECREATE");
  gFile = File;

  for (double m = 0; m < 50; m++){
    for (double n = 0; n < 50; n++){
      for (double p = 0; p < 50; p++){
        for (double q = 0; q < 50; q++){

        double xf = 0.735 + (0.07 * m);
        double yf = 0.2575 + (0.015 * n);
        double wf = 0.735 + (0.07 * p);
        double zf = 0.2575 + (0.015 * q);

        double  SignalK = ((6087.15 + (-1115.78 * xf) + (-65.9315 * xf * xf) + TMath::Exp(8.98038 + (3.60867 * xf) + (-4.52552 * xf * xf))) * (9643.84 + (-26988.8 * yf) + (22798.9 * yf * yf) + TMath::Exp(-9.44465 + (-31.0399 * yf) + (51.3565 * yf * yf))));
        double  BackgroundK = ((2945.82 + (511.272 * xf) + (-268.896 * xf * xf) + TMath::Exp(12.1195 + (-1.34962 * xf) + (-1.58848 * xf * xf))) * (13724.6 + (-35777.2 * yf) + (28770.7 * yf * yf) + TMath::Exp(-8.98010 + (-31.2435 * yf) + (50.6150 * yf * yf))));
        double  SignalPi = ((6152.03 + (-1623.24 * wf) + (38.7262 * wf * wf) + TMath::Exp(10.2967 + (0.801222 * wf) + (-2.49779 * wf * wf))) * (4902.96 + (-15915.0 * zf) + (17364.9 * zf * zf) + TMath::Exp(-5.38716 + (-22.8278 * zf) + (38.9474 * zf * zf))));
        double  BackgroundPi = ((2911.82 + (236.126 * wf) + (-221.903 * wf * wf) + TMath::Exp(10.4840 + (2.60065 * wf) + (-3.68753 * wf * wf))) * (6700.33 + (-20491.8 * zf) + (22428.1 * zf * zf) + TMath::Exp(-3.26972 + (-23.0203 * zf) + (36.6292 * zf * zf))));

        fomKPi[m][n][p][q] = ((SignalPi * SignalK) / (BackgroundPi * BackgroundK));}}}}

   c1 = new TCanvas("canvas", "Test Canvas");

}

Bool_t FitFunctionAnalysis::Process(Long64_t entry)
{
  GetEntry(entry);
   fReader.SetEntry(entry);

      //    Corrected Xi Mass
          double CorrectedXiMass = ((*Xi_MM) - (*Lambda_MM) + 1115.683);

      //    Corrected Lambda Mass
          double CorrectedLambdaMass = ((*Lc_MM) - (*Xi_MM) + (1321.71));

bool BorderCut = (
               (TMath::Log10(*PromptPi_IPCHI2_OWNPV) > 0.7)
            && (TMath::Log10(*PromptPi_IPCHI2_OWNPV) < 4.2)
            && (*PromptPi_MC15TuneV1_ProbNNpi > 0.25)
            && (TMath::Log10(*PromptK_IPCHI2_OWNPV) > 0.7)
            && (TMath::Log10(*PromptK_IPCHI2_OWNPV) < 4.2)
            && (*PromptK_MC15TuneV1_ProbNNk > 0.25)
            );

   bool AdditionalCuts = (
      (*Lc_PT > 2000.)
&&    (*Xi_PT > 1000.)
&&    (CorrectedXiMass > 1310. && CorrectedXiMass < 1330.)
   );


int   m = floor(abs(TMath::Log10(*PromptPi_IPCHI2_OWNPV) - 0.7)/(0.07));
int   n = floor(abs(*PromptPi_MC15TuneV1_ProbNNpi - 0.25)/(0.015));
int   p = floor(abs(TMath::Log10(*PromptK_IPCHI2_OWNPV) - 0.7)/(0.07));
int   q = floor(abs(*PromptK_MC15TuneV1_ProbNNk - 0.25)/(0.015));

double FOM = fomKPi[m][n][p][q];

if (FOM > 1.0 && AdditionalCuts && BorderCut){
TestHistogram->Fill(CorrectedLambdaMass);}

   return kTRUE;
}

void FitFunctionAnalysis::SlaveTerminate()
{
}

void FitFunctionAnalysis::Terminate()
{
  Double_t sigma;
  Double_t deltaSigma;
  Double_t mu;
  Double_t deltaMu;
  Double_t total;
  Double_t deltaTotal;
  TString sigmaStr;
  TString deltaSigmaStr;
  TString muStr;
  TString deltaMuStr;
  TString totalStr;
  TString deltaTotalStr;

  TF1 *SingleGaussian = new TF1("SingleGaussian",fit2MeV_Gaussian,2100.,2500.,5);
  SingleGaussian->SetParameter(0,400.);
  SingleGaussian->SetParameter(1,2286);
  SingleGaussian->SetParameter(2, 6);
  SingleGaussian->SetParLimits(2, 0.,20.);
  SingleGaussian->SetParameter(3, 0.);
  SingleGaussian->SetParameter(4, 0.);

  TF1 *DG1Mu2Tot = new TF1("DG1Mu2Tot",DGOneMuTwoTotal,2100.,2500.,7);
  DG1Mu2Tot->SetParameter(0, 400.);
  DG1Mu2Tot->SetParameter(1, 2286);
  DG1Mu2Tot->SetParameter(2, 6);
  DG1Mu2Tot->SetParLimits(2, 0., 20.);
  DG1Mu2Tot->SetParameter(3, 400.);
  DG1Mu2Tot->SetParameter(4, 6);
  DG1Mu2Tot->SetParLimits(4, 0., 20.);
  DG1Mu2Tot->SetParameter(5, 0.);
  DG1Mu2Tot->SetParameter(6, 0.);

  TF1 *DG1Mu1Tot = new TF1("DG1Mu1Tot",DGOneMuOneTotal,2100.,2500.,7);
  DG1Mu1Tot->SetParameter(0, 0.5);
  DG1Mu1Tot->SetParameter(1, 400);
  DG1Mu1Tot->SetParameter(2, 2286.);
  DG1Mu1Tot->SetParameter(3, 6);
  DG1Mu1Tot->SetParameter(4, 6);
  DG1Mu1Tot->SetParLimits(3, 0., 20.);
  DG1Mu1Tot->SetParLimits(4, 0., 20.);
  DG1Mu1Tot->SetParameter(5, 0.);
  DG1Mu1Tot->SetParameter(6, 0.);

  TF1 *DG2Mu2Tot = new TF1("DG2Mu2Tot",DGTwoMuTwoTotal,2100.,2500.,8);
  DG2Mu2Tot->SetParameter(0, 1000.);
  DG2Mu2Tot->SetParameter(1, 2286);
  DG2Mu2Tot->SetParameter(2, 6);
  DG2Mu2Tot->SetParLimits(2, 0., 20.);
  DG2Mu2Tot->SetParameter(3, 1000.);
  DG2Mu2Tot->SetParameter(4, 2287);
  DG2Mu2Tot->SetParameter(5, 6);
  DG2Mu2Tot->SetParLimits(5, 0., 20.);
  DG2Mu2Tot->SetParameter(6, 0.);
  DG2Mu2Tot->SetParameter(7, 0.);

  TF1 *DG2Mu1Tot = new TF1("DG2Mu1Tot",DGTwoMuOneTotal,2100.,2500.,8);
  DG2Mu1Tot->SetParameter(0, 0.5);
  DG2Mu1Tot->SetParameter(1, 10000);
  DG2Mu1Tot->SetParameter(2, 2286.);
  DG2Mu1Tot->SetParameter(3, 5);
  DG2Mu1Tot->SetParameter(4, 2287);
  DG2Mu1Tot->SetParameter(5, 5);
  DG2Mu1Tot->SetParLimits(3, 0., 20.);
  DG2Mu1Tot->SetParLimits(5, 0., 20.);
  DG2Mu1Tot->SetParameter(6, 0.);
  DG2Mu1Tot->SetParameter(7, 0.);

  TestHistogram->Fit("SingleGaussian");
  TestHistogram->SetMinimum(0);
  TestHistogram->SetTitle("Lc_MM - Single Gaussian");
  c1->Write("Single Gaussian Fit");

  TestHistogram->Fit("DG1Mu2Tot");
  TestHistogram->SetMinimum(0);
  TestHistogram->SetTitle("Lc_MM - Double Gaussian w/ One Mean and Two Totals");
  c1->Write("DG (One Mean and Two Totals)");

  TestHistogram->Fit("DG1Mu1Tot");
  TestHistogram->SetMinimum(0);
  TestHistogram->SetTitle("Lc_MM - Double Gaussian w/ One Mean and One Total");
  c1->Write("DG (One Mean and One Total w/ Scale)");

  TestHistogram->Fit("DG2Mu2Tot");
  TestHistogram->SetMinimum(0);
  TestHistogram->SetTitle("Lc_MM - Double Gaussian w/ Two Means and Two Totals");
  c1->Write("DG (Two Means and Two Totals)");

  TestHistogram->Fit("DG2Mu1Tot");
  TestHistogram->SetMinimum(0);
  TestHistogram->SetTitle("Lc_MM - Double Gaussian w/ Two Means and One Total");
  c1->Write("DG (Two Means and One Total w/ Scale)");

  Double_t MuAVG = DG2Mu1Tot->GetParameter(0)*DG2Mu1Tot->GetParameter(2) + (1. - DG2Mu1Tot->GetParameter(0))*DG2Mu1Tot->GetParameter(4);
  Double_t MuAVGErr = DG2Mu1Tot->GetParameter(0)*DG2Mu1Tot->GetParError(2) + (1. - DG2Mu1Tot->GetParameter(0))*DG2Mu1Tot->GetParError(4);
  Double_t SigmaAVG = DG2Mu1Tot->GetParameter(0)*DG2Mu1Tot->GetParameter(3) + (1. - DG2Mu1Tot->GetParameter(0))*DG2Mu1Tot->GetParameter(5) + DG2Mu1Tot->GetParameter(0)*(1. - DG2Mu1Tot->GetParameter(0))*(DG2Mu1Tot->GetParameter(2) - DG2Mu1Tot->GetParameter(4))*(DG2Mu1Tot->GetParameter(2) - DG2Mu1Tot->GetParameter(4));
    cout << "Average Mu Value: " << MuAVG << endl;
    cout << "Average Mu Error: " << MuAVGErr << endl;
    cout << "Average Sigma: " << SigmaAVG << endl;
    File->Close();
}
