//Testing Various Gaussian Fit Functions with Pull Plots
//Done on Best Subset of Chained4D Analysis of Chained ntuples
#define ChainedPull_cxx
#include "ChainedPull.h"


#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>

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

TH1D * SGTest = nullptr;
TH1D * DG1M1TTest = nullptr;
TH1D * DG1M2TTest = nullptr;
TH1D * DG2M1TTest = nullptr;
TH1D * DG2M2TTest = nullptr;

TFile * File = nullptr;
TCanvas * c1 = nullptr;

void ChainedPull::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void ChainedPull::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();


         SGTest = new TH1D("Mass [MeV]", "Lc_MM - Single Gaussian", 100, 2185, 2385);
         SGTest->GetXaxis()->SetTitle("MeV");
         SGTest->GetYaxis()->SetTitle("Events Per 2 MeV");

         DG1M1TTest = new TH1D("Mass [MeV]", "Lc_MM - Double Gaussian w/ One Mean and One Total", 100, 2185, 2385);
         DG1M1TTest->GetXaxis()->SetTitle("MeV");
         DG1M1TTest->GetYaxis()->SetTitle("Events Per 2 MeV");

         DG1M2TTest = new TH1D("Mass [MeV]", "Lc_MM - Double Gaussian w/ One Mean and Two Totals", 100, 2185, 2385);
         DG1M2TTest->GetXaxis()->SetTitle("MeV");
         DG1M2TTest->GetYaxis()->SetTitle("Events Per 2 MeV");

         DG2M1TTest = new TH1D("Mass [MeV]", "Lc_MM - Double Gaussian w/ Two Means and One Total", 100, 2185, 2385);
         DG2M1TTest->GetXaxis()->SetTitle("MeV");
         DG2M1TTest->GetYaxis()->SetTitle("Events Per 2 MeV");

         DG2M2TTest = new TH1D("Mass [MeV]", "Lc_MM - Double Gaussian w/ Two Means and Two Totals", 100, 2185, 2385);
         DG2M2TTest->GetXaxis()->SetTitle("MeV");
         DG2M2TTest->GetYaxis()->SetTitle("Events Per 2 MeV");

         File = new TFile("ChainPull.root", "RECREATE");
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

Bool_t ChainedPull::Process(Long64_t entry)
{
  GetEntry(entry);
   fReader.SetLocalEntry(entry);

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

   double FOM = 0;

   if (m < 50 & n < 50 & p < 50 & q < 50){
   FOM = fomKPi[m][n][p][q];
   }

   if (FOM > 1.0 && AdditionalCuts && BorderCut){
   SGTest->Fill(CorrectedLambdaMass);
   DG1M1TTest->Fill(CorrectedLambdaMass);
   DG1M2TTest->Fill(CorrectedLambdaMass);
   DG2M1TTest->Fill(CorrectedLambdaMass);
   DG2M2TTest->Fill(CorrectedLambdaMass);
   }

   return kTRUE;
}

void ChainedPull::SlaveTerminate()
{
}

void ChainedPull::Terminate()
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

  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
   pad2->SetTopMargin(0.03030303);
  pad1->Draw();
  pad2->Draw();

//////////////////////////////////////////////////////////////////////////////////////////

  TF1 *SingleGaussian = new TF1("SingleGaussian",fit2MeV_Gaussian,2100.,2500.,5);
  SingleGaussian->SetParameter(0,15000.);
  SingleGaussian->SetParameter(1,2287);
  SingleGaussian->SetParameter(2, 4);
  SingleGaussian->SetParLimits(2, 0.,20.);
  SingleGaussian->SetParameter(3, 0.);
  SingleGaussian->SetParameter(4, 0.);

pad1->cd();
  SGTest->Fit("SingleGaussian");
  SGTest->SetMinimum(0);

int BinHeight1[100];
double Pullx[100];
int FitHeight1[100];
double Pull1[100];

double count1 = 0;
double count2 = 0;
double count3 = 0;

     for (int bin = 0; bin < 100; bin++){
BinHeight1[bin] = SGTest->GetBinContent(bin + 1);
Pullx[bin] = (bin + 1);
int xvalue = 2186 + 2*(bin);
FitHeight1[bin] = round(SingleGaussian->Eval(xvalue));
Pull1[bin] = (BinHeight1[bin] - FitHeight1[bin])/TMath::Sqrt(FitHeight1[bin]);

if (Pull1[bin] > -1 && Pull1[bin] < 1){
  count1 += 1;
}

if (Pull1[bin] > -2 && Pull1[bin] < 2){
  count2 += 1;
}

if (Pull1[bin] > -3 && Pull1[bin] < 3){
  count3 += 1;
}
}
  cout << "Single Gaussian Fit" << endl;
  cout << "Pulls between -1 & 1 = " << count1 << endl;
  cout << "Pulls between -2 & 2 = " << count2 << endl;
  cout << "Pulls between -3 & 3 = " << count3 << endl;

pad2->cd();
TGraph* PullPlot1 = new TGraph(100,Pullx, Pull1);
PullPlot1->GetXaxis()->SetLimits(0,100);
PullPlot1->GetXaxis()->SetTickLength(0.);
PullPlot1->GetYaxis()->SetTickLength(0.);
PullPlot1->SetFillColor(38);
PullPlot1->GetYaxis()->SetTitle("Pull");
PullPlot1->GetYaxis()->CenterTitle();
PullPlot1->GetYaxis()->SetTitleSize(0.10);
PullPlot1->GetYaxis()->SetTitleOffset(0.2);
PullPlot1->GetXaxis()->SetLabelSize(0);
PullPlot1->GetYaxis()->SetLabelFont(42);
PullPlot1->GetYaxis()->SetLabelSize(0.06);
PullPlot1->SetTitle("");
PullPlot1->SetMinimum(-5);
PullPlot1->SetMaximum(5);
PullPlot1->Draw("AB");

 c1->Write("Single Gaussian");
////////////////////////////////////////////////////

TF1 *DG1Mu2Tot = new TF1("DG1Mu2Tot",DGOneMuTwoTotal,2100.,2500.,7);
DG1Mu2Tot->SetParameter(0, 2000.);
DG1Mu2Tot->SetParameter(1, 2286);
DG1Mu2Tot->SetParameter(2, 7);
DG1Mu2Tot->SetParLimits(2, 0., 20.);
DG1Mu2Tot->SetParameter(3, 4000.);
DG1Mu2Tot->SetParameter(4, 4);
DG1Mu2Tot->SetParLimits(4, 0., 20.);
DG1Mu2Tot->SetParameter(5, 0.);
DG1Mu2Tot->SetParameter(6, 0.);

pad1->cd();
DG1M2TTest->Fit("DG1Mu2Tot");
DG1M2TTest->SetMinimum(0);

int BinHeight2[100];
int FitHeight2[100];
double Pull2[100];

 count1 = 0;
 count2 = 0;
 count3 = 0;

for (int bin = 0; bin < 100; bin++){
BinHeight2[bin] = DG1M2TTest->GetBinContent(bin + 1);
int xvalue = 2186 + 2*(bin);
FitHeight2[bin] = round(DG1Mu2Tot->Eval(xvalue));
Pull2[bin] = (BinHeight2[bin] - FitHeight2[bin])/TMath::Sqrt(FitHeight2[bin]);

if (Pull2[bin] > -1 && Pull2[bin] < 1){
count1 += 1;
}

if (Pull2[bin] > -2 && Pull2[bin] < 2){
count2 += 1;
}

if (Pull2[bin] > -3 && Pull2[bin] < 3){
count3 += 1;
}
}

cout << "DG1Mu2Tot Fit" << endl;
cout << "Pulls between -1 & 1 = " << count1 << endl;
cout << "Pulls between -2 & 2 = " << count2 << endl;
cout << "Pulls between -3 & 3 = " << count3 << endl;

pad2->cd();
TGraph* PullPlot2 = new TGraph(100,Pullx, Pull2);
PullPlot2->GetXaxis()->SetLimits(0,100);
PullPlot2->GetXaxis()->SetTickLength(0.);
PullPlot2->GetYaxis()->SetTickLength(0.);
PullPlot2->SetFillColor(38);
PullPlot2->GetYaxis()->SetTitle("Pull");
PullPlot2->GetYaxis()->CenterTitle();
PullPlot2->GetYaxis()->SetTitleSize(0.10);
PullPlot2->GetYaxis()->SetTitleOffset(0.2);
PullPlot2->GetXaxis()->SetLabelSize(0);
PullPlot2->GetYaxis()->SetLabelFont(42);
PullPlot2->GetYaxis()->SetLabelSize(0.06);
PullPlot2->SetTitle("");
PullPlot2->SetMinimum(-5);
PullPlot2->SetMaximum(5);
PullPlot2->Draw("AB");
c1->Write("DG w/ 1 Mu & 2 Totals");

////////////////////////////////////////////////////////

TF1 *DG1Mu1Tot = new TF1("DG1Mu1Tot",DGOneMuOneTotal,2100.,2500.,7);
DG1Mu1Tot->SetParameter(0, 0.2);
DG1Mu1Tot->SetParameter(1, 20000);
DG1Mu1Tot->SetParameter(2, 2287.);
DG1Mu1Tot->SetParameter(3, 7);
DG1Mu1Tot->SetParameter(4, 4);
DG1Mu1Tot->SetParLimits(3, 0., 20.);
DG1Mu1Tot->SetParLimits(4, 0., 20.);
DG1Mu1Tot->SetParameter(5, 0.);
DG1Mu1Tot->SetParameter(6, 0.);

pad1->cd();
DG1M1TTest->Fit("DG1Mu1Tot");
DG1M1TTest->SetMinimum(0);

int BinHeight3[100];
int FitHeight3[100];
double Pull3[100];

 count1 = 0;
 count2 = 0;
 count3 = 0;

for (int bin = 0; bin < 100; bin++){
BinHeight3[bin] = DG1M1TTest->GetBinContent(bin + 1);
int xvalue = 2186 + 2*(bin);
FitHeight3[bin] = round(DG1Mu1Tot->Eval(xvalue));
Pull3[bin] = (BinHeight3[bin] - FitHeight3[bin])/TMath::Sqrt(FitHeight3[bin]);

if (Pull3[bin] > -1 && Pull3[bin] < 1){
count1 += 1;
}

if (Pull3[bin] > -2 && Pull3[bin] < 2){
count2 += 1;
}

if (Pull3[bin] > -3 && Pull3[bin] < 3){
count3 += 1;
}
}

cout << "DG1Mu1Tot Fit" << endl;
cout << "Pulls between -1 & 1 = " << count1 << endl;
cout << "Pulls between -2 & 2 = " << count2 << endl;
cout << "Pulls between -3 & 3 = " << count3 << endl;

pad2->cd();
TGraph* PullPlot3 = new TGraph(100,Pullx, Pull3);
PullPlot3->GetXaxis()->SetLimits(0,100);
PullPlot3->GetXaxis()->SetTickLength(0.);
PullPlot3->GetYaxis()->SetTickLength(0.);
PullPlot3->SetFillColor(38);
PullPlot3->GetYaxis()->SetTitle("Pull");
PullPlot3->GetYaxis()->CenterTitle();
PullPlot3->GetYaxis()->SetTitleSize(0.10);
PullPlot3->GetYaxis()->SetTitleOffset(0.2);
PullPlot3->GetXaxis()->SetLabelSize(0);
PullPlot3->GetYaxis()->SetLabelFont(42);
PullPlot3->GetYaxis()->SetLabelSize(0.06);
PullPlot3->SetTitle("");
PullPlot3->SetMinimum(-5);
PullPlot3->SetMaximum(5);
PullPlot3->Draw("AB");
c1->Write("DG w/ 1 Mu & 1 Total");

///////////////////////////////////////////////////////

TF1 *DG2Mu2Tot = new TF1("DG2Mu2Tot",DGTwoMuTwoTotal,2100.,2500.,8);
DG2Mu2Tot->SetParameter(0, 4000.);
DG2Mu2Tot->SetParameter(1, 2286);
DG2Mu2Tot->SetParameter(2, 4);
DG2Mu2Tot->SetParLimits(2, 0., 20.);
DG2Mu2Tot->SetParameter(3, 4000.);
DG2Mu2Tot->SetParameter(4, 2286);
DG2Mu2Tot->SetParameter(5, 4);
DG2Mu2Tot->SetParLimits(5, 0., 20.);
DG2Mu2Tot->SetParameter(6, 0.);
DG2Mu2Tot->SetParameter(7, 0.);

pad1->cd();
DG2M2TTest->Fit("DG2Mu2Tot");
DG2M2TTest->SetMinimum(0);

int BinHeight4[100];
int FitHeight4[100];
double Pull4[100];

 count1 = 0;
 count2 = 0;
 count3 = 0;

for (int bin = 0; bin < 100; bin++){
BinHeight4[bin] = DG2M2TTest->GetBinContent(bin + 1);
int xvalue = 2186 + 2*(bin);
FitHeight4[bin] = round(DG2Mu2Tot->Eval(xvalue));
Pull4[bin] = (BinHeight4[bin] - FitHeight4[bin])/TMath::Sqrt(FitHeight4[bin]);

if (Pull4[bin] > -1 && Pull4[bin] < 1){
count1 += 1;
}

if (Pull4[bin] > -2 && Pull4[bin] < 2){
count2 += 1;
}

if (Pull4[bin] > -3 && Pull4[bin] < 3){
count3 += 1;
}
}

cout << "DG2Mu2Tot Fit" << endl;
cout << "Pulls between -1 & 1 = " << count1 << endl;
cout << "Pulls between -2 & 2 = " << count2 << endl;
cout << "Pulls between -3 & 3 = " << count3 << endl;

pad2->cd();
TGraph* PullPlot4 = new TGraph(100,Pullx, Pull4);
PullPlot4->GetXaxis()->SetLimits(0,100);
PullPlot4->GetXaxis()->SetTickLength(0.);
PullPlot4->GetYaxis()->SetTickLength(0.);
PullPlot4->SetFillColor(38);
PullPlot4->GetYaxis()->SetTitle("Pull");
PullPlot4->GetYaxis()->CenterTitle();
PullPlot4->GetYaxis()->SetTitleSize(0.10);
PullPlot4->GetYaxis()->SetTitleOffset(0.2);
PullPlot4->GetXaxis()->SetLabelSize(0);
PullPlot4->GetYaxis()->SetLabelFont(42);
PullPlot4->GetYaxis()->SetLabelSize(0.06);
PullPlot4->SetTitle("");
PullPlot4->SetMinimum(-5);
PullPlot4->SetMaximum(5);
PullPlot4->Draw("AB");
c1->Write("DG w/ 2 Mu & 2 Totals");

///////////////////////////////////////////////////////

TF1 *DG2Mu1Tot = new TF1("DG2Mu1Tot",DGTwoMuOneTotal,2100.,2500.,8);
DG2Mu1Tot->SetParameter(0, 0.2);
DG2Mu1Tot->SetParameter(1, 20000);
DG2Mu1Tot->SetParameter(2, 2286.);
DG2Mu1Tot->SetParameter(3, 7);
DG2Mu1Tot->SetParameter(4, 2287);
DG2Mu1Tot->SetParameter(5, 4);
DG2Mu1Tot->SetParLimits(3, 0., 20.);
DG2Mu1Tot->SetParLimits(5, 0., 20.);
DG2Mu1Tot->SetParameter(6, 0.);
DG2Mu1Tot->SetParameter(7, 0.);

pad1->cd();
DG2M1TTest->Fit("DG2Mu1Tot");
DG2M1TTest->SetMinimum(0);

int BinHeight5[100];
int FitHeight5[100];
double Pull5[100];

 count1 = 0;
 count2 = 0;
 count3 = 0;

for (int bin = 0; bin < 100; bin++){
BinHeight5[bin] = DG2M1TTest->GetBinContent(bin + 1);
int xvalue = 2186 + 2*(bin);
FitHeight5[bin] = round(DG2Mu1Tot->Eval(xvalue));
Pull5[bin] = (BinHeight5[bin] - FitHeight5[bin])/TMath::Sqrt(FitHeight5[bin]);

if (Pull5[bin] > -1 && Pull5[bin] < 1){
count1 += 1;
}

if (Pull5[bin] > -2 && Pull5[bin] < 2){
count2 += 1;
}

if (Pull5[bin] > -3 && Pull5[bin] < 3){
count3 += 1;
}
}

cout << "DG2Mu1Tot Fit" << endl;
cout << "Pulls between -1 & 1 = " << count1 << endl;
cout << "Pulls between -2 & 2 = " << count2 << endl;
cout << "Pulls between -3 & 3 = " << count3 << endl;

pad2->cd();
TGraph* PullPlot5 = new TGraph(100,Pullx, Pull5);
PullPlot5->GetXaxis()->SetLimits(0,100);
PullPlot5->GetXaxis()->SetTickLength(0.);
PullPlot5->GetYaxis()->SetTickLength(0.);
PullPlot5->SetFillColor(38);
PullPlot5->GetYaxis()->SetTitle("Pull");
PullPlot5->GetYaxis()->CenterTitle();
PullPlot5->GetYaxis()->SetTitleSize(0.10);
PullPlot5->GetYaxis()->SetTitleOffset(0.2);
PullPlot5->GetXaxis()->SetLabelSize(0);
PullPlot5->GetYaxis()->SetLabelFont(42);
PullPlot5->GetYaxis()->SetLabelSize(0.06);
PullPlot5->SetTitle("");
PullPlot5->SetMinimum(-5);
PullPlot5->SetMaximum(5);
PullPlot5->Draw("AB");
c1->Write("DG w/ 2 Mu & 1 Total");

///////////////////////////////////////////////////////

File->Close();
}
