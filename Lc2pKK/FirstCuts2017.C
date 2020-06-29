//Initial Cuts with Various Fits and Pull Plots

#define FirstCuts2017_cxx
#include "FirstCuts2017.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit1MeV_Gaussian.C"
#include "fit2MeV_Gaussian.C"
#include "DGOneMuTwoTotals1MeV.C"
#include "DGOneMuTwoTotals2MeV.C"

TH1D * MassHistLoose2MeV = nullptr;
TH1D * MassHistMid2MeV = nullptr;
TH1D * MassHistTight2MeV = nullptr;

TH1D * MassHistLoose1MeV = nullptr;
TH1D * MassHistMid1MeV = nullptr;
TH1D * MassHistTight1MeV = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;
TCanvas * ex1 = nullptr;

void FirstCuts2017::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   MassHistLoose2MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistLoose2MeV->GetXaxis()->SetTitle("MeV");
   MassHistLoose2MeV->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistMid2MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistMid2MeV->GetXaxis()->SetTitle("MeV");
   MassHistMid2MeV->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistTight2MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistTight2MeV->GetXaxis()->SetTitle("MeV");
   MassHistTight2MeV->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistLoose1MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 150, 2210, 2360);
   MassHistLoose1MeV->GetXaxis()->SetTitle("MeV");
   MassHistLoose1MeV->GetYaxis()->SetTitle("Events Per 1 MeV");

   MassHistMid1MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 150, 2210, 2360);
   MassHistMid1MeV->GetXaxis()->SetTitle("MeV");
   MassHistMid1MeV->GetYaxis()->SetTitle("Events Per 1 MeV");

   MassHistTight1MeV = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 150, 2210, 2360);
   MassHistTight1MeV->GetXaxis()->SetTitle("MeV");
   MassHistTight1MeV->GetYaxis()->SetTitle("Events Per 1 MeV");

   File = new TFile("MassHist2017.root", "RECREATE");
  gFile = File;

   c1 = new TCanvas("canvas", "Test Canvas");
   ex1 = new TCanvas("ex1","Latex",500,600);
}

void FirstCuts2017::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t FirstCuts2017::Process(Long64_t entry)
{

    GetEntry(entry);
   fReader.SetLocalEntry(entry);

   bool LooseCuts = (
      (((*Kminus_ProbNNk)*(*Kplus_ProbNNk)*(*Proton_ProbNNp)) > 0.85)
   && (*Proton_PT > 850)
   && (*Kminus_PT > 700)
   && (*Kplus_PT > 700)
   );

   bool MidCuts = (
      (((*Kminus_ProbNNk)*(*Kplus_ProbNNk)*(*Proton_ProbNNp)) > 0.9)
   && (*Proton_PT > 850)
   && (*Kminus_PT > 700)
   && (*Kplus_PT > 700)
   );

   bool TightCuts = (
      (((*Kminus_ProbNNk)*(*Kplus_ProbNNk)*(*Proton_ProbNNp)) > 0.95)
   && (*Proton_PT > 850)
   && (*Kminus_PT > 700)
   && (*Kplus_PT > 700)
   );

   if (MidCuts){
     MassHistMid2MeV->Fill(*Lcplus_M);
     MassHistMid1MeV->Fill(*Lcplus_M);
   }

   if (LooseCuts){
     MassHistLoose2MeV->Fill(*Lcplus_M);
     MassHistLoose1MeV->Fill(*Lcplus_M);
   }

   if (TightCuts){
     MassHistTight2MeV->Fill(*Lcplus_M);
     MassHistTight1MeV->Fill(*Lcplus_M);
   }

   return kTRUE;
}

void FirstCuts2017::SlaveTerminate()
{
}

void FirstCuts2017::Terminate()
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

c1->cd();

double Pullx2MeV[75];
double Pullx1MeV[150];

TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
pad2->SetTopMargin(0.03030303);
pad1->Draw();
pad2->Draw();

/////////////////////////////

//TF1 *GaussianLoose2MeV = new TF1("GaussianLoose2MeV",DGOneMuTwoTotals2MeV,2200.,2400.,7);
//GaussianLoose2MeV->SetParameter(0, 2000.);
//GaussianLoose2MeV->SetParLimits(0, 0, 25000);
//GaussianLoose2MeV->SetParameter(1, 2286);
//GaussianLoose2MeV->SetParameter(2, 5);
//GaussianLoose2MeV->SetParLimits(2, 0., 20.);
//GaussianLoose2MeV->SetParameter(3, 2000.);
//GaussianLoose2MeV->SetParLimits(3, 0, 25000);
//GaussianLoose2MeV->SetParameter(4, 5);
//GaussianLoose2MeV->SetParLimits(4, 0., 20.);
//GaussianLoose2MeV->SetParameter(5, 0.);
//GaussianLoose2MeV->SetParameter(6, 0.);

TF1 *GaussianLoose2MeV = new TF1("GaussianLoose2MeV",fit2MeV_Gaussian,2100.,2500.,5);
GaussianLoose2MeV->SetParameter(0,400.);
GaussianLoose2MeV->SetParameter(1,2287);
GaussianLoose2MeV->SetParameter(2, 6);
GaussianLoose2MeV->SetParLimits(2, 0.,20.);
GaussianLoose2MeV->SetParameter(3, 0.);
GaussianLoose2MeV->SetParameter(4, 0.);

pad1->cd();
MassHistLoose2MeV->SetMinimum(0);
MassHistLoose2MeV->Fit("GaussianLoose2MeV");

int BinHeightLoose2MeV[75];
int FitHeightLoose2MeV[75];
double PullLoose2MeV[75];

double Loose2MeVcount1 = 0;
double Loose2MeVcount2 = 0;
double Loose2MeVcount3 = 0;

for (int bin = 0; bin < 75; bin++){
BinHeightLoose2MeV[bin] = MassHistLoose2MeV->GetBinContent(bin + 1);
Pullx2MeV[bin] = (bin + 1);
int xvalue2MeV = 2211 + 2*(bin);
FitHeightLoose2MeV[bin] = round(GaussianLoose2MeV->Eval(xvalue2MeV));
PullLoose2MeV[bin] = (BinHeightLoose2MeV[bin] - FitHeightLoose2MeV[bin])/TMath::Sqrt(FitHeightLoose2MeV[bin]);

if (PullLoose2MeV[bin] > -1 && PullLoose2MeV[bin] < 1){
  Loose2MeVcount1 += 1;
}

if (PullLoose2MeV[bin] > -2 && PullLoose2MeV[bin] < 2){
  Loose2MeVcount2 += 1;
}

if (PullLoose2MeV[bin] > -3 && PullLoose2MeV[bin] < 3){
  Loose2MeVcount3 += 1;
}
}

pad2->cd();
TGraph* PullPlotLoose2MeV = new TGraph(75, Pullx2MeV, PullLoose2MeV);
PullPlotLoose2MeV->GetXaxis()->SetLimits(0.5,75.5);
PullPlotLoose2MeV->GetXaxis()->SetTickLength(0.);
PullPlotLoose2MeV->GetYaxis()->SetTickLength(0.);
PullPlotLoose2MeV->SetFillColor(38);
PullPlotLoose2MeV->GetYaxis()->SetTitle("Pull");
PullPlotLoose2MeV->GetYaxis()->CenterTitle();
PullPlotLoose2MeV->GetYaxis()->SetTitleSize(0.10);
PullPlotLoose2MeV->GetYaxis()->SetTitleOffset(0.2);
PullPlotLoose2MeV->GetXaxis()->SetLabelSize(0);
PullPlotLoose2MeV->GetYaxis()->SetLabelFont(42);
PullPlotLoose2MeV->GetYaxis()->SetLabelSize(0.06);
PullPlotLoose2MeV->SetTitle("");
PullPlotLoose2MeV->SetMinimum(-5);
PullPlotLoose2MeV->SetMaximum(5);
PullPlotLoose2MeV->Draw("AB");
      c1->Write("Lc Mass - Loose/2MeV");

ex1->cd();
   TLatex Tl;
   Tl.SetTextAlign(12);
   Tl.SetTextSize(0.04);
   Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianLoose2MeV->GetParameter(0)));
   Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianLoose2MeV->GetParError(0)));
   Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianLoose2MeV->GetParameter(1)));
   Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianLoose2MeV->GetParError(1)));
   Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianLoose2MeV->GetParameter(2)));
   Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianLoose2MeV->GetParError(2)));
   Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Loose2MeVcount1));
   Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Loose2MeVcount2));
   Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Loose2MeVcount3));
   ex1->Write("Loose Fit Values/2MeV");
c1->cd();
///////////////////////////////////////////////////////

//TF1 *GaussianMid2MeV = new TF1("GaussianMid2MeV",DGOneMuTwoTotals2MeV,2200.,2400.,7);
//GaussianMid2MeV->SetParameter(0, 2000.);
//GaussianMid2MeV->SetParLimits(0, 0, 25000);
//GaussianMid2MeV->SetParameter(1, 2286);
//GaussianMid2MeV->SetParameter(2, 5);
//GaussianMid2MeV->SetParLimits(2, 0., 20.);
//GaussianMid2MeV->SetParameter(3, 2000.);
//GaussianMid2MeV->SetParLimits(3, 0, 25000);
//GaussianMid2MeV->SetParameter(4, 5);
//GaussianMid2MeV->SetParLimits(4, 0., 20.);
//GaussianMid2MeV->SetParameter(5, 0.);
//GaussianMid2MeV->SetParameter(6, 0.);

TF1 *GaussianMid2MeV = new TF1("GaussianMid2MeV",fit2MeV_Gaussian,2100.,2500.,5);
GaussianMid2MeV->SetParameter(0,400.);
GaussianMid2MeV->SetParameter(1,2287);
GaussianMid2MeV->SetParameter(2, 6);
GaussianMid2MeV->SetParLimits(2, 0.,20.);
GaussianMid2MeV->SetParameter(3, 0.);
GaussianMid2MeV->SetParameter(4, 0.);

  pad1->cd();
  MassHistMid2MeV->SetMinimum(0);
  MassHistMid2MeV->Fit("GaussianMid2MeV");

  int BinHeightMid2MeV[75];
  int FitHeightMid2MeV[75];
  double PullMid2MeV[75];

double   Mid2MeVcount1 = 0;
double   Mid2MeVcount2 = 0;
double   Mid2MeVcount3 = 0;

  for (int bin = 0; bin < 75; bin++){
  BinHeightMid2MeV[bin] = MassHistMid2MeV->GetBinContent(bin + 1);
  int xvalue2MeV = 2211 + 2*(bin);
  FitHeightMid2MeV[bin] = round(GaussianMid2MeV->Eval(xvalue2MeV));
  PullMid2MeV[bin] = (BinHeightMid2MeV[bin] - FitHeightMid2MeV[bin])/TMath::Sqrt(FitHeightMid2MeV[bin]);

  if (PullMid2MeV[bin] > -1 && PullMid2MeV[bin] < 1){
  Mid2MeVcount1 += 1;
  }

  if (PullMid2MeV[bin] > -2 && PullMid2MeV[bin] < 2){
  Mid2MeVcount2 += 1;
  }

  if (PullMid2MeV[bin] > -3 && PullMid2MeV[bin] < 3){
  Mid2MeVcount3 += 1;
  }
  }

  pad2->cd();
TGraph* PullPlotMid2MeV = new TGraph(75, Pullx2MeV, PullMid2MeV);
PullPlotMid2MeV->GetXaxis()->SetLimits(0.5,75.5);
PullPlotMid2MeV->GetXaxis()->SetTickLength(0.);
PullPlotMid2MeV->GetYaxis()->SetTickLength(0.);
PullPlotMid2MeV->SetFillColor(38);
PullPlotMid2MeV->GetYaxis()->SetTitle("Pull");
PullPlotMid2MeV->GetYaxis()->CenterTitle();
PullPlotMid2MeV->GetYaxis()->SetTitleSize(0.10);
PullPlotMid2MeV->GetYaxis()->SetTitleOffset(0.2);
PullPlotMid2MeV->GetXaxis()->SetLabelSize(0);
PullPlotMid2MeV->GetYaxis()->SetLabelFont(42);
PullPlotMid2MeV->GetYaxis()->SetLabelSize(0.06);
PullPlotMid2MeV->SetTitle("");
PullPlotMid2MeV->SetMinimum(-5);
PullPlotMid2MeV->SetMaximum(5);
PullPlotMid2MeV->Draw("AB");

      c1->Write("Lc Mass - Mid/2MeV");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianMid2MeV->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianMid2MeV->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianMid2MeV->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianMid2MeV->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianMid2MeV->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianMid2MeV->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Mid2MeVcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Mid2MeVcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Mid2MeVcount3));
   ex1->Write("Mid Fit Values/2MeV");
c1->cd();
//////////////////////////////////////////////////////

//TF1 *GaussianTight2MeV = new TF1("GaussianTight2MeV",DGOneMuTwoTotals2MeVs2MeV,2200.,2400.,7);
//GaussianTight2MeV->SetParameter(0, 2000.);
//GaussianTight2MeV->SetParLimits(0, 0, 25000);
//GaussianTight2MeV->SetParameter(1, 2286);
//GaussianTight2MeV->SetParameter(2, 5);
//GaussianTight2MeV->SetParLimits(2, 0., 20.);
//GaussianTight2MeV->SetParameter(3, 2000.);
//GaussianTight2MeV->SetParLimits(3, 0, 25000);
//GaussianTight2MeV->SetParameter(4, 5);
//GaussianTight2MeV->SetParLimits(4, 0., 20.);
//GaussianTight2MeV->SetParameter(5, 0.);
//GaussianTight2MeV->SetParameter(6, 0.);

TF1 *GaussianTight2MeV = new TF1("GaussianTight2MeV",fit2MeV_Gaussian,2100.,2500.,5);
GaussianTight2MeV->SetParameter(0,400.);
GaussianTight2MeV->SetParameter(1,2286);
GaussianTight2MeV->SetParameter(2, 6);
GaussianTight2MeV->SetParLimits(2, 0.,20.);
GaussianTight2MeV->SetParameter(3, 0.);
GaussianTight2MeV->SetParameter(4, 0.);

pad1->cd();
MassHistTight2MeV->SetMinimum(0);
MassHistTight2MeV->Fit("GaussianTight2MeV");

int BinHeightTight2MeV[75];
int FitHeightTight2MeV[75];
double PullTight2MeV[75];

double Tight2MeVcount1 = 0;
double Tight2MeVcount2 = 0;
double Tight2MeVcount3 = 0;

for (int bin = 0; bin < 75; bin++){
BinHeightTight2MeV[bin] = MassHistTight2MeV->GetBinContent(bin + 1);
int xvalue2MeV = 2211 + 2*(bin);
FitHeightTight2MeV[bin] = round(GaussianTight2MeV->Eval(xvalue2MeV));

PullTight2MeV[bin] = (BinHeightTight2MeV[bin] - FitHeightTight2MeV[bin])/TMath::Sqrt(FitHeightTight2MeV[bin]);

if (PullTight2MeV[bin] > -1 && PullTight2MeV[bin] < 1){
Tight2MeVcount1 += 1;
}

if (PullTight2MeV[bin] > -2 && PullTight2MeV[bin] < 2){
Tight2MeVcount2 += 1;
}

if (PullTight2MeV[bin] > -3 && PullTight2MeV[bin] < 3){
Tight2MeVcount3 += 1;
}
}
     pad2->cd();
     TGraph* PullPlotTight2MeV = new TGraph(75, Pullx2MeV, PullTight2MeV);
     PullPlotTight2MeV->GetXaxis()->SetLimits(0.5,75.5);
     PullPlotTight2MeV->GetXaxis()->SetTickLength(0.);
     PullPlotTight2MeV->GetYaxis()->SetTickLength(0.);
     PullPlotTight2MeV->SetFillColor(38);
     PullPlotTight2MeV->GetYaxis()->SetTitle("Pull");
     PullPlotTight2MeV->GetYaxis()->CenterTitle();
     PullPlotTight2MeV->GetYaxis()->SetTitleSize(0.10);
     PullPlotTight2MeV->GetYaxis()->SetTitleOffset(0.2);
     PullPlotTight2MeV->GetXaxis()->SetLabelSize(0);
     PullPlotTight2MeV->GetYaxis()->SetLabelFont(42);
     PullPlotTight2MeV->GetYaxis()->SetLabelSize(0.06);
     PullPlotTight2MeV->SetTitle("");
     PullPlotTight2MeV->SetMinimum(-5);
     PullPlotTight2MeV->SetMaximum(5);
     PullPlotTight2MeV->Draw("AB");
     c1->Write("Lc Mass - Tight/2MeV");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianTight2MeV->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianTight2MeV->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianTight2MeV->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianTight2MeV->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianTight2MeV->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianTight2MeV->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Tight2MeVcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Tight2MeVcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Tight2MeVcount3));
   ex1->Write("Tight Fit Values/2MeV");
c1->cd();
/////////////////////////////

//TF1 *GaussianLoose1MeV = new TF1("GaussianLoose1MeV",DGOneMuTwoTotals1MeV,2200.,2400.,7);
//GaussianLoose1MeV->SetParameter(0, 2000.);
//GaussianLoose1MeV->SetParLimits(0, 0, 25000);
//GaussianLoose1MeV->SetParameter(1, 2286);
//GaussianLoose1MeV->SetParameter(2, 5);
//GaussianLoose1MeV->SetParLimits(2, 0., 20.);
//GaussianLoose1MeV->SetParameter(3, 2000.);
//GaussianLoose1MeV->SetParLimits(3, 0, 25000);
//GaussianLoose1MeV->SetParameter(4, 5);
//GaussianLoose1MeV->SetParLimits(4, 0., 20.);
//GaussianLoose1MeV->SetParameter(5, 0.);
//GaussianLoose1MeV->SetParameter(6, 0.);

TF1 *GaussianLoose1MeV = new TF1("GaussianLoose1MeV",fit1MeV_Gaussian,2100.,2500.,5);
GaussianLoose1MeV->SetParameter(0,400.);
GaussianLoose1MeV->SetParameter(1,2287);
GaussianLoose1MeV->SetParameter(2, 6);
GaussianLoose1MeV->SetParLimits(2, 0.,20.);
GaussianLoose1MeV->SetParameter(3, 0.);
GaussianLoose1MeV->SetParameter(4, 0.);

pad1->cd();
MassHistLoose1MeV->SetMinimum(0);
MassHistLoose1MeV->Fit("GaussianLoose1MeV");

int BinHeightLoose1MeV[150];
int FitHeightLoose1MeV[150];
double PullLoose1MeV[150];

double Loose1MeVcount1 = 0;
double Loose1MeVcount2 = 0;
double Loose1MeVcount3 = 0;

for (int bin = 0; bin < 150; bin++){
BinHeightLoose1MeV[bin] = MassHistLoose1MeV->GetBinContent(bin + 1);
Pullx1MeV[bin] = (bin + 1);
int xvalue1MeV = 2210.5 + 1*(bin);
FitHeightLoose1MeV[bin] = round(GaussianLoose1MeV->Eval(xvalue1MeV));
PullLoose1MeV[bin] = (BinHeightLoose1MeV[bin] - FitHeightLoose1MeV[bin])/TMath::Sqrt(FitHeightLoose1MeV[bin]);

if (PullLoose1MeV[bin] > -1 && PullLoose1MeV[bin] < 1){
  Loose1MeVcount1 += 1;
}

if (PullLoose1MeV[bin] > -2 && PullLoose1MeV[bin] < 2){
  Loose1MeVcount2 += 1;
}

if (PullLoose1MeV[bin] > -3 && PullLoose1MeV[bin] < 3){
  Loose1MeVcount3 += 1;
}
}

pad2->cd();
TGraph* PullPlotLoose1MeV = new TGraph(150, Pullx1MeV, PullLoose1MeV);
PullPlotLoose1MeV->GetXaxis()->SetLimits(0.5,150.5);
PullPlotLoose1MeV->GetXaxis()->SetTickLength(0.);
PullPlotLoose1MeV->GetYaxis()->SetTickLength(0.);
PullPlotLoose1MeV->SetFillColor(38);
PullPlotLoose1MeV->GetYaxis()->SetTitle("Pull");
PullPlotLoose1MeV->GetYaxis()->CenterTitle();
PullPlotLoose1MeV->GetYaxis()->SetTitleSize(0.10);
PullPlotLoose1MeV->GetYaxis()->SetTitleOffset(0.2);
PullPlotLoose1MeV->GetXaxis()->SetLabelSize(0);
PullPlotLoose1MeV->GetYaxis()->SetLabelFont(42);
PullPlotLoose1MeV->GetYaxis()->SetLabelSize(0.06);
PullPlotLoose1MeV->SetTitle("");
PullPlotLoose1MeV->SetMinimum(-5);
PullPlotLoose1MeV->SetMaximum(5);
PullPlotLoose1MeV->Draw("AB");
      c1->Write("Lc Mass - Loose/1MeV");

ex1->cd();
   TLatex Tl;
   Tl.SetTextAlign(12);
   Tl.SetTextSize(0.04);
   Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianLoose1MeV->GetParameter(0)));
   Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianLoose1MeV->GetParError(0)));
   Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianLoose1MeV->GetParameter(1)));
   Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianLoose1MeV->GetParError(1)));
   Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianLoose1MeV->GetParameter(2)));
   Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianLoose1MeV->GetParError(2)));
   Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Loose1MeVcount1));
   Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Loose1MeVcount2));
   Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Loose1MeVcount3));
   ex1->Write("Loose Fit Values/1MeV");
c1->cd();
///////////////////////////////////////////////////////

//TF1 *GaussianMid1MeV = new TF1("GaussianMid1MeV",DGOneMuTwoTotals1MeV,2200.,2400.,7);
//GaussianMid1MeV->SetParameter(0, 2000.);
//GaussianMid1MeV->SetParLimits(0, 0, 25000);
//GaussianMid1MeV->SetParameter(1, 2286);
//GaussianMid1MeV->SetParameter(2, 5);
//GaussianMid1MeV->SetParLimits(2, 0., 20.);
//GaussianMid1MeV->SetParameter(3, 2000.);
//GaussianMid1MeV->SetParLimits(3, 0, 25000);
//GaussianMid1MeV->SetParameter(4, 5);
//GaussianMid1MeV->SetParLimits(4, 0., 20.);
//GaussianMid1MeV->SetParameter(5, 0.);
//GaussianMid1MeV->SetParameter(6, 0.);

TF1 *GaussianMid1MeV = new TF1("GaussianMid1MeV",fit1MeV_Gaussian,2100.,2500.,5);
GaussianMid1MeV->SetParameter(0,400.);
GaussianMid1MeV->SetParameter(1,2287);
GaussianMid1MeV->SetParameter(2, 6);
GaussianMid1MeV->SetParLimits(2, 0.,20.);
GaussianMid1MeV->SetParameter(3, 0.);
GaussianMid1MeV->SetParameter(4, 0.);

  pad1->cd();
  MassHistMid1MeV->SetMinimum(0);
  MassHistMid1MeV->Fit("GaussianMid1MeV");

  int BinHeightMid1MeV[150];
  int FitHeightMid1MeV[150];
  double PullMid1MeV[150];

double   Mid1MeVcount1 = 0;
double   Mid1MeVcount2 = 0;
double   Mid1MeVcount3 = 0;

  for (int bin = 0; bin < 150; bin++){
  BinHeightMid1MeV[bin] = MassHistMid1MeV->GetBinContent(bin + 1);
  int xvalue1MeV = 2210.5 + 1*(bin);
  FitHeightMid1MeV[bin] = round(GaussianMid1MeV->Eval(xvalue1MeV));
  PullMid1MeV[bin] = (BinHeightMid1MeV[bin] - FitHeightMid1MeV[bin])/TMath::Sqrt(FitHeightMid1MeV[bin]);

  if (PullMid1MeV[bin] > -1 && PullMid1MeV[bin] < 1){
  Mid1MeVcount1 += 1;
  }

  if (PullMid1MeV[bin] > -2 && PullMid1MeV[bin] < 2){
  Mid1MeVcount2 += 1;
  }

  if (PullMid1MeV[bin] > -3 && PullMid1MeV[bin] < 3){
  Mid1MeVcount3 += 1;
  }
  }

  pad2->cd();
TGraph* PullPlotMid1MeV = new TGraph(150, Pullx1MeV, PullMid1MeV);
PullPlotMid1MeV->GetXaxis()->SetLimits(0.5,150.5);
PullPlotMid1MeV->GetXaxis()->SetTickLength(0.);
PullPlotMid1MeV->GetYaxis()->SetTickLength(0.);
PullPlotMid1MeV->SetFillColor(38);
PullPlotMid1MeV->GetYaxis()->SetTitle("Pull");
PullPlotMid1MeV->GetYaxis()->CenterTitle();
PullPlotMid1MeV->GetYaxis()->SetTitleSize(0.10);
PullPlotMid1MeV->GetYaxis()->SetTitleOffset(0.2);
PullPlotMid1MeV->GetXaxis()->SetLabelSize(0);
PullPlotMid1MeV->GetYaxis()->SetLabelFont(42);
PullPlotMid1MeV->GetYaxis()->SetLabelSize(0.06);
PullPlotMid1MeV->SetTitle("");
PullPlotMid1MeV->SetMinimum(-5);
PullPlotMid1MeV->SetMaximum(5);
PullPlotMid1MeV->Draw("AB");

      c1->Write("Lc Mass - Mid/1MeV");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianMid1MeV->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianMid1MeV->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianMid1MeV->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianMid1MeV->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianMid1MeV->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianMid1MeV->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Mid1MeVcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Mid1MeVcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Mid1MeVcount3));
   ex1->Write("Mid Fit Values/1MeV");
c1->cd();
//////////////////////////////////////////////////////

//TF1 *GaussianTight1MeV = new TF1("GaussianTight1MeV",DGOneMuTwoTotals1MeVs1MeV,2200.,2400.,7);
//GaussianTight1MeV->SetParameter(0, 2000.);
//GaussianTight1MeV->SetParLimits(0, 0, 25000);
//GaussianTight1MeV->SetParameter(1, 2286);
//GaussianTight1MeV->SetParameter(2, 5);
//GaussianTight1MeV->SetParLimits(2, 0., 20.);
//GaussianTight1MeV->SetParameter(3, 2000.);
//GaussianTight1MeV->SetParLimits(3, 0, 25000);
//GaussianTight1MeV->SetParameter(4, 5);
//GaussianTight1MeV->SetParLimits(4, 0., 20.);
//GaussianTight1MeV->SetParameter(5, 0.);
//GaussianTight1MeV->SetParameter(6, 0.);

TF1 *GaussianTight1MeV = new TF1("GaussianTight1MeV",fit1MeV_Gaussian,2100.,2500.,5);
GaussianTight1MeV->SetParameter(0,400.);
GaussianTight1MeV->SetParameter(1,2286);
GaussianTight1MeV->SetParameter(2, 6);
GaussianTight1MeV->SetParLimits(2, 0.,20.);
GaussianTight1MeV->SetParameter(3, 0.);
GaussianTight1MeV->SetParameter(4, 0.);

pad1->cd();
MassHistTight1MeV->SetMinimum(0);
MassHistTight1MeV->Fit("GaussianTight1MeV");

int BinHeightTight1MeV[150];
int FitHeightTight1MeV[150];
double PullTight1MeV[150];

double Tight1MeVcount1 = 0;
double Tight1MeVcount2 = 0;
double Tight1MeVcount3 = 0;

for (int bin = 0; bin < 150; bin++){
BinHeightTight1MeV[bin] = MassHistTight1MeV->GetBinContent(bin + 1);
int xvalue1MeV = 2210.5 + 1*(bin);
FitHeightTight1MeV[bin] = round(GaussianTight1MeV->Eval(xvalue1MeV));

PullTight1MeV[bin] = (BinHeightTight1MeV[bin] - FitHeightTight1MeV[bin])/TMath::Sqrt(FitHeightTight1MeV[bin]);

if (PullTight1MeV[bin] > -1 && PullTight1MeV[bin] < 1){
Tight1MeVcount1 += 1;
}

if (PullTight1MeV[bin] > -2 && PullTight1MeV[bin] < 2){
Tight1MeVcount2 += 1;
}

if (PullTight1MeV[bin] > -3 && PullTight1MeV[bin] < 3){
Tight1MeVcount3 += 1;
}
}
     pad2->cd();
     TGraph* PullPlotTight1MeV = new TGraph(150, Pullx1MeV, PullTight1MeV);
     PullPlotTight1MeV->GetXaxis()->SetLimits(0.5,150.5);
     PullPlotTight1MeV->GetXaxis()->SetTickLength(0.);
     PullPlotTight1MeV->GetYaxis()->SetTickLength(0.);
     PullPlotTight1MeV->SetFillColor(38);
     PullPlotTight1MeV->GetYaxis()->SetTitle("Pull");
     PullPlotTight1MeV->GetYaxis()->CenterTitle();
     PullPlotTight1MeV->GetYaxis()->SetTitleSize(0.10);
     PullPlotTight1MeV->GetYaxis()->SetTitleOffset(0.2);
     PullPlotTight1MeV->GetXaxis()->SetLabelSize(0);
     PullPlotTight1MeV->GetYaxis()->SetLabelFont(42);
     PullPlotTight1MeV->GetYaxis()->SetLabelSize(0.06);
     PullPlotTight1MeV->SetTitle("");
     PullPlotTight1MeV->SetMinimum(-5);
     PullPlotTight1MeV->SetMaximum(5);
     PullPlotTight1MeV->Draw("AB");
     c1->Write("Lc Mass - Tight/1MeV");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianTight1MeV->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianTight1MeV->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianTight1MeV->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianTight1MeV->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianTight1MeV->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianTight1MeV->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Tight1MeVcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Tight1MeVcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Tight1MeVcount3));
   ex1->Write("Tight Fit Values/1MeV");
c1->cd();
/////////////////////////////

        File->Close();
}
