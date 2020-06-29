//Initial Cuts with Various Fits and Pull Plots

#define FirstCuts2017_cxx
#include "FirstCuts2017.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"
#include "DGOneMuTwoTotals.C"

TH1D * MassHistLoose = nullptr;
TH1D * MassHistMid = nullptr;
TH1D * MassHistTight = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;
TCanvas * ex1 = nullptr;

void FirstCuts2017::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   MassHistLoose = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistLoose->GetXaxis()->SetTitle("MeV");
   MassHistLoose->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistMid = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistMid->GetXaxis()->SetTitle("MeV");
   MassHistMid->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistTight = new TH1D("Mass [MeV]", "Lc->pKK - Lc Mass", 75, 2210, 2360);
   MassHistTight->GetXaxis()->SetTitle("MeV");
   MassHistTight->GetYaxis()->SetTitle("Events Per 2 MeV");

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
     MassHistMid->Fill(*Lcplus_M);
   }

   if (LooseCuts){
     MassHistLoose->Fill(*Lcplus_M);
   }

   if (TightCuts){
     MassHistTight->Fill(*Lcplus_M);
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
double Pullx[75];
TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
pad2->SetTopMargin(0.03030303);
pad1->Draw();
pad2->Draw();

/////////////////////////////

//TF1 *GaussianLoose = new TF1("GaussianLoose",DGOneMuTwoTotal,2200.,2400.,7);
//GaussianLoose->SetParameter(0, 2000.);
//GaussianLoose->SetParLimits(0, 0, 25000);
//GaussianLoose->SetParameter(1, 2286);
//GaussianLoose->SetParameter(2, 5);
//GaussianLoose->SetParLimits(2, 0., 20.);
//GaussianLoose->SetParameter(3, 2000.);
//GaussianLoose->SetParLimits(3, 0, 25000);
//GaussianLoose->SetParameter(4, 5);
//GaussianLoose->SetParLimits(4, 0., 20.);
//GaussianLoose->SetParameter(5, 0.);
//GaussianLoose->SetParameter(6, 0.);

TF1 *GaussianLoose = new TF1("GaussianLoose",fit2MeV_Gaussian,2100.,2500.,5);
GaussianLoose->SetParameter(0,400.);
GaussianLoose->SetParameter(1,2287);
GaussianLoose->SetParameter(2, 6);
GaussianLoose->SetParLimits(2, 0.,20.);
GaussianLoose->SetParameter(3, 0.);
GaussianLoose->SetParameter(4, 0.);

pad1->cd();
MassHistLoose->SetMinimum(0);
MassHistLoose->Fit("GaussianLoose");

int BinHeightLoose[75];
int FitHeightLoose[75];
double PullLoose[75];

double Loosecount1 = 0;
double Loosecount2 = 0;
double Loosecount3 = 0;

for (int bin = 0; bin < 75; bin++){
BinHeightLoose[bin] = MassHistLoose->GetBinContent(bin + 1);
Pullx[bin] = (bin + 1);
int xvalue = 2211 + 2*(bin);
FitHeightLoose[bin] = round(GaussianLoose->Eval(xvalue));
PullLoose[bin] = (BinHeightLoose[bin] - FitHeightLoose[bin])/TMath::Sqrt(FitHeightLoose[bin]);

if (PullLoose[bin] > -1 && PullLoose[bin] < 1){
  Loosecount1 += 1;
}

if (PullLoose[bin] > -2 && PullLoose[bin] < 2){
  Loosecount2 += 1;
}

if (PullLoose[bin] > -3 && PullLoose[bin] < 3){
  Loosecount3 += 1;
}
}

pad2->cd();
TGraph* PullPlotLoose = new TGraph(75, Pullx, PullLoose);
PullPlotLoose->GetXaxis()->SetLimits(0.5,75.5);
PullPlotLoose->GetXaxis()->SetTickLength(0.);
PullPlotLoose->GetYaxis()->SetTickLength(0.);
PullPlotLoose->SetFillColor(38);
PullPlotLoose->GetYaxis()->SetTitle("Pull");
PullPlotLoose->GetYaxis()->CenterTitle();
PullPlotLoose->GetYaxis()->SetTitleSize(0.10);
PullPlotLoose->GetYaxis()->SetTitleOffset(0.2);
PullPlotLoose->GetXaxis()->SetLabelSize(0);
PullPlotLoose->GetYaxis()->SetLabelFont(42);
PullPlotLoose->GetYaxis()->SetLabelSize(0.06);
PullPlotLoose->SetTitle("");
PullPlotLoose->SetMinimum(-5);
PullPlotLoose->SetMaximum(5);
PullPlotLoose->Draw("AB");
      c1->Write("Lc Mass - Loose");

ex1->cd();
   TLatex Tl;
   Tl.SetTextAlign(12);
   Tl.SetTextSize(0.04);
   Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianLoose->GetParameter(0)));
   Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianLoose->GetParError(0)));
   Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianLoose->GetParameter(1)));
   Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianLoose->GetParError(1)));
   Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianLoose->GetParameter(2)));
   Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianLoose->GetParError(2)));
   Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Loosecount1));
   Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Loosecount2));
   Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Loosecount3));
   ex1->Write("Loose Fit Values");
c1->cd();
///////////////////////////////////////////////////////

//TF1 *GaussianMid = new TF1("GaussianMid",DGOneMuTwoTotal,2200.,2400.,7);
//GaussianMid->SetParameter(0, 2000.);
//GaussianMid->SetParLimits(0, 0, 25000);
//GaussianMid->SetParameter(1, 2286);
//GaussianMid->SetParameter(2, 5);
//GaussianMid->SetParLimits(2, 0., 20.);
//GaussianMid->SetParameter(3, 2000.);
//GaussianMid->SetParLimits(3, 0, 25000);
//GaussianMid->SetParameter(4, 5);
//GaussianMid->SetParLimits(4, 0., 20.);
//GaussianMid->SetParameter(5, 0.);
//GaussianMid->SetParameter(6, 0.);

TF1 *GaussianMid = new TF1("GaussianMid",fit2MeV_Gaussian,2100.,2500.,5);
GaussianMid->SetParameter(0,400.);
GaussianMid->SetParameter(1,2287);
GaussianMid->SetParameter(2, 6);
GaussianMid->SetParLimits(2, 0.,20.);
GaussianMid->SetParameter(3, 0.);
GaussianMid->SetParameter(4, 0.);

  pad1->cd();
  MassHistMid->SetMinimum(0);
  MassHistMid->Fit("GaussianMid");

  int BinHeightMid[75];
  int FitHeightMid[75];
  double PullMid[75];

double   Midcount1 = 0;
double   Midcount2 = 0;
double   Midcount3 = 0;

  for (int bin = 0; bin < 75; bin++){
  BinHeightMid[bin] = MassHistMid->GetBinContent(bin + 1);
  int xvalue = 2211 + 2*(bin);
  FitHeightMid[bin] = round(GaussianMid->Eval(xvalue));
  PullMid[bin] = (BinHeightMid[bin] - FitHeightMid[bin])/TMath::Sqrt(FitHeightMid[bin]);

  if (PullMid[bin] > -1 && PullMid[bin] < 1){
  Midcount1 += 1;
  }

  if (PullMid[bin] > -2 && PullMid[bin] < 2){
  Midcount2 += 1;
  }

  if (PullMid[bin] > -3 && PullMid[bin] < 3){
  Midcount3 += 1;
  }
  }

  pad2->cd();
TGraph* PullPlotMid = new TGraph(75, Pullx, PullMid);
PullPlotMid->GetXaxis()->SetLimits(0.5,75.5);
PullPlotMid->GetXaxis()->SetTickLength(0.);
PullPlotMid->GetYaxis()->SetTickLength(0.);
PullPlotMid->SetFillColor(38);
PullPlotMid->GetYaxis()->SetTitle("Pull");
PullPlotMid->GetYaxis()->CenterTitle();
PullPlotMid->GetYaxis()->SetTitleSize(0.10);
PullPlotMid->GetYaxis()->SetTitleOffset(0.2);
PullPlotMid->GetXaxis()->SetLabelSize(0);
PullPlotMid->GetYaxis()->SetLabelFont(42);
PullPlotMid->GetYaxis()->SetLabelSize(0.06);
PullPlotMid->SetTitle("");
PullPlotMid->SetMinimum(-5);
PullPlotMid->SetMaximum(5);
PullPlotMid->Draw("AB");

      c1->Write("Lc Mass - Mid");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianMid->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianMid->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianMid->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianMid->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianMid->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianMid->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Midcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Midcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Midcount3));
   ex1->Write("Mid Fit Values");
c1->cd();
//////////////////////////////////////////////////////

//TF1 *GaussianTight = new TF1("GaussianTight",DGOneMuTwoTotal,2200.,2400.,7);
//GaussianTight->SetParameter(0, 2000.);
//GaussianTight->SetParLimits(0, 0, 25000);
//GaussianTight->SetParameter(1, 2286);
//GaussianTight->SetParameter(2, 5);
//GaussianTight->SetParLimits(2, 0., 20.);
//GaussianTight->SetParameter(3, 2000.);
//GaussianTight->SetParLimits(3, 0, 25000);
//GaussianTight->SetParameter(4, 5);
//GaussianTight->SetParLimits(4, 0., 20.);
//GaussianTight->SetParameter(5, 0.);
//GaussianTight->SetParameter(6, 0.);

TF1 *GaussianTight = new TF1("GaussianTight",fit2MeV_Gaussian,2100.,2500.,5);
GaussianTight->SetParameter(0,400.);
GaussianTight->SetParameter(1,2286);
GaussianTight->SetParameter(2, 6);
GaussianTight->SetParLimits(2, 0.,20.);
GaussianTight->SetParameter(3, 0.);
GaussianTight->SetParameter(4, 0.);

pad1->cd();
MassHistTight->SetMinimum(0);
MassHistTight->Fit("GaussianTight");

int BinHeightTight[75];
int FitHeightTight[75];
double PullTight[75];

double Tightcount1 = 0;
double Tightcount2 = 0;
double Tightcount3 = 0;

for (int bin = 0; bin < 75; bin++){
BinHeightTight[bin] = MassHistTight->GetBinContent(bin + 1);
int xvalue = 2211 + 2*(bin);
FitHeightTight[bin] = round(GaussianTight->Eval(xvalue));

PullTight[bin] = (BinHeightTight[bin] - FitHeightTight[bin])/TMath::Sqrt(FitHeightTight[bin]);

if (PullTight[bin] > -1 && PullTight[bin] < 1){
Tightcount1 += 1;
}

if (PullTight[bin] > -2 && PullTight[bin] < 2){
Tightcount2 += 1;
}

if (PullTight[bin] > -3 && PullTight[bin] < 3){
Tightcount3 += 1;
}
}
     pad2->cd();
     TGraph* PullPlotTight = new TGraph(75, Pullx, PullTight);
     PullPlotTight->GetXaxis()->SetLimits(0.5,75.5);
     PullPlotTight->GetXaxis()->SetTickLength(0.);
     PullPlotTight->GetYaxis()->SetTickLength(0.);
     PullPlotTight->SetFillColor(38);
     PullPlotTight->GetYaxis()->SetTitle("Pull");
     PullPlotTight->GetYaxis()->CenterTitle();
     PullPlotTight->GetYaxis()->SetTitleSize(0.10);
     PullPlotTight->GetYaxis()->SetTitleOffset(0.2);
     PullPlotTight->GetXaxis()->SetLabelSize(0);
     PullPlotTight->GetYaxis()->SetLabelFont(42);
     PullPlotTight->GetYaxis()->SetLabelSize(0.06);
     PullPlotTight->SetTitle("");
     PullPlotTight->SetMinimum(-5);
     PullPlotTight->SetMaximum(5);
     PullPlotTight->Draw("AB");
     c1->Write("Lc Mass - Tight");

ex1->cd();
   ex1->Clear();
Tl.DrawLatex(0.1,0.9,Form("Number of Signal Entries: %f Events", GaussianTight->GetParameter(0)));
Tl.DrawLatex(0.1,0.85,Form("Error: %f Events", GaussianTight->GetParError(0)));
Tl.DrawLatex(0.1,0.7,Form("Mean Value: %f MeV", GaussianTight->GetParameter(1)));
Tl.DrawLatex(0.1,0.65,Form("Error: %f MeV", GaussianTight->GetParError(1)));
Tl.DrawLatex(0.1,0.5,Form("Sigma: %f MeV", GaussianTight->GetParameter(2)));
Tl.DrawLatex(0.1,0.45,Form("Error: %f MeV", GaussianTight->GetParError(2)));
Tl.DrawLatex(0.1,0.3,Form("Bins Between -1 & 1 %f Bins", Tightcount1));
Tl.DrawLatex(0.1,0.2,Form("Bins Between -2 & 2 %f Bins", Tightcount2));
Tl.DrawLatex(0.1,0.1,Form("Bins Between -3 & 3 %f Bins", Tightcount3));
   ex1->Write("Tight Fit Values");
c1->cd();
/////////////////////////////

        File->Close();
}
