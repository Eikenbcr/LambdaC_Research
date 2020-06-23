#define FirstCuts2017_cxx
#include "FirstCuts2017.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"

TH1D * MassHistLoose = nullptr;
TH1D * MassHistMid = nullptr;
TH1D * MassHistTight = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;

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
}

void FirstCuts2017::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t FirstCuts2017::Process(Long64_t entry)
{

    GetEntry(entry);
   fReader.SetLocalEntry(entry);

   bool MidCuts = (
       (*Kminus_ProbNNk > 0.75)
   &&  (*Kplus_ProbNNk > 0.75)
   &&  (*Proton_ProbNNp > 0.75)
   &&  (*Kminus_PT > 500)
   &&  (*Kplus_PT > 500)
   );

   bool LooseCuts = (
       (*Kminus_ProbNNk > 0.75)
   &&  (*Kplus_ProbNNk > 0.75)
   &&  (*Proton_ProbNNp > 0.75)
   );

   bool TightCuts = (
       (*Kminus_ProbNNk > 0.85)
   &&  (*Kplus_ProbNNk > 0.85)
   &&  (*Proton_ProbNNp > 0.9)
   &&  (*Kminus_PT > 1000)
   &&  (*Kplus_PT > 1000)
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

  TF1 *Gaussian2MeV = new TF1("Gaussian2MeV",fit2MeV_Gaussian,2100.,2500.,5);
  Gaussian2MeV->SetParameter(0,400.);
  Gaussian2MeV->SetParameter(1,2286);
  Gaussian2MeV->SetParameter(2, 6);
  Gaussian2MeV->SetParLimits(2, 0.,20.);
  Gaussian2MeV->SetParameter(3, 0.);
  Gaussian2MeV->SetParameter(4, 0.);

  MassHistLoose->SetMinimum(0);
  MassHistLoose->Fit("Gaussian2MeV");

      c1->Write("Lc Mass - Loose");

  MassHistMid->SetMinimum(0);
  MassHistMid->Fit("Gaussian2MeV");

      c1->Write("Lc Mass - Original");

  MassHistTight->SetMinimum(0);
  MassHistTight->Fit("Gaussian2MeV");

     c1->Write("Lc Mass - Tight");

        File->Close();
}
