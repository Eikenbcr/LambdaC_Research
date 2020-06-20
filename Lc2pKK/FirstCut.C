#define FirstCut_cxx
#include "FirstCut.h"

//ROOT Libraries Needed for Analysis
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

void FirstCut::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   MassHistLoose = new TH1D("Mass [MeV]", "LambdaC Mass - Loose Cut", 75, 2210, 2360);
   MassHistLoose->GetXaxis()->SetTitle("MeV");
   MassHistLoose->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistMid = new TH1D("Mass [MeV]", "LambdaC Mass - Original Cut", 75, 2210, 2360);
   MassHistMid->GetXaxis()->SetTitle("MeV");
   MassHistMid->GetYaxis()->SetTitle("Events Per 2 MeV");

   MassHistTight = new TH1D("Mass [MeV]", "LambdaC Mass - Tight Cut", 75, 2210, 2360);
   MassHistTight->GetXaxis()->SetTitle("MeV");
   MassHistTight->GetYaxis()->SetTitle("Events Per 2 MeV");

   File = new TFile("MassHistogram.root", "RECREATE");
  gFile = File;

   c1 = new TCanvas("canvas", "Test Canvas");

}

void FirstCut::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t FirstCut::Process(Long64_t entry)
{

  GetEntry(entry);
  fReader.SetLocalEntry(entry);

  bool OriginalCuts = (
      (*Kminus_IPCHI2_OWNPV > 15)
  &&  (*Kplus_IPCHI2_OWNPV > 15)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.6)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.6)
  &&  (*Proton_MC15TuneV1_ProbNNp > 0.6)
  );

  bool LooseCuts = (
      (*Kminus_IPCHI2_OWNPV > 10)
  &&  (*Kplus_IPCHI2_OWNPV > 10)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.4)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.4)
  &&  (*Proton_MC15TuneV1_ProbNNp > 0.4)
  );

  bool TightCuts = (
      (*Kminus_IPCHI2_OWNPV > 18)
  &&  (*Kplus_IPCHI2_OWNPV > 18)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.65)
  &&  (*Kminus_MC15TuneV1_ProbNNk > 0.65)
  &&  (*Proton_MC15TuneV1_ProbNNp > 0.65)
  );

  if (OriginalCuts){
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

void FirstCut::SlaveTerminate()
{
}

void FirstCut::Terminate()
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
