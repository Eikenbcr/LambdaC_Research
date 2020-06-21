#define FirstCuts2017_cxx
#include "FirstCuts2017.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"

TH1D * MassHist = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;

void FirstCuts2017::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   MassHist = new TH1D("Mass [MeV]", "LambdaC Mass - 2017 Data", 75, 2210, 2360);
   MassHist->GetXaxis()->SetTitle("MeV");
   MassHist->GetYaxis()->SetTitle("Events Per 2 MeV");

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

   bool SimpleCuts = (
       (*Kminus_IPCHI2_OWNPV > 18)
   &&  (*Kplus_IPCHI2_OWNPV > 18)
   &&  (*Kminus_ProbNNk > 0.7)
   &&  (*Kplus_ProbNNk > 0.7)
   &&  (*Proton_ProbNNp > 0.7)
   );

   if (SimpleCuts){
     MassHist->Fill(*Lcplus_M);
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

MassHist->SetMinimum(0);
MassHist->Fit("Gaussian2MeV");
    c1->Write("Lc Mass");
      File->Close();
}
