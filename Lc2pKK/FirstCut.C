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

TH1D * MassHistogram = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;

void FirstCut::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   MassHistogram = new TH1D("Mass [MeV]", "Lc_MM - Single Gaussian", 75, 2210, 2360);
   MassHistogram->GetXaxis()->SetTitle("MeV");
   MassHistogram->GetYaxis()->SetTitle("Events Per 2 MeV");

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

  bool SimpleCuts = (
      (*Lc_h1_IPCHI2_OWNPV > 15)
  &&  (*Lc_h1_MC12TuneV4_ProbNNk > 0.6)
  &&  (*Lc_h2_MC12TuneV4_ProbNNk > 0.6)
  &&  (*Lc_p_MC12TuneV4_ProbNNp > 0.6)
  );

  if (SimpleCuts){
    MassHistogram->Fill(*Lc_M);
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

MassHistogram->SetMinimum(0);
MassHistogram->Fit("Gaussian2MeV");

    c1->Write("LambdaC Mass");

      File->Close();
}
