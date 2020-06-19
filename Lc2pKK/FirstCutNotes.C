#define FirstCut_cxx
#include "FirstCut.h"

//ROOT Libraries Needed for Analysis
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

//Includes Gaussian Fit File Needed for Analysis
#include "fit2MeV_Gaussian.C"

//Creating Histogram used for Mass Analysis
TH1D * MassHistogram = nullptr;

//Creates Output File for Histogram (pdf/png/root)
TFile * File = nullptr;

//Creates Canvas which the Histogram will be created on
TCanvas * c1 = nullptr;

void FirstCut::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

//Defines Parameters For Histogram such as title, units, range, axis labels
   MassHistogram = new TH1D("Mass [MeV]", "Lc_MM - Single Gaussian", 75, 2210, 2360);
   MassHistogram->GetXaxis()->SetTitle("MeV");
   MassHistogram->GetYaxis()->SetTitle("Events Per 2 MeV");

//Output Defined as "MassHistogram.root" which is a root file
   File = new TFile("MassHistogram.root", "RECREATE");
  gFile = File;

//Names the Canvas
   c1 = new TCanvas("canvas", "Test Canvas");

}

void FirstCut::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t FirstCut::Process(Long64_t entry)
{
//These Two Lines Indicate that the code will loop over all Entries/Events
  GetEntry(entry);
  fReader.SetLocalEntry(entry);

/////////
//Notes//
/////////

//Decay = (Charmed Lambda Baryon -> Proton + Kaon + AntiKaon)
//Ntuple Definitions For Particles = Charmed Lambda(Lc), Proton(Lc_p), Kaon(Lc_h2), AntiKaon(Lc_h1)

//SimpleCuts is a list of requirements that we can set on the data to increase significance of statistics
  bool SimpleCuts = (

      //Anti Kaon's (K-) Impact Parameter distance from Primary Vertex
      (*Lc_h1_IPCHI2_OWNPV > 15)

      //Kaon (K+) Impact Parameter distance from Primary Vertex
  &&  (*Lc_h2_IPCHI2_OWNPV > 15)

      //Probability that charged track is identified as a Kaon
  &&  (*Lc_h1_MC12TuneV4_ProbNNk > 0.6)

      //Probability that charged track is identified as a Kaon
  &&  (*Lc_h2_MC12TuneV4_ProbNNk > 0.6)

      //Probability that charged track is identified as a Proton
  &&  (*Lc_p_MC12TuneV4_ProbNNp > 0.6)
  );

//If data passes SimpleCuts, they will be filled into the Histogram of LambdaC Mass
  if (SimpleCuts){
    MassHistogram->Fill(*Lc_M);
  }

//Required to End Event Loop
   return kTRUE;
}

void FirstCut::SlaveTerminate()
{
}

void FirstCut::Terminate()
{
  //Necessary For Initializing Gaussian Fit (Always Include)
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

//Creates Function from Gaussian Fit
  TF1 *Gaussian2MeV = new TF1("Gaussian2MeV",fit2MeV_Gaussian,2100.,2500.,5);

//We can set parameters to certain values to improve function fit (Does not need to be exact)

//Sets Events in Signal to be 400
  Gaussian2MeV->SetParameter(0,400.);

//Sets Mean Value to be 2286 (Approximately the Mass of the LambdaC)
  Gaussian2MeV->SetParameter(1,2286);

//Sets Sigma to 6
  Gaussian2MeV->SetParameter(2, 6);

//Requires Sigma to be positive (Or between 0 and 20)
  Gaussian2MeV->SetParLimits(2, 0.,20.);

//Sets slope of linear background to 0
  Gaussian2MeV->SetParameter(3, 0.);

//Sets y intercept of linear background to 0
  Gaussian2MeV->SetParameter(4, 0.);

//Fits Data with Gaussian Fit
MassHistogram->Fit("Gaussian2MeV");

//Names the Canvas for Root File Output
    c1->Write("LambdaC Mass");

//Creates Output
      File->Close();
}
