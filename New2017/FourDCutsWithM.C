//Same as New4DAnalysis.C but with Lc_M instead of Lc_MM

#define FourDCutsWithM_cxx
#include "FourDCutsWithM.h"


#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"

using namespace std;

#include <vector>
#include <algorithm>
#include <functional>
#include <chrono>
#include <thread>
#include <array>

std::array<std::array<double, 50>, 50> fomK;
std::array<std::array<double, 50>, 50> fomPi;
std::array<std::array<std::array<std::array<double, 50>, 50>, 50>, 50> fomKPi;

TH1D * FourDimensional05 = nullptr;
TH1D * FourDimensional1 = nullptr;
TH1D * FourDimensional15 = nullptr;

TFile * File = nullptr;
TCanvas * c1 = nullptr;


void FourDCutsWithM::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void FourDCutsWithM::SlaveBegin(TTree * /*tree*/)
{

   TString option = GetOption();

         FourDimensional05 = new TH1D("Figures of Merit", "Lc->XiKpi - Lc Mass", 100, 2185, 2385);
         FourDimensional05->GetXaxis()->SetTitle("MeV");
         FourDimensional05->GetYaxis()->SetTitle("Events Per 2 MeV");

         FourDimensional1 = new TH1D("Figures of Merit", "Lc->XiKpi - Lc Mass", 100, 2185, 2385);
         FourDimensional1->GetXaxis()->SetTitle("MeV");
         FourDimensional1->GetYaxis()->SetTitle("Events Per 2 MeV");

         FourDimensional15 = new TH1D("Figures of Merit", "Lc->XiKpi - Lc Mass", 100, 2185, 2385);
         FourDimensional15->GetXaxis()->SetTitle("MeV");
         FourDimensional15->GetYaxis()->SetTitle("Events Per 2 MeV");

         File = new TFile("LcToXiKpiBestHistograms.root", "RECREATE");
        gFile = File;

         c1 = new TCanvas("canvas", "Test Canvas");

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
}

Bool_t FourDCutsWithM::Process(Long64_t entry)
{

  GetEntry(entry);
   fReader.SetLocalEntry(entry);

                     //    Corrected Xi Mass
                         double CorrectedXiMass = ((*Xi_M) - (*Lambda_M) + 1115.683);

                     //    Corrected Lambda Mass
                         double CorrectedLambdaMass = ((*Lc_M) - (*Xi_M) + (1321.71));

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

               if (FOM > 0.5 && AdditionalCuts && BorderCut){
               FourDimensional05->Fill(CorrectedLambdaMass);}

               if (FOM > 1.0 && AdditionalCuts && BorderCut){
               FourDimensional1->Fill(CorrectedLambdaMass);}

               if (FOM > 1.5 && AdditionalCuts && BorderCut){
               FourDimensional15->Fill(CorrectedLambdaMass);}

   return kTRUE;
}

void FourDCutsWithM::SlaveTerminate()
{
}

void FourDCutsWithM::Terminate()
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

    TF1 *myLambdaFit2 = new TF1("myLambdaFit2",fit2MeV_Gaussian,2100.,2500.,5);
    myLambdaFit2->SetParameter(0,400.);
    myLambdaFit2->SetParameter(1,2286);
    myLambdaFit2->SetParameter(2, 6);
    myLambdaFit2->SetParLimits(2, 0.,20.);
    myLambdaFit2->SetParameter(3, 0.);
    myLambdaFit2->SetParameter(4, 0.);

    FourDimensional05->Fit("myLambdaFit2");
    FourDimensional05->SetMinimum(0);
    c1->Write("FourDimensional 0.5");

    FourDimensional1->Fit("myLambdaFit2");
    FourDimensional1->SetMinimum(0);
    c1->Write("FourDimensional 1.0");

    FourDimensional15->Fit("myLambdaFit2");
    FourDimensional15->SetMinimum(0);
    c1->Write("FourDimensional 1.5");

  File->Close();
}
