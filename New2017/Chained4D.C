//Cretaing mass histograms using various FOM requirements on Chained ntuples

#define Chained4D_cxx
#include "Chained4D.h"

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

TH2D * FigurePi = nullptr;
TH1D * HistogramPi05 = nullptr;
TH1D * HistogramPi1 = nullptr;
TH1D * HistogramPi15 = nullptr;

TH2D * FigureK = nullptr;
TH1D * HistogramK05 = nullptr;
TH1D * HistogramK1 = nullptr;
TH1D * HistogramK15 = nullptr;

TH1D * FourDimensional05 = nullptr;
TH1D * FourDimensional1 = nullptr;
TH1D * FourDimensional15 = nullptr;

TH1D * FOMKaon = nullptr;
TH1D * FOMPion = nullptr;
TH1D * FOMKPi = nullptr;

TFile * File = nullptr;
TCanvas * c1 = nullptr;

void Chained4D::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void Chained4D::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();


         HistogramK05 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon > 0.5", 100, 2185, 2385);
         HistogramK05->GetXaxis()->SetTitle("MeV");
         HistogramK05->GetYaxis()->SetTitle("Events Per 2 MeV");

         HistogramK1 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon > 1.0", 100, 2185, 2385);
         HistogramK1->GetXaxis()->SetTitle("MeV");
         HistogramK1->GetYaxis()->SetTitle("Events Per 2 MeV");

         HistogramK15 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon > 1.5", 100, 2185, 2385);
         HistogramK15->GetXaxis()->SetTitle("MeV");
         HistogramK15->GetYaxis()->SetTitle("Events Per 2 MeV");

         HistogramPi05 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Pion > 0.5", 100, 2185, 2385);
         HistogramPi05->GetXaxis()->SetTitle("MeV");
         HistogramPi05->GetYaxis()->SetTitle("Events Per 2 MeV");

         HistogramPi1 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Pion > 1.0", 100, 2185, 2385);
         HistogramPi1->GetXaxis()->SetTitle("MeV");
         HistogramPi1->GetYaxis()->SetTitle("Events Per 2 MeV");

         HistogramPi15 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Pion > 1.5", 100, 2185, 2385);
         HistogramPi15->GetXaxis()->SetTitle("MeV");
         HistogramPi15->GetYaxis()->SetTitle("Events Per 2 MeV");

         FourDimensional05 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon & Pion > 0.5", 100, 2185, 2385);
         FourDimensional05->GetXaxis()->SetTitle("MeV");
         FourDimensional05->GetYaxis()->SetTitle("Events Per 2 MeV");

         FourDimensional1 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon & Pion > 1.0", 100, 2185, 2385);
         FourDimensional1->GetXaxis()->SetTitle("MeV");
         FourDimensional1->GetYaxis()->SetTitle("Events Per 2 MeV");

         FourDimensional15 = new TH1D("Figures of Merit", "Mass of LambdaC - FOM Kaon & Pion > 1.5", 100, 2185, 2385);
         FourDimensional15->GetXaxis()->SetTitle("MeV");
         FourDimensional15->GetYaxis()->SetTitle("Events Per 2 MeV");

         FigureK = new TH2D("Figures of Merit", "Color Plot of LambdaC Mass - One Dimensional Projections", 50, 0.7, 4.2, 50, 0.25, 1.0);
         FigureK->GetXaxis()->SetTitle("Log(PromptK_IPCHI2_OWNPV)");
         FigureK->GetYaxis()->SetTitle("PromptK_MC15TuneV1_ProbNNk");
         FigureK->GetZaxis()->SetTitle("Signal*Signal / Background*Background");

         FigurePi = new TH2D("Figures of Merit", "Color Plot of LambdaC Mass - One Dimensional Projections", 50, 0.7, 4.2, 50, 0.25, 1.0);
         FigurePi->GetXaxis()->SetTitle("Log(PromptPi_IPCHI2_OWNPV)");
         FigurePi->GetYaxis()->SetTitle("PromptPi_MC15TuneV1_ProbNNpi");
         FigurePi->GetZaxis()->SetTitle("Signal*Signal / Background*Background");

         FOMKaon = new TH1D("Figures of Merit", "Frequency of Kaon S:B FOMs", 100, 0, 2);
         FOMKaon->GetXaxis()->SetTitle("FOM");
         FOMKaon->GetYaxis()->SetTitle("Events Per 0.2 S:B");

         FOMPion = new TH1D("Figures of Merit", "Frequency of Pion S:B FOMs", 100, 0, 2);
         FOMPion->GetXaxis()->SetTitle("FOM");
         FOMPion->GetYaxis()->SetTitle("Events Per 0.2 S:B");

         FOMKPi = new TH1D("Figures of Merit", "Frequency of Kaon/Pion S:B FOMs", 100, 0, 2);
         FOMKPi->GetXaxis()->SetTitle("FOM");
         FOMKPi->GetYaxis()->SetTitle("Events Per 0.2 S:B");

            File = new TFile("New4DFast.root", "RECREATE");
           gFile = File;

            c1 = new TCanvas("canvas", "Test Canvas");


            for (double i = 0; i < 50; i++){
              for (double j = 0; j < 50; j++){

                      double x = 0.735 + (0.07 * i);
                      double y = 0.2575 + (0.015 * j);
                      double  SignalK = ((6087.15 + (-1115.78 * x) + (-65.9315 * x * x) + TMath::Exp(8.98038 + (3.60867 * x) + (-4.52552 * x * x))) * (9643.84 + (-26988.8 * y) + (22798.9 * y * y) + TMath::Exp(-9.44465 + (-31.0399 * y) + (51.3565 * y * y))));
                      double  BackgroundK = ((2945.82 + (511.272 * x) + (-268.896 * x * x) + TMath::Exp(12.1195 + (-1.34962 * x) + (-1.58848 * x * x))) * (13724.6 + (-35777.2 * y) + (28770.7 * y * y) + TMath::Exp(-8.98010 + (-31.2435 * y) + (50.6150 * y * y))));
                      fomK[i][j] = (SignalK / BackgroundK);}}

           for (double k = 0; k < 50; k++){
             for (double l = 0; l < 50; l++){

                     double w = 0.735 + (0.07 * k);
                     double z = 0.2575 + (0.015 * l);
                     double  SignalPi = ((6152.03 + (-1623.24 * w) + (38.7262 * w * w) + TMath::Exp(10.2967 + (0.801222 * w) + (-2.49779 * w * w))) * (4902.96 + (-15915.0 * z) + (17364.9 * z * z) + TMath::Exp(-5.38716 + (-22.8278 * z) + (38.9474 * z * z))));
                     double  BackgroundPi = ((2911.82 + (236.126 * w) + (-221.903 * w * w) + TMath::Exp(10.4840 + (2.60065 * w) + (-3.68753 * w * w))) * (6700.33 + (-20491.8 * z) + (22428.1 * z * z) + TMath::Exp(-3.26972 + (-23.0203 * z) + (36.6292 * z * z))));
                     fomPi[k][l] = (SignalPi / BackgroundPi);}}

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

Bool_t Chained4D::Process(Long64_t entry)
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

double FOMK = 0;
double FOMPi = 0;
double FOM = 0;

if (m < 50 & n < 50 & p < 50 & q < 50){
FOMK = fomK[p][q];
FOMPi = fomK[m][n];
FOM = fomKPi[m][n][p][q];
}

if (FOMK > 0.5 && AdditionalCuts && BorderCut){
HistogramK05->Fill(CorrectedLambdaMass);}

if (FOMK > 1.0 && AdditionalCuts && BorderCut){
HistogramK1->Fill(CorrectedLambdaMass);}

if (FOMK > 1.5 && AdditionalCuts && BorderCut){
HistogramK15->Fill(CorrectedLambdaMass);}

if (FOMPi > 0.5 && AdditionalCuts && BorderCut){
HistogramPi05->Fill(CorrectedLambdaMass);}

if (FOMPi > 1.0 && AdditionalCuts && BorderCut){
HistogramPi1->Fill(CorrectedLambdaMass);}

if (FOMPi > 1.5 && AdditionalCuts && BorderCut){
HistogramPi15->Fill(CorrectedLambdaMass);}

if (FOM > 0.5 && AdditionalCuts && BorderCut){
FourDimensional05->Fill(CorrectedLambdaMass);}

if (FOM > 1.0 && AdditionalCuts && BorderCut){
FourDimensional1->Fill(CorrectedLambdaMass);}

if (FOM > 1.5 && AdditionalCuts && BorderCut){
FourDimensional15->Fill(CorrectedLambdaMass);}

   return kTRUE;
}

void Chained4D::SlaveTerminate()
{
}

void Chained4D::Terminate()
{

      for (double i = 0; i < 50; i++){
        for (double j = 0; j < 50; j++){

          double x = 0.735 + (0.07 * i);
          double y = 0.2575 + (0.015 * j);

                 FOMKaon->Fill(fomK[i][j]);
                 FigureK->Fill(x, y, fomK[i][j]);}}

                 for (double k = 0; k < 50; k++){
                   for (double l = 0; l < 50; l++){

                     double w = 0.735 + (0.07 * k);
                     double z = 0.2575 + (0.015 * l);

                           FOMPion->Fill(fomPi[k][l]);
                           FigurePi->Fill(w, z, fomPi[k][l]);}}

      for (double m = 0; m < 50; m++){
        for (double n = 0; n < 50; n++){
          for (double p = 0; p < 50; p++){
            for (double q = 0; q < 50; q++){

    FOMKPi->Fill(fomKPi[m][n][p][q]);
    }}}}


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
    TF1 *myLambdaFit = new TF1("myLambdaFit",fit2MeV_Gaussian,2100.,2500.,5);
    myLambdaFit->SetParameter(0,400.);
    myLambdaFit->SetParameter(1,2286);
    myLambdaFit->SetParameter(2, 6);
    myLambdaFit->SetParLimits(2, 0.,20.);
    myLambdaFit->SetParameter(3, 0.);
    myLambdaFit->SetParameter(4, 0.);

    TF1 *myLambdaFit2 = new TF1("myLambdaFit2",fit2MeV_Gaussian,2100.,2500.,5);
    myLambdaFit2->SetParameter(0,400.);
    myLambdaFit2->SetParameter(1,2286);
    myLambdaFit2->SetParameter(2, 6);
    myLambdaFit2->SetParLimits(2, 0.,20.);
    myLambdaFit2->SetParameter(3, 0.);
    myLambdaFit2->SetParameter(4, 0.);

    HistogramK05->Fit("myLambdaFit");
    HistogramK05->SetMinimum(0);
    c1->Write("Kaon 0.5");

    HistogramK1->Fit("myLambdaFit");
    HistogramK1->SetMinimum(0);
    c1->Write("Kaon 1.0");

    HistogramK15->Fit("myLambdaFit");
    HistogramK15->SetMinimum(0);
    c1->Write("Kaon 1.5");

    HistogramPi05->Fit("myLambdaFit");
    HistogramPi05->SetMinimum(0);
    c1->Write("Pion 0.5");

    HistogramPi1->Fit("myLambdaFit2");
    HistogramPi1->SetMinimum(0);
    c1->Write("Pion 1.0");

    HistogramPi15->Fit("myLambdaFit2");
    HistogramPi15->SetMinimum(0);
    c1->Write("Pion 1.5");

    FourDimensional05->Fit("myLambdaFit2");
    FourDimensional05->SetMinimum(0);
    c1->Write("FourDimensional 0.5");

    FourDimensional1->Fit("myLambdaFit2");
    FourDimensional1->SetMinimum(0);
    c1->Write("FourDimensional 1.0");

    FourDimensional15->Fit("myLambdaFit2");
    FourDimensional15->SetMinimum(0);
    c1->Write("FourDimensional 1.5");

    FigureK->Draw("COLZ");
     c1->Write("Kaon Color Plot");

     FigurePi->Draw("COLZ");
      c1->Write("Pion Color Plot");

      FOMKaon->Draw();
      c1->Write("Kaon FOMs");

      FOMPion->Draw();
      c1->Write("Pion FOMs");

    FOMKPi->Draw();
    c1->Write("Kaon/Pion FOMs");

    File->Close();

}
