#define ChainedStatistics_cxx
#include "ChainedStatistics.h"


#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include "fit2MeV_Gaussian.C"
#include "DGOneMuOneTotal.C"

using namespace std;

#include <vector>
#include <algorithm>
#include <functional>
#include <chrono>
#include <thread>
#include <array>

std::array<std::array<std::array<std::array<double, 50>, 50>, 50>, 50> fomKPi;

TH1D * OptimalCut = nullptr;

TH1D * LcLowP = nullptr;
TH1D * LcMidLowP = nullptr;
TH1D * LcMidHighP = nullptr;
TH1D * LcHighP = nullptr;
TH1D * LcPDistribution = nullptr;

TH1D * XiLowP = nullptr;
TH1D * XiMidLowP = nullptr;
TH1D * XiMidHighP = nullptr;
TH1D * XiHighP = nullptr;
TH1D * XiPDistribution = nullptr;

TH1D * PolarityMagUp = nullptr;
TH1D * PolarityMagDown = nullptr;

TH1D * Particle = nullptr;
TH1D * AntiParticle = nullptr;

TGraphErrors * gLcP = nullptr;
TGraphErrors * gXiP = nullptr;
TGraphErrors * gPolarity = nullptr;
TGraphErrors * gPID = nullptr;

TGraphErrors * gLcPDG = nullptr;
TGraphErrors * gXiPDG = nullptr;
TGraphErrors * gPolarityDG = nullptr;
TGraphErrors * gPIDDG = nullptr;

TFile * File = nullptr;
TCanvas * c1 = nullptr;

void ChainedStatistics::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void ChainedStatistics::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

            OptimalCut = new TH1D("Figures of Merit", "Optimal Mass Cut on LambdaC Baryon", 100, 2185, 2385);
            OptimalCut->GetXaxis()->SetTitle("MeV");
            OptimalCut->GetYaxis()->SetTitle("Events Per 2 MeV");

            LcLowP = new TH1D("Figures of Merit", "LambdaC Mass - Low LcP", 100, 2185, 2385);
            LcLowP->GetXaxis()->SetTitle("MeV");
            LcLowP->GetYaxis()->SetTitle("Events Per 2 MeV");

            LcMidLowP = new TH1D("Figures of Merit", "LambdaC Mass - MidLow LcP", 100, 2185, 2385);
            LcMidLowP->GetXaxis()->SetTitle("MeV");
            LcMidLowP->GetYaxis()->SetTitle("Events Per 2 MeV");

            LcMidHighP = new TH1D("Figures of Merit", "LambdaC Mass - MidHigh LcP", 100, 2185, 2385);
            LcMidHighP->GetXaxis()->SetTitle("MeV");
            LcMidHighP->GetYaxis()->SetTitle("Events Per 2 MeV");

            LcHighP = new TH1D("Figures of Merit", "LambdaC Mass - High LcP", 100, 2185, 2385);
            LcHighP->GetXaxis()->SetTitle("MeV");
            LcHighP->GetYaxis()->SetTitle("Events Per 2 MeV");

            LcPDistribution = new TH1D("Figures of Merit", "LcP Distribution After Cuts", 100, 0, 250000);
            LcPDistribution->GetXaxis()->SetTitle("MeV");
            LcPDistribution->GetYaxis()->SetTitle("Events Per 150 MeV");

            XiLowP = new TH1D("Figures of Merit", "LambdaC Mass - Low XiP", 100, 2185, 2385);
            XiLowP->GetXaxis()->SetTitle("MeV");
            XiLowP->GetYaxis()->SetTitle("Events Per 2 MeV");

            XiMidLowP = new TH1D("Figures of Merit", "LambdaC Mass - MidLow XiP", 100, 2185, 2385);
            XiMidLowP->GetXaxis()->SetTitle("MeV");
            XiMidLowP->GetYaxis()->SetTitle("Events Per 2 MeV");

            XiMidHighP = new TH1D("Figures of Merit", "LambdaC Mass - MidHigh XiP", 100, 2185, 2385);
            XiMidHighP->GetXaxis()->SetTitle("MeV");
            XiMidHighP->GetYaxis()->SetTitle("Events Per 2 MeV");

            XiHighP = new TH1D("Figures of Merit", "LambdaC Mass - High XiP", 100, 2185, 2385);
            XiHighP->GetXaxis()->SetTitle("MeV");
            XiHighP->GetYaxis()->SetTitle("Events Per 2 MeV");

            XiPDistribution = new TH1D("Figures of Merit", "XiP Distribution After Cuts", 100, 0, 250000);
            XiPDistribution->GetXaxis()->SetTitle("MeV");
            XiPDistribution->GetYaxis()->SetTitle("Events Per 150 MeV");

            PolarityMagDown = new TH1D("Figures of Merit", "LambdaC Mass - MagnetDown", 100, 2185, 2385);
            PolarityMagDown->GetXaxis()->SetTitle("MeV");
            PolarityMagDown->GetYaxis()->SetTitle("Events Per 2 MeV");

            PolarityMagUp = new TH1D("Figures of Merit", "LambdaC Mass - MagnetUp", 100, 2185, 2385);
            PolarityMagUp->GetXaxis()->SetTitle("MeV");
            PolarityMagUp->GetYaxis()->SetTitle("Events Per 2 MeV");

            Particle = new TH1D("Figures of Merit", "LambdaC Mass - LambdaC Baryon", 100, 2185, 2385);
            Particle->GetXaxis()->SetTitle("MeV");
            Particle->GetYaxis()->SetTitle("Events Per 2 MeV");

            AntiParticle = new TH1D("Figures of Merit", "LambdaC Mass - LambdaC AntiBaryon", 100, 2185, 2385);
            AntiParticle->GetXaxis()->SetTitle("MeV");
            AntiParticle->GetYaxis()->SetTitle("Events Per 2 MeV");

      File = new TFile("ChainStats.root", "RECREATE");
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

Bool_t ChainedStatistics::Process(Long64_t entry)
{
    GetEntry(entry);
   fReader.SetLocalEntry(entry);

   //    Corrected Xi Mass
   double CorrectedXiMass = ((*Xi_MM) - (*Lambda_MM) + 1115.683);

   //    Corrected Lambda Mass
   double CorrectedLambdaMass = ((*Lc_MM) - (*Xi_MM) + (1321.71));

   bool LcLowPCut = (
     (*Lc_P < 50200.)
   );

   bool LcMidLowPCut = (
     (*Lc_P > 50200. && *Lc_P < 65700.)
   );

   bool LcMidHighPCut = (
     (*Lc_P > 65700. && *Lc_P < 87900.)
   );

   bool LcHighPCut = (
     (*Lc_P > 87900.)
   );

   bool XiLowPCut = (
     (*Xi_P < 27150.)
   );

   bool XiMidLowPCut = (
      (*Xi_P > 27150. && *Xi_P < 36550.)
   );

   bool XiMidHighPCut = (
      (*Xi_P > 36550. && *Xi_P < 50100.)
   );

   bool XiHighPCut = (
      (*Xi_P > 50100.)
   );

   bool MagDown = (
      (*Polarity < 0.)
   );

   bool MagUp = (
      (*Polarity > 0.)
   );

   bool ParticleCut = (
      (*Lc_ID > 0.)
   );

   bool AntiParticleCut = (
      (*Lc_ID < 0.)
   );


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
   OptimalCut->Fill(CorrectedLambdaMass);
   LcPDistribution->Fill(*Lc_P);
   XiPDistribution->Fill(*Xi_P);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && LcLowPCut){
   LcLowP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut &&LcMidLowPCut){
   LcMidLowP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && LcMidHighPCut){
   LcMidHighP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && LcHighPCut){
   LcHighP->Fill(CorrectedLambdaMass);}

     //////////////

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && XiLowPCut){
   XiLowP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && XiMidLowPCut){
   XiMidLowP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && XiMidHighPCut){
   XiMidHighP->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && XiHighPCut){
   XiHighP->Fill(CorrectedLambdaMass);}

     //////////////

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && MagDown){
   PolarityMagDown->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && MagUp){
   PolarityMagUp->Fill(CorrectedLambdaMass);}

     /////////////

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && ParticleCut){
   Particle->Fill(CorrectedLambdaMass);}

                 if (FOM > 1.0 && AdditionalCuts && BorderCut && AntiParticleCut){
   AntiParticle->Fill(CorrectedLambdaMass);}


   return kTRUE;
}

void ChainedStatistics::SlaveTerminate()
{

}

void ChainedStatistics::Terminate()
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

      TF1 *myLambdaFitSG = new TF1("myLambdaFitSG",fit2MeV_Gaussian,2100.,2500.,5);
      myLambdaFitSG->SetParameter(0,400.);
      myLambdaFitSG->SetParameter(1,2286.5);
      myLambdaFitSG->SetParameter(2, 6);
      myLambdaFitSG->SetParLimits(2, 0.,20.);
      myLambdaFitSG->SetParameter(3, 0.);
      myLambdaFitSG->SetParameter(4, 0.);

      TF1 *LcFitSG1 = new TF1("LcFitSG1",fit2MeV_Gaussian,2100.,2500.,5);
      LcFitSG1->SetParameter(0,400.);
      LcFitSG1->SetParameter(1,2286.5);
      LcFitSG1->SetParameter(2, 6);
      LcFitSG1->SetParLimits(2, 0.,20.);
      LcFitSG1->SetParameter(3, 0.);
      LcFitSG1->SetParameter(4, 0.);

      TF1 *LcFitSG2 = new TF1("LcFitSG2",fit2MeV_Gaussian,2100.,2500.,5);
      LcFitSG2->SetParameter(0,400.);
      LcFitSG2->SetParameter(1,2286.5);
      LcFitSG2->SetParameter(2, 6);
      LcFitSG2->SetParLimits(2, 0.,20.);
      LcFitSG2->SetParameter(3, 0.);
      LcFitSG2->SetParameter(4, 0.);

      TF1 *LcFitSG3 = new TF1("LcFitSG3",fit2MeV_Gaussian,2100.,2500.,5);
      LcFitSG3->SetParameter(0,400.);
      LcFitSG3->SetParameter(1,2286.5);
      LcFitSG3->SetParameter(2, 6);
      LcFitSG3->SetParLimits(2, 0.,20.);
      LcFitSG3->SetParameter(3, 0.);
      LcFitSG3->SetParameter(4, 0.);

      TF1 *LcFitSG4 = new TF1("LcFitSG4",fit2MeV_Gaussian,2100.,2500.,5);
      LcFitSG4->SetParameter(0,400.);
      LcFitSG4->SetParameter(1,2286.5);
      LcFitSG4->SetParameter(2, 6);
      LcFitSG4->SetParLimits(2, 0.,20.);
      LcFitSG4->SetParameter(3, 0.);
      LcFitSG4->SetParameter(4, 0.);

      TF1 *XiFitSG1 = new TF1("XiFitSG1",fit2MeV_Gaussian,2100.,2500.,5);
      XiFitSG1->SetParameter(0,400.);
      XiFitSG1->SetParameter(1,2286.5);
      XiFitSG1->SetParameter(2, 6);
      XiFitSG1->SetParLimits(2, 0.,20.);
      XiFitSG1->SetParameter(3, 0.);
      XiFitSG1->SetParameter(4, 0.);

      TF1 *XiFitSG2 = new TF1("XiFitSG2",fit2MeV_Gaussian,2100.,2500.,5);
      XiFitSG2->SetParameter(0,400.);
      XiFitSG2->SetParameter(1,2286.5);
      XiFitSG2->SetParameter(2, 6);
      XiFitSG2->SetParLimits(2, 0.,20.);
      XiFitSG2->SetParameter(3, 0.);
      XiFitSG2->SetParameter(4, 0.);

      TF1 *XiFitSG3 = new TF1("XiFitSG3",fit2MeV_Gaussian,2100.,2500.,5);
      XiFitSG3->SetParameter(0,400.);
      XiFitSG3->SetParameter(1,2286.5);
      XiFitSG3->SetParameter(2, 6);
      XiFitSG3->SetParLimits(2, 0.,20.);
      XiFitSG3->SetParameter(3, 0.);
      XiFitSG3->SetParameter(4, 0.);

      TF1 *XiFitSG4 = new TF1("XiFitSG4",fit2MeV_Gaussian,2100.,2500.,5);
      XiFitSG4->SetParameter(0,400.);
      XiFitSG4->SetParameter(1,2286.5);
      XiFitSG4->SetParameter(2, 6);
      XiFitSG4->SetParLimits(2, 0.,20.);
      XiFitSG4->SetParameter(3, 0.);
      XiFitSG4->SetParameter(4, 0.);

      TF1 *MagDownFitSG = new TF1("MagDownFitSG",fit2MeV_Gaussian,2100.,2500.,5);
      MagDownFitSG->SetParameter(0,400.);
      MagDownFitSG->SetParameter(1,2286.5);
      MagDownFitSG->SetParameter(2, 6);
      MagDownFitSG->SetParLimits(2, 0.,20.);
      MagDownFitSG->SetParameter(3, 0.);
      MagDownFitSG->SetParameter(4, 0.);

      TF1 *MagUpFitSG = new TF1("MagUpFitSG",fit2MeV_Gaussian,2100.,2500.,5);
      MagUpFitSG->SetParameter(0,400.);
      MagUpFitSG->SetParameter(1,2286.5);
      MagUpFitSG->SetParameter(2, 6);
      MagUpFitSG->SetParLimits(2, 0.,20.);
      MagUpFitSG->SetParameter(3, 0.);
      MagUpFitSG->SetParameter(4, 0.);

      TF1 *ParticleFitSG = new TF1("ParticleFitSG",fit2MeV_Gaussian,2100.,2500.,5);
      ParticleFitSG->SetParameter(0,400.);
      ParticleFitSG->SetParameter(1,2286.5);
      ParticleFitSG->SetParameter(2, 6);
      ParticleFitSG->SetParLimits(2, 0.,20.);
      ParticleFitSG->SetParameter(3, 0.);
      ParticleFitSG->SetParameter(4, 0.);

      TF1 *AntiParticleFitSG = new TF1("AntiParticleFitSG",fit2MeV_Gaussian,2100.,2500.,5);
      AntiParticleFitSG->SetParameter(0,400.);
      AntiParticleFitSG->SetParameter(1,2286.5);
      AntiParticleFitSG->SetParameter(2, 6);
      AntiParticleFitSG->SetParLimits(2, 0.,20.);
      AntiParticleFitSG->SetParameter(3, 0.);
      AntiParticleFitSG->SetParameter(4, 0.);

      OptimalCut->Fit("myLambdaFitSG");
      OptimalCut->SetMinimum(0);
      c1->Write("OptimalCut");

      LcLowP->Fit("LcFitSG1");
      LcLowP->SetMinimum(0);
      c1->Write("LcPLow");

      LcMidLowP->Fit("LcFitSG2");
      LcMidLowP->SetMinimum(0);
      c1->Write("LcPMidLow");

      LcMidHighP->Fit("LcFitSG3");
      LcMidHighP->SetMinimum(0);
      c1->Write("LcPMidHigh");

      LcHighP->Fit("LcFitSG4");
      LcHighP->SetMinimum(0);
      c1->Write("LcPHigh");

      LcPDistribution->Draw();
      LcPDistribution->SetMinimum(0);
      c1->Write("LcPDistribution");

      double y1 = LcFitSG1->GetParameter(1);
      double y2 = LcFitSG2->GetParameter(1);
      double y3 = LcFitSG3->GetParameter(1);
      double y4 = LcFitSG4->GetParameter(1);

      double yerr1 = LcFitSG1->GetParError(1);
      double yerr2 = LcFitSG2->GetParError(1);
      double yerr3 = LcFitSG3->GetParError(1);
      double yerr4 = LcFitSG4->GetParError(1);

      const Int_t n = 4;
      Double_t x[n] = {1,2,3,4};
      Double_t xerr[n] = {0,0,0,0};

      double LcPAvG = ((y1 + y2 + y3 + y4)/n);
      Double_t y[n] = {y1 - LcPAvG,y2 - LcPAvG,y3 - LcPAvG,y4 - LcPAvG};
      Double_t yerr[n] = {yerr1,yerr2,yerr3,yerr4};
  //    Double_t chi2LcP = (((y1 - LcPAvG)*(y1 - LcPAvG))+((y2 - LcPAvG)*(y2 - LcPAvG))+((y3 - LcPAvG)*(y3 - LcPAvG))+((y4 - LcPAvG)*(y4 - LcPAvG)))/LcPAvG;
      gLcP = new TGraphErrors(n,x,y,xerr,yerr);
      gLcP->SetMarkerColor(2);
      gLcP->SetMarkerStyle(20);
      gLcP->GetXaxis()->SetNdivisions(4);
      gLcP->GetXaxis()->SetTitle("LambdaC_P Regions");
      gLcP->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
      gLcP->SetTitle("Deviations in LambdaC_MM of Different LambdaC_P Regions");
      gLcP->Draw("ALP");
      c1->Write("LcPDivision");

      XiLowP->Fit("XiFitSG1");
      XiLowP->SetMinimum(0);
      c1->Write("XiPLow");

      XiMidLowP->Fit("XiFitSG2");
      XiMidLowP->SetMinimum(0);
      c1->Write("XiPMidLow");

      XiMidHighP->Fit("XiFitSG3");
      XiMidHighP->SetMinimum(0);
      c1->Write("XiPMidHigh");

      XiHighP->Fit("XiFitSG4");
      XiHighP->SetMinimum(0);
      c1->Write("XiPHigh");

      XiPDistribution->Draw();
      XiPDistribution->SetMinimum(0);
      c1->Write("XiPDistribution");

    double yy1 = XiFitSG1->GetParameter(1);
    double yy2 = XiFitSG2->GetParameter(1);
    double yy3 = XiFitSG3->GetParameter(1);
    double yy4 = XiFitSG4->GetParameter(1);

    double yyerr1 = XiFitSG1->GetParError(1);
    double yyerr2 = XiFitSG2->GetParError(1);
    double yyerr3 = XiFitSG3->GetParError(1);
    double yyerr4 = XiFitSG4->GetParError(1);

    double XiPAvG = ((yy1 + yy2 + yy3 + yy4)/n);
    Double_t yy[n] = {yy1 - XiPAvG,yy2 - XiPAvG,yy3 - XiPAvG,yy4 - XiPAvG};
    Double_t yyerr[n] = {yyerr1,yyerr2,yyerr3,yyerr4};

    gXiP = new TGraphErrors(n,x,yy,xerr,yyerr);
    gXiP->SetMarkerColor(2);
    gXiP->SetMarkerStyle(20);
    gXiP->GetXaxis()->SetNdivisions(4);
    gXiP->GetXaxis()->SetTitle("Xi_P Regions");
    gXiP->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gXiP->SetTitle("Deviations in LambdaC_MM of Different Xi_P Regions");
    gXiP->Draw("ALP");
    c1->Write("XiPDivision");

    PolarityMagDown->Fit("MagDownFitSG");
    PolarityMagDown->SetMinimum(0);
    c1->Write("PolarityMagDown");

    PolarityMagUp->Fit("MagUpFitSG");
    PolarityMagUp->SetMinimum(0);
    c1->Write("PolarityMagUp");

    double p1 = MagDownFitSG->GetParameter(1);
    double p2 = MagUpFitSG->GetParameter(1);

    double perr1 = MagDownFitSG->GetParError(1);
    double perr2 = MagUpFitSG->GetParError(1);

      const Int_t m = 2;

    double PolarityAvG = ((p1 + p2)/m);
    Double_t p[m] = {p1 - PolarityAvG,p2 - PolarityAvG};
    Double_t perr[m] = {perr1,perr2};

    gPolarity = new TGraphErrors(m,x,p,xerr,perr);
    gPolarity->SetMarkerColor(2);
    gPolarity->SetMarkerStyle(20);
    gPolarity->GetXaxis()->SetNdivisions(2);
    gPolarity->GetXaxis()->SetTitle("Magnet Orientations");
    gPolarity->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gPolarity->SetTitle("Deviations in LambdaC_MM of Different Magnet Orientations");
    gPolarity->Draw("ALP");
    c1->Write("PolarityDivision");

    Particle->Fit("ParticleFitSG");
    Particle->SetMinimum(0);
    c1->Write("Baryon");

    AntiParticle->Fit("AntiParticleFitSG");
    AntiParticle->SetMinimum(0);
    c1->Write("AntiBaryon");

    double pp1 = ParticleFitSG->GetParameter(1);
    double pp2 = AntiParticleFitSG->GetParameter(1);

    double pperr1 = ParticleFitSG->GetParError(1);
    double pperr2 = AntiParticleFitSG->GetParError(1);

    double ParticleAvG = ((pp1 + pp2)/m);
    Double_t pp[m] = {pp1 - ParticleAvG,pp2 - ParticleAvG};
    Double_t pperr[m] = {pperr1,pperr2};

    gPID = new TGraphErrors(m,x,pp,xerr,pperr);
    gPID->SetMarkerColor(2);
    gPID->SetMarkerStyle(20);
    gPID->GetXaxis()->SetNdivisions(2);
    gPID->GetXaxis()->SetTitle("Particle ID");
    gPID->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gPID->SetTitle("Deviations in LambdaC_MM of Baryon or AntiBaryon");
    gPID->Draw("ALP");
    c1->Write("ParticleDivision");

////////////////////////////////////////////////////////////////////////


      TF1 *myLambdaFitDG = new TF1("myLambdaFitDG",DGOneMuOneTotal,2100.,2500.,7);
      myLambdaFitDG->SetParameter(0, 0.2);
      myLambdaFitDG->SetParameter(1, 20000);
      myLambdaFitDG->SetParameter(2, 2287.);
      myLambdaFitDG->SetParameter(3, 7);
      myLambdaFitDG->SetParameter(4, 4);
      myLambdaFitDG->SetParLimits(3, 0., 20.);
      myLambdaFitDG->SetParLimits(4, 0., 20.);
      myLambdaFitDG->SetParameter(5, 0.);
      myLambdaFitDG->SetParameter(6, 0.);

      TF1 *LcFitDG1 = new TF1("LcFitDG1",DGOneMuOneTotal,2100.,2500.,7);
      LcFitDG1->SetParameter(0, 0.2);
      LcFitDG1->SetParameter(1, 20000);
      LcFitDG1->SetParameter(2, 2287.);
      LcFitDG1->SetParameter(3, 7);
      LcFitDG1->SetParameter(4, 4);
      LcFitDG1->SetParLimits(3, 0., 20.);
      LcFitDG1->SetParLimits(4, 0., 20.);
      LcFitDG1->SetParameter(5, 0.);
      LcFitDG1->SetParameter(6, 0.);

      TF1 *LcFitDG2 = new TF1("LcFitDG2",DGOneMuOneTotal,2100.,2500.,7);
      LcFitDG2->SetParameter(0, 0.2);
      LcFitDG2->SetParameter(1, 20000);
      LcFitDG2->SetParameter(2, 2287.);
      LcFitDG2->SetParameter(3, 7);
      LcFitDG2->SetParameter(4, 4);
      LcFitDG2->SetParLimits(3, 0., 20.);
      LcFitDG2->SetParLimits(4, 0., 20.);
      LcFitDG2->SetParameter(5, 0.);
      LcFitDG2->SetParameter(6, 0.);

      TF1 *LcFitDG3 = new TF1("LcFitDG3",DGOneMuOneTotal,2100.,2500.,7);
      LcFitDG3->SetParameter(0, 0.2);
      LcFitDG3->SetParameter(1, 20000);
      LcFitDG3->SetParameter(2, 2287.);
      LcFitDG3->SetParameter(3, 7);
      LcFitDG3->SetParameter(4, 4);
      LcFitDG3->SetParLimits(3, 0., 20.);
      LcFitDG3->SetParLimits(4, 0., 20.);
      LcFitDG3->SetParameter(5, 0.);
      LcFitDG3->SetParameter(6, 0.);

      TF1 *LcFitDG4 = new TF1("LcFitDG4",DGOneMuOneTotal,2100.,2500.,7);
      LcFitDG4->SetParameter(0, 0.2);
      LcFitDG4->SetParameter(1, 20000);
      LcFitDG4->SetParameter(2, 2287.);
      LcFitDG4->SetParameter(3, 7);
      LcFitDG4->SetParameter(4, 4);
      LcFitDG4->SetParLimits(3, 0., 20.);
      LcFitDG4->SetParLimits(4, 0., 20.);
      LcFitDG4->SetParameter(5, 0.);
      LcFitDG4->SetParameter(6, 0.);

      TF1 *XiFitDG1 = new TF1("XiFitDG1",DGOneMuOneTotal,2100.,2500.,7);
      XiFitDG1->SetParameter(0, 0.2);
      XiFitDG1->SetParameter(1, 20000);
      XiFitDG1->SetParameter(2, 2287.);
      XiFitDG1->SetParameter(3, 7);
      XiFitDG1->SetParameter(4, 4);
      XiFitDG1->SetParLimits(3, 0., 20.);
      XiFitDG1->SetParLimits(4, 0., 20.);
      XiFitDG1->SetParameter(5, 0.);
      XiFitDG1->SetParameter(6, 0.);

      TF1 *XiFitDG2 = new TF1("XiFitDG2",DGOneMuOneTotal,2100.,2500.,7);
      XiFitDG2->SetParameter(0, 0.2);
      XiFitDG2->SetParameter(1, 20000);
      XiFitDG2->SetParameter(2, 2287.);
      XiFitDG2->SetParameter(3, 7);
      XiFitDG2->SetParameter(4, 4);
      XiFitDG2->SetParLimits(3, 0., 20.);
      XiFitDG2->SetParLimits(4, 0., 20.);
      XiFitDG2->SetParameter(5, 0.);
      XiFitDG2->SetParameter(6, 0.);

      TF1 *XiFitDG3 = new TF1("XiFitDG3",DGOneMuOneTotal,2100.,2500.,7);
      XiFitDG3->SetParameter(0, 0.2);
      XiFitDG3->SetParameter(1, 20000);
      XiFitDG3->SetParameter(2, 2287.);
      XiFitDG3->SetParameter(3, 7);
      XiFitDG3->SetParameter(4, 4);
      XiFitDG3->SetParLimits(3, 0., 20.);
      XiFitDG3->SetParLimits(4, 0., 20.);
      XiFitDG3->SetParameter(5, 0.);
      XiFitDG3->SetParameter(6, 0.);

      TF1 *XiFitDG4 = new TF1("XiFitDG4",DGOneMuOneTotal,2100.,2500.,7);
      XiFitDG4->SetParameter(0, 0.2);
      XiFitDG4->SetParameter(1, 20000);
      XiFitDG4->SetParameter(2, 2287.);
      XiFitDG4->SetParameter(3, 7);
      XiFitDG4->SetParameter(4, 4);
      XiFitDG4->SetParLimits(3, 0., 20.);
      XiFitDG4->SetParLimits(4, 0., 20.);
      XiFitDG4->SetParameter(5, 0.);
      XiFitDG4->SetParameter(6, 0.);

      TF1 *MagDownFitDG = new TF1("MagDownFitDG",DGOneMuOneTotal,2100.,2500.,7);
      MagDownFitDG->SetParameter(0, 0.2);
      MagDownFitDG->SetParameter(1, 20000);
      MagDownFitDG->SetParameter(2, 2287.);
      MagDownFitDG->SetParameter(3, 7);
      MagDownFitDG->SetParameter(4, 4);
      MagDownFitDG->SetParLimits(3, 0., 20.);
      MagDownFitDG->SetParLimits(4, 0., 20.);
      MagDownFitDG->SetParameter(5, 0.);
      MagDownFitDG->SetParameter(6, 0.);

      TF1 *MagUpFitDG = new TF1("MagUpFitDG",DGOneMuOneTotal,2100.,2500.,7);
      MagUpFitDG->SetParameter(0, 0.2);
      MagUpFitDG->SetParameter(1, 20000);
      MagUpFitDG->SetParameter(2, 2287.);
      MagUpFitDG->SetParameter(3, 7);
      MagUpFitDG->SetParameter(4, 4);
      MagUpFitDG->SetParLimits(3, 0., 20.);
      MagUpFitDG->SetParLimits(4, 0., 20.);
      MagUpFitDG->SetParameter(5, 0.);
      MagUpFitDG->SetParameter(6, 0.);

      TF1 *ParticleFitDG = new TF1("ParticleFitDG",DGOneMuOneTotal,2100.,2500.,7);
      ParticleFitDG->SetParameter(0, 0.2);
      ParticleFitDG->SetParameter(1, 20000);
      ParticleFitDG->SetParameter(2, 2287.);
      ParticleFitDG->SetParameter(3, 7);
      ParticleFitDG->SetParameter(4, 4);
      ParticleFitDG->SetParLimits(3, 0., 20.);
      ParticleFitDG->SetParLimits(4, 0., 20.);
      ParticleFitDG->SetParameter(5, 0.);
      ParticleFitDG->SetParameter(6, 0.);

      TF1 *AntiParticleFitDG = new TF1("AntiParticleFitDG",DGOneMuOneTotal,2100.,2500.,7);
      AntiParticleFitDG->SetParameter(0, 0.2);
      AntiParticleFitDG->SetParameter(1, 20000);
      AntiParticleFitDG->SetParameter(2, 2287.);
      AntiParticleFitDG->SetParameter(3, 7);
      AntiParticleFitDG->SetParameter(4, 4);
      AntiParticleFitDG->SetParLimits(3, 0., 20.);
      AntiParticleFitDG->SetParLimits(4, 0., 20.);
      AntiParticleFitDG->SetParameter(5, 0.);
      AntiParticleFitDG->SetParameter(6, 0.);

      OptimalCut->Fit("myLambdaFitDG");
      OptimalCut->SetMinimum(0);
      c1->Write("OptimalCut - DG");

      LcLowP->Fit("LcFitDG1");
      LcLowP->SetMinimum(0);
      c1->Write("LcPLow - DG");

      LcMidLowP->Fit("LcFitDG2");
      LcMidLowP->SetMinimum(0);
      c1->Write("LcPMidLow - DG");

      LcMidHighP->Fit("LcFitDG3");
      LcMidHighP->SetMinimum(0);
      c1->Write("LcPMidHigh - DG");

      LcHighP->Fit("LcFitDG4");
      LcHighP->SetMinimum(0);
      c1->Write("LcPHigh - DG");

      LcPDistribution->Draw();
      LcPDistribution->SetMinimum(0);
      c1->Write("LcPDistribution - DG");

      double z1 = LcFitDG1->GetParameter(2);
      double z2 = LcFitDG2->GetParameter(2);
      double z3 = LcFitDG3->GetParameter(2);
      double z4 = LcFitDG4->GetParameter(2);

      double zerr1 = LcFitDG1->GetParError(2);
      double zerr2 = LcFitDG2->GetParError(2);
      double zerr3 = LcFitDG3->GetParError(2);
      double zerr4 = LcFitDG4->GetParError(2);

      double LcPAvGDG = ((z1 + z2 + z3 + z4)/n);
      Double_t z[n] = {z1 - LcPAvGDG,z2 - LcPAvGDG,z3 - LcPAvGDG,z4 - LcPAvGDG};
      Double_t zerr[n] = {zerr1,zerr2,zerr3,zerr4};
      gLcPDG = new TGraphErrors(n,x,z,xerr,zerr);
      gLcPDG->SetMarkerColor(2);
      gLcPDG->SetMarkerStyle(20);
      gLcPDG->GetXaxis()->SetNdivisions(4);
      gLcPDG->GetXaxis()->SetTitle("LambdaC_P Regions");
      gLcPDG->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
      gLcPDG->SetTitle("Deviations in LambdaC_MM of Different LambdaC_P Regions");
      gLcPDG->Draw("ALP");
      c1->Write("LcPDivision - DG");

      XiLowP->Fit("XiFitDG1");
      XiLowP->SetMinimum(0);
      c1->Write("XiPLow - DG");

      XiMidLowP->Fit("XiFitDG2");
      XiMidLowP->SetMinimum(0);
      c1->Write("XiPMidLow - DG");

      XiMidHighP->Fit("XiFitDG3");
      XiMidHighP->SetMinimum(0);
      c1->Write("XiPMidHigh - DG");

      XiHighP->Fit("XiFitDG4");
      XiHighP->SetMinimum(0);
      c1->Write("XiPHigh - DG");

      XiPDistribution->Draw();
      XiPDistribution->SetMinimum(0);
      c1->Write("XiPDistribution - DG");

    double zz1 = XiFitDG1->GetParameter(2);
    double zz2 = XiFitDG2->GetParameter(2);
    double zz3 = XiFitDG3->GetParameter(2);
    double zz4 = XiFitDG4->GetParameter(2);

    double zzerr1 = XiFitDG1->GetParError(2);
    double zzerr2 = XiFitDG2->GetParError(2);
    double zzerr3 = XiFitDG3->GetParError(2);
    double zzerr4 = XiFitDG4->GetParError(2);

    double XiPAvGDG = ((zz1 + zz2 + zz3 + zz4)/n);
    Double_t zz[n] = {zz1 - XiPAvGDG,zz2 - XiPAvGDG,zz3 - XiPAvGDG,zz4 - XiPAvGDG};
    Double_t zzerr[n] = {zzerr1,zzerr2,zzerr3,zzerr4};

    gXiPDG = new TGraphErrors(n,x,zz,xerr,zzerr);
    gXiPDG->SetMarkerColor(2);
    gXiPDG->SetMarkerStyle(20);
    gXiPDG->GetXaxis()->SetNdivisions(4);
    gXiPDG->GetXaxis()->SetTitle("Xi_P Regions");
    gXiPDG->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gXiPDG->SetTitle("Deviations in LambdaC_MM of Different Xi_P Regions");
    gXiPDG->Draw("ALP");
    c1->Write("XiPDivision - DG");

    PolarityMagDown->Fit("MagDownFitDG");
    PolarityMagDown->SetMinimum(0);
    c1->Write("PolarityMagDown - DG");

    PolarityMagUp->Fit("MagUpFitDG");
    PolarityMagUp->SetMinimum(0);
    c1->Write("PolarityMagUp - DG");

    double p1DG = MagDownFitDG->GetParameter(2);
    double p2DG = MagUpFitDG->GetParameter(2);

    double perrDG1 = MagDownFitDG->GetParError(2);
    double perrDG2 = MagUpFitDG->GetParError(2);

    double PolarityAvGDG = ((p1DG + p2DG)/m);
    Double_t pDG[m] = {p1DG - PolarityAvGDG,p2DG - PolarityAvGDG};
    Double_t perrDG[m] = {perrDG1,perrDG2};

    gPolarityDG = new TGraphErrors(m,x,pDG,xerr,perrDG);
    gPolarityDG->SetMarkerColor(2);
    gPolarityDG->SetMarkerStyle(20);
    gPolarityDG->GetXaxis()->SetNdivisions(2);
    gPolarityDG->GetXaxis()->SetTitle("Magnet Orientations");
    gPolarityDG->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gPolarityDG->SetTitle("Deviations in LambdaC_MM of Different Magnet Orientations");
    gPolarityDG->Draw("ALP");
    c1->Write("PolarityDivision - DG");

    Particle->Fit("ParticleFitDG");
    Particle->SetMinimum(0);
    c1->Write("Baryon - DG");

    AntiParticle->Fit("AntiParticleFitDG");
    AntiParticle->SetMinimum(0);
    c1->Write("AntiBaryon - DG");

    double pp1DG = ParticleFitDG->GetParameter(2);
    double pp2DG = AntiParticleFitDG->GetParameter(2);

    double pperrDG1 = ParticleFitDG->GetParError(2);
    double pperrDG2 = AntiParticleFitDG->GetParError(2);

    double ParticleAvGDG = ((pp1DG + pp2DG)/m);
    Double_t ppDG[m] = {pp1DG - ParticleAvGDG,pp2DG - ParticleAvGDG};
    Double_t pperrDG[m] = {pperrDG1,pperrDG2};

    gPIDDG = new TGraphErrors(m,x,ppDG,xerr,pperrDG);
    gPIDDG->SetMarkerColor(2);
    gPIDDG->SetMarkerStyle(20);
    gPIDDG->GetXaxis()->SetNdivisions(2);
    gPIDDG->GetXaxis()->SetTitle("Particle ID");
    gPIDDG->GetYaxis()->SetTitle("Measured Mass - Average Mass [MeV]");
    gPIDDG->SetTitle("Deviations in LambdaC_MM of Baryon or AntiBaryon");
    gPIDDG->Draw("ALP");
    c1->Write("ParticleDivision - DG");

      File->Close();
}
