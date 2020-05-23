#define ChainedStatistics_cxx
#include "ChainedStatistics.h"


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

      TF1 *myLambdaFit = new TF1("myLambdaFit",fit2MeV_Gaussian,2100.,2500.,5);
      myLambdaFit->SetParameter(0,400.);
      myLambdaFit->SetParameter(1,2286.5);
      myLambdaFit->SetParameter(2, 6);
      myLambdaFit->SetParLimits(2, 0.,20.);
      myLambdaFit->SetParameter(3, 0.);
      myLambdaFit->SetParameter(4, 0.);

      TF1 *LcFit1 = new TF1("LcFit1",fit2MeV_Gaussian,2100.,2500.,5);
      LcFit1->SetParameter(0,400.);
      LcFit1->SetParameter(1,2286.5);
      LcFit1->SetParameter(2, 6);
      LcFit1->SetParLimits(2, 0.,20.);
      LcFit1->SetParameter(3, 0.);
      LcFit1->SetParameter(4, 0.);

      TF1 *LcFit2 = new TF1("LcFit2",fit2MeV_Gaussian,2100.,2500.,5);
      LcFit2->SetParameter(0,400.);
      LcFit2->SetParameter(1,2286.5);
      LcFit2->SetParameter(2, 6);
      LcFit2->SetParLimits(2, 0.,20.);
      LcFit2->SetParameter(3, 0.);
      LcFit2->SetParameter(4, 0.);

      TF1 *LcFit3 = new TF1("LcFit3",fit2MeV_Gaussian,2100.,2500.,5);
      LcFit3->SetParameter(0,400.);
      LcFit3->SetParameter(1,2286.5);
      LcFit3->SetParameter(2, 6);
      LcFit3->SetParLimits(2, 0.,20.);
      LcFit3->SetParameter(3, 0.);
      LcFit3->SetParameter(4, 0.);

      TF1 *LcFit4 = new TF1("LcFit4",fit2MeV_Gaussian,2100.,2500.,5);
      LcFit4->SetParameter(0,400.);
      LcFit4->SetParameter(1,2286.5);
      LcFit4->SetParameter(2, 6);
      LcFit4->SetParLimits(2, 0.,20.);
      LcFit4->SetParameter(3, 0.);
      LcFit4->SetParameter(4, 0.);

      TF1 *XiFit1 = new TF1("XiFit1",fit2MeV_Gaussian,2100.,2500.,5);
      XiFit1->SetParameter(0,400.);
      XiFit1->SetParameter(1,2286.5);
      XiFit1->SetParameter(2, 6);
      XiFit1->SetParLimits(2, 0.,20.);
      XiFit1->SetParameter(3, 0.);
      XiFit1->SetParameter(4, 0.);

      TF1 *XiFit2 = new TF1("XiFit2",fit2MeV_Gaussian,2100.,2500.,5);
      XiFit2->SetParameter(0,400.);
      XiFit2->SetParameter(1,2286.5);
      XiFit2->SetParameter(2, 6);
      XiFit2->SetParLimits(2, 0.,20.);
      XiFit2->SetParameter(3, 0.);
      XiFit2->SetParameter(4, 0.);

      TF1 *XiFit3 = new TF1("XiFit3",fit2MeV_Gaussian,2100.,2500.,5);
      XiFit3->SetParameter(0,400.);
      XiFit3->SetParameter(1,2286.5);
      XiFit3->SetParameter(2, 6);
      XiFit3->SetParLimits(2, 0.,20.);
      XiFit3->SetParameter(3, 0.);
      XiFit3->SetParameter(4, 0.);

      TF1 *XiFit4 = new TF1("XiFit4",fit2MeV_Gaussian,2100.,2500.,5);
      XiFit4->SetParameter(0,400.);
      XiFit4->SetParameter(1,2286.5);
      XiFit4->SetParameter(2, 6);
      XiFit4->SetParLimits(2, 0.,20.);
      XiFit4->SetParameter(3, 0.);
      XiFit4->SetParameter(4, 0.);

      TF1 *MagDownFit = new TF1("MagDownFit",fit2MeV_Gaussian,2100.,2500.,5);
      MagDownFit->SetParameter(0,400.);
      MagDownFit->SetParameter(1,2286.5);
      MagDownFit->SetParameter(2, 6);
      MagDownFit->SetParLimits(2, 0.,20.);
      MagDownFit->SetParameter(3, 0.);
      MagDownFit->SetParameter(4, 0.);

      TF1 *MagUpFit = new TF1("MagUpFit",fit2MeV_Gaussian,2100.,2500.,5);
      MagUpFit->SetParameter(0,400.);
      MagUpFit->SetParameter(1,2286.5);
      MagUpFit->SetParameter(2, 6);
      MagUpFit->SetParLimits(2, 0.,20.);
      MagUpFit->SetParameter(3, 0.);
      MagUpFit->SetParameter(4, 0.);

      TF1 *ParticleFit = new TF1("ParticleFit",fit2MeV_Gaussian,2100.,2500.,5);
      ParticleFit->SetParameter(0,400.);
      ParticleFit->SetParameter(1,2286.5);
      ParticleFit->SetParameter(2, 6);
      ParticleFit->SetParLimits(2, 0.,20.);
      ParticleFit->SetParameter(3, 0.);
      ParticleFit->SetParameter(4, 0.);

      TF1 *AntiParticleFit = new TF1("AntiParticleFit",fit2MeV_Gaussian,2100.,2500.,5);
      AntiParticleFit->SetParameter(0,400.);
      AntiParticleFit->SetParameter(1,2286.5);
      AntiParticleFit->SetParameter(2, 6);
      AntiParticleFit->SetParLimits(2, 0.,20.);
      AntiParticleFit->SetParameter(3, 0.);
      AntiParticleFit->SetParameter(4, 0.);

      OptimalCut->Fit("myLambdaFit");
      OptimalCut->SetMinimum(0);
      c1->Write("OptimalCut");

      LcLowP->Fit("LcFit1");
      LcLowP->SetMinimum(0);
      c1->Write("LcPLow");

      LcMidLowP->Fit("LcFit2");
      LcMidLowP->SetMinimum(0);
      c1->Write("LcPMidLow");

      LcMidHighP->Fit("LcFit3");
      LcMidHighP->SetMinimum(0);
      c1->Write("LcPMidHigh");

      LcHighP->Fit("LcFit4");
      LcHighP->SetMinimum(0);
      c1->Write("LcPHigh");

      LcPDistribution->Draw();
      LcPDistribution->SetMinimum(0);
      c1->Write("LcPDistribution");

      double y1 = LcFit1->GetParameter(1);
      double y2 = LcFit2->GetParameter(1);
      double y3 = LcFit3->GetParameter(1);
      double y4 = LcFit4->GetParameter(1);

      double yerr1 = LcFit1->GetParError(1);
      double yerr2 = LcFit2->GetParError(1);
      double yerr3 = LcFit3->GetParError(1);
      double yerr4 = LcFit4->GetParError(1);

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

      XiLowP->Fit("XiFit1");
      XiLowP->SetMinimum(0);
      c1->Write("XiPLow");

      XiMidLowP->Fit("XiFit2");
      XiMidLowP->SetMinimum(0);
      c1->Write("XiPMidLow");

      XiMidHighP->Fit("XiFit3");
      XiMidHighP->SetMinimum(0);
      c1->Write("XiPMidHigh");

      XiHighP->Fit("XiFit4");
      XiHighP->SetMinimum(0);
      c1->Write("XiPHigh");

      XiPDistribution->Draw();
      XiPDistribution->SetMinimum(0);
      c1->Write("XiPDistribution");

    double yy1 = XiFit1->GetParameter(1);
    double yy2 = XiFit2->GetParameter(1);
    double yy3 = XiFit3->GetParameter(1);
    double yy4 = XiFit4->GetParameter(1);

    double yyerr1 = XiFit1->GetParError(1);
    double yyerr2 = XiFit2->GetParError(1);
    double yyerr3 = XiFit3->GetParError(1);
    double yyerr4 = XiFit4->GetParError(1);

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

    PolarityMagDown->Fit("MagDownFit");
    PolarityMagDown->SetMinimum(0);
    c1->Write("PolarityMagDown");

    PolarityMagUp->Fit("MagUpFit");
    PolarityMagUp->SetMinimum(0);
    c1->Write("PolarityMagUp");

    double p1 = MagDownFit->GetParameter(1);
    double p2 = MagUpFit->GetParameter(1);

    double perr1 = MagDownFit->GetParError(1);
    double perr2 = MagUpFit->GetParError(1);

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

    Particle->Fit("ParticleFit");
    Particle->SetMinimum(0);
    c1->Write("Baryon");

    AntiParticle->Fit("AntiParticleFit");
    AntiParticle->SetMinimum(0);
    c1->Write("AntiBaryon");

    double pp1 = ParticleFit->GetParameter(1);
    double pp2 = AntiParticleFit->GetParameter(1);

    double pperr1 = ParticleFit->GetParError(1);
    double pperr2 = AntiParticleFit->GetParError(1);

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

      File->Close();
}
