#define Chained1D_cxx
#include "Chained1D.h"

#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>

#include "fit2MeV_Gaussian.C"

TCanvas *c1 = new TCanvas("Background Subtraction","Histograms",1000,500);
TFile * File = nullptr;

TH1 *PreliminaryMass = nullptr;

TH1 *PionIPCHI2Signal = nullptr;
TH1 *PionIPCHI2bkgd = nullptr;
TH1 *PionIPCHI2SignalEstimate = nullptr;

TH1 *PionProbSignal = nullptr;
TH1 *PionProbbkgd = nullptr;
TH1 *PionProbSignalEstimate = nullptr;

TH1 *KaonIPCHI2Signal = nullptr;
TH1 *KaonIPCHI2bkgd = nullptr;
TH1 *KaonIPCHI2SignalEstimate = nullptr;

TH1 *KaonProbSignal = nullptr;
TH1 *KaonProbbkgd = nullptr;
TH1 *KaonProbSignalEstimate = nullptr;

double UB1 = 1.;
double LB1 = 0.25;
const int Bin1 = 50;

double UB2 = 4.1;
double LB2 = 0.7;
const int Bin2 = 50;

void Chained1D::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void Chained1D::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

      File = new TFile("New1DProjection.root", "RECREATE");
     gFile = File;

   TH1::SetDefaultSumw2(kTRUE);

   PreliminaryMass = new TH1D("Mass", "Histogram of Preliminary Mass Cuts", 100, 2185, 2385);

   PionProbSignal = new TH1D("Probability", "Signal Region", Bin1, LB1, UB1);
   PionProbbkgd = new TH1D("Probability", "Background Region", Bin1, LB1, UB1);
   PionProbSignalEstimate = new TH1D("Probability", "Signal Estimation", Bin1, LB1, UB1);

   PionIPCHI2Signal = new TH1D("IPCHI2", "Signal Region", Bin2, LB2, UB2);
   PionIPCHI2bkgd = new TH1D("IPCHI2", "Background Region", Bin2, LB2, UB2);
   PionIPCHI2SignalEstimate = new TH1D("IPCHI2", "Signal Estimation", Bin2, LB2, UB2);

   KaonProbSignal = new TH1D("Probability", "Signal Region", Bin1, LB1, UB1);
   KaonProbbkgd = new TH1D("Probability", "Background Region", Bin1, LB1, UB1);
   KaonProbSignalEstimate = new TH1D("Probability", "Signal Estimation", Bin1, LB1, UB1);

   KaonIPCHI2Signal = new TH1D("IPCHI2", "Signal Region", Bin2, LB2, UB2);
   KaonIPCHI2bkgd = new TH1D("IPCHI2", "Background Region", Bin2, LB2, UB2);
   KaonIPCHI2SignalEstimate = new TH1D("IPCHI2", "Signal Estimation", Bin2, LB2, UB2);

   PionProbSignal->SetLineColor(kBlue);
   PionProbbkgd->SetLineColor(kRed);
   PionProbSignalEstimate->SetLineColor(kGreen+3);
   PionIPCHI2Signal->SetLineColor(kBlue);
   PionIPCHI2bkgd->SetLineColor(kRed);
   PionIPCHI2SignalEstimate->SetLineColor(kGreen+3);

   KaonProbSignal->SetLineColor(kBlue);
   KaonProbbkgd->SetLineColor(kRed);
   KaonProbSignalEstimate->SetLineColor(kGreen+3);
   KaonIPCHI2Signal->SetLineColor(kBlue);
   KaonIPCHI2bkgd->SetLineColor(kRed);
   KaonIPCHI2SignalEstimate->SetLineColor(kGreen+3);
}

Bool_t Chained1D::Process(Long64_t entry)
{
   GetEntry(entry);
   fReader.SetLocalEntry(entry);

      //Corrected Xi Mass
      double CorrectedXiMass = ((*Xi_MM) - (*Lambda_MM) + 1115.683);

      //Corrected Lambda Mass
      double CorrectedLambdaMass = ((*Lc_MM) - (*Xi_MM) + (1321.71));

      bool MassCuts = (
         (*Lc_MM > 2185. && *Lc_MM < 2385.)
 //       && (CorrectedXiMass > 1310. && CorrectedXiMass < 1330.)
   );

   bool PreliminaryCuts = (
   (*Lc_PT > 2000.)
   &&  (TMath::Log10(*PromptK_IPCHI2_OWNPV) > 0.7)
   &&	(TMath::Log10(*PromptPi_IPCHI2_OWNPV) > 0.7)
   &&  (TMath::Log10(*PromptK_IPCHI2_OWNPV) < 4.1)
   &&	(TMath::Log10(*PromptPi_IPCHI2_OWNPV) < 4.1)
   &&	(*PromptK_MC15TuneV1_ProbNNk > 0.25)
   &&	(*PromptPi_MC15TuneV1_ProbNNpi > 0.25)
   );

 if (PreliminaryCuts && MassCuts)
 PreliminaryMass->Fill(CorrectedLambdaMass);

 //Defining Signal Region & Background Region//

   bool SignalRegion = CorrectedLambdaMass > 2274. && CorrectedLambdaMass < 2300.;
//   bool BackgroundRegion = (CorrectedLambdaMass > 2215. && CorrectedLambdaMass < 2228.) || (CorrectedLambdaMass > 2342. && CorrectedLambdaMass < 2355.);
   bool BackgroundRegion = (CorrectedLambdaMass > 2185. && CorrectedLambdaMass < 2211.) || (CorrectedLambdaMass > 2359. && CorrectedLambdaMass < 2385.);

 ////////
 //Pion//
 ////////

   if (SignalRegion && PreliminaryCuts && MassCuts)
      PionIPCHI2Signal->Fill(TMath::Log10(*PromptPi_IPCHI2_OWNPV));

   if (BackgroundRegion && PreliminaryCuts && MassCuts)
     PionIPCHI2bkgd->Fill(TMath::Log10(*PromptPi_IPCHI2_OWNPV));

   if (SignalRegion && PreliminaryCuts && MassCuts)
         PionProbSignal->Fill(*PromptPi_MC15TuneV1_ProbNNpi);

   if (BackgroundRegion && PreliminaryCuts && MassCuts)
        PionProbbkgd->Fill(*PromptPi_MC15TuneV1_ProbNNpi);

//        PionProbSignalEstimate->Add(PionProbSignal,PionProbbkgd,1.0,-1.0);
//        PionIPCHI2SignalEstimate->Add(PionIPCHI2Signal,PionIPCHI2bkgd,1.0,-1.0);

        PionProbSignalEstimate->Add(PionProbSignal,PionProbbkgd,1.0,-0.5);
        PionIPCHI2SignalEstimate->Add(PionIPCHI2Signal,PionIPCHI2bkgd,1.0,-0.5);

 ////////
 //Kaon//
 ////////

   if (SignalRegion && PreliminaryCuts && MassCuts)
      KaonIPCHI2Signal->Fill(TMath::Log10(*PromptK_IPCHI2_OWNPV));

   if (BackgroundRegion && PreliminaryCuts && MassCuts)
     KaonIPCHI2bkgd->Fill(TMath::Log10(*PromptK_IPCHI2_OWNPV));

   if (SignalRegion && PreliminaryCuts && MassCuts)
         KaonProbSignal->Fill(*PromptK_MC15TuneV1_ProbNNk);

   if (BackgroundRegion && PreliminaryCuts && MassCuts)
        KaonProbbkgd->Fill(*PromptK_MC15TuneV1_ProbNNk);

//   KaonProbSignalEstimate->Add(KaonProbSignal,KaonProbbkgd,1.0,-1.0);
//   KaonIPCHI2SignalEstimate->Add(KaonIPCHI2Signal,KaonIPCHI2bkgd,1.0,-1.0);

   KaonProbSignalEstimate->Add(KaonProbSignal,KaonProbbkgd,1.0,-0.5);
   KaonIPCHI2SignalEstimate->Add(KaonIPCHI2Signal,KaonIPCHI2bkgd,1.0,-0.5);

   return kTRUE;
}

void Chained1D::SlaveTerminate()
{
}

void Chained1D::Terminate()
{

    PreliminaryMass->Draw();
    TPaveText *t0 = new TPaveText(0.3, 0.91, 0.7, 1.0, "brNDC");
    t0->AddText("Preliminary Mass Histogram");
    PreliminaryMass->SetMinimum(0);
    PreliminaryMass->GetYaxis()->SetTitle("Events per 2 MeV");
    PreliminaryMass->GetXaxis()->SetTitle("MeV");
    c1->Write("Preliminary Mass Histogram");

    TF1 *PionProbFit = new TF1("f1","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.25, 1.);

    TF1 *PionIPCHI2BackgroundFit = new TF1("f2","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.7, 4.1);

    TF1 *PionIPCHI2SignalFit = new TF1("f3","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.7, 4.1);


    TF1 *KaonProbFit = new TF1("f4","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.25, 1.);

    TF1 *KaonIPCHI2BackgroundFit = new TF1("f5","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.7, 4.1);

    TF1 *KaonIPCHI2SignalFit = new TF1("f6","[0] + [1]*x + [2]*x*x + exp([3] + [4]*x + [5]*x*x)", 0.7, 4.1);

  cout << "Pion IPCHI2 Signal" << endl;
    PionIPCHI2SignalEstimate->Fit("f3");
    PionIPCHI2SignalEstimate->GetFunction("f3")->SetLineColor(kGreen+3);
  cout << "Pion IPCHI2 Background" << endl;
    PionIPCHI2bkgd->Fit("f2");
    PionIPCHI2bkgd->GetFunction("f2")->SetLineColor(kRed);
  cout << "Pion Prob Signal" << endl;
    PionProbSignalEstimate->Fit("f1");
    PionProbSignalEstimate->GetFunction("f1")->SetLineColor(kGreen+3);
  cout << "Pion Prob Background" << endl;
    PionProbbkgd->Fit("f1");
    PionProbbkgd->GetFunction("f1")->SetLineColor(kRed);

  cout << "Kaon IPCHI2 Signal" << endl;
    KaonIPCHI2SignalEstimate->Fit("f6");
    KaonIPCHI2SignalEstimate->GetFunction("f6")->SetLineColor(kGreen+3);
  cout << "Kaon IPCHI2 Background" << endl;
    KaonIPCHI2bkgd->Fit("f5");
    KaonIPCHI2bkgd->GetFunction("f5")->SetLineColor(kRed);
  cout << "Kaon Prob Signal" << endl;
    KaonProbSignalEstimate->Fit("f4");
    KaonProbSignalEstimate->GetFunction("f4")->SetLineColor(kGreen+3);
  cout << "Kaon Prob Background" << endl;
    KaonProbbkgd->Fit("f4");
    KaonProbbkgd->GetFunction("f4")->SetLineColor(kRed);

  gStyle->SetOptTitle(0);
  TPaveText *t1 = new TPaveText(0.3, 0.91, 0.7, 1.0, "brNDC");
  t1->AddText("Signal and Background Estimation for PromptPi_IPCHI2_OWNPV");
  TPaveText *t2 = new TPaveText(0.3, 0.91, 0.7, 1.0, "brNDC");
  t2->AddText("Signal and Background Estimation for PromptPi_MC15TuneV1_ProbNNpi");

  TPaveText *t3 = new TPaveText(0.3, 0.91, 0.7, 1.0, "brNDC");
  t3->AddText("Signal and Background Estimation for PromptK_IPCHI2_OWNPV");
  TPaveText *t4 = new TPaveText(0.3, 0.91, 0.7, 1.0, "brNDC");
  t4->AddText("Signal and Background Estimation for PromptK_MC15TuneV1_ProbNNk");

  PionIPCHI2SignalEstimate->SetMaximum(22000);
  PionIPCHI2SignalEstimate->SetMinimum(0);
  KaonIPCHI2SignalEstimate->SetMaximum(22000);
  KaonIPCHI2SignalEstimate->SetMinimum(0);


    PionIPCHI2SignalEstimate->Draw();
    PionIPCHI2Signal->Draw("SAME");
    PionIPCHI2bkgd->Draw("SAME");
    t1->Draw("SAME");
    PionIPCHI2SignalEstimate->GetYaxis()->SetTitle("Events per 1 mm");
    PionIPCHI2SignalEstimate->GetXaxis()->SetTitle("IPCHI2");
    gPad->BuildLegend(0.78,0.75,0.98,0.95);
    c1->Write("Pion IPCHI2 Estimations");

    c1->Clear();
    PionProbSignalEstimate->Draw();
    PionProbSignal->Draw("SAME");
    PionProbbkgd->Draw("SAME");
    t2->Draw("SAME");
    PionProbSignalEstimate->GetYaxis()->SetTitle("Events per 1 mm");
    PionProbSignalEstimate->GetXaxis()->SetTitle("Probability");
    gPad->BuildLegend(0.78,0.75,0.98,0.95);
    c1->Write("Pion Prob Estimations");


    KaonIPCHI2SignalEstimate->Draw();
    KaonIPCHI2Signal->Draw("SAME");
    KaonIPCHI2bkgd->Draw("SAME");
    t3->Draw("SAME");
    KaonIPCHI2SignalEstimate->GetYaxis()->SetTitle("Events per 1 mm");
    KaonIPCHI2SignalEstimate->GetXaxis()->SetTitle("IPCHI2");
    gPad->BuildLegend(0.78,0.75,0.98,0.95);
    c1->Write("Kaon IPCHI2 Estimations");

    c1->Clear();
    KaonProbSignalEstimate->Draw();
    KaonProbSignal->Draw("SAME");
    KaonProbbkgd->Draw("SAME");
    t4->Draw("SAME");
    KaonProbSignalEstimate->GetYaxis()->SetTitle("Events per 1 mm");
    KaonProbSignalEstimate->GetXaxis()->SetTitle("Probability");
    gPad->BuildLegend(0.78,0.75,0.98,0.95);
    c1->Write("Kaon Prob Estimations");
}
