#define DalitzPlot_cxx
#include "DalitzPlot.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

TH2D * DalitzPlotDs = nullptr;
TH1D * KpKmMassHist = nullptr;
TH1D * PipKmMassHist = nullptr;
TH1D * PipKpMassHist = nullptr;

TFile * File = nullptr;

TCanvas * c1 = nullptr;

void DalitzPlot::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
     DalitzPlotDs = new TH2D("Dalitz Plot", "Dalitz Plot", 100, 0.8, 2.2, 100, 0.5, 2.0);
         DalitzPlotDs->GetXaxis()->SetTitle("m^{2}(K^{-}K^{+})[GeV^{2}/c^{4}]");
         DalitzPlotDs->GetYaxis()->SetTitle("m^{2}(pi^{+}K^{-})[GeV^{2}/c^{4}]");
         DalitzPlotDs->GetZaxis()->SetTitle("Events");
                                          
         KpKmMassHist = new TH1D("M^{2} [GeV^{2}/c^{4}]", "Kplus & Kplus Invariant Mass Combination", 100, 0.95, 2);
         KpKmMassHist->GetXaxis()->SetTitle("m^{2}(K^{-}K^{+})[GeV^{2}/c^{4}]");                    
         KpKmMassHist->GetYaxis()->SetTitle("Events");
 
         PipKmMassHist = new TH1D("M^{2} [GeV^{2}/c^{4}]", "Piplus & Kminus Invariant Mass Combination", 100, 0.5, 2);
         PipKmMassHist->GetXaxis()->SetTitle("m^{2}(pi^{+}K^{-})[GeV^{2}/c^{4}]");                   
         PipKmMassHist->GetYaxis()->SetTitle("Events");
   
         PipKpMassHist = new TH1D("M^{2} [GeV^{2}/c^{4}]", "Piplus & Kplus Invariant Mass Combination", 100, 0.5, 2);
         PipKpMassHist->GetXaxis()->SetTitle("m^{2}(pi^{+}K^{+})[GeV^{2}/c^{4}]");                   
         PipKpMassHist->GetYaxis()->SetTitle("Events");  
   
       File = new TFile("DalitzDs.root", "RECREATE");
  gFile = File;

   c1 = new TCanvas("canvas", "Test Canvas");
}

void DalitzPlot::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t DalitzPlot::Process(Long64_t entry)
{
  GetEntry(entry);
   fReader.SetLocalEntry(entry);
 
  double P_Pip  = *Piplus_P;
  double P_Kp = *Kplus_P;
  double P_Km = *Kminus_P;
 
  double M_Pip  = *Piplus_M;
  double M_Kp = *Kplus_M;
  double M_Km = *Kminus_M;
     
  double E_Pip  = TMath::Sqrt(((P_Pip)*(P_Pip))+((M_Pip)*(M_Pip)));
  double E_Kp = TMath::Sqrt(((P_Kp)*(P_Kp))+((M_Kp)*(M_Kp)));
  double E_Km = TMath::Sqrt(((P_Km)*(P_Km))+((M_Km)*(M_Km)));
   
double M2_KpKm = ((((E_Kp)+(E_Km))*((E_Kp)+(E_Km))) - (((P_Kp)+(P_Km))*((P_Kp)+(P_Km))))/(1000*1000);
double M2_PipKm  = ((((E_Pip)+(E_Km))*((E_Pip)+(E_Km))) - (((P_Pip)+(P_Km))*((P_Pip)+(P_Km))))/(1000*1000);
double M2_PipKp  = ((((E_Pip)+(E_Kp))*((E_Pip)+(E_Kp))) - (((P_Pip)+(P_Kp))*((P_Pip)+(P_Kp))))/(1000*1000);

 KpKmMassHist->Fill(M2_KpKm);
 PipKmMassHist->Fill(M2_PipKm);
 PipKpMassHist->Fill(M2_PipKp);
 DalitzPlotDs->Fill(M2_KpKm, M2_PipKm);
   return kTRUE;
}

void DalitzPlot::SlaveTerminate()
{
}

void DalitzPlot::Terminate()
{
c1->cd();
DalitzPlotDs->Draw();
 c1->Write("Dalitz Plot");
DalitzPlotDs->Draw("COLZ");
 c1->Write("Dalitz Plot - COLZ");
DalitzPlotDs->Draw("CONTZ");
 c1->Write("Dalitz Plot - CONTZ");
   
KpKmMassHist->Draw();
 c1->Write("Kp & Km Mass"); 
   
 PipKmMassHist->Draw();
 c1->Write("Pip & Km Mass");
   
 PipKpMassHist->Draw();
 c1->Write("Pip & Kp Mass");
 File->Close();
}
