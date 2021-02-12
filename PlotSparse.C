#include "TLegend.h"
#include "THnSparse.h"

void ProcessGraph(TGraph *h, const TString gname="gn", Int_t sty=20, Int_t col=1, Double_t size=1.4)
{
    h->SetTitle(gname);
    h->SetMarkerStyle(sty);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetMarkerSize(size);
}
void ProcessHisto(TH1 *h, Double_t size=1.4, Int_t col=1, Int_t style=20)
{
    gPad->SetTickx();
    gPad->SetTicky();
    h->SetMarkerSize(size);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetMarkerStyle(style);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetTitleFont(42);
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.045);
}
void ProcessHisto2D(TH2 *h)
{
    //  h->SetLogz();
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetTitleFont(42);
    h->GetXaxis()->SetTitleFont(42);
    h->GetZaxis()->SetLabelFont(42);
    h->GetZaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
}
void ProcessLegend(TLegend *leg)
{
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();
}
TH1D *projectInPtbin(THnSparse *sparse, Int_t Axis,const Double_t Ptmin=0,const Double_t Ptmax=1000, const TString hname = "hnamePro")
{
    sparse->GetAxis(0)->SetRangeUser(Ptmin, Ptmax);
    TH1D *h = sparse->Projection(Axis,"");
    h->SetName(hname.Data());
    //h->Sumw2();
    return h;
}

//---------------------------------------
void PlotSparse()
{
    gROOT->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    
    TFile *f1 = new TFile("AnalysisResults.root");
    TDirectory *QA_emc = (TDirectory*)f1->Get("PWGHF_HFEBESpectraEMC_CentCentralcaloClusters");
    TList *QA_emc1 = (TList*)QA_emc->Get("HFEBESpectraEMC_Cent_wTend_0_10CentralcaloClusters");
    TH1F *NEvents = (TH1F*)QA_emc1->FindObject("fNevents");
    TH1F *Centality = (TH1F*)QA_emc1->FindObject("fCent");
    
    THnSparse *electron_EMC = (THnSparse *) QA_emc1->FindObject("Electron");
    
    cout << "N events Cent 0-10%: " << Centality->GetEntries()<< ", " << NEvents->GetBinContent(3) << endl;
    
    Double_t pTbin[5] = {3,4,5,6,8};
    TString PtBinNames[4] = {"3 < p_{T} < 4 GeV/c", "4 < p_{T} < 5 GeV/c", "5 < p_{T} < 6 GeV/c", "6 < p_{T} < 8 GeV/c"};
    
    electron_EMC->GetAxis(1)->SetRangeUser(-1,3); //nsig
    electron_EMC->GetAxis(3)->SetRangeUser(0.,20000); //M20
    electron_EMC->GetAxis(4)->SetRangeUser(0.05,0.9); //M02
    
    TH1 *Eovp_EMC[4];
    for(int i=0; i<4; i++){
        Eovp_EMC[i] = projectInPtbin(electron_EMC, 2, pTbin[i]+0.01, pTbin[i+1]-0.05); //E/p
        Eovp_EMC[i]->SetName(Form("Eovp_EMC_%i",i));
        Eovp_EMC[i]->SetTitle(PtBinNames[i]);
    }
    
    electron_EMC->GetAxis(1)->SetRangeUser(-10,-4); //nsig
    TH1 *HEovp_EMC[12];
    for(int i=0; i<4; i++){
        HEovp_EMC[i] = projectInPtbin(electron_EMC, 2, pTbin[i]+0.01, pTbin[i+1]-0.05); //hadron E/p
        HEovp_EMC[i]->SetName(Form("HEovp_EMC_%i",i));
        HEovp_EMC[i]->SetTitle(PtBinNames[i]);
    }

    Int_t a=0,b=0;
     TH1D *ElePure[12];
     Double_t Yield_Ele[12], Yield_Ele_Err[12], Yield_Ele_Denom[12];
     Double_t Yield_Ele_CC[12];
     TCanvas *c3 = new TCanvas("EovP_Ele", "E/p distribution with hadron contamination",50,50,1000,600);
     c3->Divide(2,2);
     for(int i=0;i<4;i++){
         c3->cd(i+1);
         // gPad->SetLogy();
         gPad->SetLeftMargin(0.11);
         //gPad->SetRightMargin(0.05);
         gPad->SetBottomMargin(0.10);

         Eovp_EMC[i]->Sumw2();
         HEovp_EMC[i]->Sumw2();
         //HEovp_EMC[i]->Rebin();

         a = Eovp_EMC[1]->FindBin(0.35);
         b = Eovp_EMC[1]->FindBin(0.45);

         Double_t E=0, H=0;
         Int_t j=a;
         while(j<b+1){
             E = E + Eovp_EMC[i]->GetBinContent(j);
             H = H + HEovp_EMC[i]->GetBinContent(j);
             j=j+1;
         }
         double EH = E/H;
         HEovp_EMC[i]->Scale(EH);
         
         ProcessHisto(Eovp_EMC[i],1.4,1,20);
         ProcessHisto(HEovp_EMC[i],1.4,2,22);

         ElePure[i] = (TH1D *)Eovp_EMC[i]->Clone();
         ElePure[i]->Reset();
         ElePure[i]->Add(Eovp_EMC[i],HEovp_EMC[i],1,-1);
         ProcessHisto(ElePure[i],1.4,4,24);

         Eovp_EMC[i]->GetXaxis()->SetTitle("E/p");
         Eovp_EMC[i]->GetYaxis()->SetTitle("dN");
         Eovp_EMC[i]->GetXaxis()->SetRangeUser(0.2,1.8);
         Eovp_EMC[i]->Draw();
         HEovp_EMC[i]->Draw("same");
         ElePure[i]->Draw("same");

         if(i==0) cout << ElePure[i]->FindBin(0.81) << ", " << ElePure[i]->FindBin(1.19) << endl;
         Yield_Ele[i] = ElePure[i]->IntegralAndError(ElePure[i]->FindBin(0.81), ElePure[i]->FindBin(1.19), Yield_Ele_Err[i]);

         Eovp_EMC[i]->GetXaxis()->SetRangeUser(0.8,1.2);
         HEovp_EMC[i]->GetXaxis()->SetRangeUser(0.8,1.2);
         Yield_Ele_CC[i] = Eovp_EMC[i]->Integral() - HEovp_EMC[i]->Integral();
         cout << "yield check and err: " << Yield_Ele[i] << ", " << Yield_Ele_CC[i] << ", " << Yield_Ele_Err[i]/Yield_Ele[i] << endl;
         Yield_Ele_Denom[i] = Eovp_EMC[i]->Integral();
         cout << "Electron purity : " << Yield_Ele[i]/(Eovp_EMC[i]->Integral()) << endl;

         Eovp_EMC[i]->GetXaxis()->SetRangeUser(0.2,1.8);
         HEovp_EMC[i]->GetXaxis()->SetRangeUser(0.2,1.8);
         ElePure[i]->GetXaxis()->SetRangeUser(0.2,1.8);

         if(i==1){
             TLegend *leg = new TLegend(0.157,0.736,0.407,0.855);
             leg->AddEntry(Eovp_EMC[1],"Elec candidate","p");
             leg->AddEntry(HEovp_EMC[1],"Hadrons","p");
             leg->AddEntry(ElePure[1],"Elec Cand - Hadrons","p");
             ProcessLegend(leg);
             leg->SetTextSize(0.055);
         }
     }
}
