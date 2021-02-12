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
void InclusiveEleDCA()
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
    TH1F *InclsElecPt = (TH1F*)QA_emc1->FindObject("fInclsElecPt");
    TH2F *InclElecDCAPt = (TH2F*)QA_emc1->FindObject("fInclElecDCA");

    cout << "N events Cent 0-10%: " << Centality->GetEntries()<< ", " << NEvents->GetBinContent(3) << endl;
    
    Double_t PtBin[5] = {3,4,5,6,8};
    TString PtBinNames[4] = {"3 < p_{T} < 4 GeV/c", "4 < p_{T} < 5 GeV/c", "5 < p_{T} < 6 GeV/c", "6 < p_{T} < 8 GeV/c"};

    TH1D *IncEDCA[10];
    for(int i=0;i<4; i++){
        
        //DCA for inclusive electrons//
        IncEDCA[i] = InclElecDCAPt->ProjectionY(Form("IncEDCA_%i",i),InclElecDCAPt->GetXaxis()->FindBin(PtBin[i]+0.1), InclElecDCAPt->GetXaxis()->FindBin(PtBin[i+1] - 0.1));
        IncEDCA[i]->SetTitle(PtBinNames[i]);
    }
    
    TCanvas *c3 = new TCanvas("DCA_IncwHad", "DCA-IncE, Pt bins, with HadCont",50,50,1200,700);
    c3->Divide(2,2);
    for(int i=0;i<4;i++){
        c3->cd(i+1);
        gPad->SetLeftMargin(0.13);
        //gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.13);
        gPad->SetLogy();
        
        IncEDCA[i]->Sumw2();
        ProcessHisto(IncEDCA[i],1.2,kBlack,20);
        
        IncEDCA[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
        IncEDCA[i]->Draw();
        
        c3->cd(2);
        TLegend *leg1 = new TLegend(0.33,0.72,0.58,0.89);
        leg1->AddEntry(IncEDCA[0],"Inclusive Elec","pl");
        ProcessLegend(leg1);
    }
}

