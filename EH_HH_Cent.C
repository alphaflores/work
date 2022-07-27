#include "TLegend.h"
#include "THnSparse.h"
#include "TMath.h"

/* Plots electron to hadron and hadron to hadron correlations in SE and ME Pb-Pb collisions 
 *
 * This macro is a copy of E_Had_Corr1 but with Mixed Event Correction techniques applied
 */

void ProcessGraph(TGraph *h, const TString gname="gn", Int_t sty=20, Int_t col=1, Double_t size=1.4)
{
    h->SetTitle(gname);
    h->SetMarkerStyle(sty);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetMarkerSize(size);
}
void ProcessHisto(TH1 *h, Double_t size=1.4, Int_t col=1, Int_t style=20)//Function to set histograms into one style
{
    //gPad->SetTickx();
    //gPad->SetTicky();
    //h->SetMarkerSize(size);
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


TH2 *Delta_2DSelectEH(THnSparse *sparse3, Int_t Axis2D1, Int_t Axis2D2, const Double_t PtTrigmin3=0, const Double_t PtTrigmax3=1000, const Double_t PtAssomin3=0, const Double_t PtAssomax3=1000, const TString hname3 = "DP2DEHname")
{
    sparse3->GetAxis(0)->SetRangeUser(PtTrigmin3, PtTrigmax3);//Electron trigger
    sparse3->GetAxis(1)->SetRangeUser(PtAssomin3, PtAssomax3);//Associate hadron
    TH2 *h = sparse3->Projection(Axis2D2,Axis2D1,"");//Project Axis 2 (delta phi) and Axis 3 (delta eta)
    h->SetName(hname3.Data());
    //h->Sumw2();
    return h;
}


void ProcessTPave(TPaveText* t2, Int_t font=42, Double_t size=0.06)
{
    t2->SetFillStyle(0);
    t2->SetBorderSize(0);
    t2->SetTextColor(kBlack);
    t2->SetTextSize(size);
    t2->SetTextFont(font);
    t2->Draw();
}

Double_t GetScaleFactor(TH2 *h)
{
  Int_t Etabin[2], Phibin[2];
  Phibin[0] = h->GetXaxis()->FindBin(-0.01);
  Phibin[1] = h->GetXaxis()->FindBin(0.01);
  Etabin[0] = h->GetYaxis()->FindBin(-0.01);
  Etabin[1] = h->GetYaxis()->FindBin(0.01);

  Double_t GlobalScale=0;
  Int_t Nbins = 0;
  for(int k=0;k<2;k++){
    for(int l=0;l<2;l++){
      Int_t bin = h->GetBin(Phibin[k],Etabin[l]);
      GlobalScale = GlobalScale + h->GetBinContent(bin);
      Nbins++;
    }
  }
  GlobalScale = GlobalScale/Nbins;

  return GlobalScale;
}

//---------------------------------------
void EH_HH_Cent() {
    gROOT->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
 
    TFile *file = new TFile("AnalysisResults_PbPblhc18.root","READ");

    TDirectoryFile *EH_020 = (TDirectoryFile*)file->Get("EH_020_EH_CentTrig"); //subfile with INT7 trigger results
    TH1F *NEvents = (TH1F*)EH_020->FindObject("fNevents");//Get number of events from central clusters
    int events = 0;
    int eventselect = 0;
    events = NEvents->GetEntries();
    eventselect = NEvents->GetBinContent(3);
    //cout << "Number of events (Centrality trigger): " << events << endl;
    //cout << "Number of events selected (Centrality trigger): " << eventselect << endl;
    TH1F *Centrality = (TH1F*)EH_020->FindObject("fCentralityPass");//Get centrality passed
    int centrality_events = 0;
    centrality_events = Centrality->GetEntries();
    cout << "Number of events in 0 - 10% centrality: " << centrality_events << endl;

    //%%%%%%%%%
    //THnSparse with electron ID information
    //axes: 0=pTe, 1=pTh, 2=deltaphi, 3=deltaeta   
    //%%%%%%%%%
    //
    //H-H
    THnSparse *SprsHadHCorrl = (THnSparse *) EH_020->FindObject("fSprsHadHCorrl");//Getting EMCal electron THnSparse object
    THnSparse *SprsMixHadHCorrl = (THnSparse *) EH_020->FindObject("fSprsMixHadHCorrl");//Getting EMCal electron THnSparse object

    //E-H 
    THnSparse *SprsInclusiveEHCorrl = (THnSparse *) EH_020->FindObject("fSprsInclusiveEHCorrl");//Getting EMCal electron THnSparse object (inclusive ele-h)
    THnSparse *SprsMixInclusiveEHCorrl = (THnSparse *) EH_020->FindObject("fSprsMixInclusiveEHCorrl");//Getting EMCal electron THnSparse object (Mixed event; inclusive ele-h)

    Double_t pTbinTrigger1[2] = {4,12};
    Double_t ptbinAsso1[6] = {1,2,3,4,5,7};
    TString PtBinNames[5] = {"1 < Hp_{T} < 2 GeV/c", "2 < Hp_{T} < 3 GeV/c", "3 < Hp_{T} < 4 GeV/c", "4 < Hp_{T} < 5 GeV/c", "5 < Hp_{T} < 7 GeV/c"};

    //1D Delta Phi (-1 < Delta Eta < 1)
    //
    SprsHadHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta
    SprsMixHadHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta
    SprsInclusiveEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta
    SprsMixInclusiveEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta

    //Obtain EH Correlation Histograms

    //InclusiveE-H 2D histograms
    TH2 *DPhiDEtaEH_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaEH_1[i] = Delta_2DSelectEH(SprsInclusiveEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaEH_1[i]->SetName(Form("Delta Phi Delta Eta EH SE_%i",i));
        DPhiDEtaEH_1[i]->SetTitle(PtBinNames[i]);
    }
    //InclusiveE-H Mixed event
    TH2 *DPhiDEtaEHMix_1[6];//2D Delta phi delta eta histograms mixed event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaEHMix_1[i] = Delta_2DSelectEH(SprsMixInclusiveEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaEHMix_1[i]->SetName(Form("Delta Phi Delta Eta EH ME_%i",i));
        DPhiDEtaEHMix_1[i]->SetTitle(PtBinNames[i]);
    }

    //H-H 2D histograms
    TH2 *DPhiDEtaHH_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaHH_1[i] = Delta_2DSelectEH(SprsHadHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaHH_1[i]->SetName(Form("Delta Phi Delta Eta HH SE_%i",i));
        DPhiDEtaHH_1[i]->SetTitle(PtBinNames[i]);
    }
    //H-H Mixed event
    TH2 *DPhiDEtaHHMix_1[6];//2D Delta phi delta eta histograms mixed event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaHHMix_1[i] = Delta_2DSelectEH(SprsMixHadHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaHHMix_1[i]->SetName(Form("Delta Phi Delta Eta HH ME_%i",i));
        DPhiDEtaHHMix_1[i]->SetTitle(PtBinNames[i]);
    }
    
    //Histograms obtained from ThNSparses
    //DPhiDEtaEH_1
    //DPhiDEtaEHMix_1

    //Show different scale factor binning
    TH1 *DPhiScaleEHMix1[6];
    TH1 *DPhiScaleEHMix2[6];

/*for(int i=0;i<5;i++){
    DPhiScaleEHMix1[i] = DPhiDEtaEHMix_1[i]->ProjectionY(Form("IncDp%i",i),25,26); //25, 26 is the bin nos of DeltaEta(-0.01,0.01)
  }
  TCanvas *cChk = new TCanvas("MEScaleFactorDPhiIncl_pPb", "Scale factors,EH Dphi mixed event",50,50,1000,800);
  cChk->Divide(3,2);
  for(int i=0;i<6;i++){
    cChk->cd(i+1);
    if(i<5){
      DPhiScaleEHMix1[i]->Sumw2();
      ProcessHisto(DPhiScaleEHMix1[i]);
      //DPhiScaleEHMix1[i]->Rebin(2);

      DPhiScaleEHMix2[i] = (TH1 *) DPhiScaleEHMix1[i]->Clone();
      ProcessHisto(DPhiScaleEHMix2[i],1.4,2);

      Double_t scale = DPhiScaleEHMix1[i]->GetBinContent(DPhiScaleEHMix1[i]->FindBin(-0.01)) + DPhiScaleEHMix1[i]->GetBinContent(DPhiScaleEHMix1[i]->FindBin(0.01));
      scale = scale/2.0;
      DPhiScaleEHMix1[i]->Scale(1/scale);

      Double_t scaleClone = DPhiScaleEHMix2[i]->GetEntries();
      scaleClone = scaleClone/(DPhiScaleEHMix2[i]->GetNbinsX());
      DPhiScaleEHMix2[i]->Scale(1/scaleClone);

      DPhiScaleEHMix1[i]->GetYaxis()->SetRangeUser(0,1.4*DPhiScaleEHMix1[i]->GetBinContent(DPhiScaleEHMix1[i]->GetMaximumBin()));
      DPhiScaleEHMix1[i]->Draw();
      DPhiScaleEHMix2[i]->Draw("same");
    }
    if(i==5){
      TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
      ProcessLegend(leg);
      leg->SetTextSize(0.05);
      leg->AddEntry(DPhiScaleEHMix1[0],"Avg (#Delta#eta, #Delta#phi #approx 0,0)","pl");
      leg->AddEntry(DPhiScaleEHMix2[0],"Avg (#Delta#eta #approx 0)","pl");
      leg->Draw();

      TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
      t1->AddText("Scale factor calculation");
      t1->AddText("4 < p_{T}^{trigE} < 12 GeV/c");
      ProcessTPave(t1);
    }
  }
*/
  //Scale EH ME
  //Plot EH ME 2D
  TH2 *DPhiDEtaEHMix_1Clone[6];
  TH2 *DPhiDEtaEHMix_1Clone2[6];

  TCanvas *cNS = new TCanvas("DphiEtaEh_ME_PbPb_NoScale", "Dphi Deta mixed event (no scaling), IncE-H",50,50,1000,800);
  cNS->Divide(3,2);

    for(int i=0; i<5; i++){
    //Clone DPhiDEtaEHMix_1
        DPhiDEtaEHMix_1Clone[i] = (TH2 *)DPhiDEtaEHMix_1[i]->Clone();//Not scaled
        DPhiDEtaEHMix_1Clone2[i] = (TH2 *)DPhiDEtaEHMix_1[i]->Clone();//Scaled
    }
    for(int i=0; i<6; i++){
        //E-H ME Scale Factor 
                
        cNS->cd(i+1);

        if(i<5){
        DPhiDEtaEHMix_1Clone[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        DPhiDEtaEHMix_1Clone[i]->GetYaxis()->SetTitle("#Delta#eta");
        DPhiDEtaEHMix_1Clone[i]->Sumw2();

        DPhiDEtaEHMix_1Clone[i]->Draw("SURF1");
        }

        if (i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaEHMix_1Clone[0],"#Delta#eta#Delta#phi InclusiveEH Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Mixed Event");
            t1->AddText("4 < p_{T}^{trig e} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }

    TCanvas *cS = new TCanvas("DphiEtaEh_ME_PbPb_Scaled", "Dphi Deta mixed event (scaled), IncE-H",50,50,1000,800);
    cS->Divide(3,2);
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
        double ScaleFac = GetScaleFactor(DPhiDEtaEHMix_1[i]);
        DPhiDEtaEHMix_1[i]->Scale(1/ScaleFac);

        double ScaleFac2 = GetScaleFactor(DPhiDEtaEHMix_1Clone2[i]);
        DPhiDEtaEHMix_1Clone2[i]->Scale(1/ScaleFac2);
    }
    for(int i=0; i<6; i++){//Changing i<5 to i<6 causes a seg fault and I don't know whyyyyy
        
        cS->cd(i+1);

        if(i<5){
        DPhiDEtaEHMix_1Clone2[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        DPhiDEtaEHMix_1Clone2[i]->GetYaxis()->SetTitle("#Delta#eta");
        DPhiDEtaEHMix_1Clone2[i]->Sumw2();
        
        DPhiDEtaEHMix_1Clone2[i]->Draw("SURF1");
        }

        if (i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaEHMix_1Clone2[0],"#Delta#eta#Delta#phi InclusiveE-H Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Scaled Mixed Event");
            t1->AddText("4 < p_{T}^{trig e} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }

    //Divide Same event by mixed event 2D (HH) with scaling
    //Plot SE/ME 2D
    TCanvas *c222 = new TCanvas("DphiEtaEh_SE/ME_PbPb", "Dphi Deta same event/mixed event (scaled), IncE-H",50,50,1000,800);
    c222->Divide(3,2);
    TH2* DPhiDEtaEHDiv_1[6];
    for(int i=0; i<6; i++){
         c222->cd(i+1);
        if(i<5) {
         DPhiDEtaEHDiv_1[i] = (TH2 *)DPhiDEtaEH_1[i]->Clone();
         DPhiDEtaEHDiv_1[i]->Reset();
         DPhiDEtaEHDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaEHDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaEHDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaEHDiv_1[i]->Divide(DPhiDEtaEH_1[i],DPhiDEtaEHMix_1[i],1,1);

         DPhiDEtaEHDiv_1[i]->Draw("SURF1");
         }
        
        if(i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaEHDiv_1[0],"InclusiveE-H Same Event/Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Same Event/(Scaled)Mixed Event");
            t1->AddText("4 < p_{T}^{trig e} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }
         
    //Plot EH 2D SE
TCanvas *c23 = new TCanvas("DphiEtaEh_SE_PbPb", "Dphi Deta same event, IncE-H",50,50,1000,800);
    c23->Divide(3,2);
    for(int i=0; i<6; i++){
        c23->cd(i+1);
        if(i<5){
            DPhiDEtaEH_1[i]->Sumw2();
            ProcessHisto(DPhiDEtaEH_1[i],1.2,1);
            DPhiDEtaEH_1[i]->Draw("SURF1");
        }

        if(i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaEH_1[0],"InclusiveE-H Same Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Same Event");
            t1->AddText("4 < p_{T}^{trig e} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
        
    }
   
    //Generate (Projection X) 1D SE/ME Histogram
    TH1* DPhiEHDiv1[6];
    TH1* DPhiEHProj[6];
    for(int i=0; i<5; i++){
    Int_t bin1 = DPhiDEtaEHDiv_1[i]->GetYaxis()->FindBin(-1);
    Int_t bin2 = DPhiDEtaEHDiv_1[i]->GetYaxis()->FindBin(1);
    DPhiEHDiv1[i] = DPhiDEtaEHDiv_1[i]->ProjectionX(Form("IncDPhiCor%i",i),bin1,bin2);
    DPhiEHProj[i] = DPhiDEtaEH_1[i]->ProjectionX(Form("IncDPhiSE%i",i),bin1,bin2);
    }

    //Plot DPhi before and after ME correction
    TCanvas *DphiChk = new TCanvas("MECorrDPhiIncl_pPb","Dphi before and after ME correction, IncE",50,50,1000,800);
    DphiChk->Divide(3,2);
    for(int i=0;i<6;i++){
        DphiChk->cd(i+1);
        if(i<5){
        DPhiEHDiv1[i]->Sumw2();
        ProcessHisto(DPhiEHDiv1[i],1.2,1);
        DPhiEHProj[i]->Sumw2();
        ProcessHisto(DPhiEHProj[i],1.2,2);
        DPhiEHDiv1[i]->GetYaxis()->SetTitle("1/N dN");
        DPhiEHProj[i]->GetYaxis()->SetTitle("1/N dN");

      
      DPhiEHDiv1[i]->Scale(1/DPhiEHDiv1[i]->Integral());
      DPhiEHProj[i]->Scale(1/DPhiEHProj[i]->Integral());
        
      DPhiEHDiv1[i]->GetYaxis()->SetRangeUser(0.8*DPhiEHDiv1[i]->GetBinContent(DPhiEHDiv1[i]->GetMinimumBin()),1.1*DPhiEHDiv1[i]->GetBinContent(DPhiEHDiv1[i]->GetMaximumBin()));
      DPhiEHProj[i]->GetYaxis()->SetRangeUser(0.8*DPhiEHProj[i]->GetBinContent(DPhiEHProj[i]->GetMinimumBin()),1.1*DPhiEHProj[i]->GetBinContent(DPhiEHProj[i]->GetMaximumBin()));

      DPhiEHProj[i]->Draw();
      DPhiEHDiv1[i]->Draw("same");

    
    }

    if(i==5){
      TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
      ProcessLegend(leg);
      leg->SetTextSize(0.05);
      leg->AddEntry(DPhiEHProj[0],"#Delta#varphi InclE-H Before ME correct","pl");
      leg->AddEntry(DPhiEHDiv1[0],"#Delta#varphi InclE-H After ME correct","pl");
      leg->Draw();
      TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
      t1->AddText("4 < p_{T}^{trig e} < 12 GeV/c");
      t1->AddText("Centrality trigger");
      t1->AddText("0 - 10%");
      ProcessTPave(t1);
    }
  }

    //%%%%%%%%%%%
    //H-H Correlations
    //%%%%%%%%%%%

    //Scale HH ME
  //Plot HH ME 2D
  TH2 *DPhiDEtaHHMix_1Clone[6];
  TH2 *DPhiDEtaHHMix_1Clone2[6];

  TCanvas *cNSHH = new TCanvas("DphiEtaHh_ME_PbPb_NoScale", "Dphi Deta mixed event (no scaling), H-H",50,50,1000,800);
  cNSHH->Divide(3,2);

    for(int i=0; i<5; i++){
    //Clone DPhiDEtaHHMix_1
        DPhiDEtaHHMix_1Clone[i] = (TH2 *)DPhiDEtaHHMix_1[i]->Clone();//Not scaled
        DPhiDEtaHHMix_1Clone2[i] = (TH2 *)DPhiDEtaHHMix_1[i]->Clone();//Scaled
    }
    for(int i=0; i<6; i++){
        //E-H ME Scale Factor 
                
        cNSHH->cd(i+1);

        if(i<5){
        DPhiDEtaHHMix_1Clone[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        DPhiDEtaHHMix_1Clone[i]->GetYaxis()->SetTitle("#Delta#eta");
        DPhiDEtaHHMix_1Clone[i]->Sumw2();

        DPhiDEtaHHMix_1Clone[i]->Draw("SURF1");
        }

        if (i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaHHMix_1Clone[0],"#Delta#eta#Delta#phi H-H Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Mixed Event");
            t1->AddText("4 < p_{T}^{trig h} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }

    TCanvas *cSHH = new TCanvas("DphiEtaHh_ME_PbPb_Scaled", "Dphi Deta mixed event (scaled), H-H",50,50,1000,800);
    cSHH->Divide(3,2);
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
        double ScaleFac = GetScaleFactor(DPhiDEtaHHMix_1[i]);
        DPhiDEtaHHMix_1[i]->Scale(1/ScaleFac);

        double ScaleFac2 = GetScaleFactor(DPhiDEtaHHMix_1Clone2[i]);
        DPhiDEtaHHMix_1Clone2[i]->Scale(1/ScaleFac2);
    }
    for(int i=0; i<6; i++){//Changing i<5 to i<6 causes a seg fault and I don't know whyyyyy
        
        cSHH->cd(i+1);

        if(i<5){
        DPhiDEtaHHMix_1Clone2[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        DPhiDEtaHHMix_1Clone2[i]->GetYaxis()->SetTitle("#Delta#eta");
        DPhiDEtaHHMix_1Clone2[i]->Sumw2();
        
        DPhiDEtaHHMix_1Clone2[i]->Draw("SURF1");
        }

        if (i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaHHMix_1Clone2[0],"#Delta#eta#Delta#phi H-H Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Scaled Mixed Event");
            t1->AddText("4 < p_{T}^{trig h} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }

    //Divide Same event by mixed event 2D (HH) with scaling
    //Plot SE/ME 2D
    TCanvas *c222HH = new TCanvas("DphiEtaHh_SE/ME_PbPb", "Dphi Deta same event/mixed event (scaled), H-H",50,50,1000,800);
    c222HH->Divide(3,2);
    TH2* DPhiDEtaHHDiv_1[6];
    for(int i=0; i<6; i++){
         c222HH->cd(i+1);
        if(i<5) {
         DPhiDEtaHHDiv_1[i] = (TH2 *)DPhiDEtaHH_1[i]->Clone();
         DPhiDEtaHHDiv_1[i]->Reset();
         DPhiDEtaHHDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaHHDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaHHDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaHHDiv_1[i]->Divide(DPhiDEtaHH_1[i],DPhiDEtaHHMix_1[i],1,1);

         DPhiDEtaHHDiv_1[i]->Draw("SURF1");
         }
        
        if(i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaHHDiv_1[0],"H-H Same Event/Mixed Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Same Event/(Scaled)Mixed Event");
            t1->AddText("4 < p_{T}^{trig h} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
    }
         
    //Plot EH 2D SE
TCanvas *c23HH = new TCanvas("DphiEtaHh_SE_PbPb", "Dphi Deta same event, H-H",50,50,1000,800);
    c23HH->Divide(3,2);
    for(int i=0; i<6; i++){
        c23HH->cd(i+1);
        if(i<5){
            DPhiDEtaHH_1[i]->Sumw2();
            ProcessHisto(DPhiDEtaHH_1[i],1.2,1);
            DPhiDEtaHH_1[i]->Draw("SURF1");
        }

        if(i==5){
            TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
            ProcessLegend(leg);
            leg->SetTextSize(0.05);
            leg->AddEntry(DPhiDEtaHH_1[0],"H-H Same Event","pl");
            leg->Draw();
            TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
            t1->AddText("Same Event");
            t1->AddText("4 < p_{T}^{trig h} < 12 GeV/c");
            t1->AddText("Centrality trigger");
            t1->AddText("0 - 10%");
            ProcessTPave(t1);
        }
        
    }
   
    //Generate (Projection X) 1D SE/ME Histogram
    TH1* DPhiHHDiv1[6];
    TH1* DPhiHHProj[6];
    for(int i=0; i<5; i++){
    Int_t bin1 = DPhiDEtaHHDiv_1[i]->GetYaxis()->FindBin(-1);
    Int_t bin2 = DPhiDEtaHHDiv_1[i]->GetYaxis()->FindBin(1);
    DPhiHHDiv1[i] = DPhiDEtaHHDiv_1[i]->ProjectionX(Form("HHDPhiCor%i",i),bin1,bin2);
    DPhiHHProj[i] = DPhiDEtaHH_1[i]->ProjectionX(Form("HHDPhiSE%i",i),bin1,bin2);
    }

    //Plot DPhi before and after ME correction
    TCanvas *DphiChkHH = new TCanvas("MECorrDPhiH_pPb","Dphi before and after ME correction, H-H",50,50,1000,800);
    DphiChkHH->Divide(3,2);
    for(int i=0;i<6;i++){
        DphiChkHH->cd(i+1);
        if(i<5){
        DPhiHHDiv1[i]->Sumw2();
        ProcessHisto(DPhiHHDiv1[i],1.2,1);
        DPhiHHProj[i]->Sumw2();
        ProcessHisto(DPhiHHProj[i],1.2,2);
        DPhiHHDiv1[i]->GetYaxis()->SetTitle("1/N dN");
        DPhiHHProj[i]->GetYaxis()->SetTitle("1/N dN");

      
      DPhiHHDiv1[i]->Scale(1/DPhiHHDiv1[i]->Integral());
      DPhiHHProj[i]->Scale(1/DPhiHHProj[i]->Integral());
        
      DPhiHHDiv1[i]->GetYaxis()->SetRangeUser(0.8*DPhiHHDiv1[i]->GetBinContent(DPhiHHDiv1[i]->GetMinimumBin()),1.1*DPhiHHDiv1[i]->GetBinContent(DPhiHHDiv1[i]->GetMaximumBin()));
      DPhiHHProj[i]->GetYaxis()->SetRangeUser(0.8*DPhiHHProj[i]->GetBinContent(DPhiHHProj[i]->GetMinimumBin()),1.1*DPhiHHProj[i]->GetBinContent(DPhiHHProj[i]->GetMaximumBin()));

      DPhiHHProj[i]->Draw();
      DPhiHHDiv1[i]->Draw("same");

    
    }

    if(i==5){
      TLegend *leg = new TLegend(0.259,0.235,0.509,0.377);
      ProcessLegend(leg);
      leg->SetTextSize(0.05);
      leg->AddEntry(DPhiHHProj[0],"#Delta#varphi H-H Before ME correct","pl");
      leg->AddEntry(DPhiHHDiv1[0],"#Delta#varphi H-H After ME correct","pl");
      leg->Draw();
      TPaveText *t1 = new TPaveText(0.371,0.431,0.621,0.693);
      t1->AddText("4 < p_{T}^{trig h} < 12 GeV/c");
      t1->AddText("Centrality trigger");
      t1->AddText("0 - 10%");
      ProcessTPave(t1);
    }
  }

    //%%%%%%%%
    //Obtain DPhiDeta with ME correction of remaining histograms
    //%%%%%%%%
    //ULSE

    THnSparse *SprsULSEHCorrl = (THnSparse *) EH_020->FindObject("fSprsULSEHCorrl");//Getting EMCal electron THnSparse object (ULSE)
    THnSparse *SprsMixULSEHCorrl = (THnSparse *) EH_020->FindObject("fSprsMixULSEHCorrl");//Getting EMCal electron THnSparse object (Mixed event; inclusive ele-h)

    SprsULSEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta
    SprsMixULSEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta

    //Obtain EH Correlation Histograms

    //ULS 2D histograms
    TH2 *DPhiDEtaULS_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaULS_1[i] = Delta_2DSelectEH(SprsULSEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaULS_1[i]->SetName(Form("Delta Phi Delta Eta ULS SE_%i",i));
        DPhiDEtaULS_1[i]->SetTitle(PtBinNames[i]);
    }
    //ULS Mixed event
    TH2 *DPhiDEtaULSMix_1[6];//2D Delta phi delta eta histograms mixed event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaULSMix_1[i] = Delta_2DSelectEH(SprsMixULSEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaULSMix_1[i]->SetName(Form("Delta Phi Delta Eta ULS ME_%i",i));
        DPhiDEtaULSMix_1[i]->SetTitle(PtBinNames[i]);
    }

    //Scale ME then SE/ME
    TH2 *DPhiDEtaULSDiv_1[6];
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
        double ScaleFac = GetScaleFactor(DPhiDEtaULSMix_1[i]);
        DPhiDEtaULSMix_1[i]->Scale(1/ScaleFac);

         DPhiDEtaULSDiv_1[i] = (TH2 *)DPhiDEtaULS_1[i]->Clone();
         DPhiDEtaULSDiv_1[i]->Reset();
         DPhiDEtaULSDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaULSDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaULSDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaULSDiv_1[i]->Divide(DPhiDEtaULS_1[i],DPhiDEtaULSMix_1[i],1,1);
   
    }

    //LSE

    THnSparse *SprsLSEHCorrl = (THnSparse *) EH_020->FindObject("fSprsLSEHCorrl");//Getting EMCal electron THnSparse object (ULSE)
    THnSparse *SprsMixLSEHCorrl = (THnSparse *) EH_020->FindObject("fSprsMixLSEHCorrl");//Getting EMCal electron THnSparse object (Mixed event; inclusive ele-h)

    SprsLSEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta
    SprsMixLSEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta

    //Obtain EH Correlation Histograms

    //LS 2D histograms
    TH2 *DPhiDEtaLS_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaLS_1[i] = Delta_2DSelectEH(SprsLSEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaLS_1[i]->SetName(Form("Delta Phi Delta Eta LS SE_%i",i));
        DPhiDEtaLS_1[i]->SetTitle(PtBinNames[i]);
    }
    //LS Mixed event
    TH2 *DPhiDEtaLSMix_1[6];//2D Delta phi delta eta histograms mixed event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaLSMix_1[i] = Delta_2DSelectEH(SprsMixLSEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaLSMix_1[i]->SetName(Form("Delta Phi Delta Eta LS ME_%i",i));
        DPhiDEtaLSMix_1[i]->SetTitle(PtBinNames[i]);
    }

    //Scale ME then SE/ME
    TH2 *DPhiDEtaLSDiv_1[6];
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
        double ScaleFac = GetScaleFactor(DPhiDEtaLSMix_1[i]);
        DPhiDEtaLSMix_1[i]->Scale(1/ScaleFac);

         DPhiDEtaLSDiv_1[i] = (TH2 *)DPhiDEtaULS_1[i]->Clone();
         DPhiDEtaLSDiv_1[i]->Reset();
         DPhiDEtaLSDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaLSDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaLSDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaLSDiv_1[i]->Divide(DPhiDEtaLS_1[i],DPhiDEtaLSMix_1[i],1,1);
   
    }
    //ULS NoPE
    
    THnSparse *SprsULSNoPartnerEHCorrl = (THnSparse *) EH_020->FindObject("fSprsULSNoPartnerEHCorrl");//Getting EMCal electron THnSparse object (ULSE)

    SprsULSNoPartnerEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta

    //Obtain EH Correlation Histograms

    //ULS 2D histograms
    TH2 *DPhiDEtaULSNoPE_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaULSNoPE_1[i] = Delta_2DSelectEH(SprsULSNoPartnerEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaULSNoPE_1[i]->SetName(Form("Delta Phi Delta Eta ULS NoPE SE_%i",i));
        DPhiDEtaULSNoPE_1[i]->SetTitle(PtBinNames[i]);
    }

    //Scale ME then SE/ME
    TH2 *DPhiDEtaULSNoPEDiv_1[6];
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
         DPhiDEtaULSNoPEDiv_1[i] = (TH2 *)DPhiDEtaULSNoPE_1[i]->Clone();
         DPhiDEtaULSNoPEDiv_1[i]->Reset();
         DPhiDEtaULSNoPEDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaULSNoPEDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaULSNoPEDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaULSNoPEDiv_1[i]->Divide(DPhiDEtaULSNoPE_1[i],DPhiDEtaULSMix_1[i],1,1);
   
    }
    //LS NoPE

    THnSparse *SprsLSNoPartnerEHCorrl = (THnSparse *) EH_020->FindObject("fSprsLSNoPartnerEHCorrl");//Getting EMCal electron THnSparse object (LSE)

    SprsLSNoPartnerEHCorrl->GetAxis(3)->SetRangeUser(-1, 1); //Setting delta eta

    //Obtain EH Correlation Histograms

    //ULS 2D histograms
    TH2 *DPhiDEtaLSNoPE_1[6];//2D Delta phi delta eta histograms same event (E-H)
    for(int i=0; i<5; i++){ //Fill histogram 
        DPhiDEtaLSNoPE_1[i] = Delta_2DSelectEH(SprsLSNoPartnerEHCorrl, 2, 3, pTbinTrigger1[0]+0.01, pTbinTrigger1[1]-0.05, ptbinAsso1[i]+0.01, ptbinAsso1[i+1]-0.05); //E/p  Projects the cut of E/p
        DPhiDEtaLSNoPE_1[i]->SetName(Form("Delta Phi Delta Eta LS NoPE SE_%i",i));
        DPhiDEtaLSNoPE_1[i]->SetTitle(PtBinNames[i]);
    }

    //Scale ME then SE/ME
    TH2 *DPhiDEtaLSNoPEDiv_1[6];
    for(int i=0; i<5; i++){ //Apply 4 bin scale factor to HH ME
         DPhiDEtaLSNoPEDiv_1[i] = (TH2 *)DPhiDEtaULSNoPE_1[i]->Clone();
         DPhiDEtaLSNoPEDiv_1[i]->Reset();
         DPhiDEtaLSNoPEDiv_1[i]->SetTitle(PtBinNames[i]);
         DPhiDEtaLSNoPEDiv_1[i]->SetXTitle("#Delta #varphi (rad)");
         DPhiDEtaLSNoPEDiv_1[i]->SetYTitle("#Delta#eta");
         DPhiDEtaLSNoPEDiv_1[i]->Divide(DPhiDEtaLSNoPE_1[i],DPhiDEtaLSMix_1[i],1,1);
   
    }
    //Output root file

    cout << "Writing output root file" << endl;
    TFile fout("MECorrections_Cent020.root","RECREATE");
    Centrality->Write("NEvents_FromCent");
    for(int i=0; i<5; i++){
        DPhiDEtaEHDiv_1[i]->Write(Form("DPhiDEtaIncEH_Div_%i",i));
        DPhiDEtaHHDiv_1[i]->Write(Form("DPhiDEtaHH_Div_%i",i));
        DPhiDEtaULSDiv_1[i]->Write(Form("DPhiDEtaULS_Div_%i",i));
        DPhiDEtaLSDiv_1[i]->Write(Form("DPhiDEtaLS_Div_%i",i));
        DPhiDEtaULSNoPEDiv_1[i]->Write(Form("DPhiDEtaULSNoPE_Div_%i",i)); 
        DPhiDEtaLSNoPEDiv_1[i]->Write(Form("DPhiDEtaLSNoPE_Div_%i",i));
    }
    fout.Close();


}//End of Macro
