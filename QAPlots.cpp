//
// QAPlots.cpp
// Opens a root file from the electron QA task and extracts histograms.
//
// by Erin Gauger
// Oct. 29th, 2019
//

void QAPlots(){
    
    void makePretty(TH1 *histo, Int_t style, Int_t color);
    gStyle->SetOptStat(0); //gets rid of statistics box when drawing
    
    ///////////////////////////
    //Open existing root file//
    ///////////////////////////
    TFile *file = new TFile("QATrain891_Runlist3_AnalysisResults.root","READ");
    /*Small note on file options
     NEW or CREATE: Create a new file and open it for writing, if the file already exists the file is not opened
     RECREATE: Create a new file, if the file already exists it will be overwritten
     UPDATE: Open an existing file for writing. If no file exists, it is created.
     READ: Open an existing file for reading (default).*/
    
    /////////////////////////////
    //Get objects from the file//
    /////////////////////////////
    TDirectoryFile *int7Dir = (TDirectoryFile*)file->Get("PWGHF_hfeHFEemcQAINT7GblTrk_EMC"); //subfile with INT7 trigger results
    TList *list010 = (TList*)int7Dir->Get("HFEemcQAINT7_wTender_0_10GblTrk_EMC"); //list with 0-10% centrality results
    
    //Number of events
    TH1F *nEventsHisto = (TH1F*)list010->FindObject("fNevents");
    Double_t nEvents = nEventsHisto->GetBinContent(3);
    cout<<"Number of INT7 events in 0-10% central collisions = "<<nEvents<<endl;
    
    //Primary vertex position along z-axis
    TH1F *vtxZ = (TH1F*)list010->FindObject("fVtxZ");
    
    //Transverse momentum of all tracks
    TH1F *trkPt = (TH1F*)list010->FindObject("fTrkPt");
    makePretty(trkPt,20,kRed+1);
    
    //Eta and phi distribution of all tracks
    TH1F *trkEta = (TH1F*)list010->FindObject("fTrketa");
    TH1F *trkPhi = (TH1F*)list010->FindObject("fTrkphi");
    
    //TPC nSigma vs. p
    TH2F *nSigP = (TH2F*)list010->FindObject("fTPCnsig");

    //THnSparse with electron ID information
    //axes: 0=trigger info, 1=pT, 2=TPC nSigma, 3=E/p, 4=M20, 5=M02, 6=MC eID, 7=supermod no., 8=centrality
    THnSparse *elecSparse = (THnSparse*)list010->FindObject("Electron");
    //Note on sparses: sparses are like histograms with many axes. You can make a "cut" in one axis and project on another.
    //Later in this macro, I make a cut on TPC nSigma and show how the E/p distribution changes.
    
    //////////////////////
    //Project the sparse//
    //////////////////////
    
    //TPC nSigma vs. pT
    TH2F *nSigPt = (TH2F*)elecSparse->Projection(2,1,"");
    nSigPt->SetTitle("TPC nSigma distribution from Sparse");
    nSigPt->SetName("nSigPt");
    
    //E/p histo before cuts
    TH1D *eopAll = (TH1D*)elecSparse->Projection(3,"");
    eopAll->SetName("eopAll");
    makePretty(eopAll,20,kBlack);
    
    //Make -1<TPC nSigma<3 cut (for electrons)
    elecSparse->GetAxis(2)->SetRangeUser(-1.,3.);
    
    //E/p histo after electron cut
    TH1D *eopElec = (TH1D*)elecSparse->Projection(3,"");
    eopElec->SetName("eopElec");
    makePretty(eopElec,20,kBlue+1);
    
    //Make -10<TPC nSigma<-4 cut (for hadrons)
    elecSparse->GetAxis(2)->SetRangeUser(-10.,-4.);
    
    //E/p histo after hadron cut
    TH1D *eopHad = (TH1D*)elecSparse->Projection(3,"");
    eopHad->SetName("eopHad");
    makePretty(eopHad,20,kRed+1);
    
    ///////////////////
    //Draw histograms//
    ///////////////////
    TCanvas *canvas1 = new TCanvas("canvas1","Events",20,20,900,500);
    canvas1->Divide(2,1);
    TCanvas *canvas2 = new TCanvas("canvas2","Eta and Phi",40,40,900,500);
    canvas2->Divide(2,1);
    TCanvas *canvas3 = new TCanvas("canvas3","nSig vs. p/pT",60,60,900,500);
    canvas3->Divide(2,1);
    TCanvas *canvas5 = new TCanvas("canvas5","E/p",100,100,900,500);
    canvas5->Divide(2,1);
    
    canvas1->cd(1);
    nEventsHisto->Draw();
    canvas1->cd(2);
    vtxZ->Draw();
    
    canvas2->cd(1);
    trkEta->Draw();
    canvas2->cd(2);
    trkPhi->Draw();
    
    canvas3->cd(1);
    gPad->SetLogz(); //makes z-axis a log scale
    nSigP->GetXaxis()->SetRangeUser(2,30); //set axes range to be same as the other
    nSigP->GetYaxis()->SetRangeUser(-8,8);
    nSigP->Draw("colz");
    canvas3->cd(2);
    gPad->SetLogz(); //makes z-axis a log scale
    nSigPt->Draw("colz");
    
    canvas5->cd(1);
    eopAll->Draw("p");
    eopHad->Draw("psame");
    
    TLegend *leg5 = new TLegend(0.375,0.641,0.622,0.840);
    leg5->SetFillColor(kWhite);
    leg5->SetBorderSize(0);
    leg5->SetTextSize(0.05);
    leg5->SetTextFont(42);
    leg5->AddEntry(eopAll,"No cuts","pl");
    leg5->AddEntry(eopHad,"-10<nSig<-4","pl");
    leg5->AddEntry(eopElec,"-1<nSig<3","pl");
    leg5->Draw();
    
    canvas5->cd(2);
    eopElec->Draw("p");
    
}
///////////////
// FUNCTIONS //
///////////////
void makePretty(TH1 *histo, Int_t style, Int_t color){
    
    histo->SetMarkerStyle(style);
    histo->SetMarkerColor(color);
    histo->SetLineColor(color);
    
}

