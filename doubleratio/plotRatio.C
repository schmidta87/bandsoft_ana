void plotRatio(TString inFileDatName , TString inFileSimBacName , char* qualifier){
	
	// Define some functions
	void PrettyTH1D(TH1D*hist,TString xtit);


	TFile * inFileDat = new TFile(inFileDatName);
	TFile * inFileSimBac = new TFile(inFileSimBacName);

	// Get X' Distributions
	TH1D * hXp_full_dat = (TH1D*) inFileDat->Get("hXp");
	hXp_full_dat->Rebin(4);
	PrettyTH1D(hXp_full_dat,"x'");
	TH1D * hXp_full_sim = (TH1D*) inFileSimBac->Get("hXp_sim");
	hXp_full_sim->Rebin(4);
	TH1D * hXp_full_bac = (TH1D*) inFileSimBac->Get("hXp_bac");
	hXp_full_bac->Rebin(4);
	// Get X' vs Al Distributions
	TH1D * hXp_dat_0 = (TH1D*) inFileDat->Get("hXpBins_0");
	TH1D * hXp_dat_1 = (TH1D*) inFileDat->Get("hXpBins_1");
	TH1D * hXp_dat_2 = (TH1D*) inFileDat->Get("hXpBins_2");
	TH1D * hXp_dat_3 = (TH1D*) inFileDat->Get("hXpBins_3");
	TH1D * hXp_dat_4 = (TH1D*) inFileDat->Get("hXpBins_4");
	TH1D * hXp_dat_5 = (TH1D*) inFileDat->Get("hXpBins_5");
	TH1D * hXp_dat_6 = (TH1D*) inFileDat->Get("hXpBins_6");
	TH1D * hXp_dat_7 = (TH1D*) inFileDat->Get("hXpBins_7");
	PrettyTH1D(hXp_dat_0,"x'");
	PrettyTH1D(hXp_dat_1,"x'");
	PrettyTH1D(hXp_dat_2,"x'");
	PrettyTH1D(hXp_dat_3,"x'");
	PrettyTH1D(hXp_dat_4,"x'");
	PrettyTH1D(hXp_dat_5,"x'");
	PrettyTH1D(hXp_dat_6,"x'");
	PrettyTH1D(hXp_dat_7,"x'");
	TH1D * hXp_bac_0 = (TH1D*) inFileSimBac->Get("hXpBins_0_bac");
	TH1D * hXp_bac_1 = (TH1D*) inFileSimBac->Get("hXpBins_1_bac");
	TH1D * hXp_bac_2 = (TH1D*) inFileSimBac->Get("hXpBins_2_bac");
	TH1D * hXp_bac_3 = (TH1D*) inFileSimBac->Get("hXpBins_3_bac");
	TH1D * hXp_bac_4 = (TH1D*) inFileSimBac->Get("hXpBins_4_bac");
	TH1D * hXp_bac_5 = (TH1D*) inFileSimBac->Get("hXpBins_5_bac");
	TH1D * hXp_bac_6 = (TH1D*) inFileSimBac->Get("hXpBins_6_bac");
	TH1D * hXp_bac_7 = (TH1D*) inFileSimBac->Get("hXpBins_7_bac");
	TH1D * hXp_sim_0 = (TH1D*) inFileSimBac->Get("hXpBins_0_sim");
	TH1D * hXp_sim_1 = (TH1D*) inFileSimBac->Get("hXpBins_1_sim");
	TH1D * hXp_sim_2 = (TH1D*) inFileSimBac->Get("hXpBins_2_sim");
	TH1D * hXp_sim_3 = (TH1D*) inFileSimBac->Get("hXpBins_3_sim");
	TH1D * hXp_sim_4 = (TH1D*) inFileSimBac->Get("hXpBins_4_sim");
	TH1D * hXp_sim_5 = (TH1D*) inFileSimBac->Get("hXpBins_5_sim");
	TH1D * hXp_sim_6 = (TH1D*) inFileSimBac->Get("hXpBins_6_sim");
	TH1D * hXp_sim_7 = (TH1D*) inFileSimBac->Get("hXpBins_7_sim");
	// Get the As distributions for Low and Hi x'
	TH1D * hAs_lo_dat = (TH1D*) inFileDat->Get("hAs_lo");
	PrettyTH1D(hAs_lo_dat,"Alpha_s");
	TH1D * hAs_lo_bac = (TH1D*) inFileSimBac->Get("hAs_lo_bac");
	TH1D * hAs_lo_sim = (TH1D*) inFileSimBac->Get("hAs_lo_sim");
	TH1D * hAs_hi_dat = (TH1D*) inFileDat->Get("hAs_hi");
	PrettyTH1D(hAs_hi_dat,"Alpha_s");
	TH1D * hAs_hi_bac = (TH1D*) inFileSimBac->Get("hAs_hi_bac");
	TH1D * hAs_hi_sim = (TH1D*) inFileSimBac->Get("hAs_hi_sim");
	
	// X' Distribution Data to Background
	TCanvas * c_xp = new TCanvas("c_xp","X' Distribution");
	c_xp->cd();
	hXp_full_dat->Draw();
	hXp_full_bac->Draw("h,same");
	c_xp->Modified();
	c_xp->Update();
	c_xp->SaveAs(Form("xp_hist_%s.pdf",qualifier));

	// X' Distribution to Simulation
	TCanvas * c_xp_sim = new TCanvas("c_xp_sim","X' Distribution to Sim");
	c_xp_sim->cd();
	hXp_full_dat->Add( hXp_full_bac , -1 );
	hXp_full_dat->Draw();
	double sim_scaling = hXp_full_dat->Integral() / hXp_full_sim->Integral();
	hXp_full_sim->Scale( sim_scaling );
	hXp_full_sim->Draw("h,same");
	c_xp_sim->Modified();
	c_xp_sim->Update();
	c_xp_sim->SaveAs(Form("xp_sim_hist_%s.pdf",qualifier));


	// X' Distributions in Bins of Alpha S to Background
	TCanvas * c_xp_al_1 = new TCanvas("c_xp_al_1","X' Distributions in Bins of Alpha S",1200,900);
	c_xp_al_1->Divide(2,2);
	c_xp_al_1->cd(1);
  	hXp_dat_0->Draw();
	hXp_bac_0->Draw("h,same");
  	c_xp_al_1->cd(2);
  	hXp_dat_1->Draw();
	hXp_bac_1->Draw("h,same");
  	c_xp_al_1->cd(3);
  	hXp_dat_2->Draw();
	hXp_bac_2->Draw("h,same");
  	c_xp_al_1->cd(4);
  	hXp_dat_3->Draw();
	hXp_bac_3->Draw("h,same");
	c_xp_al_1->SaveAs(Form("xp_al_bins_1_%s.pdf",qualifier));
	TCanvas * c_xp_al_2 = new TCanvas("c_xp_al_2","X' Distributions in Bins of Alpha S",1200,900);
	c_xp_al_2->Divide(2,2);
	c_xp_al_2->cd(1);
  	hXp_dat_4->Draw();
	hXp_bac_4->Draw("h,same");
  	c_xp_al_2->cd(2);
  	hXp_dat_5->Draw();
	hXp_bac_5->Draw("h,same");
  	c_xp_al_2->cd(3);
  	hXp_dat_6->Draw();
	hXp_bac_6->Draw("h,same");
  	c_xp_al_2->cd(4);
  	hXp_dat_7->Draw();
	hXp_bac_7->Draw("h,same");
	c_xp_al_2->SaveAs(Form("xp_al_bins_2_%s.pdf",qualifier));

	hXp_dat_0->Add( hXp_bac_0 , -1 );
	hXp_dat_1->Add( hXp_bac_1 , -1 );
	hXp_dat_2->Add( hXp_bac_2 , -1 );
	hXp_dat_3->Add( hXp_bac_3 , -1 );
	hXp_dat_4->Add( hXp_bac_4 , -1 );
	hXp_dat_5->Add( hXp_bac_5 , -1 );
	hXp_dat_6->Add( hXp_bac_6 , -1 );
	hXp_dat_7->Add( hXp_bac_7 , -1 );
  	hXp_sim_0->Scale( sim_scaling );
  	hXp_sim_1->Scale( sim_scaling );
  	hXp_sim_2->Scale( sim_scaling );
  	hXp_sim_3->Scale( sim_scaling );
  	hXp_sim_4->Scale( sim_scaling );
  	hXp_sim_5->Scale( sim_scaling );
  	hXp_sim_6->Scale( sim_scaling );
  	hXp_sim_7->Scale( sim_scaling );

	TCanvas * c_xp_al_3 = new TCanvas("c_xp_al_3","X' Distributions in Bins of Alpha S",1200,900);
	c_xp_al_3->Divide(2,2);
	c_xp_al_3->cd(1);
  	hXp_dat_0->Draw();
	hXp_sim_0->Draw("h,same");
  	c_xp_al_3->cd(2);
  	hXp_dat_1->Draw();
	hXp_sim_1->Draw("h,same");
  	c_xp_al_3->cd(3);
  	hXp_dat_2->Draw();
	hXp_sim_2->Draw("h,same");
  	c_xp_al_3->cd(4);
  	hXp_dat_3->Draw();
	hXp_sim_3->Draw("h,same");
	c_xp_al_3->SaveAs(Form("xp_al_bins_3_%s.pdf",qualifier));
	TCanvas * c_xp_al_4 = new TCanvas("c_xp_al_4","X' Distributions in Bins of Alpha S",1200,900);
	c_xp_al_4->Divide(2,2);
	c_xp_al_4->cd(1);
  	hXp_dat_4->Draw();
	hXp_sim_4->Draw("h,same");
  	c_xp_al_4->cd(2);
  	hXp_dat_5->Draw();
	hXp_sim_5->Draw("h,same");
  	c_xp_al_4->cd(3);
  	hXp_dat_6->Draw();
	hXp_sim_6->Draw("h,same");
  	c_xp_al_4->cd(4);
  	hXp_dat_7->Draw();
	hXp_sim_7->Draw("h,same");
	c_xp_al_4->SaveAs(Form("xp_al_bins_4_%s.pdf",qualifier));

	// Plot As distributions for low and high x'	
  	TCanvas * c_as = new TCanvas("c_as","Alpha S Distribution",900,900);
	c_as->Divide(1,2);
	c_as->cd(1);
	hAs_lo_dat->Draw();
	hAs_lo_bac->Draw("h,same");
	c_as->cd(2);
	hAs_hi_dat->Draw();
	hAs_hi_bac->Draw("h,same");
	c_as->SaveAs(Form("as_dist_%s.pdf",qualifier));

	hAs_lo_dat->Add( hAs_lo_bac , -1 );
	hAs_hi_dat->Add( hAs_hi_bac , -1 );
	hAs_lo_sim->Scale( sim_scaling );
	hAs_hi_sim->Scale( sim_scaling );
  	TCanvas * c_as_sim = new TCanvas("c_as_sim","Alpha S Distribution",900,900);
	c_as_sim->Divide(1,2);
	c_as_sim->cd(1);
	hAs_lo_sim->Draw("h");
	hAs_lo_dat->Draw("same");
	c_as_sim->cd(2);
	hAs_hi_dat->SetMinimum(0);
	hAs_hi_dat->Draw();
	hAs_hi_sim->Draw("h,same");
	c_as_sim->SaveAs(Form("as_dist_sim_%s.pdf",qualifier));

	// Now plot the ratios:
	TCanvas * c_ratios = new TCanvas("c_ratios","Ratios",900,900);
	c_ratios->Divide(1,2);
	c_ratios->cd(1);
	hAs_lo_dat->Divide( hAs_lo_sim );
	if( hAs_lo_dat->GetMaximum() > 20 ) hAs_lo_dat->SetMaximum(10);
	hAs_lo_dat->Draw();
	c_ratios->cd(2);
	hAs_hi_dat->Divide( hAs_hi_sim );
	if( hAs_hi_dat->GetMaximum() > 20 ) hAs_hi_dat->SetMaximum(10);
	hAs_hi_dat->Draw();
	c_ratios->SaveAs(Form("singleratios_%s.pdf",qualifier));

	TCanvas * c_doublerat = new TCanvas("c_doublerat","Double Ratio");
	c_doublerat->cd();
	hAs_hi_dat->Divide( hAs_lo_dat );
	hAs_hi_dat->SetMinimum(0);
	if( hAs_hi_dat->GetMaximum() > 20 ) hAs_hi_dat->SetMaximum(10);
	//if( hAs_hi_dat->GetBinContent(5) < 10 ) hAs_hi_dat->SetMaximum(2);
	hAs_hi_dat->Draw();
	TLine * line = new TLine(1.2,1,1.6,1);
	line->SetLineWidth(3);
	line->Draw("same");
	c_doublerat->SaveAs(Form("doubleratio_%s.pdf",qualifier));

	TFile * outFile = new TFile(Form("out_doubleratios_%s.root",qualifier),"RECREATE");
	outFile->cd();
	hAs_hi_dat->Write();
	outFile->Close();
  
	return;
}

void PrettyTH1D( TH1D * hist , TString xtit ){
	hist->SetStats(0);
	hist->SetLineColor(4);
	hist->SetLineWidth(2);
	hist->Scale(1);
	hist->SetTitle("");
	hist->SetMinimum(0);
	hist->GetXaxis()->SetTitle(xtit);


	return;
}
