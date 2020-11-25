void wp_xpas_plot(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1F* data, TH1F* sim, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileBac = new TFile(inBac);
	TFile * inFileSim = new TFile(inSim);

	// Get background normalization
	TVector3 * bacvector = (TVector3*)inFileDat->Get("bacnorm");
	double bacnorm = bacvector->X();
	TH1F * Wp_bac = (TH1F*)inFileBac->Get("Wp");
	bacnorm = bacnorm/Wp_bac->Integral() ;
	Wp_bac->Scale( bacnorm );

	// Get simulation normalization
	TH1F * Wp_sim = (TH1F*)inFileSim->Get("Wp");
	TH1F * Wp_dat = (TH1F*)inFileDat->Get("Wp");
	Wp_dat->Add( Wp_bac, -1 );
	double full_simnorm = Wp_dat->Integral() / Wp_sim->Integral();
	Wp_sim->Scale( full_simnorm );

	TCanvas * c0 = new TCanvas("c0","",800,600);
	c0->cd();
	Wp_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(Wp_dat,Wp_sim,"W' [GeV]","Counts");
	c0->SaveAs("full_Wp.pdf");

	
	// Grab intereseted histograms
	const int xb_bins = 4;
	const int xp_bins = 6;
	TH1F ** Wp_xb_bins_dat = new TH1F*[xb_bins];
	TH1F ** Wp_xp_bins_dat = new TH1F*[xp_bins];
	TH1F ** Wp_xb_bins_bac = new TH1F*[xb_bins];
	TH1F ** Wp_xp_bins_bac = new TH1F*[xp_bins];
	TH1F ** Wp_xb_bins_sim = new TH1F*[xb_bins];
	TH1F ** Wp_xp_bins_sim = new TH1F*[xp_bins];

	TCanvas * c1 = new TCanvas("c1","",800,600);
	c1->Divide(2,2);
	for( int i = 0 ; i < xb_bins ; i++ ){
		Wp_xb_bins_dat[i] = (TH1F*) inFileDat->Get(Form("Wp_aS_bin_%i",i));
		Wp_xb_bins_bac[i] = (TH1F*) inFileBac->Get(Form("Wp_aS_bin_%i",i));
		Wp_xb_bins_sim[i] = (TH1F*) inFileSim->Get(Form("Wp_aS_bin_%i",i));


		if( Wp_xb_bins_sim[i]->Integral() == 0 ) continue;
		
		Wp_xb_bins_bac[i] -> Scale( bacnorm );
		Wp_xb_bins_dat[i] -> Add( Wp_xb_bins_bac[i] , -1 );

		Wp_xb_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = Wp_xb_bins_dat[i]->Integral() / Wp_xb_bins_sim[i]->Integral();
		//pn_xb_bins_sim[i] -> Scale( simnorm );

		c1->cd(i+1);
		Wp_xb_bins_dat[i]->Rebin(2);
		Wp_xb_bins_sim[i]->Rebin(2);
		TString current_title = Wp_xb_bins_dat[i]->GetTitle();
		Wp_xb_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(Wp_xb_bins_dat[i],Wp_xb_bins_sim[i],"W' [GeV]","Counts");
	}
	c1->SaveAs("wp_as_bins.pdf");
	
	TCanvas * c2 = new TCanvas("c2","",800,600);
	c2->Divide(3,2);
	for( int i = 0 ; i < xp_bins ; i++ ){
		Wp_xp_bins_dat[i] = (TH1F*) inFileDat->Get(Form("Wp_xP_bin_%i",i));
		Wp_xp_bins_bac[i] = (TH1F*) inFileBac->Get(Form("Wp_xP_bin_%i",i));
		Wp_xp_bins_sim[i] = (TH1F*) inFileSim->Get(Form("Wp_xP_bin_%i",i));
		if( Wp_xp_bins_sim[i]->Integral() == 0 ) continue;
		
		Wp_xp_bins_bac[i] -> Scale( bacnorm );
		Wp_xp_bins_dat[i] -> Add( Wp_xp_bins_bac[i] , -1 );

		Wp_xp_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = Wp_xp_bins_dat[i]->Integral() / Wp_xp_bins_sim[i]->Integral();
		//pn_xp_bins_sim[i] -> Scale( simnorm );

		c2->cd(i+1);
		Wp_xp_bins_dat[i]->Rebin(2);
		Wp_xp_bins_sim[i]->Rebin(2);
		TString current_title = Wp_xp_bins_dat[i]->GetTitle();
		Wp_xp_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(Wp_xp_bins_dat[i],Wp_xp_bins_sim[i],"W' [GeV]","Counts");
	}
	c2->SaveAs("wp_xp_bins.pdf");
	
	return;
}

void label1D(TH1F* data, TH1F* sim, TString xlabel, TString ylabel){
	data->SetLineColor(4);
	data->SetMarkerColor(4);
	data->SetMarkerStyle(8);
	data->SetMarkerSize(1);
	data->SetStats(0);

	sim->SetLineColor(2);
	sim->SetLineWidth(1);
	sim->SetStats(0);

	data->Draw("p");
	sim->Draw("hist,same");

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	data->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}
