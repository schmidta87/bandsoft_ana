void pn_xbxp_plot(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1F* data, TH1F* sim, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileBac = new TFile(inBac);
	TFile * inFileSim = new TFile(inSim);

	// Get background normalization
	TVector3 * bacvector = (TVector3*)inFileDat->Get("bacnorm");
	double bacnorm = bacvector->X();
	TH1F * pn_bac = (TH1F*)inFileBac->Get("pn");
	bacnorm = bacnorm/pn_bac->Integral() ;
	pn_bac->Scale( bacnorm );

	// Get simulation normalization
	TH1F * pn_sim = (TH1F*)inFileSim->Get("pn");
	TH1F * pn_dat = (TH1F*)inFileDat->Get("pn");
	pn_dat->Add( pn_bac, -1 );
	double full_simnorm = pn_dat->Integral() / pn_sim->Integral();
	pn_sim->Scale( full_simnorm );

	TCanvas * c0 = new TCanvas("c0","",800,600);
	c0->cd();
	pn_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(pn_dat,pn_sim,"|p_{n}| [GeV/c]","Counts");
	c0->SaveAs("full_pn.pdf");

	// Grab intereseted histograms
	const int xb_bins = 4;
	const int xp_bins = 6;
	TH1F ** pn_xb_bins_dat = new TH1F*[xb_bins];
	TH1F ** pn_xp_bins_dat = new TH1F*[xp_bins];
	TH1F ** pn_xb_bins_bac = new TH1F*[xb_bins];
	TH1F ** pn_xp_bins_bac = new TH1F*[xp_bins];
	TH1F ** pn_xb_bins_sim = new TH1F*[xb_bins];
	TH1F ** pn_xp_bins_sim = new TH1F*[xp_bins];

	TCanvas * c1 = new TCanvas("c1","",800,600);
	c1->Divide(2,2);
	for( int i = 0 ; i < xb_bins ; i++ ){
		pn_xb_bins_dat[i] = (TH1F*) inFileDat->Get(Form("pn_xB_bin_%i",i));
		pn_xb_bins_bac[i] = (TH1F*) inFileBac->Get(Form("pn_xB_bin_%i",i));
		pn_xb_bins_sim[i] = (TH1F*) inFileSim->Get(Form("pn_xB_bin_%i",i));


		if( pn_xb_bins_sim[i]->Integral() == 0 ) continue;
		
		pn_xb_bins_bac[i] -> Scale( bacnorm );
		pn_xb_bins_dat[i] -> Add( pn_xb_bins_bac[i] , -1 );

		pn_xb_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = pn_xb_bins_dat[i]->Integral() / pn_xb_bins_sim[i]->Integral();
		//pn_xb_bins_sim[i] -> Scale( simnorm );

		c1->cd(i+1);
		pn_xb_bins_dat[i]->Rebin(2);
		pn_xb_bins_sim[i]->Rebin(2);
		TString current_title = pn_xb_bins_dat[i]->GetTitle();
		pn_xb_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(pn_xb_bins_dat[i],pn_xb_bins_sim[i],"|p_{n}| [GeV/c]","Counts");
	}
	c1->SaveAs("pn_xb_bins.pdf");

	TCanvas * c2 = new TCanvas("c2","",800,600);
	c2->Divide(3,2);
	for( int i = 0 ; i < xp_bins ; i++ ){
		pn_xp_bins_dat[i] = (TH1F*) inFileDat->Get(Form("pn_xP_bin_%i",i));
		pn_xp_bins_bac[i] = (TH1F*) inFileBac->Get(Form("pn_xP_bin_%i",i));
		pn_xp_bins_sim[i] = (TH1F*) inFileSim->Get(Form("pn_xP_bin_%i",i));
		if( pn_xp_bins_sim[i]->Integral() == 0 ) continue;
		
		pn_xp_bins_bac[i] -> Scale( bacnorm );
		pn_xp_bins_dat[i] -> Add( pn_xp_bins_bac[i] , -1 );

		pn_xp_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = pn_xp_bins_dat[i]->Integral() / pn_xp_bins_sim[i]->Integral();
		//pn_xp_bins_sim[i] -> Scale( simnorm );

		c2->cd(i+1);
		pn_xp_bins_dat[i]->Rebin(2);
		pn_xp_bins_sim[i]->Rebin(2);
		TString current_title = pn_xp_bins_dat[i]->GetTitle();
		pn_xp_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(pn_xp_bins_dat[i],pn_xp_bins_sim[i],"|p_{n}| [GeV/c]","Counts");
	}
	c2->SaveAs("pn_xp_bins.pdf");

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
