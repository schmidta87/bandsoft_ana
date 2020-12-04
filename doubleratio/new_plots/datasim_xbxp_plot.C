void datasim_xbxp_plot(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1F* data, TH1F* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1F* data, TH1F* sim, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileBac = new TFile(inBac);
	TFile * inFileSim = new TFile(inSim);

	// Get background normalization
	TVector3 * bacvector = (TVector3*)inFileDat->Get("bacnorm");
	double bacnorm = bacvector->X();
	TH1F * xP_bac = (TH1F*)inFileBac->Get("xP");
	TH1F * xB_bac = (TH1F*)inFileBac->Get("xB");
	bacnorm = bacnorm/xP_bac->Integral() ;
	xP_bac->Scale( bacnorm );
	xB_bac->Scale( bacnorm );

	// Get simulation normalization
	TH1F * xP_sim = (TH1F*)inFileSim->Get("xP");
	TH1F * xP_dat = (TH1F*)inFileDat->Get("xP");
	TH1F * xB_sim = (TH1F*)inFileSim->Get("xB");
	TH1F * xB_dat = (TH1F*)inFileDat->Get("xB");
	xP_dat->Add( xP_bac, -1 );
	xB_dat->Add( xB_bac, -1 );
	double full_simnorm = xP_dat->Integral() / xP_sim->Integral();
	xP_sim->Scale( full_simnorm );
	xB_sim->Scale( full_simnorm );


	// x' absolute and ratio
	TCanvas * c0 = new TCanvas("c0","",800,600);
	xP_dat->Rebin(4);
	xP_sim->Rebin(4);
	c0->Divide(1,2);
	c0->cd(1);
	xP_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(xP_dat,xP_sim,"x' ","Counts");
	c0->cd(2);
	xB_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D_ratio(xP_dat,xP_sim,"x' ","Data/Sim");
	c0->SaveAs("full_xp.pdf");

	// xB absolute and ratio
	TCanvas * c1 = new TCanvas("c1","",800,600);
	xB_dat->Rebin(4);
	xB_sim->Rebin(4);
	c1->Divide(1,2);
	c1->cd(1);
	xB_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(xB_dat,xB_sim,"x_{B} ","Counts");
	c1->cd(2);
	xB_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D_ratio(xB_dat,xB_sim,"x_{B} ","Data/Sim");
	c1->SaveAs("full_xB.pdf");

	// Grab intereseted histograms
	const int as_bins = 4;
	TH1F ** xb_bins_dat = new TH1F*[as_bins];
	TH1F ** xp_bins_dat = new TH1F*[as_bins];
	TH1F ** xb_bins_bac = new TH1F*[as_bins];
	TH1F ** xp_bins_bac = new TH1F*[as_bins];
	TH1F ** xb_bins_sim = new TH1F*[as_bins];
	TH1F ** xp_bins_sim = new TH1F*[as_bins];

	TCanvas * c2 = new TCanvas("c2","",800,600);
	c2->Divide(3,2);
	for( int i = 0 ; i < as_bins ; i++ ){
		xb_bins_dat[i] = (TH1F*) inFileDat->Get(Form("xB_aS_bin_%i",i));
		xb_bins_bac[i] = (TH1F*) inFileBac->Get(Form("xB_aS_bin_%i",i));
		xb_bins_sim[i] = (TH1F*) inFileSim->Get(Form("xB_aS_bin_%i",i));


		if( xb_bins_sim[i]->Integral() == 0 ) continue;
		
		xb_bins_bac[i] -> Scale( bacnorm );
		xb_bins_dat[i] -> Add( xb_bins_bac[i] , -1 );

		xb_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = xb_bins_dat[i]->Integral() / xb_bins_sim[i]->Integral();
		//pn_xb_bins_sim[i] -> Scale( simnorm );

		c2->cd(i);
		xb_bins_dat[i]->Rebin(4);
		xb_bins_sim[i]->Rebin(4);
		TString current_title = xb_bins_dat[i]->GetTitle();
		xb_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));
		label1D(xb_bins_dat[i],xb_bins_sim[i],"x_{B}","Counts");

		c2->cd(i+3);
		label1D_ratio(xb_bins_dat[i],xb_bins_sim[i],"x_{B}","Data/Sim");
	}
	c2->SaveAs("ratioxB_as_bins.pdf");

	TCanvas * c3 = new TCanvas("c3","",800,600);
	c3->Divide(3,2);
	for( int i = 0 ; i < as_bins ; i++ ){
		xp_bins_dat[i] = (TH1F*) inFileDat->Get(Form("xP_aS_bin_%i",i));
		xp_bins_bac[i] = (TH1F*) inFileBac->Get(Form("xP_aS_bin_%i",i));
		xp_bins_sim[i] = (TH1F*) inFileSim->Get(Form("xP_aS_bin_%i",i));


		if( xp_bins_sim[i]->Integral() == 0 ) continue;
		
		xp_bins_bac[i] -> Scale( bacnorm );
		xp_bins_dat[i] -> Add( xp_bins_bac[i] , -1 );

		xp_bins_sim[i] -> Scale( full_simnorm );

		// Rescale each plot to compare shape
		double simnorm = xp_bins_dat[i]->Integral() / xp_bins_sim[i]->Integral();
		//pn_xb_bins_sim[i] -> Scale( simnorm );

		c3->cd(i);
		xp_bins_dat[i]->Rebin(4);
		xp_bins_sim[i]->Rebin(4);
		TString current_title = xp_bins_dat[i]->GetTitle();
		xp_bins_dat[i]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));
		label1D(xp_bins_dat[i],xp_bins_sim[i],"xB","Counts");

		c3->cd(i+3);
		label1D_ratio(xp_bins_dat[i],xp_bins_sim[i],"x'","Data/sim");
	}
	c3->SaveAs("ratioxP_as_bins.pdf");

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

void label1D_ratio(TH1F* data, TH1F* sim, TString xlabel, TString ylabel){
	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0.5,1.5);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
