void asBins_wpxp_ratios(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label2D(TH2D* data, TString xlabel, TString ylabel);
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax );

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileBac = new TFile(inBac);
	TFile * inFileSim = new TFile(inSim);

	// Get the TTrees
	TTree * inTreeDat = (TTree*) inFileDat->Get("tagged");
	TTree * inTreeBac = (TTree*) inFileBac->Get("tagged");
	TTree * inTreeSim = (TTree*) inFileSim->Get("tagged");

	// Get and set the background normalization
	TVector3 * datnorm = (TVector3*)inFileDat->Get("bacnorm");
	TVector3 * bacnorm = (TVector3*)inFileBac->Get("bacnorm");
	inTreeBac->SetWeight( datnorm->X() / bacnorm->X() );


	// Define histograms we want to plot:
	const int nAs_bins = 6;
	const double As_min = 1.35;
	const double As_max = 1.65;
	const int nXp_bins = 10;
	const int nXp_coarse_bins = 4;
	const double Xp_min = 0.25;
	const double Xp_max = 0.65;
	const int nWp_bins = 10;
	const int nWp_coarse_bins = 4;
	const double Wp_min = 1.8;
	const double Wp_max = 2.6;
	TH1D ** xp_asBin_dat = new TH1D*[nAs_bins];
	TH1D ** xp_asBin_bac = new TH1D*[nAs_bins];
	TH1D ** xp_asBin_sim = new TH1D*[nAs_bins];
	TH1D ** wp_asBin_dat = new TH1D*[nAs_bins];
	TH1D ** wp_asBin_bac = new TH1D*[nAs_bins];
	TH1D ** wp_asBin_sim = new TH1D*[nAs_bins];
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		xp_asBin_dat[bin] = new TH1D(Form("h1_xp_asBin_dat_%i",bin),"",nXp_bins,Xp_min,Xp_max);
		xp_asBin_bac[bin] = new TH1D(Form("h1_xp_asBin_bac_%i",bin),"",nXp_bins,Xp_min,Xp_max);
		xp_asBin_sim[bin] = new TH1D(Form("h1_xp_asBin_sim_%i",bin),"",nXp_bins,Xp_min,Xp_max);
		wp_asBin_dat[bin] = new TH1D(Form("h1_wp_asBin_dat_%i",bin),"",nWp_bins,Wp_min,Wp_max);
		wp_asBin_bac[bin] = new TH1D(Form("h1_wp_asBin_bac_%i",bin),"",nWp_bins,Wp_min,Wp_max);
		wp_asBin_sim[bin] = new TH1D(Form("h1_wp_asBin_sim_%i",bin),"",nWp_bins,Wp_min,Wp_max);
	}

	// Define any global cuts we want:
	TCut pN_cut = "tag[nleadindex]->getMomentumN().Mag() > 0.3";
	TCut pt_cut = "tag[nleadindex]->getPt().Mag() < 0.1";
	TCut glob_cuts = pN_cut && pt_cut;

	// Figure out sim normalization by integral in As
	TH1D * h1_as_dat = new TH1D("h1_as_dat","",50,1,2);
	TH1D * h1_as_bac = new TH1D("h1_as_bac","",50,1,2);
	TH1D * h1_as_sim = new TH1D("h1_as_sim","",50,1,2);
	inTreeDat->Draw("tag[nleadindex]->getAs() >> h1_as_dat", glob_cuts );
	inTreeBac->Draw("tag[nleadindex]->getAs() >> h1_as_bac", glob_cuts );
	inTreeSim->Draw("tag[nleadindex]->getAs() >> h1_as_sim", glob_cuts );
	// 	background subtraction
	h1_as_dat->Add( h1_as_bac, -1 );
	//	normalize simulation
	double dat_int = h1_as_dat->Integral();
	double sim_int = h1_as_sim->Integral();
	double simnorm = dat_int / sim_int;
	
	TCanvas * c1_as_wpxp_tag = new TCanvas("c1_as_wpxp_tag","",800,600);
	c1_as_wpxp_tag->Divide(nAs_bins,2);
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		// Define this alphaS bin min/max cut
		double this_min_as = As_min + bin*(As_max - As_min)/nAs_bins;
		double this_max_as = As_min + (bin+1)*(As_max - As_min)/nAs_bins;
		TCut this_as_cut = Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",this_min_as,this_max_as);
		//	set title of histograms
		TString current_title = Form("%.2f < As < %.2f",this_min_as,this_max_as);
		wp_asBin_dat[bin]->SetTitle( current_title );
		wp_asBin_sim[bin]->SetTitle( current_title );
		xp_asBin_dat[bin]->SetTitle( current_title );
		xp_asBin_sim[bin]->SetTitle( current_title );


		// 	fill the histograms of W' for given As bin
		c1_as_wpxp_tag->cd(bin+1);
		inTreeDat->Draw(Form("tag[nleadindex]->getWp() >> h1_wp_asBin_dat_%i",bin),this_as_cut && glob_cuts );
		inTreeBac->Draw(Form("tag[nleadindex]->getWp() >> h1_wp_asBin_bac_%i",bin),this_as_cut && glob_cuts );
		inTreeSim->Draw(Form("tag[nleadindex]->getWp() >> h1_wp_asBin_sim_%i",bin),this_as_cut && glob_cuts );
		// 	fill the histograms of X' for given As bin
		c1_as_wpxp_tag->cd(bin+1 + nAs_bins);
		inTreeDat->Draw(Form("tag[nleadindex]->getXp() >> h1_xp_asBin_dat_%i",bin),this_as_cut && glob_cuts );
		inTreeBac->Draw(Form("tag[nleadindex]->getXp() >> h1_xp_asBin_bac_%i",bin),this_as_cut && glob_cuts );
		inTreeSim->Draw(Form("tag[nleadindex]->getXp() >> h1_xp_asBin_sim_%i",bin),this_as_cut && glob_cuts );

		// 	do background subtraction for both W',X':
		wp_asBin_dat[bin]->Add( wp_asBin_bac[bin] , -1 );
		xp_asBin_dat[bin]->Add( xp_asBin_bac[bin] , -1 );
		// 	do simulation scaling for both
		wp_asBin_sim[bin]->Scale( simnorm );
		xp_asBin_sim[bin]->Scale( simnorm );

		c1_as_wpxp_tag->cd(bin+1);
		label1D_ratio( wp_asBin_dat[bin] , wp_asBin_sim[bin], "W'","Data/Sim",0,2);
		c1_as_wpxp_tag->cd(bin+1 + nAs_bins);
		label1D_ratio( xp_asBin_dat[bin] , xp_asBin_sim[bin], "x'","Data/Sim",0,2);
	
	}
	c1_as_wpxp_tag->SaveAs("asBins_wpxp_ratios.pdf");


	return;
}

void label2D(TH2D* data, TString xlabel, TString ylabel){
	//data->SetLineColor(4);
	//data->SetMarkerColor(4);
	//data->SetMarkerStyle(8);
	//data->SetMarkerSize(1);
	data->SetStats(0);

	data->Draw("colz");

	
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);


	return;
}
void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(4);
	data->SetMarkerColor(4);
	data->SetMarkerStyle(8);
	data->SetMarkerSize(1);
	data->SetStats(0);

	sim->SetLineColor(2);
	sim->SetLineWidth(1);
	sim->SetStats(0);

	sim->Draw("hist");
	data->Draw("p,same");

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	sim->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	sim->GetXaxis()->SetTitle(xlabel);
	sim->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}

void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax ){
	gStyle->SetOptFit(1);
	
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
	data_copy->SetTitle( sim_copy->GetTitle() );
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	data_copy->GetYaxis()->SetRangeUser(ymin,ymax);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
