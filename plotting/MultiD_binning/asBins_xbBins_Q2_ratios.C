void asBins_xbBins_Q2_ratios(TString inDat, TString inBac, TString inSim){

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
	const int nAs_bins = 5;
	const double As_min = 1.35;
	const double As_max = 1.6;
	const int nXb_bins = 10;
	const int nXb_coarse_bins = 4;
	const double Xb_min = 0.15;
	const double Xb_max = 0.55;
	const int nQ2_bins = 16;
	const int nQ2_coarse_bins = 4;
	const double Q2_min = 2.0;
	const double Q2_max = 6.0;
	TH1D *** Q2_dat = new TH1D**[nAs_bins];
	TH1D *** Q2_bac = new TH1D**[nAs_bins];
	TH1D *** Q2_sim = new TH1D**[nAs_bins];
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		Q2_dat[bin] = new TH1D*[nXb_coarse_bins];
		Q2_bac[bin] = new TH1D*[nXb_coarse_bins];
		Q2_sim[bin] = new TH1D*[nXb_coarse_bins];
		for( int xb_bin = 0 ; xb_bin < nXb_coarse_bins ; xb_bin++ ){
			Q2_dat[bin][xb_bin] = new TH1D(Form("h1_Q2_dat_%i_%i",bin,xb_bin),"",nQ2_bins,Q2_min,Q2_max);
			Q2_bac[bin][xb_bin] = new TH1D(Form("h1_Q2_bac_%i_%i",bin,xb_bin),"",nQ2_bins,Q2_min,Q2_max);
			Q2_sim[bin][xb_bin] = new TH1D(Form("h1_Q2_sim_%i_%i",bin,xb_bin),"",nQ2_bins,Q2_min,Q2_max);
		}
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
	
	TCanvas * c1_as_xb_Q2_tag = new TCanvas("c1_as_xb_Q2_tag","",800,600);
	TCanvas * c1_as_xb_Q2_ratio = new TCanvas("c1_as_xb_Q2_ratio","",800,600);
	c1_as_xb_Q2_tag->Divide(nXb_coarse_bins,	nAs_bins);
	c1_as_xb_Q2_ratio->Divide(nXb_coarse_bins,	nAs_bins);
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		// Define this alphaS bin min/max cut
		double this_min_as = As_min + bin*(As_max - As_min)/nAs_bins;
		double this_max_as = As_min + (bin+1)*(As_max - As_min)/nAs_bins;
		TCut this_as_cut = Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",this_min_as,this_max_as);
		
		// Now loop over the X' bins we want:
		for( int xb_bin = 0 ; xb_bin < nXb_coarse_bins ; xb_bin++ ){
			double this_min_xb = Xb_min + xb_bin*(Xb_max - Xb_min)/nXb_coarse_bins;
			double this_max_xb = Xb_min + (xb_bin+1)*(Xb_max - Xb_min)/nXb_coarse_bins;
			TCut this_xb_cut = Form("eHit->getXb() > %f && eHit->getXb() < %f",this_min_xb,this_max_xb);
		
			//	set title of histograms
			TString current_title = Form("%.2f < As < %.2f, %.2f < x_{B} < %.2f",this_min_as,this_max_as,this_min_xb,this_max_xb);
			Q2_dat[bin][xb_bin]->SetTitle( current_title );
			Q2_sim[bin][xb_bin]->SetTitle( current_title );


			// 	fill the histograms of X' for given As bin, W' bin
			c1_as_xb_Q2_ratio->cd(xb_bin+1 + bin*nXb_coarse_bins);
			c1_as_xb_Q2_tag->cd(xb_bin+1 + bin*nXb_coarse_bins);
			inTreeDat->Draw(Form("eHit->getQ2() >> h1_Q2_dat_%i_%i",bin,xb_bin), this_as_cut && this_xb_cut && glob_cuts );
			inTreeBac->Draw(Form("eHit->getQ2() >> h1_Q2_bac_%i_%i",bin,xb_bin), this_as_cut && this_xb_cut && glob_cuts );
			inTreeSim->Draw(Form("eHit->getQ2() >> h1_Q2_sim_%i_%i",bin,xb_bin), this_as_cut && this_xb_cut && glob_cuts );

			if( Q2_dat[bin][xb_bin]->Integral() == 0 || Q2_bac[bin][xb_bin]->Integral() == 0 || Q2_sim[bin][xb_bin]->Integral() == 0 ){ 
				c1_as_xb_Q2_tag->cd(xb_bin+1 + bin*nXb_coarse_bins);
				c1_as_xb_Q2_ratio->cd(xb_bin+1 + bin*nXb_coarse_bins);
				continue;
			}

			// 	do background subtraction for X':
			Q2_dat[bin][xb_bin]->Add( Q2_bac[bin][xb_bin] , -1 );
			//	do simulation scaling
			Q2_sim[bin][xb_bin]->Scale( simnorm );

			
			c1_as_xb_Q2_tag->cd(xb_bin+1 + bin*nXb_coarse_bins);
			label1D( Q2_dat[bin][xb_bin] , Q2_sim[bin][xb_bin], "Q^{2}","Counts [a.u.]");

			c1_as_xb_Q2_ratio->cd(xb_bin+1 + bin*nXb_coarse_bins);
			label1D_ratio( Q2_dat[bin][xb_bin] , Q2_sim[bin][xb_bin], "Q^{2}","Data/Sim",0,2);
		}
	}
	c1_as_xb_Q2_tag		->SaveAs("asBins_xbBins_Q2.pdf");
	c1_as_xb_Q2_ratio	->SaveAs("asBins_xbBins_Q2_ratio.pdf");


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
