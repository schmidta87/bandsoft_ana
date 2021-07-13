void must_plot(TString inFileDatTagName, TString inFileBacTagName,
	       	 TString inFileSimTagName, TString inFileDatIncName, 
		 TString inFileSimIncName ){
	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, int color);



	// Load TFiles
	TFile * inFileDatTag = new TFile(inFileDatTagName);
	TFile * inFileBacTag = new TFile(inFileBacTagName);
	TFile * inFileSimTag = new TFile(inFileSimTagName);
	TFile * inFileDatInc = new TFile(inFileDatIncName);
	TFile * inFileSimInc = new TFile(inFileSimIncName);

	// Get the TTrees
	TTree * inTreeDatTag = (TTree*) inFileDatTag->Get("tagged");
	TTree * inTreeBacTag = (TTree*) inFileBacTag->Get("tagged");
	TTree * inTreeSimTag = (TTree*) inFileSimTag->Get("tagged");
	TTree * inTreeDatInc = (TTree*) inFileDatInc->Get("electrons");
	TTree * inTreeSimInc = (TTree*) inFileSimInc->Get("electrons");

	// Get and set the background normalization for tagged
	TVector3 * datnorm = (TVector3*)inFileDatTag->Get("bacnorm");
	TVector3 * bacnorm = (TVector3*)inFileBacTag->Get("bacnorm");
	inTreeBacTag->SetWeight( datnorm->X() / bacnorm->X() );


	// Get simulation normalization
	double 		L_tag = 2.28e10 / (130.1);
	double 		L_inc = 9.61e6 / (130.1);
	double 		Q_inc = 535885;
	double		Q_tag = 3.72e+07;

	// Define the histograms for the ratio:
	const int nAs_bins = 3;
	const double As_min = 1.3;
	const double As_max = 1.6;
	const int nXp_bins = 4;
	const double Xp_min = 0.2;
	const double Xp_max = 0.6;
	const int nXb_bins = 4;
	const double Xb_min = 0.2;
	const double Xb_max = 0.6;
	// What we want:
	//			[ N_tag,data (Q2,theta_nq,x';alphaS) / Q_tag] / [ N_tag,sim(Q2,theta_nq,x';alphaS) / L_tag,sim ]
	//	R(tag/inc) = 	------------------------------------------------------------------------------------------------
	//				     [ N_inc,data (Q2, x=x') / Q_inc] / [ N_inc,sim (Q2, x=x') / L_inc,sim ]
	//
	//			[ N_tag,data / Q_tag ]	 [ N_inc,sim / L_inc,sim ]
	//	R(tag/inc) = 	---------------------- * -------------------------
	//			[ N_inc,data / Q_inc ]   [ N_tag,sim / L_tag,sim ]
	//
	//			 [ N_tag,data / Q_tag ]	     [ N_inc,sim / L_inc,sim ]
	//	R(tag/inc) = 	------------------------   * -------------------------
	//			[ N_tag,sim / L_tag,sim ]      [ N_inc,dat / Q_inc ]
	//
	// These histograms all share a Q2 bin of [2,10]
	// The tag histograms have cos_theta_nq bin of [-1,-0.8]
	TH1D ** h1_as_xp_dat = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for dat
	TH1D ** h1_as_xp_bac = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for dat
	TH1D ** h1_as_xp_sim = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for sim
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		h1_as_xp_dat[bin]	= new TH1D(Form("h1_as_xp_dat_%i",bin),"",nAs_bins,As_min,As_max);
		h1_as_xp_bac[bin]	= new TH1D(Form("h1_as_xp_bac_%i",bin),"",nAs_bins,As_min,As_max);
		h1_as_xp_sim[bin]	= new TH1D(Form("h1_as_xp_sim_%i",bin),"",nAs_bins,As_min,As_max);
	}
	TH1D ** h1_q2_xb_dat = new TH1D*[10];
	TH1D ** h1_q2_xb_sim = new TH1D*[10];
	for( int bin = 0 ; bin < nXb_bins ; bin++ ){
		h1_q2_xb_dat[bin]	= new TH1D(Form("h1_q2_xb_dat_%i",bin),"",10,0,10);
		h1_q2_xb_sim[bin]	= new TH1D(Form("h1_q2_xb_sim_%i",bin),"",10,0,10);
	}


	// First let's get:
	//	R(tag) = 	N_tag,data (Q2,theta_nq,x';alphaS) / [ N_tag,sim(Q2,theta_nq,x';alphaS) ]
	//	as just individual histograms for numerator and denominator
	TCanvas * c1_as_xp_tag = new TCanvas("c1_as_xp_tag","",800,600);
	c1_as_xp_tag->Divide(3,2);
	double TAGSUM = 0;
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_as_xp_tag->cd(bin+1);
		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[nleadindex]->getXp() > %f && tag[nleadindex]->getXp() < %f",this_min_xp,this_max_xp);

		inTreeDatTag->Draw(Form("tag[nleadindex]->getAs() >> h1_as_xp_dat_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeBacTag->Draw(Form("tag[nleadindex]->getAs() >> h1_as_xp_bac_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeSimTag->Draw(Form("tag[nleadindex]->getAs() >> h1_as_xp_sim_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");

		// Do the background subtraction:
		h1_as_xp_dat[bin]->Add(h1_as_xp_bac[bin],-1);

		// Do the normalization for simulation and data:
		h1_as_xp_dat[bin]->Scale(1./Q_tag);
		h1_as_xp_sim[bin]->Scale(1./L_tag);

		TAGSUM += h1_as_xp_dat[bin]->Integral();
		// Plot data and simulation:
		double new_simnorm = h1_as_xp_dat[bin]->Integral() / h1_as_xp_sim[bin]->Integral();
		TString current_title = Form("%f < x' < %f",this_min_xp,this_max_xp);
		h1_as_xp_dat[bin]->SetTitle( current_title );
		label1D(h1_as_xp_dat[bin],h1_as_xp_sim[bin],"Alpha_{S}","Counts");
	}
	c1_as_xp_tag->SaveAs("as_xpBins_tag.pdf");



	// Now we can go get:
	//	R(inc) = 	N_inc,data (Q2,x=x') / N_inc,sim(Q2,x=x') 
	//	by saving each number individiually
	TCanvas * c1_q2_xb_inc = new TCanvas("c1_q2_xb_inc","",800,600);
	c1_q2_xb_inc->Divide(3,2);
	double DatIncCounts[nXb_bins] = {0};
	double SimIncCounts[nXb_bins] = {0};
	double INCSUM = 0;
	for( int bin = 0 ; bin < nXb_bins ; bin++ ){
		c1_q2_xb_inc->cd(bin+1);
		double this_min_xb = Xb_min+0.05 + bin*(Xb_max - Xb_min)/nXb_bins;
		double this_max_xb = Xb_min+0.05 + (bin+1)*(Xb_max - Xb_min)/nXb_bins;
		TCut this_xb_cut = Form("eHit->getXb() > %f && eHit->getXb() < %f",this_min_xb,this_max_xb);

		inTreeDatInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_dat_%i",bin),this_xb_cut);
		inTreeSimInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_sim_%i",bin),this_xb_cut);

		// Do the normalization for simulation and data:
		h1_q2_xb_dat[bin]->Scale(1./Q_inc);
		h1_q2_xb_sim[bin]->Scale(1./L_inc);

		DatIncCounts[bin] = h1_q2_xb_dat[bin]->Integral();
		SimIncCounts[bin] = h1_q2_xb_sim[bin]->Integral();
		INCSUM += DatIncCounts[bin];

		// Plot data and simulation:
		h1_q2_xb_dat[bin]->SetTitle(Form("%f < x_{B} < %f",this_min_xb,this_max_xb));
		label1D(h1_q2_xb_dat[bin],h1_q2_xb_sim[bin],"Q^{2}","Counts");


	}
	c1_q2_xb_inc->SaveAs("q2_xbBins_inc.pdf");


	// Now do the R_tag/R_inc ratio
	TCanvas * c1_must = new TCanvas("c1_must","",800,600);
	c1_must->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_must->cd(bin+1);


		//h1_as_xp_dat[bin] -> Scale( 1./Q_tag );  // N_data,tag
		//h1_as_xp_sim[bin] -> Scale( 1./L_tag );	 // N_sim,tag

		double N_inc_dat = (DatIncCounts[bin]);
		double N_inc_sim = (SimIncCounts[bin]);

		h1_as_xp_dat[bin]->Scale( N_inc_sim / N_inc_dat );	// N_dat,tag / N_dat,inc

		
		label1D_ratio(h1_as_xp_dat[bin],h1_as_xp_sim[bin],"Alpha_{S}","R=R_{tag}/R_{inc}",2);

	}
	c1_must->SaveAs("ratio_tag_inc.pdf");

	TCanvas * c1_must_xp03 = new TCanvas("c1_must_xp03","",800,600);
	c1_must_xp03->Divide(2,2);

	int colors[nXp_bins] = {2,1,8,9};
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_must_xp03->cd(bin+1);

		h1_as_xp_dat[bin]->Divide( h1_as_xp_sim[bin] );
		label1D_ratio(h1_as_xp_dat[bin],h1_as_xp_dat[0],"Alpha_{S}","R_{normed}",colors[bin]);

	}
	c1_must_xp03->SaveAs("ratio_tag_inc_xp03.pdf");

	// Overlay all the plots on top of each other
	TCanvas * c1_must_xp03_overlaid = new TCanvas("c1_must_xp03_overlaid","",800,600);
	c1_must_xp03_overlaid->cd();
	
	//h1_as_xp_dat[0]->Divide(h1_as_xp_dat[0]);
	h1_as_xp_dat[1]->Divide(h1_as_xp_dat[0]);
	h1_as_xp_dat[2]->Divide(h1_as_xp_dat[0]);
	h1_as_xp_dat[3]->Divide(h1_as_xp_dat[0]);
	
	h1_as_xp_dat[1]->SetLineColor(colors[1]);
	h1_as_xp_dat[1]->SetMarkerColor(colors[1]);
	h1_as_xp_dat[1]->SetLineWidth(3);
	h1_as_xp_dat[1]->SetStats(0);

	h1_as_xp_dat[2]->SetLineColor(colors[2]);
	h1_as_xp_dat[2]->SetMarkerColor(colors[2]);
	h1_as_xp_dat[2]->SetLineWidth(3);
	h1_as_xp_dat[2]->SetStats(0);

	h1_as_xp_dat[3]->SetLineColor(colors[3]);
	h1_as_xp_dat[3]->SetMarkerColor(colors[3]);
	h1_as_xp_dat[3]->SetLineWidth(3);
	h1_as_xp_dat[3]->SetStats(0);

	h1_as_xp_dat[1]->Draw("ep");
	h1_as_xp_dat[2]->Draw("ep,same");
	h1_as_xp_dat[3]->Draw("ep,same");

	TLine* line = new TLine(1.3, 1.,1.6, 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	h1_as_xp_dat[1]->GetYaxis()->SetRangeUser(0.7,3);
	h1_as_xp_dat[1]->GetXaxis()->SetTitle("Alpha_{S}");
	h1_as_xp_dat[1]->GetYaxis()->SetTitle("R_{normed}");
	h1_as_xp_dat[1]->SetTitle("");
	
	c1_must_xp03_overlaid->SaveAs("ratio_tag_inc_xp03_overlaid.pdf");

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

void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, int color){
	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(color);
	data_copy->SetMarkerColor(color);
	data_copy->SetLineWidth(3);
	data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);


	data_copy->Draw("ep");
	TLine* line = new TLine(1.3, 1.,1.6, 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0.7,3);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	//data_copy->Fit("pol1","ESR","",1.3,1.6);
	return;
}
