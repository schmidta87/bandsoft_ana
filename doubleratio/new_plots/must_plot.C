void must_plot(TString inFileDatTagName, TString inFileBacTagName,
	       	 TString inFileSimTagName, TString inFileDatIncName, 
		 TString inFileSimIncName ){
	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);

	TFile * inFileDatTag = new TFile(inFileDatTagName);
	TFile * inFileBacTag = new TFile(inFileBacTagName);
	TFile * inFileSimTag = new TFile(inFileSimTagName);
	TFile * inFileDatInc = new TFile(inFileDatIncName);
	TFile * inFileSimInc = new TFile(inFileSimIncName);

	TTree * inTreeDatTag 	= (TTree*) inFileDatTag->Get("tagged");
	TTree * inTreeBacTag 	= (TTree*) inFileBacTag->Get("tagged");
	TTree * inTreeSimTag 	= (TTree*) inFileSimTag->Get("tagged");
	TTree * inTreeDatInc	= (TTree*) inFileDatInc->Get("electrons");
	TTree * inTreeSimInc	= (TTree*) inFileSimInc->Get("electrons");

	// Number of events in tree to read:
	const int nEvents = 100000;

	// Get background normalization
	TVector3 * bacvector = (TVector3*)inFileDatTag->Get("bacnorm");
	double bacnorm = bacvector->X();

	// Get simulation normalization
	double 		L_tag = 2.28e10;
	double 		L_inc = 9.61e6;
	double 		Q_inc = 535885;
	double		Q_tag = 4.37663e+07;

	// Define the histograms for the ratio:
	const int nAs_bins = 5;
	const double As_min = 1.2;
	const double As_max = 1.7;
	const int nXp_bins = 6;
	const double Xp_min = 0.1;
	const double Xp_max = 0.7;
	const int nXb_bins = 6;
	const double Xb_min = 0.1;
	const double Xb_max = 0.7;
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
	// These histograms all share a Q2 bin of [2,8]
	// The tag histograms have theta_nq bin of [150,180]
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

	// Define the bins of theta_nq and Q2:
	const double min_Q2 = 2;
	const double max_Q2 = 8;
	const double min_CosThetaNQ = -1;
	const double max_CosThetaNQ = -0.8;
	TCut Q2_bin = Form("eHit->getQ2() > %f && eHit->getQ2() < %f",min_Q2,max_Q2);
	TCut ThetaNQ_bin = Form("cos(tag[0]->getThetaNQ()) > %f && cos(tag[0]->getThetaNQ()) < %f",min_CosThetaNQ,max_CosThetaNQ);
	TCut tagged = Q2_bin && ThetaNQ_bin;
	TCut inclusive = Q2_bin;

	// Full histograms for normalization of background:
	TH1D * h1_as_dat = new TH1D("h1_as_dat","",nAs_bins,As_min,As_max);
	TH1D * h1_as_bac = new TH1D("h1_as_bac","",nAs_bins,As_min,As_max);
	TH1D * h1_as_sim = new TH1D("h1_as_sim","",nAs_bins,As_min,As_max);
	inTreeDatTag->Draw("tag[0]->getAs() >> h1_as_dat",tagged,"",nEvents);
	inTreeBacTag->Draw("tag[0]->getAs() >> h1_as_bac",tagged,"",nEvents);
	inTreeSimTag->Draw("tag[0]->getAs() >> h1_as_sim",tagged,"",nEvents);
	bacnorm = bacnorm / h1_as_bac->Integral();

	TCanvas * c1_as_full = new TCanvas("c1_as_full","",800,600);
	// Do the background subtraction:
	h1_as_bac->Scale(bacnorm);
	h1_as_dat->Add(h1_as_bac,-1);
	//h1_as_dat->Scale(Q_inc/Q_tag);
	// Do the simulation normalization:
	double simnorm_tag = h1_as_dat->Integral() / h1_as_sim->Integral();
	simnorm_tag = 1./L_tag;
	//h1_as_sim->Scale(simnorm_tag);
	// Plot the data and simulation:
	h1_as_dat->SetTitle(Form("Full distribution, C_{bac} = %f, C_{sim} = %f",bacnorm,simnorm_tag));
	label1D(h1_as_dat,h1_as_sim,"Alpha_{S}","Counts");
	c1_as_full->SaveAs("full_as_tag.pdf");

	// First let's get:
	//	R(tag) = 	N_tag,data (Q2,theta_nq,x';alphaS) / [ N_tag,sim(Q2,theta_nq,x';alphaS) ]
	//	as just individual histograms for numerator and denominator
	TCanvas * c1_as_xp_tag = new TCanvas("c1_as_xp_tag","",800,600);
	c1_as_xp_tag->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_as_xp_tag->cd(bin+1);
		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[0]->getXp() > %f && tag[0]->getXp() < %f",this_min_xp,this_max_xp);
		TCut this_tag = tagged && this_xp_cut;

		inTreeDatTag->Draw(Form("tag[0]->getAs() >> h1_as_xp_dat_%i",bin),this_tag,"",nEvents);
		inTreeBacTag->Draw(Form("tag[0]->getAs() >> h1_as_xp_bac_%i",bin),this_tag,"",nEvents);
		inTreeSimTag->Draw(Form("tag[0]->getAs() >> h1_as_xp_sim_%i",bin),this_tag,"",nEvents);

		// Do the background subtraction:
		h1_as_xp_bac[bin]->Scale(bacnorm);
		h1_as_xp_dat[bin]->Add(h1_as_xp_bac[bin],-1);
		// Plot data and simulation:
		h1_as_xp_dat[bin]->SetTitle(Form("%f < x' < %f",this_min_xp,this_max_xp));

		label1D(h1_as_xp_dat[bin],h1_as_xp_sim[bin],"Alpha_{S}","Counts");
	}
	c1_as_xp_tag->SaveAs("as_xpBins_tag.pdf");
	

	TCanvas * c1_q2_full = new TCanvas("c1_q2_full","",800,600);
	// Now that we have the tagged histograms, we just need to get the numbers of inclusive at xB=x'
	TH1D * h1_xb_dat = new TH1D("h1_xb_dat","",nXb_bins,Xb_min,Xb_max);
	TH1D * h1_xb_sim = new TH1D("h1_xb_sim","",nXb_bins,Xb_min,Xb_max);
	inTreeDatInc->Draw("eHit->getXb() >> h1_xb_dat",inclusive,"",nEvents);
	inTreeSimInc->Draw("eHit->getXb() >> h1_xb_sim",inclusive,"",nEvents);
	// Do the simulation normalization:
	double simnorm_inc = h1_xb_dat->Integral() / h1_xb_sim->Integral();
	simnorm_inc = 1./L_inc;
	//h1_xb_sim->Scale(simnorm_inc);
	// Plot data and simulation:
	h1_as_dat->SetTitle(Form("Full distribution, C_{sim} = %f",simnorm_inc));
	label1D(h1_xb_dat,h1_xb_sim,"x_{B}","Counts");

	c1_q2_full->SaveAs("q2_full_inc.pdf");

	// Now we can go get:
	//	R(inc) = 	N_inc,data (Q2,x=x') / N_inc,sim(Q2,x=x') 
	//	by saving each number individiually
	TCanvas * c1_q2_xb_inc = new TCanvas("c1_q2_xb_inc","",800,600);
	c1_q2_xb_inc->Divide(3,2);
	double DatIncCounts[nXb_bins] = {0};
	double SimIncCounts[nXb_bins] = {0};
	for( int bin = 0 ; bin < nXb_bins ; bin++ ){
		c1_q2_xb_inc->cd(bin+1);
		double this_min_xb = Xb_min+0.05 + bin*(Xb_max - Xb_min)/nXb_bins;
		double this_max_xb = Xb_min+0.05 + (bin+1)*(Xb_max - Xb_min)/nXb_bins;
		TCut this_xb_cut = Form("eHit->getXb() > %f && eHit->getXb() < %f",this_min_xb,this_max_xb);
		TCut this_inc = inclusive && this_xb_cut;

		inTreeDatInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_dat_%i",bin),this_inc,"",nEvents);
		inTreeSimInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_sim_%i",bin),this_inc,"",nEvents);

		DatIncCounts[bin] = h1_q2_xb_dat[bin]->Integral();
		SimIncCounts[bin] = h1_q2_xb_sim[bin]->Integral();

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

		cout << "x' bin: " << bin << "\n";
	        cout << "Y_tag data integral: " 	<< h1_as_xp_dat[bin]->Integral() << "\n";
	        cout << "N_tag sim integral: " 		<< h1_as_xp_sim[bin]->Integral() << "\n";
		cout << "Y_inc data integral: " 	<< DatIncCounts[bin] << "\n";
		cout << "N_inc sim integral: " 		<< SimIncCounts[bin] << "\n";
		cout << "Q_tag: " 			<< Q_tag << "\n";
		cout << "Q_inc: " 			<< Q_inc << "\n";
		cout << "L_tag: " 			<< L_tag << "\n";
		cout << "L_inc: " 			<< L_inc << "\n";
		cout << "\n";


		h1_as_xp_dat[bin] -> Scale( 1./Q_tag );  // N_data,tag
		h1_as_xp_sim[bin] -> Scale( 1./L_tag );	 // N_sim,tag

		double N_inc_dat = (DatIncCounts[bin] / Q_inc);
		double N_inc_sim = (SimIncCounts[bin] / L_inc);

		h1_as_xp_dat[bin]->Scale( N_inc_sim / N_inc_dat );	// N_dat,tag / N_dat,inc

		
		label1D_ratio(h1_as_xp_dat[bin],h1_as_xp_sim[bin],"Alpha_{S}","R=R_{tag}/R_{sim}");

	}
	c1_must->SaveAs("ratio_tag_inc.pdf");

	TCanvas * ctest = new TCanvas("ctest","",800,600);
	ctest->Divide(3,2);

	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		ctest->cd(bin+1);

		h1_as_xp_dat[bin]->Divide( h1_as_xp_sim[bin] );
		label1D_ratio(h1_as_xp_dat[bin],h1_as_xp_dat[1],"Alpha_{S}","R=R_{tag}/R_{sim} / R_(x'=0.3)");
	}
	ctest->SaveAs("ratio_tag_inc_xp03.pdf");

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

void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
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
	data_copy->GetYaxis()->SetRangeUser(0,2);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
