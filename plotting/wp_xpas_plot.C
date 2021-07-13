void wp_xpas_plot(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_sims(TH1D* sim2, TH1D* sim3, double label1, double label2);

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
	TH1D * wp_dat = new TH1D("wp_dat","wp_dat",13,1.7,3);
	TH1D * wp_bac = new TH1D("wp_bac","wp_bac",13,1.7,3);
	TH1D * wp_sim = new TH1D("wp_sim","wp_sim",13,1.7,3);
	//TH1D * wp_sim12 = new TH1D("wp_sim12","wp_sim12",13,1.7,3);
	//TH1D * wp_sim19 = new TH1D("wp_sim19","wp_sim19",13,1.7,3);
	const int nXp_bins = 6;
	const double Xp_min = 0.1;
	const double Xp_max = 0.7;
	TH1D ** wp_xp_bins_dat = new TH1D*[nXp_bins];
	TH1D ** wp_xp_bins_bac = new TH1D*[nXp_bins];
	TH1D ** wp_xp_bins_sim = new TH1D*[nXp_bins];
	//TH1D ** wp_xp_bins_sim12 = new TH1D*[nXp_bins];
	//TH1D ** wp_xp_bins_sim19 = new TH1D*[nXp_bins];
	
	// Draw the full wp distribution
	TCanvas * c_wp = new TCanvas("c_wp","",800,600);
	//c_wp->Divide(2,2);
	//c_wp->cd(1);
	inTreeDat->Draw("tag[nleadindex]->getWp() >> wp_dat","tag[nleadindex]->getMomentumN().Mag()> 0.3");
	inTreeBac->Draw("tag[nleadindex]->getWp() >> wp_bac","tag[nleadindex]->getMomentumN().Mag()> 0.3");
	inTreeSim->Draw("tag[nleadindex]->getWp() >> wp_sim","tag[nleadindex]->getMomentumN().Mag()> 0.3");
	//inTreeSim->Draw("tag->getWp() >> wp_sim12","nHits[nleadindex]->getEdep()/1E4 > 12");
	//inTreeSim->Draw("tag->getWp() >> wp_sim19","nHits[nleadindex]->getEdep()/1E4 > 19");

	wp_dat->Add(wp_bac,-1);
	double full_simnorm = (double)wp_dat->Integral() / wp_sim->Integral();
	//double full_simnorm12 = (double)wp_dat->Integral() / wp_sim12->Integral();
	//double full_simnorm19 = (double)wp_dat->Integral() / wp_sim19->Integral();
	wp_sim->Scale( full_simnorm );
	//wp_sim12->Scale( full_simnorm12 );
	//wp_sim19->Scale( full_simnorm19 );
	
	c_wp->cd(1);
	wp_sim->SetTitle(Form("Full distribution, E_{sim} > 5 , C_{sim} = %f",full_simnorm));
	label1D(wp_dat,wp_sim,"W' [GeV]","Counts");

	//c_wp->cd(2);
	//wp_sim12->SetTitle(Form("Full distribution, E_{sim} > 12 , C_{sim} = %f",full_simnorm12));
	//label1D(wp_dat,wp_sim12,"W' [GeV]","Counts");

	//c_wp->cd(3);
	//wp_sim19->SetTitle(Form("Full distribution, E_{sim} > 19 , C_{sim} = %f",full_simnorm19));
	//label1D(wp_dat,wp_sim19,"W' [GeV]","Counts");

	c_wp->SaveAs("full_wp.pdf");


	// Grab intereseted histograms
	TCanvas * c_wp_xp = new TCanvas("c_wp_xp","",800,600);
	c_wp_xp->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c_wp_xp->cd(bin+1);
		wp_xp_bins_dat[bin] = new TH1D(Form("wp_xP_bin_%i_dat",bin),"",13,1.7,3);
		wp_xp_bins_bac[bin] = new TH1D(Form("wp_xP_bin_%i_bac",bin),"",13,1.7,3);
		wp_xp_bins_sim[bin] = new TH1D(Form("wp_xP_bin_%i_sim",bin),"",13,1.7,3);
		//wp_xp_bins_sim12[bin] = new TH1D(Form("wp_xP_bin_%i_sim12",bin),"",13,1.7,3);
		//wp_xp_bins_sim19[bin] = new TH1D(Form("wp_xP_bin_%i_sim19",bin),"",13,1.7,3);

		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[nleadindex]->getXp() > %f && tag[nleadindex]->getXp() < %f",this_min_xp,this_max_xp);

		inTreeDat->Draw(Form("tag[nleadindex]->getWp() >> wp_xP_bin_%i_dat",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeBac->Draw(Form("tag[nleadindex]->getWp() >> wp_xP_bin_%i_bac",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeSim->Draw(Form("tag[nleadindex]->getWp() >> wp_xP_bin_%i_sim",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		//inTreeSim->Draw(Form("tag[nleadindex]->getWp() >> wp_xP_bin_%i_sim12",bin),this_xp_cut && "nHits[nleadindex]->getEdep()/1E4 > 12");
		//inTreeSim->Draw(Form("tag[nleadindex]->getWp() >> wp_xP_bin_%i_sim19",bin),this_xp_cut && "nHits[nleadindex]->getEdep()/1E4 > 19");
		
		if(	wp_xp_bins_dat[bin]->Integral() == 0 ||
			wp_xp_bins_bac[bin]->Integral() == 0 ||
			wp_xp_bins_sim[bin]->Integral() == 0 ) continue;
			//wp_xp_bins_sim12[bin]->Integral() == 0 ||
			//wp_xp_bins_sim19[bin]->Integral() == 0 ) continue;

		// Background subtraction
		wp_xp_bins_dat[bin] -> Add( wp_xp_bins_bac[bin] , -1 );
		// Scale simulation
		wp_xp_bins_sim[bin] -> Scale( full_simnorm );
		//wp_xp_bins_sim12[bin] -> Scale( full_simnorm12 );
		//wp_xp_bins_sim19[bin] -> Scale( full_simnorm19 );
		// 	calculate a re-normalization
		double simnorm = wp_xp_bins_dat[bin]->Integral() / wp_xp_bins_sim[bin]->Integral();
		//double simnorm12 = wp_xp_bins_dat[bin]->Integral() / wp_xp_bins_sim12[bin]->Integral();
		//double simnorm19 = wp_xp_bins_dat[bin]->Integral() / wp_xp_bins_sim19[bin]->Integral();


		TString current_title = Form("%f < xP < %f",this_min_xp,this_max_xp);
		wp_xp_bins_sim[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));
		//wp_xp_bins_sim12[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm12));
		//wp_xp_bins_sim19[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm19));

		label1D(wp_xp_bins_dat[bin],wp_xp_bins_sim[bin],"W' [GeV]","Counts");
		//label1D_sims(wp_xp_bins_sim12[bin],wp_xp_bins_sim19[bin],simnorm12,simnorm19);
	}
	c_wp_xp->SaveAs("wp_xp_bins.pdf");



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

void label1D_sims(TH1D* sim2, TH1D* sim3, double label1, double label2){
	sim2->SetLineColor(3);
	sim2->SetLineWidth(1);
	sim2->SetStats(0);

	sim3->SetLineColor(1);
	sim3->SetLineWidth(1);
	sim3->SetStats(0);

	sim2->Draw("hist,same");
	sim3->Draw("hist,same");

	TLegend * legend = new TLegend(0.6,0.7,0.9,0.8);
	TString title1 = Form("E_{sim} > 12, C_{new} = %f",label1);
        TString title2 = Form("E_{sim} > 19, C_{new} = %f",label2);

	legend->AddEntry(sim2,title1,"f");
	legend->AddEntry(sim3,title2,"f");
	legend->Draw();

	return;
}
