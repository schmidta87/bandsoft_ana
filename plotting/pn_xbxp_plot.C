void pn_xbxp_plot(TString inDat, TString inBac, TString inSim){

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
	TH1D * pn_dat = new TH1D("pn_dat","pn_dat",40,0.2,0.6);
	TH1D * pn_bac = new TH1D("pn_bac","pn_bac",40,0.2,0.6);
	TH1D * pn_sim = new TH1D("pn_sim","pn_sim",40,0.2,0.6);
	TH1D * pn_sim12 = new TH1D("pn_sim12","pn_sim12",40,0.2,0.6);
	TH1D * pn_sim19 = new TH1D("pn_sim19","pn_sim19",40,0.2,0.6);
	const int nXp_bins = 6;
	const double Xp_min = 0.1;
	const double Xp_max = 0.7;
	const int nXb_bins = 6;
	const double Xb_min = 0.1;
	const double Xb_max = 0.7;
	TH1D ** pn_xb_bins_dat = new TH1D*[nXb_bins];
	TH1D ** pn_xp_bins_dat = new TH1D*[nXp_bins];
	TH1D ** pn_xb_bins_bac = new TH1D*[nXb_bins];
	TH1D ** pn_xp_bins_bac = new TH1D*[nXp_bins];
	TH1D ** pn_xb_bins_sim = new TH1D*[nXb_bins];
	TH1D ** pn_xp_bins_sim = new TH1D*[nXp_bins];
	TH1D ** pn_xb_bins_sim12 = new TH1D*[nXb_bins];
	TH1D ** pn_xp_bins_sim12 = new TH1D*[nXp_bins];
	TH1D ** pn_xb_bins_sim19 = new TH1D*[nXb_bins];
	TH1D ** pn_xp_bins_sim19 = new TH1D*[nXp_bins];
	
	// Draw the full pn distribution
	TCanvas * c_pn = new TCanvas("c_pn","",800,600);
	c_pn->Divide(2,2);
	c_pn->cd(1);
	inTreeDat->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat");
	inTreeBac->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim12","nHits[nleadindex]->getEdep()/1E4 > 12");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim19","nHits[nleadindex]->getEdep()/1E4 > 19");

	pn_dat->Add(pn_bac,-1);
	double full_simnorm = (double)pn_dat->Integral() / pn_sim->Integral();
	pn_sim->Scale( full_simnorm );
	double full_simnorm12 = (double)pn_dat->Integral() / pn_sim12->Integral();
	pn_sim12->Scale( full_simnorm12 );
	double full_simnorm19 = (double)pn_dat->Integral() / pn_sim19->Integral();
	pn_sim19->Scale( full_simnorm19 );
	
	cout << full_simnorm << " " << full_simnorm12 << " " << full_simnorm19 << "\n";
	pn_sim->SetTitle(Form("Full distribution, E_{sim} > 5 , C_{sim} = %f",full_simnorm));
	label1D(pn_dat,pn_sim,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(2);
	pn_sim12->SetTitle(Form("Full distribution, E_{sim} > 12 , C_{sim} = %f",full_simnorm12));
	label1D(pn_dat,pn_sim12,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(3);
	pn_sim19->SetTitle(Form("Full distribution, E_{sim} > 19 , C_{sim} = %f",full_simnorm19));
	label1D(pn_dat,pn_sim19,"|p_{n}| [GeV/c]","Counts");

	c_pn->SaveAs("full_pn.pdf");


	// Grab intereseted histograms
	TCanvas * c_pn_xb = new TCanvas("c_pn_xb","",800,600);
	c_pn_xb->Divide(2,2);
	for( int bin = 0 ; bin < nXb_bins ; bin++ ){
		c_pn_xb->cd(bin+1);
		pn_xb_bins_dat[bin] = new TH1D(Form("pn_xB_bin_%i_dat",bin),"",40,0.2,0.6);
		pn_xb_bins_bac[bin] = new TH1D(Form("pn_xB_bin_%i_bac",bin),"",40,0.2,0.6);
		pn_xb_bins_sim[bin] = new TH1D(Form("pn_xB_bin_%i_sim",bin),"",40,0.2,0.6);
		pn_xb_bins_sim12[bin] = new TH1D(Form("pn_xB_bin_%i_sim12",bin),"",40,0.2,0.6);
		pn_xb_bins_sim19[bin] = new TH1D(Form("pn_xB_bin_%i_sim19",bin),"",40,0.2,0.6);

		double this_min_xb = Xb_min+0.05 + bin*(Xb_max - Xb_min)/nXb_bins;
		double this_max_xb = Xb_min+0.05 + (bin+1)*(Xb_max - Xb_min)/nXb_bins;
		TCut this_xb_cut = Form("eHit->getXb() > %f && eHit->getXb() < %f",this_min_xb,this_max_xb);

		inTreeDat->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xB_bin_%i_dat",bin),this_xb_cut);
		inTreeBac->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xB_bin_%i_bac",bin),this_xb_cut);
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xB_bin_%i_sim",bin),this_xb_cut);
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xB_bin_%i_sim12",bin),this_xb_cut && "nHits[nleadindex]->getEdep()/1E4 > 12");
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xB_bin_%i_sim19",bin),this_xb_cut && "nHits[nleadindex]->getEdep()/1E4 > 19");
		
		if(	pn_xb_bins_dat[bin]->Integral() == 0 ||
			pn_xb_bins_bac[bin]->Integral() == 0 ||
			pn_xb_bins_sim[bin]->Integral() == 0 || 
			pn_xb_bins_sim12[bin]->Integral() == 0 || 
			pn_xb_bins_sim19[bin]->Integral() == 0 ) continue;

		// Background subtraction
		pn_xb_bins_dat[bin] -> Add( pn_xb_bins_bac[bin] , -1 );
		// Scale simulation
		pn_xb_bins_sim[bin] -> Scale( full_simnorm );
		pn_xb_bins_sim12[bin] -> Scale( full_simnorm12 );
		pn_xb_bins_sim19[bin] -> Scale( full_simnorm19 );
		// 	calculate a re-normalization
		double simnorm = pn_xb_bins_dat[bin]->Integral() / pn_xb_bins_sim[bin]->Integral();
		double simnorm12 = pn_xb_bins_dat[bin]->Integral() / pn_xb_bins_sim12[bin]->Integral();
		double simnorm19 = pn_xb_bins_dat[bin]->Integral() / pn_xb_bins_sim19[bin]->Integral();


		TString current_title = Form("%f < xB < %f",this_min_xb,this_max_xb);
		pn_xb_bins_sim[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(pn_xb_bins_dat[bin],pn_xb_bins_sim[bin],"|p_{n}| [GeV/c]","Counts");
		label1D_sims(pn_xb_bins_sim12[bin],pn_xb_bins_sim19[bin],simnorm12,simnorm19);
	}
	c_pn_xb->SaveAs("pn_xb_bins.pdf");

	// Grab intereseted histograms
	TCanvas * c_pn_xp = new TCanvas("c_pn_xp","",800,600);
	c_pn_xp->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c_pn_xp->cd(bin+1);
		pn_xp_bins_dat[bin] = new TH1D(Form("pn_xP_bin_%i_dat",bin),"",40,0.2,0.6);
		pn_xp_bins_bac[bin] = new TH1D(Form("pn_xP_bin_%i_bac",bin),"",40,0.2,0.6);
		pn_xp_bins_sim[bin] = new TH1D(Form("pn_xP_bin_%i_sim",bin),"",40,0.2,0.6);
		pn_xp_bins_sim12[bin] = new TH1D(Form("pn_xP_bin_%i_sim12",bin),"",40,0.2,0.6);
		pn_xp_bins_sim19[bin] = new TH1D(Form("pn_xP_bin_%i_sim19",bin),"",40,0.2,0.6);

		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[nleadindex]->getXp() > %f && tag[nleadindex]->getXp() < %f",this_min_xp,this_max_xp);

		inTreeDat->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xP_bin_%i_dat",bin),this_xp_cut);
		inTreeBac->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xP_bin_%i_bac",bin),this_xp_cut);
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xP_bin_%i_sim",bin),this_xp_cut);
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xP_bin_%i_sim12",bin),this_xp_cut && "nHits[nleadindex]->getEdep()/1E4 > 12");
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumN().Mag() >> pn_xP_bin_%i_sim19",bin),this_xp_cut && "nHits[nleadindex]->getEdep()/1E4 > 19");
		
		if(	pn_xp_bins_dat[bin]->Integral() == 0 ||
			pn_xp_bins_bac[bin]->Integral() == 0 ||
			pn_xp_bins_sim[bin]->Integral() == 0 ||
			pn_xp_bins_sim12[bin]->Integral() == 0 || 
			pn_xp_bins_sim19[bin]->Integral() == 0 ) continue;

		// Background subtraction
		pn_xp_bins_dat[bin] -> Add( pn_xp_bins_bac[bin] , -1 );
		// Scale simulation
		pn_xp_bins_sim[bin] -> Scale( full_simnorm );
		pn_xp_bins_sim12[bin] -> Scale( full_simnorm12 );
		pn_xp_bins_sim19[bin] -> Scale( full_simnorm19 );
		// 	calculate a re-normalization
		double simnorm = pn_xp_bins_dat[bin]->Integral() / pn_xp_bins_sim[bin]->Integral();
		double simnorm12 = pn_xp_bins_dat[bin]->Integral() / pn_xp_bins_sim12[bin]->Integral();
		double simnorm19 = pn_xp_bins_dat[bin]->Integral() / pn_xp_bins_sim19[bin]->Integral();


		TString current_title = Form("%f < xP < %f",this_min_xp,this_max_xp);
		pn_xp_bins_sim[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));
		pn_xp_bins_sim12[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm12));
		pn_xp_bins_sim19[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm19));

		label1D(pn_xp_bins_dat[bin],pn_xp_bins_sim[bin],"|p_{n}| [GeV/c]","Counts");
		label1D_sims(pn_xp_bins_sim12[bin],pn_xp_bins_sim19[bin],simnorm12,simnorm19);
	}
	c_pn_xp->SaveAs("pn_xp_bins.pdf");



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
	legend->AddEntry(sim,"E_{sim}","f");
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
