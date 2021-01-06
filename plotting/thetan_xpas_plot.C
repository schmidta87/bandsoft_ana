void thetan_xpas_plot(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);

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
	TH1D * thetan_dat = new TH1D("thetan_dat","thetan_dat",20,-1,-0.9);
	TH1D * thetan_bac = new TH1D("thetan_bac","thetan_bac",20,-1,-0.9);
	TH1D * thetan_sim = new TH1D("thetan_sim","thetan_sim",20,-1,-0.9);
	const int nXp_bins = 6;
	const double Xp_min = 0.1;
	const double Xp_max = 0.7;
	TH1D ** thetan_xp_as12_bins_dat = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as12_bins_bac = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as12_bins_sim = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as13_bins_dat = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as13_bins_bac = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as13_bins_sim = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as14_bins_dat = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as14_bins_bac = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as14_bins_sim = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as15_bins_dat = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as15_bins_bac = new TH1D*[nXp_bins];
	TH1D ** thetan_xp_as15_bins_sim = new TH1D*[nXp_bins];
	
	// Draw the full thetan distribution
	TCanvas * c_thetan = new TCanvas("c_thetan","",800,600);
	c_thetan->cd();
	inTreeDat->Draw("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_dat");
	inTreeBac->Draw("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_bac");
	inTreeSim->Draw("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_sim");

	thetan_dat->Add(thetan_bac,-1);
	double full_simnorm = (double)thetan_dat->Integral() / thetan_sim->Integral();
	thetan_sim->Scale( full_simnorm );
	
	thetan_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(thetan_dat,thetan_sim,"CosTheta_{n}","Counts");
	c_thetan->SaveAs("full_thetan.pdf");

	// Now draw cos(theta_n) for fixed alphaS and in bins of x':
	TCanvas * c_thetan_xp_as12 = new TCanvas("c_thetan_xp_as12","",800,600);
	TCanvas * c_thetan_xp_as13 = new TCanvas("c_thetan_xp_as13","",800,600);
	TCanvas * c_thetan_xp_as14 = new TCanvas("c_thetan_xp_as14","",800,600);
	TCanvas * c_thetan_xp_as15 = new TCanvas("c_thetan_xp_as15","",800,600);

	c_thetan_xp_as12->Divide(3,2);
	c_thetan_xp_as13->Divide(3,2);
	c_thetan_xp_as14->Divide(3,2);
	c_thetan_xp_as15->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){

		// Define histograms for this x' bin
		thetan_xp_as12_bins_dat[bin] = new TH1D(Form("thetan_xp_as12_bins_%i_dat",bin),"",20,-1,-0.9);
		thetan_xp_as12_bins_bac[bin] = new TH1D(Form("thetan_xp_as12_bins_%i_bac",bin),"",20,-1,-0.9);
		thetan_xp_as12_bins_sim[bin] = new TH1D(Form("thetan_xp_as12_bins_%i_sim",bin),"",20,-1,-0.9);
		thetan_xp_as13_bins_dat[bin] = new TH1D(Form("thetan_xp_as13_bins_%i_dat",bin),"",20,-1,-0.9);
		thetan_xp_as13_bins_bac[bin] = new TH1D(Form("thetan_xp_as13_bins_%i_bac",bin),"",20,-1,-0.9);
		thetan_xp_as13_bins_sim[bin] = new TH1D(Form("thetan_xp_as13_bins_%i_sim",bin),"",20,-1,-0.9);
		thetan_xp_as14_bins_dat[bin] = new TH1D(Form("thetan_xp_as14_bins_%i_dat",bin),"",20,-1,-0.9);
		thetan_xp_as14_bins_bac[bin] = new TH1D(Form("thetan_xp_as14_bins_%i_bac",bin),"",20,-1,-0.9);
		thetan_xp_as14_bins_sim[bin] = new TH1D(Form("thetan_xp_as14_bins_%i_sim",bin),"",20,-1,-0.9);
		thetan_xp_as15_bins_dat[bin] = new TH1D(Form("thetan_xp_as15_bins_%i_dat",bin),"",20,-1,-0.9);
		thetan_xp_as15_bins_bac[bin] = new TH1D(Form("thetan_xp_as15_bins_%i_bac",bin),"",20,-1,-0.9);
		thetan_xp_as15_bins_sim[bin] = new TH1D(Form("thetan_xp_as15_bins_%i_sim",bin),"",20,-1,-0.9);

		// Define the cut for this x' bin
		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[nleadindex]->getXp() > %f && tag[nleadindex]->getXp() < %f",this_min_xp,this_max_xp);

		TCut as12_cut = "tag[nleadindex]->getAs() > 1.2 && tag[nleadindex]->getAs() < 1.3";
		TCut as13_cut = "tag[nleadindex]->getAs() > 1.3 && tag[nleadindex]->getAs() < 1.4";
		TCut as14_cut = "tag[nleadindex]->getAs() > 1.4 && tag[nleadindex]->getAs() < 1.5";
		TCut as15_cut = "tag[nleadindex]->getAs() > 1.5 && tag[nleadindex]->getAs() < 1.6";

		// AlphaS 1.2-1.3
		c_thetan_xp_as12->cd(bin+1);
		inTreeDat->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as12_bins_%i_dat",bin),this_xp_cut && as12_cut);
		inTreeBac->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as12_bins_%i_bac",bin),this_xp_cut && as12_cut);
		inTreeSim->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as12_bins_%i_sim",bin),this_xp_cut && as12_cut);

		// Background subtraction
		thetan_xp_as12_bins_dat[bin] -> Add( thetan_xp_as12_bins_bac[bin] , -1 );
		// Scale simulation
		thetan_xp_as12_bins_sim[bin] -> Scale( full_simnorm );
		// 	calculate a re-normalization
		double simnorm = thetan_xp_as12_bins_dat[bin]->Integral() / thetan_xp_as12_bins_sim[bin]->Integral();

		TString current_title = Form("1.2<Alpha_{S}<1.3 , %f < xP < %f",this_min_xp,this_max_xp);
		thetan_xp_as12_bins_dat[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(thetan_xp_as12_bins_dat[bin],thetan_xp_as12_bins_sim[bin],"CosTheta_{n}' [GeV]","Counts");



		// AlphaS 1.3-1.4
		c_thetan_xp_as13->cd(bin+1);
		inTreeDat->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as13_bins_%i_dat",bin),this_xp_cut && as13_cut);
		inTreeBac->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as13_bins_%i_bac",bin),this_xp_cut && as13_cut);
		inTreeSim->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as13_bins_%i_sim",bin),this_xp_cut && as13_cut);

		// Background subtraction
		thetan_xp_as13_bins_dat[bin] -> Add( thetan_xp_as13_bins_bac[bin] , -1 );
		// Scale simulation
		thetan_xp_as13_bins_sim[bin] -> Scale( full_simnorm );
		// 	calculate a re-normalization
		simnorm = thetan_xp_as13_bins_dat[bin]->Integral() / thetan_xp_as13_bins_sim[bin]->Integral();

		current_title = Form("1.3<Alpha_{S}<1.4 , %f < xP < %f",this_min_xp,this_max_xp);
		thetan_xp_as13_bins_dat[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(thetan_xp_as13_bins_dat[bin],thetan_xp_as13_bins_sim[bin],"CosTheta_{n}' [GeV]","Counts");


		// AlphaS 1.4-1.5
		c_thetan_xp_as14->cd(bin+1);
		inTreeDat->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as14_bins_%i_dat",bin),this_xp_cut && as14_cut);
		inTreeBac->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as14_bins_%i_bac",bin),this_xp_cut && as14_cut);
		inTreeSim->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as14_bins_%i_sim",bin),this_xp_cut && as14_cut);

		// Background subtraction
		thetan_xp_as14_bins_dat[bin] -> Add( thetan_xp_as14_bins_bac[bin] , -1 );
		// Scale simulation
		thetan_xp_as14_bins_sim[bin] -> Scale( full_simnorm );
		// 	calculate a re-normalization
		simnorm = thetan_xp_as14_bins_dat[bin]->Integral() / thetan_xp_as14_bins_sim[bin]->Integral();

		current_title = Form("1.4<Alpha_{S}<1.5 , %f < xP < %f",this_min_xp,this_max_xp);
		thetan_xp_as14_bins_dat[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(thetan_xp_as14_bins_dat[bin],thetan_xp_as14_bins_sim[bin],"CosTheta_{n}' [GeV]","Counts");




		// AlphaS 1.5-1.6
		c_thetan_xp_as15->cd(bin+1);
		inTreeDat->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as15_bins_%i_dat",bin),this_xp_cut && as15_cut);
		inTreeBac->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as15_bins_%i_bac",bin),this_xp_cut && as15_cut);
		inTreeSim->Draw(Form("cos(tag[nleadindex]->getMomentumN().Theta()) >> thetan_xp_as15_bins_%i_sim",bin),this_xp_cut && as15_cut);

		// Background subtraction
		thetan_xp_as15_bins_dat[bin] -> Add( thetan_xp_as15_bins_bac[bin] , -1 );
		// Scale simulation
		thetan_xp_as15_bins_sim[bin] -> Scale( full_simnorm );
		// 	calculate a re-normalization
		simnorm = thetan_xp_as15_bins_dat[bin]->Integral() / thetan_xp_as15_bins_sim[bin]->Integral();

		current_title = Form("1.5<Alpha_{S}<1.6 , %f < xP < %f",this_min_xp,this_max_xp);
		thetan_xp_as15_bins_dat[bin]->SetTitle(current_title + Form(", C_{new} = %f",simnorm));

		label1D(thetan_xp_as15_bins_dat[bin],thetan_xp_as15_bins_sim[bin],"CosTheta_{n}' [GeV]","Counts");


	}
	c_thetan_xp_as12->SaveAs("thetan_xpbins_as12.pdf");
	c_thetan_xp_as13->SaveAs("thetan_xpbins_as13.pdf");
	c_thetan_xp_as14->SaveAs("thetan_xpbins_as14.pdf");
	c_thetan_xp_as15->SaveAs("thetan_xpbins_as15.pdf");
	



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
