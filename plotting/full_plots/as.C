void as(TString inDat, TString inBac, TString inSim){

	TCut pNcut = "tag[nleadindex]->getMomentumN().Mag() < 0.32 && tag[nleadindex]->getMomentumN().Mag() > 0.25 && !(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1) && nHits[nleadindex]->getEdep()>10";
	TCut pNcut_sim = "tag_smeared[nleadindex]->getMomentumN().Mag() < 0.32 && tag_smeared[nleadindex]->getMomentumN().Mag() > 0.25 && !(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1) && nHits[nleadindex]->getEdep()>10";
	TCut pTcut[3] = {"tag[nleadindex]->getPt().Mag() < 0.2","tag[nleadindex]->getPt().Mag()< 0.1","tag[nleadindex]->getPt().Mag()>0.1 && tag[nleadindex]->getPt().Mag()<0.2"};
	TCut pTcut_sim[3] = {"tag_smeared[nleadindex]->getPt().Mag() < 0.2",
	     			"tag_smeared[nleadindex]->getPt().Mag()< 0.1",
				"tag_smeared[nleadindex]->getPt().Mag()>0.1 && tag_smeared[nleadindex]->getPt().Mag()<0.2"};

	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax );
	void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim );
	void simulation_weighting(TH1D* sim, double Ndata, double Nsim );

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
	//inTreeBac->SetWeight( datnorm->X() / bacnorm->X() );
	// normalization uncertainty of the background:
	TVector3 * fullnorm = (TVector3*)inFileDat->Get("bacnorm");
	double Cscale = fullnorm->Z(); // this is the same as bacnorm->X() for the data file
	double Sigma_Cscale = (fullnorm->Y() - fullnorm->X())/2.;  // this is the before-time and after-time levels
	double NB_sim = bacnorm->X();
	double Sigma_NB_sim = sqrt(NB_sim);
		// sigma_Cscale / Cscale ~ 7%
	double Ndata = 0;
	double Nsim = 0;
	double sim_scaling = 0;

	// Define histograms we want to plot:
	TH1D ** as_dat = new TH1D*[3];
	TH1D ** as_bac = new TH1D*[3];
	TH1D ** as_sim = new TH1D*[3];
	for(int i = 0 ; i < 3 ; i++){
		as_dat[i] = new TH1D(Form("as_dat_%i",i),"",20,1.2125,1.7125);
		as_bac[i] = new TH1D(Form("as_bac_%i",i),"",20,1.2125,1.7125);
		as_sim[i] = new TH1D(Form("as_sim_%i",i),"",20,1.2125,1.7125);
	}

	// Draw the full as distribution
	TCanvas * c_as = new TCanvas("c_as","",800,600);
	c_as->Divide(3,2);
	for( int i = 0 ; i < 3 ; i++){
		TString pTtitle = "Full pT";
		if( i == 1 ){
			pTtitle = "Low pT";
		}
		if( i == 2 ){
			pTtitle = "High pT";
		}

		c_as->cd(i+1);
		inTreeDat->Draw(Form("tag[nleadindex]->getAs() >> as_dat_%i",i),pNcut && pTcut[i]);
		inTreeBac->Draw(Form("tag[nleadindex]->getAs() >> as_bac_%i",i),pNcut && pTcut[i]);
		inTreeSim->Draw(Form("tag_smeared[nleadindex]->getAs() >> as_sim_%i",i),pNcut_sim && pTcut_sim[i]);

		// Background subraction
		background_subtraction( as_dat[i] , as_bac[i] , Cscale, NB_sim, Sigma_Cscale, Sigma_NB_sim );

		if( i == 0 ){
		}
		Ndata = as_dat[i]->Integral(); 		// total data we have
		Nsim = as_sim[i]->Integral(); 		// total simulation we have
		sim_scaling = Ndata/Nsim;		// scale of the simulation bin
		// Simulation weighting
		simulation_weighting( as_sim[i], Ndata, Nsim );
		
		cout << "Pt bin: " << i << " " << pTtitle << "\n";
		as_sim[i]->SetTitle(pTtitle+Form(", C_{sim} = %f, ",sim_scaling));
		label1D(as_dat[i],as_sim[i],"Alpha_{S}","Counts");

		c_as->cd(4+i);
		label1D_ratio(as_dat[i],as_sim[i],"Alpha_{S}","Data/Sim",0,2);
	}


	c_as->SaveAs("full_as.pdf");

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

	for( int bin = 1 ; bin < data->GetXaxis()->GetNbins() ; ++bin ){
		cerr << "Data: " << data->GetBinCenter(bin) << " " << data->GetBinContent(bin) << " " << data->GetBinError(bin) << "\n";
		cerr << "Sim: " << sim->GetBinCenter(bin) << " " << sim->GetBinContent(bin) << " " << sim->GetBinError(bin) << "\n";
	}

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
	//data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	//sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	for( int bin = 1 ; bin < data_copy->GetXaxis()->GetNbins() ; ++bin ){
		double xval = data_copy->GetBinCenter(bin);

		double D = data_copy->GetBinContent(bin);
		double S = sim_copy->GetBinContent(bin);

		double Derr = data_copy->GetBinError(bin);
		double Serr = data_copy->GetBinError(bin);

		double ratio = D/S;
		double error = sqrt( pow(Derr/S,2) + pow(D/(S*S)*Serr,2) );

		if( S == 0 ){
			data_copy->SetBinContent(bin,-1);
			data_copy->SetBinError(bin,0);
			cerr << "Data/Sim: " << xval << " " << -1 << " " << 0 << "\n";
		}
		else{
			data_copy->SetBinContent(bin,ratio);
			data_copy->SetBinError(bin,error);
			cerr << "Data/Sim: " << xval << " " << ratio << " " << error << "\n";
		}
	}

	//data_copy->Divide(sim_copy);
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


void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim ){
	for( int bin = 1 ; bin < dat->GetXaxis()->GetNbins(); bin++ ){
		double SpB_bin = dat->GetBinContent(bin);
		double NB_bin = bac->GetBinContent(bin); // before any re-weighting

		double Sigma_NB_bin = sqrt(NB_bin);
		//double Sigma_Cscale already defined above
		//double Sigma_NB_sim already defined above

		double i1 = (Cscale/NB_sim)*Sigma_NB_bin;
		double i2 = (NB_bin/NB_sim)*Sigma_Cscale;
		double i3 = (NB_bin*Cscale)/pow(NB_sim,2)*Sigma_NB_sim;
		double Sigma_B = sqrt(i1*i1 + i2*i2 + i3*i3);

		double Sigma = sqrt( SpB_bin + pow(Sigma_B,2) ); // total uncertainty per bin

		dat->SetBinContent(bin, SpB_bin - NB_bin*Cscale/NB_sim );
		dat->SetBinError(bin, Sigma );
	}



	return;
}
void simulation_weighting(TH1D* sim, double Ndata, double Nsim ){
	for( int bin = 1 ; bin < sim->GetXaxis()->GetNbins(); bin++ ){
		double Nsim_bin = sim->GetBinContent(bin); // simulation stats in bin before re-weight

		double Sigma_Nsim_bin = sqrt(Nsim_bin);
		double Sigma_Ndata = sqrt(Ndata);
		double Sigma_Nsim = sqrt(Nsim);
		double i1 = (Ndata/Nsim)*Sigma_Nsim_bin;
		double i2 = (Nsim_bin/Nsim)*Sigma_Ndata;
		double i3 = (Nsim_bin*Ndata)/pow(Nsim,2)*Sigma_Nsim;

		double Sigma = sqrt(i1*i1 + i2*i2 + i3*i3);

		sim->SetBinContent(bin, Nsim_bin * Ndata/Nsim );
		sim->SetBinError(bin, Sigma );
	}


	return;
}
