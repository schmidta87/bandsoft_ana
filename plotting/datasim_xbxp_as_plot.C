void datasim_xbxp_as_plot(TString inDat, TString inBac, TString inSim){


	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax );
	void label1D_sims(TH1D* sim2, TH1D* sim3, double label1, double label2);
	void label1D_ratio_sims(TH1D* data, TH1D* sim2, TH1D* sim3 );

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
	// 	x' and xB:
	TH1D * xp_dat = new TH1D("xp_dat","xp_dat",50,0,1);
	TH1D * xp_bac = new TH1D("xp_bac","xp_bac",50,0,1);
	TH1D * xp_sim = new TH1D("xp_sim","xp_sim",50,0,1);
	TH1D * xp_sim12 = new TH1D("xp_sim12","xp_sim12",50,0,1);
	TH1D * xp_sim19 = new TH1D("xp_sim19","xp_sim19",50,0,1);

	TH1D * xb_dat = new TH1D("xb_dat","xb_dat",25,0,0.5);
	TH1D * xb_bac = new TH1D("xb_bac","xb_bac",25,0,0.5);
	TH1D * xb_sim = new TH1D("xb_sim","xb_sim",25,0,0.5);
	TH1D * xb_sim12 = new TH1D("xb_sim12","xb_sim12",25,0,0.5);
	TH1D * xb_sim19 = new TH1D("xb_sim19","xb_sim19",25,0,0.5);
	// 	x' and xB bins of As:
	const int nAs_bins = 4;
	const double As_min = 1.2;
	const double As_max = 1.6;
	TH1D ** xb_as_bins_dat = new TH1D*[nAs_bins];
	TH1D ** xb_as_bins_bac = new TH1D*[nAs_bins];
	TH1D ** xb_as_bins_sim = new TH1D*[nAs_bins];
	TH1D ** xb_as_bins_sim12 = new TH1D*[nAs_bins];
	TH1D ** xb_as_bins_sim19 = new TH1D*[nAs_bins];
	TH1D ** xp_as_lpt_bins_dat = new TH1D*[nAs_bins];
	TH1D ** xp_as_lpt_bins_bac = new TH1D*[nAs_bins];
	TH1D ** xp_as_lpt_bins_sim = new TH1D*[nAs_bins];
	TH1D ** xp_as_lpt_bins_sim12 = new TH1D*[nAs_bins];
	TH1D ** xp_as_lpt_bins_sim19 = new TH1D*[nAs_bins];
	TH1D ** xp_as_hpt_bins_dat = new TH1D*[nAs_bins];
	TH1D ** xp_as_hpt_bins_bac = new TH1D*[nAs_bins];
	TH1D ** xp_as_hpt_bins_sim = new TH1D*[nAs_bins];
	TH1D ** xp_as_hpt_bins_sim12 = new TH1D*[nAs_bins];
	TH1D ** xp_as_hpt_bins_sim19 = new TH1D*[nAs_bins];
	
	// Draw the full xp distribution
	TCanvas * c_xp = new TCanvas("c_xp","",800,600);
	c_xp->Divide(2,2);
	c_xp->cd(1);
	inTreeDat->Draw("tag[nleadindex]->getXp2() >> xp_dat");
	inTreeBac->Draw("tag[nleadindex]->getXp2() >> xp_bac");
	inTreeSim->Draw("tag[nleadindex]->getXp2() >> xp_sim");
	inTreeSim->Draw("tag[nleadindex]->getXp2() >> xp_sim12","nHits[nleadindex]->getEdep()/1E4 > 12");
	inTreeSim->Draw("tag[nleadindex]->getXp2() >> xp_sim19","nHits[nleadindex]->getEdep()/1E4 > 19");

	xp_dat->Add(xp_bac,-1);
	double full_simnorm = (double)xp_dat->Integral() / xp_sim->Integral();
	double full_simnorm12 = (double)xp_dat->Integral() / xp_sim12->Integral();
	double full_simnorm19 = (double)xp_dat->Integral() / xp_sim19->Integral();
	xp_sim->Scale( full_simnorm );
	xp_sim12->Scale( full_simnorm12 );
	xp_sim19->Scale( full_simnorm19 );
	
	c_xp->cd(1);
	xp_sim->SetTitle(Form("Full distribution, E_{sim} > 5 , C_{sim} = %f",full_simnorm));
	label1D(xp_dat,xp_sim,"x'","Counts");
	c_xp->cd(2);
	xp_sim12->SetTitle(Form("Full distribution, E_{sim} > 12 , C_{sim} = %f",full_simnorm12));
	label1D(xp_dat,xp_sim12,"x'","Counts");
	c_xp->cd(3);
	xp_sim19->SetTitle(Form("Full distribution, E_{sim} > 19 , C_{sim} = %f",full_simnorm19));
	label1D(xp_dat,xp_sim19,"x'","Counts");

	c_xp->SaveAs("full_xp.pdf");

	// Draw the full xp distribution RATIOS
	TCanvas * c_xp_ratio = new TCanvas("c_xp_ratio","",800,600);
	c_xp_ratio->Divide(2,2);

	c_xp_ratio->cd(1);
	label1D_ratio(xp_dat,xp_sim,"x'","Data/Sim",0,2);
	c_xp_ratio->cd(2);
	label1D_ratio(xp_dat,xp_sim12,"x'","Data/Sim",0,2);
	c_xp_ratio->cd(3);
	label1D_ratio(xp_dat,xp_sim19,"x'","Data/Sim",0,2);
	
	c_xp_ratio->SaveAs("full_xp_ratio.pdf");

	// Draw the full xb distribution
	TCanvas * c_xb = new TCanvas("c_xb","",800,600);
	c_xb->Divide(1,2);
	c_xb->cd(1);
	inTreeDat->Draw("eHit->getXb() >> xb_dat");
	inTreeBac->Draw("eHit->getXb() >> xb_bac");
	inTreeSim->Draw("eHit->getXb() >> xb_sim");
	inTreeSim->Draw("eHit->getXb() >> xb_sim");

	xb_dat->Add(xb_bac,-1);
	xb_sim->Scale( full_simnorm );
	
	xb_sim->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(xb_dat,xb_sim,"x_{B}","Counts");

	c_xb->cd(2);
	label1D_ratio(xb_dat,xb_sim,"x_{B}","Data/Sim",0,2);

	c_xb->SaveAs("full_xb.pdf");


	// Now draw the data/sim of x' distributions as a function of alphaS
	TCanvas * c_xp_as_lpt = new TCanvas("c_xp_as_lpt","",800,600);
	TCanvas * c_xp_as_hpt = new TCanvas("c_xp_as_hpt","",800,600);
	c_xp_as_lpt->Divide(4,2);
	c_xp_as_hpt->Divide(4,2);
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		xp_as_lpt_bins_dat[bin] = new TH1D(Form("xp_as_lpt_bins_%i_dat",bin),"",25,0,1);
		xp_as_lpt_bins_bac[bin] = new TH1D(Form("xp_as_lpt_bins_%i_bac",bin),"",25,0,1);
		xp_as_lpt_bins_sim[bin] = new TH1D(Form("xp_as_lpt_bins_%i_sim",bin),"",25,0,1);
		xp_as_lpt_bins_sim12[bin] = new TH1D(Form("xp_as_lpt_bins_%i_sim12",bin),"",25,0,1);
		xp_as_lpt_bins_sim19[bin] = new TH1D(Form("xp_as_lpt_bins_%i_sim19",bin),"",25,0,1);

		xp_as_hpt_bins_dat[bin] = new TH1D(Form("xp_as_hpt_bins_%i_dat",bin),"",25,0,1);
		xp_as_hpt_bins_bac[bin] = new TH1D(Form("xp_as_hpt_bins_%i_bac",bin),"",25,0,1);
		xp_as_hpt_bins_sim[bin] = new TH1D(Form("xp_as_hpt_bins_%i_sim",bin),"",25,0,1);
		xp_as_hpt_bins_sim12[bin] = new TH1D(Form("xp_as_hpt_bins_%i_sim12",bin),"",25,0,1);
		xp_as_hpt_bins_sim19[bin] = new TH1D(Form("xp_as_hpt_bins_%i_sim19",bin),"",25,0,1);


		double this_min_as = As_min + bin*(As_max - As_min)/nAs_bins;
		double this_max_as = As_min + (bin+1)*(As_max - As_min)/nAs_bins;
		TCut this_as_cut = Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",this_min_as,this_max_as);

		TCut lpt = "tag[nleadindex]->getPt().Mag() <  0.1 ";
		TCut hpt = "tag[nleadindex]->getPt().Mag() >= 0.1 ";

		c_xp_as_lpt->cd(bin+1);
		inTreeDat->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_lpt_bins_%i_dat",bin),this_as_cut && lpt);
		inTreeBac->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_lpt_bins_%i_bac",bin),this_as_cut && lpt);
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_lpt_bins_%i_sim",bin),this_as_cut && lpt);
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_lpt_bins_%i_sim12",bin),this_as_cut && lpt && "nHits[nleadindex]->getEdep()/1E4 > 12");
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_lpt_bins_%i_sim19",bin),this_as_cut && lpt && "nHits[nleadindex]->getEdep()/1E4 > 19");
		c_xp_as_hpt->cd(bin+1);
		inTreeDat->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_hpt_bins_%i_dat",bin),this_as_cut && hpt);
		inTreeBac->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_hpt_bins_%i_bac",bin),this_as_cut && hpt);
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_hpt_bins_%i_sim",bin),this_as_cut && hpt);
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_hpt_bins_%i_sim12",bin),this_as_cut && hpt && "nHits[nleadindex]->getEdep()/1E4 > 12");
		inTreeSim->Draw(Form("tag[nleadindex]->getXp2() >> xp_as_hpt_bins_%i_sim19",bin),this_as_cut && hpt && "nHits[nleadindex]->getEdep()/1E4 > 19");

		if(	xp_as_lpt_bins_dat[bin]->Integral() == 0 ||
			xp_as_lpt_bins_bac[bin]->Integral() == 0 ||
			xp_as_lpt_bins_sim[bin]->Integral() == 0 ||
			xp_as_hpt_bins_dat[bin]->Integral() == 0 ||
			xp_as_hpt_bins_bac[bin]->Integral() == 0 ||
			xp_as_hpt_bins_sim[bin]->Integral() == 0 ) continue;

		// Background subtraction
		xp_as_lpt_bins_dat[bin] -> Add( xp_as_lpt_bins_bac[bin] , -1 );
		xp_as_hpt_bins_dat[bin] -> Add( xp_as_hpt_bins_bac[bin] , -1 );
		// Scale simulation
		xp_as_lpt_bins_sim[bin] -> Scale( full_simnorm );
		xp_as_lpt_bins_sim12[bin] -> Scale( full_simnorm12 );
		xp_as_lpt_bins_sim19[bin] -> Scale( full_simnorm19 );
		xp_as_hpt_bins_sim[bin] -> Scale( full_simnorm );
		xp_as_hpt_bins_sim12[bin] -> Scale( full_simnorm12 );
		xp_as_hpt_bins_sim19[bin] -> Scale( full_simnorm19 );
		// 	calculate a re-normalization
		double simnorm_lpt = xp_as_lpt_bins_dat[bin]->Integral() / xp_as_lpt_bins_sim[bin]->Integral();
		double simnorm_lpt12 = xp_as_lpt_bins_dat[bin]->Integral() / xp_as_lpt_bins_sim12[bin]->Integral();
		double simnorm_lpt19 = xp_as_lpt_bins_dat[bin]->Integral() / xp_as_lpt_bins_sim19[bin]->Integral();
		double simnorm_hpt = xp_as_hpt_bins_dat[bin]->Integral() / xp_as_hpt_bins_sim[bin]->Integral();
		double simnorm_hpt12 = xp_as_hpt_bins_dat[bin]->Integral() / xp_as_hpt_bins_sim12[bin]->Integral();
		double simnorm_hpt19 = xp_as_hpt_bins_dat[bin]->Integral() / xp_as_hpt_bins_sim19[bin]->Integral();


		TString current_title_lpt = Form("p_{T}<0.1 ,%f < Alpha_{S} < %f",this_min_as,this_max_as);
		TString current_title_hpt = Form("p_{T}>=0.1 ,%f < Alpha_{S} < %f",this_min_as,this_max_as);
		xp_as_lpt_bins_sim[bin]->SetTitle(current_title_lpt + Form(", C_{new} = %f",simnorm_lpt));
		xp_as_hpt_bins_sim[bin]->SetTitle(current_title_hpt + Form(", C_{new} = %f",simnorm_hpt));

		c_xp_as_lpt->cd(bin+1);
		label1D(xp_as_lpt_bins_dat[bin],xp_as_lpt_bins_sim[bin],"x' ","Counts");
		label1D_sims(xp_as_lpt_bins_sim12[bin],xp_as_lpt_bins_sim19[bin],simnorm_lpt12,simnorm_lpt19);

		c_xp_as_lpt->cd(bin+5);
		label1D_ratio(xp_as_lpt_bins_dat[bin],xp_as_lpt_bins_sim[bin],"x'","Data/Sim",0,5);
		label1D_ratio_sims(xp_as_lpt_bins_dat[bin],xp_as_lpt_bins_sim12[bin],xp_as_lpt_bins_sim19[bin]);

		c_xp_as_hpt->cd(bin+1);
		label1D(xp_as_hpt_bins_dat[bin],xp_as_hpt_bins_sim[bin],"x' ","Counts");
		label1D_sims(xp_as_hpt_bins_sim12[bin],xp_as_hpt_bins_sim19[bin],simnorm_hpt12,simnorm_hpt19);

		c_xp_as_hpt->cd(bin+5);
		label1D_ratio(xp_as_hpt_bins_dat[bin],xp_as_hpt_bins_sim[bin],"x'","Data/Sim",0,5);
		label1D_ratio_sims(xp_as_hpt_bins_dat[bin],xp_as_hpt_bins_sim12[bin],xp_as_hpt_bins_sim19[bin]);
	}
	c_xp_as_lpt->SaveAs("xp_as_lpt_bins.pdf");
	c_xp_as_hpt->SaveAs("xp_as_hpt_bins.pdf");

	
	TCanvas * c_xp_as_lpt_xp03 = new TCanvas("c_xp_as_lpt_xp03","",800,600);
	TCanvas * c_xp_as_hpt_xp03 = new TCanvas("c_xp_as_hpt_xp03","",800,600);
	double lpt_pars_p0[4] = {0.8651,0.872,0.2845,0};
	double lpt_pars_p1[4] = {0.3974,0.600,2.661,0};
	c_xp_as_lpt_xp03->Divide(4,1);
	c_xp_as_hpt_xp03->Divide(4,1);
	for( int bin = 0 ; bin < nAs_bins ; bin++ ){
		c_xp_as_lpt_xp03->cd(bin+1);
	
		// Get counts at x'=0.3 for data and simulation
		double xp03_data_lpt = xp_as_lpt_bins_dat[bin]->GetBinContent( xp_as_lpt_bins_dat[bin]->GetXaxis()->FindBin(0.3)  );
		double xp03_simu_lpt = xp_as_lpt_bins_sim[bin]->GetBinContent( xp_as_lpt_bins_sim[bin]->GetXaxis()->FindBin(0.3)  );
		if( xp03_data_lpt == 0 || xp03_simu_lpt == 0 ) continue;
		xp_as_lpt_bins_dat[bin]->Scale( 1./xp03_data_lpt );
		xp_as_lpt_bins_sim[bin]->Scale( 1./xp03_simu_lpt );

		// Plot [ data/data(x'=0.3) ] / [ sim/sim(x'=0.3) ]
		label1D_ratio(xp_as_lpt_bins_dat[bin],xp_as_lpt_bins_sim[bin],"x'","[Data/Data(x'=0.3)]/[Sim/Sim(x'=0.3)]",0.8,3);


		c_xp_as_hpt_xp03->cd(bin+1);
	
		// Get counts at x'=0.3 for data and simulation
		double xp03_data = xp_as_hpt_bins_dat[bin]->GetBinContent( xp_as_hpt_bins_dat[bin]->GetXaxis()->FindBin(0.3)  );
		double xp03_simu = xp_as_hpt_bins_sim[bin]->GetBinContent( xp_as_hpt_bins_sim[bin]->GetXaxis()->FindBin(0.3)  );
		if( xp03_data == 0 || xp03_simu == 0 ) continue;
		xp_as_hpt_bins_dat[bin]->Scale( 1./xp03_data );
		xp_as_hpt_bins_sim[bin]->Scale( 1./xp03_simu );

		// Plot [ data/data(x'=0.3) ] / [ sim/sim(x'=0.3) ]
		label1D_ratio(xp_as_hpt_bins_dat[bin],xp_as_hpt_bins_sim[bin],"x'","[Data/Data(x'=0.3)]/[Sim/Sim(x'=0.3)]",0.8,3);
		TF1 * lpt_fit = new TF1(Form("lpt_fit_%i",bin),"pol1",0,1);
		lpt_fit->SetParameters(lpt_pars_p0[bin],lpt_pars_p1[bin]);
		lpt_fit->SetLineStyle(2);
		lpt_fit->SetLineColor(2);
		lpt_fit->Draw("same");

	}
	c_xp_as_lpt_xp03->SaveAs("xp_as_lpt_bins_normed.pdf");
	c_xp_as_hpt_xp03->SaveAs("xp_as_hpt_bins_normed.pdf");
	
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

	data_copy->Fit("pol1","QESR","",0.25,0.6);

	return;
}

void label1D_ratio_sims(TH1D* data, TH1D* sim2, TH1D* sim3 ){
	
	TH1D * data2_copy = (TH1D*) data->Clone();
	TH1D * data3_copy = (TH1D*) data->Clone();
	TH1D * sim2_copy = (TH1D*) sim2->Clone();
	TH1D * sim3_copy = (TH1D*) sim3->Clone();
	

	data2_copy->Divide(sim2_copy);
	data2_copy->SetLineColor(3);
	data2_copy->Draw("ep,same");

	data3_copy->Divide(sim3_copy);
	data3_copy->SetLineColor(1);
	data3_copy->Draw("ep,same");

	return;
}
