void tofpm_edep_study(TString inDat, TString inBac, TString inSim){

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
	TH1D * tofpm_dat 	= new TH1D("tofpm_dat","tofpm_dat"	,40,0,20);
	TH1D * tofpm_bac 	= new TH1D("tofpm_bac","tofpm_bac"	,40,0,20);
	TH1D * tofpm_sim5 	= new TH1D("tofpm_sim5","tofpm_sim5"	,40,0,20);
	TH1D * tofpm_sim8 	= new TH1D("tofpm_sim8","tofpm_sim8"	,40,0,20);
	TH1D * tofpm_sim10 	= new TH1D("tofpm_sim10","tofpm_sim10"	,40,0,20);
	TH1D * tofpm_sim12 	= new TH1D("tofpm_sim12","tofpm_sim12"	,40,0,20);
	TH1D * tofpm_sim15 	= new TH1D("tofpm_sim15","tofpm_sim15"	,40,0,20);
	TH1D * tofpm_sim20 	= new TH1D("tofpm_sim20","tofpm_sim20"	,40,0,20);

	// Draw the full tofpm distribution
	TCanvas * c_tofpm_prebs = new TCanvas("c_tofpm_prebs","",800,600);
	c_tofpm_prebs->cd(1);
	inTreeDat->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_dat");
	inTreeBac->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_bac");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim5",	"nHits[nleadindex]->getEdep()/1E4 > 5");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim8",	"nHits[nleadindex]->getEdep()/1E4 > 8");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim10",	"nHits[nleadindex]->getEdep()/1E4 > 10");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim12",	"nHits[nleadindex]->getEdep()/1E4 > 12");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim15",	"nHits[nleadindex]->getEdep()/1E4 > 15");
	inTreeSim->Draw("nHits[nleadindex]->getTof()/(nHits[nleadindex]->getDL().Mag()/100) >> tofpm_sim20",	"nHits[nleadindex]->getEdep()/1E4 > 20");
	
	label1D(tofpm_dat,tofpm_bac,"ToF/dL [ns/m]","Counts");
	c_tofpm_prebs->SaveAs("full_tofpm_prebs.pdf");

	// Draw the full tofpm distribution
	TCanvas * c_tofpm = new TCanvas("c_tofpm","",800,600);
	c_tofpm->Divide(2,3);

	c_tofpm->cd(1);
	tofpm_dat->Add(tofpm_bac,-1);
	double data_integral 	= tofpm_dat->Integral( 	 );
	double sim5_integral	= tofpm_sim5->Integral(  );
	double sim8_integral	= tofpm_sim8->Integral(  );
	double sim10_integral	= tofpm_sim10->Integral(  );
	double sim12_integral	= tofpm_sim12->Integral(  );
	double sim15_integral	= tofpm_sim15->Integral(  );
	double sim20_integral	= tofpm_sim20->Integral(  );

	double scale5 	= data_integral / sim5_integral ;
	double scale8 	= data_integral / sim8_integral ;
	double scale10 	= data_integral / sim10_integral;
	double scale12 	= data_integral / sim12_integral;
	double scale15 	= data_integral / sim15_integral;
	double scale20 	= data_integral / sim20_integral;

	tofpm_sim5->Scale(  scale5	);
	tofpm_sim8->Scale(  scale8	);
	tofpm_sim10->Scale( scale10	);
	tofpm_sim12->Scale( scale12	);
	tofpm_sim15->Scale( scale15	);
	tofpm_sim20->Scale( scale20	);
	
	c_tofpm->cd(1);
	tofpm_sim5	->SetTitle(Form("Full distribution, E_{sim} > 5 , C_{sim} = %f",scale5	));
	label1D(tofpm_dat,tofpm_sim5,"ToF/dL [ns/m]","Counts");

	c_tofpm->cd(2);
	tofpm_sim8	->SetTitle(Form("Full distribution, E_{sim} > 8 , C_{sim} = %f",scale8	));
	label1D(tofpm_dat,tofpm_sim8,"ToF/dL [ns/m]","Counts");

	c_tofpm->cd(3);
	tofpm_sim10->SetTitle(Form("Full distribution, E_{sim} > 10 , C_{sim} = %f",scale10	));
	label1D(tofpm_dat,tofpm_sim10,"ToF/dL [ns/m]","Counts");

	c_tofpm->cd(4);
	tofpm_sim12->SetTitle(Form("Full distribution, E_{sim} > 12 , C_{sim} = %f",scale12	));
	label1D(tofpm_dat,tofpm_sim12,"ToF/dL [ns/m]","Counts");

	c_tofpm->cd(5);
	tofpm_sim15->SetTitle(Form("Full distribution, E_{sim} > 15 , C_{sim} = %f",scale15	));
	label1D(tofpm_dat,tofpm_sim15,"ToF/dL [ns/m]","Counts");

	c_tofpm->cd(6);
	tofpm_sim20->SetTitle(Form("Full distribution, E_{sim} > 20 , C_{sim} = %f",scale20	));
	label1D(tofpm_dat,tofpm_sim20,"ToF/dL| [ns/m]","Counts");



	c_tofpm->SaveAs("full_tofpm_edep.pdf");




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
	legend->AddEntry(sim,"Background","f");
	legend->Draw();

	return;
}
