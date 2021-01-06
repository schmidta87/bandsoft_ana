void pn_edep_study(TString inDat, TString inBac, TString inSim){

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
	TH1D * pn_dat = new TH1D("pn_dat","pn_dat",40,0.2,0.6);
	TH1D * pn_bac = new TH1D("pn_bac","pn_bac",40,0.2,0.6);
	TH1D * pn_sim5 	= new TH1D("pn_sim5","pn_sim5",40,0.2,0.6);
	TH1D * pn_sim8 	= new TH1D("pn_sim8","pn_sim8",40,0.2,0.6);
	TH1D * pn_sim10 = new TH1D("pn_sim10","pn_sim10",40,0.2,0.6);
	TH1D * pn_sim12 = new TH1D("pn_sim12","pn_sim12",40,0.2,0.6);
	TH1D * pn_sim15 = new TH1D("pn_sim15","pn_sim15",40,0.2,0.6);
	TH1D * pn_sim20 = new TH1D("pn_sim20","pn_sim20",40,0.2,0.6);
	
	// Draw the full pn distribution
	TCanvas * c_pn = new TCanvas("c_pn","",800,600);
	c_pn->Divide(2,3);
	c_pn->cd(1);
	inTreeDat->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat");
	inTreeBac->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim5",	"nHits[nleadindex]->getEdep()/1E4 > 5");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim8",	"nHits[nleadindex]->getEdep()/1E4 > 8");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim10",	"nHits[nleadindex]->getEdep()/1E4 > 10");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim12",	"nHits[nleadindex]->getEdep()/1E4 > 12");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim15",	"nHits[nleadindex]->getEdep()/1E4 > 15");
	inTreeSim->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim20",	"nHits[nleadindex]->getEdep()/1E4 > 20");

	pn_dat->Add(pn_bac,-1);
	double data_integral = pn_dat->Integral( pn_dat->GetXaxis()->FindBin(0.3) , pn_dat->GetXaxis()->FindBin(0.7) );
	double sim5_integral	= pn_sim5->Integral( pn_sim5->GetXaxis()->FindBin(0.3) , pn_sim5->GetXaxis()->FindBin(0.7) );
	double sim8_integral	= pn_sim8->Integral( pn_sim8->GetXaxis()->FindBin(0.3) , pn_sim8->GetXaxis()->FindBin(0.7) );
	double sim10_integral	= pn_sim10->Integral( pn_sim10->GetXaxis()->FindBin(0.3) , pn_sim10->GetXaxis()->FindBin(0.7) );
	double sim12_integral	= pn_sim12->Integral( pn_sim12->GetXaxis()->FindBin(0.3) , pn_sim12->GetXaxis()->FindBin(0.7) );
	double sim15_integral	= pn_sim15->Integral( pn_sim15->GetXaxis()->FindBin(0.3) , pn_sim15->GetXaxis()->FindBin(0.7) );
	double sim20_integral	= pn_sim20->Integral( pn_sim20->GetXaxis()->FindBin(0.3) , pn_sim20->GetXaxis()->FindBin(0.7) );

	double scale5 	= data_integral / sim5_integral ;
	double scale8 	= data_integral / sim8_integral ;
	double scale10 	= data_integral / sim10_integral;
	double scale12 	= data_integral / sim12_integral;
	double scale15 	= data_integral / sim15_integral;
	double scale20 	= data_integral / sim20_integral;

	pn_sim5->Scale(  scale5	);
	pn_sim8->Scale(  scale8	);
	pn_sim10->Scale( scale10	);
	pn_sim12->Scale( scale12	);
	pn_sim15->Scale( scale15	);
	pn_sim20->Scale( scale20	);
	
	c_pn->cd(1);
	pn_sim5	->SetTitle(Form("Full distribution, E_{sim} > 5 , C_{sim} = %f",scale5	));
	label1D(pn_dat,pn_sim5,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(2);
	pn_sim8	->SetTitle(Form("Full distribution, E_{sim} > 8 , C_{sim} = %f",scale8	));
	label1D(pn_dat,pn_sim8,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(3);
	pn_sim10->SetTitle(Form("Full distribution, E_{sim} > 10 , C_{sim} = %f",scale10	));
	label1D(pn_dat,pn_sim10,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(4);
	pn_sim12->SetTitle(Form("Full distribution, E_{sim} > 12 , C_{sim} = %f",scale12	));
	label1D(pn_dat,pn_sim12,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(5);
	pn_sim15->SetTitle(Form("Full distribution, E_{sim} > 15 , C_{sim} = %f",scale15	));
	label1D(pn_dat,pn_sim15,"|p_{n}| [GeV/c]","Counts");

	c_pn->cd(6);
	pn_sim20->SetTitle(Form("Full distribution, E_{sim} > 20 , C_{sim} = %f",scale20	));
	label1D(pn_dat,pn_sim20,"|p_{n}| [GeV/c]","Counts");



	c_pn->SaveAs("full_pn_edep.pdf");




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
