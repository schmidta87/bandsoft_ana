void must_plot_virt(TString inFileDatTagName, TString inFileBacTagName,
	       	 TString inFileSimTagName, TString inFileDatIncName, 
		 TString inFileSimIncName ){
	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);



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
	const int nVirt_bins = 5;
	const double Virt_min = -0.45;
	const double Virt_max = -0.2;
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
	// These histograms all share a Q2 bin of [2,10]
	// The tag histograms have cos_theta_nq bin of [-1,-0.8]
	TH1D ** h1_virt_xp_dat = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for dat
	TH1D ** h1_virt_xp_bac = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for dat
	TH1D ** h1_virt_xp_sim = new TH1D*[nXp_bins]; // this is an alphaS histogram in a bin of X' for sim
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		h1_virt_xp_dat[bin]	= new TH1D(Form("h1_virt_xp_dat_%i",bin),"",nVirt_bins,Virt_min,Virt_max);
		h1_virt_xp_bac[bin]	= new TH1D(Form("h1_virt_xp_bac_%i",bin),"",nVirt_bins,Virt_min,Virt_max);
		h1_virt_xp_sim[bin]	= new TH1D(Form("h1_virt_xp_sim_%i",bin),"",nVirt_bins,Virt_min,Virt_max);
	}
	TH1D ** h1_q2_xb_dat = new TH1D*[10];
	TH1D ** h1_q2_xb_sim = new TH1D*[10];
	for( int bin = 0 ; bin < nXb_bins ; bin++ ){
		h1_q2_xb_dat[bin]	= new TH1D(Form("h1_q2_xb_dat_%i",bin),"",10,0,10);
		h1_q2_xb_sim[bin]	= new TH1D(Form("h1_q2_xb_sim_%i",bin),"",10,0,10);
	}

	// Draw alphaS distribution
	TH1D * virt_dat = new TH1D("virt_dat","",nVirt_bins,Virt_min,Virt_max);
	TH1D * virt_bac = new TH1D("virt_bac","",nVirt_bins,Virt_min,Virt_max);
	TH1D * virt_sim = new TH1D("virt_sim","",nVirt_bins,Virt_min,Virt_max);

	TCanvas * c_as = new TCanvas("c_as","",800,600);
	c_as->cd();

const double mP = 0.93827231;
	// virt^2 = [ [ (mD - En)^2 - (p_n)^2 ] - mP^2 ]/mP^2
	inTreeDatTag->Draw("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_dat",	"tag[nleadindex]->getMomentumN().Mag() > 0.3");
	inTreeBacTag->Draw("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_bac",	"tag[nleadindex]->getMomentumN().Mag() > 0.3");
	inTreeSimTag->Draw("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_sim",	"tag[nleadindex]->getMomentumN().Mag() > 0.3");

	virt_dat->Add(virt_bac,-1);
	double full_simnorm = (double)virt_dat->Integral() / virt_sim->Integral();
	virt_sim->Scale( full_simnorm );
	
	virt_dat->SetTitle(Form("Full distribution, C_{sim} = %f",full_simnorm));
	label1D(virt_dat,virt_sim,"Virt","Counts");
	c_as->SaveAs("full_virt.pdf");






	// First let's get:
	//	R(tag) = 	N_tag,data (Q2,theta_nq,x';alphaS) / [ N_tag,sim(Q2,theta_nq,x';alphaS) ]
	//	as just individual histograms for numerator and denominator
	TCanvas * c1_virt_xp_tag = new TCanvas("c1_virt_xp_tag","",800,600);
	c1_virt_xp_tag->Divide(3,2);
	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_virt_xp_tag->cd(bin+1);
		double this_min_xp = Xp_min+0.05 + bin*(Xp_max - Xp_min)/nXp_bins;
		double this_max_xp = Xp_min+0.05 + (bin+1)*(Xp_max - Xp_min)/nXp_bins;
		TCut this_xp_cut = Form("tag[nleadindex]->getXp() > %f && tag[nleadindex]->getXp() < %f",this_min_xp,this_max_xp);

		inTreeDatTag->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2)>> h1_virt_xp_dat_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeBacTag->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2)>> h1_virt_xp_bac_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");
		inTreeSimTag->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2)>> h1_virt_xp_sim_%i",bin),this_xp_cut && "tag[nleadindex]->getMomentumN().Mag() > 0.3");

		// Do the background subtraction:
		h1_virt_xp_dat[bin]->Add(h1_virt_xp_bac[bin],-1);

		// Do the normalization for simulation and data:
		h1_virt_xp_dat[bin]->Scale(1./Q_tag);
		h1_virt_xp_sim[bin]->Scale(1./L_tag);

		// Plot data and simulation:
		double new_simnorm = h1_virt_xp_dat[bin]->Integral() / h1_virt_xp_sim[bin]->Integral();
		TString current_title = Form("%f < x' < %f",this_min_xp,this_max_xp);
		h1_virt_xp_dat[bin]->SetTitle( current_title + Form(", C_{new} = %f",new_simnorm)  );
		label1D(h1_virt_xp_dat[bin],h1_virt_xp_sim[bin],"Virt","Counts");
	}
	c1_virt_xp_tag->SaveAs("virt_xpBins_tag.pdf");




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

		inTreeDatInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_dat_%i",bin),this_xb_cut);
		inTreeSimInc->Draw(Form("eHit->getQ2() >> h1_q2_xb_sim_%i",bin),this_xb_cut);

		// Do the normalization for simulation and data:
		h1_q2_xb_dat[bin]->Scale(1./Q_inc);
		h1_q2_xb_sim[bin]->Scale(1./L_inc);

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


		//h1_virt_xp_dat[bin] -> Scale( 1./Q_tag );  // N_data,tag
		//h1_virt_xp_sim[bin] -> Scale( 1./L_tag );	 // N_sim,tag

		double N_inc_dat = (DatIncCounts[bin]);
		double N_inc_sim = (SimIncCounts[bin]);

		h1_virt_xp_dat[bin]->Scale( N_inc_sim / N_inc_dat );	// N_dat,tag / N_dat,inc

		
		label1D_ratio(h1_virt_xp_dat[bin],h1_virt_xp_sim[bin],"Virt","R=R_{tag}/R_{inc}");

	}
	c1_must->SaveAs("ratio_tag_inc.pdf");

	TCanvas * c1_must_xp03 = new TCanvas("c1_must_xp03","",800,600);
	c1_must_xp03->Divide(3,2);

	for( int bin = 0 ; bin < nXp_bins ; bin++ ){
		c1_must_xp03->cd(bin+1);

		h1_virt_xp_dat[bin]->Divide( h1_virt_xp_sim[bin] );
		label1D_ratio(h1_virt_xp_dat[bin],h1_virt_xp_dat[1],"Virt","R=R_{tag}/R_{inc} / R_(x'=0.3)");

	}
	c1_must_xp03->SaveAs("ratio_tag_inc_xp03.pdf");

	
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


	for( int bin = 1 ; bin < data_copy->GetXaxis()->GetNbins()+1 ; bin++){
		cout << data_copy->GetXaxis()->GetBinCenter(bin) << " " << data_copy->GetBinContent(bin) << "\n";
	}

	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0.7,3);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
