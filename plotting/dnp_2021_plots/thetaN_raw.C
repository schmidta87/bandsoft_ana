void thetaN_raw(TString inputFileNameDat, TString inputFileNameMix){

	
	TFile inFileDat = TFile(inputFileNameDat);
	TFile inFileMix = TFile(inputFileNameMix);

	TH1D *ThetaN_dat = (TH1D*) inFileDat.Get("hThetaN");
	TH1D *ThetaN_mix = (TH1D*) inFileMix.Get("hThetaN");

	TVector3 *datnorm = (TVector3*) inFileDat.Get("bacnorm");
	TVector3 *bacnorm = (TVector3*) inFileMix.Get("bacnorm");

	ThetaN_mix->Scale( ThetaN_dat->Integral() / ThetaN_mix->Integral() );


	TCanvas * c = new TCanvas("c","",800,600);
	c->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);

	ThetaN_dat->Rebin(32);
	ThetaN_mix->Rebin(32);

	ThetaN_dat->SetMinimum(0);
	ThetaN_mix->SetMinimum(0);

	ThetaN_dat->GetXaxis()->SetLabelSize(0.05);
	ThetaN_dat->GetYaxis()->SetLabelSize(0.05);


	ThetaN_dat->SetTitle("");
	ThetaN_dat->GetXaxis()->SetTitle("ThetaN [ns]");
	ThetaN_dat->GetXaxis()->SetTitleSize(0.65);
	ThetaN_dat->GetYaxis()->SetTitle("Counts");
	ThetaN_dat->GetYaxis()->SetTitleSize(0.65);

	//ThetaN_dat->GetXaxis()->SetRangeUser(12,100);
	//ThetaN_mix->GetXaxis()->SetRangeUser(12,100);

	ThetaN_dat->SetLineWidth(3);
	ThetaN_mix->SetLineWidth(3);
	ThetaN_mix->SetLineColor(2);

	ThetaN_dat->SetStats(0);
	ThetaN_mix->SetStats(0);

	ThetaN_dat->Draw();
	ThetaN_mix->Draw("H,SAME");

	for(int i = 1 ; i <= ThetaN_dat->GetXaxis()->GetNbins(); ++i ){
		cout << ThetaN_dat->GetXaxis()->GetBinCenter(i) << " " << ThetaN_dat->GetBinContent(i) << " " << ThetaN_mix->GetBinContent(i) << "\n";
	}

	c->SaveAs("ThetaNHist.pdf");


	inFileDat.Close();
	inFileMix.Close();
	return;
}
