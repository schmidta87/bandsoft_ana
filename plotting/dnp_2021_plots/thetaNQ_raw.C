void thetaNQ_raw(TString inputFileNameDat, TString inputFileNameMix){

	
	TFile inFileDat = TFile(inputFileNameDat);
	TFile inFileMix = TFile(inputFileNameMix);

	TH1D *ThetaNQ_dat = (TH1D*) inFileDat.Get("hThetaNQ");
	TH1D *ThetaNQ_mix = (TH1D*) inFileMix.Get("hThetaNQ");

	TVector3 *datnorm = (TVector3*) inFileDat.Get("bacnorm");
	TVector3 *bacnorm = (TVector3*) inFileMix.Get("bacnorm");

	ThetaNQ_mix->Scale( ThetaNQ_dat->Integral() / ThetaNQ_mix->Integral() );


	TCanvas * c = new TCanvas("c","",800,600);
	c->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);

	ThetaNQ_dat->Rebin(64);
	ThetaNQ_mix->Rebin(64);

	ThetaNQ_dat->SetMinimum(0);
	ThetaNQ_mix->SetMinimum(0);

	ThetaNQ_dat->GetXaxis()->SetLabelSize(0.05);
	ThetaNQ_dat->GetYaxis()->SetLabelSize(0.05);


	ThetaNQ_dat->SetTitle("");
	ThetaNQ_dat->GetXaxis()->SetTitle("ThetaNQ [ns]");
	ThetaNQ_dat->GetXaxis()->SetTitleSize(0.65);
	ThetaNQ_dat->GetYaxis()->SetTitle("Counts");
	ThetaNQ_dat->GetYaxis()->SetTitleSize(0.65);

	//ThetaNQ_dat->GetXaxis()->SetRangeUser(12,100);
	//ThetaNQ_mix->GetXaxis()->SetRangeUser(12,100);

	ThetaNQ_dat->SetLineWidth(3);
	ThetaNQ_mix->SetLineWidth(3);
	ThetaNQ_mix->SetLineColor(2);

	ThetaNQ_dat->SetStats(0);
	ThetaNQ_mix->SetStats(0);

	ThetaNQ_dat->Draw();
	ThetaNQ_mix->Draw("H,SAME");

	for(int i = 1 ; i <= ThetaNQ_dat->GetXaxis()->GetNbins(); ++i ){
		cout << ThetaNQ_dat->GetXaxis()->GetBinCenter(i) << " " << ThetaNQ_dat->GetBinContent(i) << " " << ThetaNQ_mix->GetBinContent(i) << "\n";
	}

	c->SaveAs("ThetaNQHist.pdf");


	inFileDat.Close();
	inFileMix.Close();
	return;
}
