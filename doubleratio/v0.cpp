#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TClonesArray.h"

#include "constants.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"
//#include "clas12fiducial.h"

// For processing data

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outfile.root] [apply fiducial (1, 0)] [inputDatafiles]\n";
		return -1;
	}


	// Or wants to see:
	// 	Pn data, sim in bins of x and x'
	// 	Let's look at kinematics in bins of x and x' (Q2, W') data and sim etc (1D and 2D plots?)
	// 	Limit Q2/W' of double ratio plot and see what happens


	const double AdcToMeVee = 2300;
	const double MeVee_cut = 5;
	const double Q2_bin_min = 2;
	const double Q2_bin_max = 10;
	const double CosThetaNQ_bin_min = -1;
	const double CosThetaNQ_bin_max = -0.8;
	const double NMomentum_min = 0.3;
	const double NMomentum_max = 0.6;
	const double Wp_min = 1.8;
	const double Wp_max = 4.5;
	const double Al_min = 1.2;
	const double Al_max = 1.6;
	const double Al_bin_width = 0.05;
	int NAl_bins = (Al_max-Al_min + Al_bin_width/2.)/Al_bin_width;
	const double pt_min = 0;
	const double pt_max = 0.5;
	int pt_bins = 20;

	const double Xp_min = 0.1;
	const double Xp_max = 0.7;
	const double Xp_bin_width = 0.1;
	int NXp_bins = (Xp_max - Xp_min + Xp_bin_width/2.)/Xp_bin_width;

	// Output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");

	// ToF histogram for background normalization
	TH1D * hToF_bac = new TH1D("hToF_bac","hToF_bac",1000,-25,75);
	TH1D * hXp = new TH1D("hXp","hXp",1000,0,1);
	TH1D * hXb = new TH1D("hXb","hXb",1000,0,1);
	TH1D * hQ2 = new TH1D("hQ2","hQ2",1000,0,10);
	TH1D * hPn = new TH1D("hPn","hPn",1000,0.0,0.65);
	TH2D * h2XpQ2 = new TH2D("h2XpQ2","h2XpQ2",100,0,10,100,0,1);
	TH2D * h2Q2Wp_hi = new TH2D("h2Q2Wp_hi","h2Q2Wp_hi",100,0,10,100,0,4);
	TH2D * h2Q2Wp_lo = new TH2D("h2Q2Wp_lo","h2Q2Wp_lo",100,0,10,100,0,4);
	TH2D * h2XpAs = new TH2D("h2XpAs","h2XpAs",6,0.2,0.8,5,1.2,1.6);

	TH1D * hWp = new TH1D("hWp","hWp",1000,0,10);
	TH1D * hAs = new TH1D("hAs","hAs",600,1.2,1.6);
	TH1D ** hXpBins = new TH1D*[NAl_bins];
	for( int i = 0 ; i < NAl_bins ; i++){
		hXpBins[i] = new TH1D(Form("hXpBins_%i",i),Form("hXpBins_%i",i),100,0,1);
	}

	TH1D * hAs_hi = new TH1D("hAs_hi","hAs_hi",NAl_bins,Al_min,Al_max);
	TH1D * hAs_lo = new TH1D("hAs_lo","hAs_lo",NAl_bins,Al_min,Al_max);

	// AlphaS histograms for bins in X'
	TH1D ** hAs_bins = new TH1D*[NXp_bins];
	for( int i = 0 ; i < NXp_bins ; i++){
		hAs_bins[i] = new TH1D(Form("hAs_bins_%i",i),Form("hAs_bins_%i",i),16,1.2,1.6);
	}
	// Pn histograms for bins in X' and Xb
	TH1D ** hPn_xp = new TH1D*[NXp_bins];
	TH1D ** hPn_xb = new TH1D*[NXp_bins];
	for( int i = 0 ; i < NXp_bins ; i++){
		hPn_xp[i] = new TH1D(Form("hPn_xp_%i",i),Form("hPn_xp_%i",i),80,0.25,0.65);
		hPn_xb[i] = new TH1D(Form("hPn_xb_%i",i),Form("hPn_xb_%i",i),80,0.25,0.65);
	}


	// Neutron momentum binned histograms
	const double NMomentum_bin_min = 0.20;
	const double NMomentum_bin_max = 0.50;
	const double NMomentum_bin_width = 0.1;
	const int NMomentum_bins = (NMomentum_bin_max-NMomentum_bin_min + NMomentum_bin_width/2.)/NMomentum_bin_width;

	TH1D ** hXb_mom_bins = new TH1D*[NMomentum_bins];
	TH1D ** hXp_mom_bins = new TH1D*[NMomentum_bins];

	for( int i = 0 ; i < NMomentum_bins ; i++){
		hXb_mom_bins[i] = new TH1D(Form("hXb_mom_bins_%i",i),Form("hXb_mom_bins_%i",i),100,0,1);
		hXp_mom_bins[i] = new TH1D(Form("hXp_mom_bins_%i",i),Form("hXp_mom_bins_%i",i),100,0,1);
	}

	// AlphaS binned histograms
	const double alphaS_bin_min = 1.2;
	const double alphaS_bin_max = 1.6;
	const double alphaS_bin_width = 0.1;
	const int NalphaS_bins = (alphaS_bin_max - alphaS_bin_min + alphaS_bin_width/2.)/alphaS_bin_width;

	TH1D ** hXb_aS_bins = new TH1D*[NalphaS_bins];
	TH1D ** hXp_aS_bins = new TH1D*[NalphaS_bins];

	for( int i = 0 ; i < NalphaS_bins ; i++){
		hXb_aS_bins[i] = new TH1D(Form("hXb_aS_bins_%i",i),Form("hXb_aS_bins_%i",i),100,0,1);
		hXp_aS_bins[i] = new TH1D(Form("hXp_aS_bins_%i",i),Form("hXp_aS_bins_%i",i),100,0,1);
	}

	TH2D * h2AsVi = new TH2D("h2AsVi","h2AsVi",50,1,2,50,-1,0);
	TH2D * h2AsPn = new TH2D("h2AsPn","h2AsPn",50,1,2,50,0,1);
	TH2D * h2ViPn = new TH2D("h2ViPn","h2ViPn",50,-1,0,50,0,1);
	// AlphaS versus pt
	TH2D * h2Aspt = new TH2D("h2Aspt","h2Aspt",NAl_bins,Al_min,Al_max,pt_bins,pt_min,pt_max);


	//Store data for 4D phase space Q2,W,alpha_s and pt (pt is one bin for now so 3D histogram)
	const int datastore_Q2_bin = 1;
	double datastore_Q2_limits[datastore_Q2_bin+1] = {2,10} ;
	const int datastore_Wp_bin = 2;
	double datastore_Wp_limits[datastore_Wp_bin+1] = {2, 3, 4.5} ;
	const int datastore_as_bin = 2;
	double datastore_as_limits[datastore_as_bin+1] = {1.3, 1.4, 1.6} ;

	TH3D * h3_datastore_alphaS_Wp_Q = new TH3D("h3_datastore_alphaS_Wp_Q","h3_datastore_alphaS_Wp_Q",datastore_Q2_bin,datastore_Q2_limits,datastore_Wp_bin,datastore_Wp_limits,datastore_as_bin,datastore_as_limits);

/*
	// virtuality binned histograms
	const double virt_bin_min = 1.2;
	const double virt_bin_max = 1.6;
	const double virt_bin_width = 0.1;
	const double Nvirt_bins = (virt_bin_max - virt_bin_min + virt_bin_width/2.)/virt_bin_width;

	TH1D ** hXb_virt_bins = new TH1D*[Nvirt_bins];
	TH1D ** hXp_virt_bins = new TH1D*[Nvirt_bins];

	for( int i = 0 ; i < Nvirt_bins ; i++){
		hXb_virt_bins[i] = new TH1D(Form("hXb_virt_bins_%i",i),Form("hXb_virt_bins_%i",i),100,0,1);
		hXp_virt_bins[i] = new TH1D(Form("hXp_virt_bins_%i",i),Form("hXp_virt_bins_%i",i),100,0,1);
	}
*/


	//int doFiducial = atoi(argv[2]);
	//clas12fiducial* fid = new clas12fiducial();

	// Loop over all the files that are given to me to get the best statistics per bar
	for( int i = 3 ; i < argc ; i++ ){

		TFile * inFile = new TFile(argv[i]);
		TTree * inTree;

		if (inFile->GetListOfKeys()->Contains("T")){
			inTree = (TTree*)inFile->Get("T");
		} else if (inFile->GetListOfKeys()->Contains("tagged")) {
			inTree = (TTree*)inFile->Get("tagged");
		} else {
			cerr << "File has no entries\n";
			return -2;
		}

		// Initialize the input branches
		int Runno		= 0;
		double Ebeam		= 0;
		double gated_charge	= 0;
		double livetime		= 0;
		double starttime	= 0;
		double current		= 0;
		int nMult		= 0;
		TClonesArray* nHit = new TClonesArray("bandhit");
		clashit* eHit = new clashit;
		TClonesArray* tag  = new TClonesArray("taghit");
		inTree->SetBranchAddress("Runno"		,&Runno			);
		inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
		inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
		inTree->SetBranchAddress("livetime"		,&livetime		);
		inTree->SetBranchAddress("starttime"		,&starttime		);
		inTree->SetBranchAddress("current"		,&current		);
		//	Neutron branches:
		inTree->SetBranchAddress("nMult"		,&nMult			);
		inTree->SetBranchAddress("nHits"		,&nHit			);
		// 	Electron branches:
		inTree->SetBranchAddress("eHit"			,&eHit			);
		//	Tagg branches:
		inTree->SetBranchAddress("tag"			,&tag			);



		// Start working on one of the files, looping over all of the events
		cout << "Working on file: " << argv[i] << "\n";
		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 100000 == 0 ) cout << "\ton event " << ev << "\n";

			// Clear all branches before getting the entry from tree
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			nMult		= 0;

			// Not necessary after switch to TClonesArray
//			for( int thishit = 0; thishit < maxNeutrons ; thishit++){
//				nHit[thishit].Clear();
//				tag[thishit].Clear();
//			}

			eHit->Clear();

			inTree->GetEntry(ev);

			// Check neutron multiplicity
			if( nMult != 1 ) continue;

			// Get band and tag hit from clones array
			bandhit* this_nHit = (bandhit*)nHit->At(0);
			taghit* this_tag = (taghit*)tag->At(0);

			if( this_nHit->getStatus() != 0 ) continue;
			if( this_nHit->getTofFadc() == 0 ) continue;
			if( this_nHit->getEdep() < AdcToMeVee*MeVee_cut ) continue;

			// Check electron information
			if( eHit->getPID() != 11 ) continue;
			if( eHit->getCharge() != -1 ) continue;
			if( eHit->getEoP() < 0.17 ) continue;
			if( eHit->getEoP() > 0.3 ) continue;
			if( eHit->getEpcal() < 0.07 ) continue;
			if( eHit->getV() < 9 ) continue;
			if( eHit->getW() < 9 ) continue;
			if( eHit->getVtz() < -8 ) continue;
			if( eHit->getVtz() > 3 ) continue;
			if( eHit->getMomentum() < 2. ) continue;

			//if (doFiducial) {
			//	int eSect = fid->GetElectronAcceptance(eHit->getTheta()*TMath::RadToDeg(), eHit->getPhi()*TMath::RadToDeg(), eHit->getMomentum());
			//	if( eSect < 0 ) continue;
			//}

			// Define our 2D bins in Q2 and Theta_nq
			if( eHit->getQ2() < Q2_bin_min ) continue;
			if( eHit->getQ2() > Q2_bin_max ) continue;
			if( cos(this_tag->getThetaNQ()) < CosThetaNQ_bin_min ) continue;
			if( cos(this_tag->getThetaNQ()) > CosThetaNQ_bin_max ) continue;
			if( eHit->getW2() < 2*2 ) continue;

			// Fill the TOF histogram to extract the background normalization:
			hToF_bac -> Fill( (this_nHit->getTofFadc())/(this_nHit->getDL().Mag()/100.) );

			// Now only look at neutrons in our signal region:
			if( this_nHit->getTofFadc() < 0. ) continue;
			if( this_tag->getMomentumN().Mag() > NMomentum_max ) continue;
			if( this_tag->getMomentumN().Mag() < NMomentum_min ) continue;
			if( this_tag->getMomentumN().Mag() != this_tag->getMomentumN().Mag() ) continue; // check if NaN
			if( this_tag->getWp() < Wp_min ) continue;
			if( this_tag->getWp() > Wp_max ) continue;
			if( this_tag->getAs() < Al_min ) continue;
			if( this_tag->getAs() > Al_max ) continue;


			// Fill full Xp,As distributions
			h2XpAs->Fill( this_tag->getXp() , this_tag->getAs() );


			int binXp = ( this_tag->getXp() - Xp_min )/Xp_bin_width;
			if( binXp > -1 && this_tag->getXp() < Xp_max){
				hAs_bins[binXp]->Fill( this_tag->getAs() );
				hPn_xp[binXp]->Fill(this_tag->getMomentumN().Mag());
			}
			int binXb = ( eHit->getXb() - Xp_min)/Xp_bin_width;
			if( binXb > -1 && eHit->getXb() < Xp_max){
				hPn_xb[binXb]->Fill(this_tag->getMomentumN().Mag());
			}


			// Fill Hi and Lo Xp phase space in Q2,W'
			if( this_tag->getXp() < 0.35 && this_tag->getXp() > 0.25)
				h2Q2Wp_lo->Fill( eHit->getQ2() , this_tag->getWp() );
			else if( this_tag->getXp() > 0.5 )
				h2Q2Wp_hi->Fill( eHit->getQ2() , this_tag->getWp() );


			hXp->Fill( this_tag->getXp() );
			hXb->Fill( eHit->getXb() );
			hAs->Fill( this_tag->getAs() );
			hWp->Fill( this_tag->getWp() );
			hQ2->Fill( eHit->getQ2()  );
			hPn->Fill( this_tag->getMomentumN().Mag());
			h2XpQ2->Fill( eHit->getQ2() , this_tag->getXp() );

			// Fill Xp distribution for bin in alpha_s
			int binAl = ( this_tag->getAs() - Al_min )/Al_bin_width ;
			hXpBins[binAl]->Fill( this_tag->getXp() );


			double thisPn = this_tag->getMomentumN().Mag();
			if( thisPn > NMomentum_bin_min && thisPn < NMomentum_bin_max) {
				int binPn = (this_tag->getMomentumN().Mag() - NMomentum_bin_min)/NMomentum_bin_width;
				hXb_mom_bins[binPn]->Fill(eHit->getXb());
				hXp_mom_bins[binPn]->Fill(this_tag->getXp());
			}

			double thisaS = this_tag->getAs();
			if( thisaS > alphaS_bin_min && thisaS < alphaS_bin_max) {
				int binaS = (this_tag->getAs() - alphaS_bin_min)/alphaS_bin_width;
				hXb_aS_bins[binaS]->Fill(eHit->getXb());
				hXp_aS_bins[binaS]->Fill(this_tag->getXp());
			}

			// Fill As distribution for lo and hi Xp
			if( this_tag->getXp() > 0.25 && this_tag->getXp() < 0.35 )
				hAs_lo -> Fill( this_tag->getAs() );
			else if( this_tag->getXp() > 0.5 )
				hAs_hi -> Fill( this_tag->getAs() );

			// (E_i,p_i) = (mD,0) - (E_n,p_n)
			// virt = (E_i^2-p_i^2 - mP^2)/mP^2
			double En = sqrt(thisPn*thisPn + mN*mN);
			double E_i = mD - En;
			double p_i = thisPn;
			double thisVirt = (E_i*E_i - p_i*p_i - mP*mP)/(mP*mP);
			double thisptmag = this_tag->getPt().Mag();
			h2AsVi->Fill( thisaS , thisVirt );
			h2AsPn->Fill( thisaS , thisPn );
			h2ViPn->Fill( thisVirt, thisPn );
			h2Aspt->Fill( thisaS , thisptmag);

			h3_datastore_alphaS_Wp_Q->Fill(eHit->getQ2(), this_tag->getWp(), thisaS);


		} // end loop over events

		inFile->Close();
	}// end loop over files

	// Get the normalization for the background:
	TFitResultPtr fit = (TFitResultPtr)hToF_bac->Fit("pol0","QESR","",-20,0);
	double norm_per_bin = fit->Parameter(0);
		// Given our momentum max and min, solve for bins in ToF/m
	double beta_min = 1./sqrt(1.+ pow(mN/NMomentum_min,2));
	double beta_max = 1./sqrt(1.+ pow(mN/NMomentum_max,2));
	// max beta = min ToF and vice versa
	double TofpM_max = 1./(cAir*beta_min)*100;
	double TofpM_min = 1./(cAir*beta_max)*100;
	int TofpM_min_bin = hToF_bac->FindBin( TofpM_min );
	int TofpM_max_bin = hToF_bac->FindBin( TofpM_max );
	int nBins = (TofpM_max_bin - TofpM_min_bin);
	double background_counts = norm_per_bin * nBins;
	TVector3 bacnorm(background_counts,0,0);

	outFile->cd();
	bacnorm.Write("bacnorm");
	hToF_bac->Write();
	hXp->Write();
	hXb->Write();
	hAs->Write();
	hQ2->Write();
	h2XpQ2->Write();
	h2Q2Wp_hi->Write();
	h2Q2Wp_lo->Write();
	h2XpAs->Write();
	hWp->Write();
	hPn->Write();
	for( int i = 0 ; i < NAl_bins ; i++){
		hXpBins[i]->Write();
	}
	for( int i = 0 ; i < NXp_bins ; i++){
		hAs_bins[i]->Write();
	}
	for( int i = 0 ; i < NMomentum_bins ; i++){
		hXb_mom_bins[i]->Write();
		hXp_mom_bins[i]->Write();
	}
	for( int i = 0 ; i < NalphaS_bins ; i++){
		hXb_aS_bins[i]->Write();
		hXp_aS_bins[i]->Write();
	}
	for( int i = 0 ; i < NXp_bins ; i++){
		hPn_xp[i]->Write();
		hPn_xb[i]->Write();
	}
	hAs_lo->Write();
	hAs_hi->Write();

	h2AsVi->Write();
	h2AsPn->Write();
	h2ViPn->Write();
	h2Aspt->Write();

	h3_datastore_alphaS_Wp_Q->Write();


	outFile->Close();
	return 0;
}
