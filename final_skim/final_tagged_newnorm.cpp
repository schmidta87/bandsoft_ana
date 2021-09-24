#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TEventList.h"
#include "TCut.h"
#include "TChain.h"
#include "TVector3.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "constants.h"
#include "kinematic_cuts.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

// For processing data

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	if (argc < 3){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFiles] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[<MC,DATA,MIXED> = <0, 1, 2> \n";
		cerr << "\t\t[inputFile] = ____.root ____.root ____.root ...\n\n";
		return -1;
	}
	int MC_DATA_OPT = atoi(argv[2]);
	TCut ePID 		= "";
        TCut eCharge 		= "";
        TCut eEoP		= "";
        TCut eEpcal		= "";
        TCut eVW		= "";
        TCut eVtx		= "";
        TCut eMom		= "";
        TCut eQ2		= "";
        TCut eW			= "";

	// Conditions for a final accepted event electron
	if( MC_DATA_OPT == 1 || MC_DATA_OPT == 2 ){
		ePID 	= Form("eHit->getPID() == %i",					ECUT_PID);
		eCharge = Form("eHit->getCharge() == %i",				ECUT_charge); 
		eEoP	= Form("eHit->getEoP() > %f && eHit->getEoP() < %f",		ECUT_EoP_min,		ECUT_EoP_max);
		eEpcal	= Form("eHit->getEpcal() > %f",					ECUT_Epcal_min);
		eVW	= Form("eHit->getV() > %f && eHit->getW() > %f",		ECUT_V_min,		ECUT_W_min);
		eVtx	= Form("eHit->getVtz() > %f && eHit->getVtz() < %f",		ECUT_vtx_min,		ECUT_vtx_max);
		eMom	= Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f",	ECUT_pE_min,		ECUT_pE_max);
		eQ2	= Form("eHit->getQ2() > %f && eHit->getQ2() < %f",		ECUT_Q2_min,		ECUT_Q2_max);
		eW	= Form("eHit->getW2() > %f",					ECUT_W2_min);
	}
	if( MC_DATA_OPT == 0 ){
		ePID 	= Form("eHit_smeared->getPID() == %i",						ECUT_PID);
		eCharge = Form("eHit_smeared->getCharge() == %i",					ECUT_charge); 
		eEoP	= Form("eHit_smeared->getEoP() > %f && eHit_smeared->getEoP() < %f",		ECUT_EoP_min,		ECUT_EoP_max);
		eEpcal	= Form("eHit_smeared->getEpcal() > %f",						ECUT_Epcal_min);
		eVW	= Form("eHit_smeared->getV() > %f && eHit_smeared->getW() > %f",		ECUT_V_min,		ECUT_W_min);
		eVtx	= Form("eHit_smeared->getVtz() > %f && eHit_smeared->getVtz() < %f",		ECUT_vtx_min,		ECUT_vtx_max);
		eMom	= Form("eHit_smeared->getMomentum() > %f && eHit_smeared->getMomentum() < %f",	ECUT_pE_min,		ECUT_pE_max);
		eQ2	= Form("eHit_smeared->getQ2() > %f && eHit_smeared->getQ2() < %f",		ECUT_Q2_min,		ECUT_Q2_max);
		eW	= Form("eHit_smeared->getW2() > %f",						ECUT_W2_min);
	}
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eMom && eQ2 && eW;

	// Conditions for a final accepted event neutron - signal or background
	TCut nGood	= Form("goodneutron == %i",					NCUT_goodneutron);
	TCut nLeadIdx	= Form("nleadindex != %i",					NCUT_leadindex);
	TCut nEdep	= Form("nHits[nleadindex]->getEdep() > %f",			NCUT_Edep);
	TCut nThetaNQ 		= "";
		// make a fiducial cut around the edge of BAND:
	TCut fiducial1	= Form("!(nHits[nleadindex]->getSector()==1 && (nHits[nleadindex]->getX() >  70 || nHits[nleadindex]->getX() < -70))");
	TCut fiducial2	= Form("!(nHits[nleadindex]->getSector()==2 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() < -90))");
	TCut fiducial3	= Form("!(nHits[nleadindex]->getSector()==3 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() <  60))");
	TCut fiducial4	= Form("!(nHits[nleadindex]->getSector()==4 && (nHits[nleadindex]->getX() < -90 || nHits[nleadindex]->getX() > -60))");
	TCut fiducial5	= Form("!(nHits[nleadindex]->getSector()==5 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() < -90))");
		// to this fiducial, implement a ThetaN cut for beam pipe issue:
	TCut fiducialTheta 	= "";
	if( MC_DATA_OPT == 1 || MC_DATA_OPT == 2 ){
		fiducialTheta	= Form("tag[nleadindex]->getMomentumN().Theta() < 168.5*TMath::Pi()/180.");
		nThetaNQ	= Form("tag[nleadindex]->getThetaNQ() > %f && tag[nleadindex]->getThetaNQ() < %f",NCUT_THETANQ_min,NCUT_THETANQ_max);
	}
	if( MC_DATA_OPT == 0 ){
		fiducialTheta	= Form("tag_smeared[nleadindex]->getMomentumN().Theta() < 168.5*TMath::Pi()/180.");
		nThetaNQ	= Form("tag_smeared[nleadindex]->getThetaNQ() > %f && tag_smeared[nleadindex]->getThetaNQ() < %f",NCUT_THETANQ_min,NCUT_THETANQ_max);
	}
	TCut fiducial = fiducial1 && fiducial2 && fiducial3 && fiducial4 && fiducial5 && fiducialTheta;

	// Bad bars for v3.1, consistent across run periods:
	//TCut nBad_111	= Form("!(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_121	= Form("!(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_255	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==5)");
	//TCut nBad_352	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==2)");

	TCut tagged 	= inclusive && nGood && nLeadIdx && nEdep && nThetaNQ  && fiducial;	// basic cuts 
	//tagged		= tagged && nBad_111 && nBad_121 && nBad_255 && nBad_352;		// bad bar removal


	// Conditions for a final accepted event neutron in signal region
	TCut nToF	= Form("nHits[nleadindex]->getTof() > %f && nHits[nleadindex]->getTof() < %f"	,NCUT_Tof_min, NCUT_Tof_max);
	TCut nPn	= "";
	TCut nPnNaN	= "";
	TCut nWp	= "";
	TCut nAs	= "";
	if( MC_DATA_OPT == 1 || MC_DATA_OPT == 2 ){
		nPn	= Form("tag[nleadindex]->getMomentumN().Mag() > %f && tag[nleadindex]->getMomentumN().Mag() < %f", NCUT_Pn_min, NCUT_Pn_max );
		nPnNaN	= Form("tag[nleadindex]->getMomentumN().Mag() == tag[nleadindex]->getMomentumN().Mag()");
		nWp	= Form("tag[nleadindex]->getWp() > %f && tag[nleadindex]->getWp() < %f",	NCUT_Wp_min,	NCUT_Wp_max);
		nAs	= Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",	NCUT_As_min,	NCUT_As_max);
	}
	if( MC_DATA_OPT == 0 ){
		nPn	= Form("tag_smeared[nleadindex]->getMomentumN().Mag() > %f && tag_smeared[nleadindex]->getMomentumN().Mag() < %f", NCUT_Pn_min, NCUT_Pn_max );
		nPnNaN	= Form("tag_smeared[nleadindex]->getMomentumN().Mag() == tag_smeared[nleadindex]->getMomentumN().Mag()");
		nWp	= Form("tag_smeared[nleadindex]->getWp() > %f && tag_smeared[nleadindex]->getWp() < %f",	NCUT_Wp_min,	NCUT_Wp_max);
		nAs	= Form("tag_smeared[nleadindex]->getAs() > %f && tag_smeared[nleadindex]->getAs() < %f",	NCUT_As_min,	NCUT_As_max);
	}
	TCut tagged_signal = tagged && nToF && nPn && nPnNaN && nWp && nAs;
	
	// Final TCut:
	TString cut = Form("%s",tagged_signal.GetTitle());

	// Load input files
	TChain* inTree = new TChain("tagged");
	for( int i = 3 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree->Add(argv[i]);
	}

	// Create an output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");
	inTree->LoadTree(0);
	TTree * outTree = inTree->CloneTree(0);


	// get background normalization level for this neutron PID
	TH1D * hToF_bac = NULL;
	TH1D * hToF_far = NULL;
	TH1D * hToF_full = new TH1D("hToF_full","hToF_full",8000,-100,300);
	TH1D * hThetaNQ = new TH1D("hThetaNQ","hThetaNQ",8000,140,180);
	TH1D * hThetaN 	= new TH1D("hThetaN","hThetaN",8000,140,180);
	TH1D * hPe = new TH1D("hPe","hPe",4000,3,7);
	TVector3 bacnorm;
	if( MC_DATA_OPT == 1 ){ // if this is data file
		
		// Use the data file to get the number of counts in the background region
		// and rescale it by the region sizes of Signal and Background.
		hToF_bac = new TH1D("hToF_bac","hToF_bac",8000,-100,300);
		hToF_far = new TH1D("hToF_far","hToF_far",8000,-100,300);
		TCut background_time = Form("nHits[nleadindex]->getTof() > %f && nHits[nleadindex]->getTof() < %f",
						NCUT_BACK_Tof_min , NCUT_BACK_Tof_max );
		TCut far_time = Form("nHits[nleadindex]->getTof() > %f && nHits[nleadindex]->getTof() < %f",
						NCUT_FARTIME_Tof_min , NCUT_FARTIME_Tof_max );
		//inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_bac",tagged && (background_time || far_time) );
		inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_bac",tagged && background_time );
		inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_far",tagged && far_time );
		inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_full",tagged );

		inTree->Draw("tag[nleadindex]->getThetaNQ() * 180./TMath::Pi() >> hThetaNQ",tagged && background_time );
		inTree->Draw("tag[nleadindex]->getMomentumN().Theta() * 180./TMath::Pi() >> hThetaN",tagged && background_time );
		inTree->Draw("eHit->getMomentum() >> hPe", tagged && background_time);

		// Now we need to rescale this by the amount of bunches in the signal region:
		double nSignalBunches 	= (NCUT_Tof_max - NCUT_Tof_min)/BEAM_BUNCH;
		double bac_nBkgrdBunches 	= (NCUT_BACK_Tof_max - NCUT_BACK_Tof_min)/BEAM_BUNCH;
		double far_nBkgrdBunches 	= (NCUT_FARTIME_Tof_max - NCUT_FARTIME_Tof_min)/BEAM_BUNCH;

		bacnorm.SetXYZ( hToF_bac->Integral() * nSignalBunches/bac_nBkgrdBunches , 
			     	hToF_far->Integral() * nSignalBunches/far_nBkgrdBunches ,
				0 );
		bacnorm.SetZ( (bacnorm.X() + bacnorm.Y())/2. );

		//TFitResultPtr fit = (TFitResultPtr)hToF_bac->Fit("pol0","QESR","",-20,0);
		//double norm_per_bin = fit->Parameter(0);
		//	// Given our momentum max and min, solve for bins in ToF/m
		//double beta_min = 1./sqrt(1.+ pow(mN/NCUT_Pn_min,2));
		//double beta_max = 1./sqrt(1.+ pow(mN/NCUT_Pn_max,2));
		//	// max beta = min ToF and vice versa
		//double TofpM_max = 1./(cAir*beta_min)*100;
		//double TofpM_min = 1./(cAir*beta_max)*100;
		//int TofpM_min_bin = hToF_bac->FindBin( TofpM_min );
		//int TofpM_max_bin = hToF_bac->FindBin( TofpM_max );
		//int nBins = (TofpM_max_bin - TofpM_min_bin); 	
		//double background_counts = norm_per_bin * nBins;
	
		//bacnorm.SetXYZ(background_counts,0,0);
	}
	if( MC_DATA_OPT == 2 ){ // if this is mixed file for background
		// Grab the number of mixed events we have created
		double N_mixed		= inTree->GetEntries(); // N_mixed is different from the correct normalization since
								// now we have fiducial and bad bar cuts.
		inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_full",tagged );
		inTree->Draw("tag[nleadindex]->getThetaNQ() * 180./TMath::Pi() >> hThetaNQ",tagged );
		inTree->Draw("tag[nleadindex]->getMomentumN().Theta() * 180./TMath::Pi() >> hThetaN",tagged );
		inTree->Draw("eHit->getMomentum() >> hPe", tagged );
		N_mixed = hToF_full->Integral();
		
		bacnorm.SetXYZ( N_mixed , 0 , 0 );

	}
	if( MC_DATA_OPT == 0 ){ // if this is a simulation file
		inTree->Draw("nHits[nleadindex]->getTof()  >> hToF_full",tagged );
		inTree->Draw("tag_smeared[nleadindex]->getThetaNQ() * 180./TMath::Pi() >> hThetaNQ",tagged );
		inTree->Draw("tag_smeared[nleadindex]->getMomentumN().Theta() * 180./TMath::Pi() >> hThetaN",tagged );
		inTree->Draw("eHit_smeared->getMomentum() >> hPe", tagged );

	}

	// Create the good event list from our cuts defined above
	inTree->Draw(">>goodEvents",cut);
	TEventList * goodEvents = (TEventList*) gDirectory->Get("goodEvents");
	int nEvents = goodEvents->GetN();


	// Loop over all the good events and write the output tree
	for( int ev = 0 ; ev < nEvents ; ev++ ){
		if( ev % 100000 == 0 ) cout << "\ton event " << ev << "\n";

		int entry = goodEvents->GetEntry(ev);
		inTree->GetEntry(entry);

		outTree->Fill();

	} // end loop over events

	cout << "writing stuff to file:\n";
	// Write and close
	outFile->cd();
	hToF_full->Write();
	hThetaNQ->Write();
	hThetaN->Write();
	hPe->Write();
	if( MC_DATA_OPT == 1 ){ // if this is data file
		hToF_bac->Write();
		hToF_far->Write();
	}
	if( MC_DATA_OPT == 2 ){
	}
	if( MC_DATA_OPT == 0 ){
	}
	outTree->Write();
	bacnorm.Write("bacnorm");	
	outFile->Close();
	return 0;
}
