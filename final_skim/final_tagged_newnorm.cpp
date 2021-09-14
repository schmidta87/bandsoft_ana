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

	// Conditions for a final accepted event electron
	TCut ePID 	= Form("eHit->getPID() == %i",					ECUT_PID);
	TCut eCharge 	= Form("eHit->getCharge() == %i",				ECUT_charge); 
	TCut eEoP	= Form("eHit->getEoP() > %f && eHit->getEoP() < %f",		ECUT_EoP_min,		ECUT_EoP_max);
	TCut eEpcal	= Form("eHit->getEpcal() > %f",					ECUT_Epcal_min);
	TCut eVW	= Form("eHit->getV() > %f && eHit->getW() > %f",		ECUT_V_min,		ECUT_W_min);
	TCut eVtx	= Form("eHit->getVtz() > %f && eHit->getVtz() < %f",		ECUT_vtx_min,		ECUT_vtx_max);
	TCut eMom	= Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f",	ECUT_pE_min,		ECUT_pE_max);
	TCut eQ2	= Form("eHit->getQ2() > %f && eHit->getQ2() < %f",		ECUT_Q2_min,		ECUT_Q2_max);
	TCut eW		= Form("eHit->getW2() > %f",					ECUT_W2_min);
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eMom && eQ2 && eW;

	// Conditions for a final accepted event neutron - signal or background
	TCut nGood	= Form("goodneutron == %i",					NCUT_goodneutron);
	TCut nLeadIdx	= Form("nleadindex != %i",					NCUT_leadindex);
	TCut nEdep;
		// data or background:
	if( MC_DATA_OPT == 0 || MC_DATA_OPT == 2)
		nEdep	= Form("nHits[nleadindex]->getEdep() > %f",			NCUT_Edep);
		// sim:
	else if( MC_DATA_OPT == 1)
		nEdep	= Form("nHits[nleadindex]->getEdep() > %f",			NCUT_Edep);
	TCut nThetaNQ	= Form("tag[nleadindex]->getThetaNQ() > %f && tag[nleadindex]->getThetaNQ() < %f",NCUT_THETANQ_min,NCUT_THETANQ_max);
		// make a fiducial cut around the edge of BAND:
	TCut fiducial1	= Form("!(nHits[nleadindex]->getSector()==1 && (nHits[nleadindex]->getX() >  70 || nHits[nleadindex]->getX() < -70))");
	TCut fiducial2	= Form("!(nHits[nleadindex]->getSector()==2 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() < -90))");
	TCut fiducial3	= Form("!(nHits[nleadindex]->getSector()==3 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() <  60))");
	TCut fiducial4	= Form("!(nHits[nleadindex]->getSector()==4 && (nHits[nleadindex]->getX() < -90 || nHits[nleadindex]->getX() > -60))");
	TCut fiducial5	= Form("!(nHits[nleadindex]->getSector()==5 && (nHits[nleadindex]->getX() >  90 || nHits[nleadindex]->getX() < -90))");
		// to this fiducial, implement a ThetaN cut for beam pipe issue:
	TCut fiducialTheta	= Form("tag[nleadindex]->getMomentumN().Theta() < 168.5*TMath::Pi()/180.");
	TCut fiducial = fiducial1 && fiducial2 && fiducial3 && fiducial4 && fiducial5 && fiducialTheta;
	// kill any bad bars:
		// Spring 2019 bad bars that will always be killed:
	//TCut nBad_341	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_432	= Form("!(nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==2)");
	//TCut nBad_425	= Form("!(nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==5)");
	//	// This is a hole in BAND due to a TDC shift for 10.2 GeV data set:
	//TCut nHole_1 	= Form("!(nHits[nleadindex]->getSector()==2 && (nHits[nleadindex]->getComponent()==4 || nHits[nleadindex]->getComponent()==5 || nHits[nleadindex]->getComponent()==6 || nHits[nleadindex]->getComponent()==7) && (nHits[nleadindex]->getX()>90 || nHits[nleadindex]->getX()<-110))");
	//TCut nHole_2	= Form("!(nHits[nleadindex]->getSector()==3 && (nHits[nleadindex]->getComponent()==1 || nHits[nleadindex]->getComponent()==2 ) && (nHits[nleadindex]->getX()>80 || nHits[nleadindex]->getX() < 45) )");
	//	// This are killed due to discrepancy between MC and simulation:
	//	// WITH theta cut v2.0 04-21-2021
	//TCut nBad_111	= Form("!(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_121	= Form("!(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_225	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==5)");
	//TCut nBad_227	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==7)");
	//TCut nBad_424	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==4)");
	//TCut nBad_426	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_521	= Form("!(nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_236	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_431	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_241	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_242	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==2)");
	//TCut nBad_244	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==4)");
	//TCut nBad_346	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_353	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==3)");
		// WITHOUT theta cut v2.0 04-21-2021
	//TCut nBad_111	= Form("!(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_217	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getComponent()==7)");
	//TCut nBad_225	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==5)");
	//TCut nBad_226	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_424	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==4)");
	//TCut nBad_426	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_521	= Form("!(nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_236	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_336	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_431	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_532	= Form("!(nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getComponent()==2)");
	//TCut nBad_241	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_242	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==2)");
	//TCut nBad_244	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==4)");
	//TCut nBad_346	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_445	= Form("!(nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==5)");
	//TCut nBad_541	= Form("!(nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==1)");
	//TCut nBad_542	= Form("!(nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getComponent()==2)");
	//TCut nBad_255	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==5)");
	//TCut nBad_256	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==6)");
	//TCut nBad_257	= Form("!(nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==7)");
	//TCut nBad_351	= Form("!(nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getComponent()==1)");

	
	
	TCut tagged 	= inclusive && nGood && nLeadIdx && nEdep && nThetaNQ  && fiducial;	// basic cuts 
				//&& nBad_341 	&& nBad_432 && nBad_425				// bad bars from calibration for Spring2019
				//&& nHole_1 	&& nHole_2					// hole for Spring2019
				//&& nBad_111	&& nBad_121	&& nBad_225	&& nBad_227	// bad bars in MC vs data v2.0 04-21-2021 with theta cut 
				//&& nBad_424	&& nBad_426	&& nBad_521	&& nBad_236	
				//&& nBad_431	&& nBad_241	&& nBad_242	&& nBad_244	
				//&& nBad_346	&& nBad_353;					
				//&& nBad_111	&& nBad_217	&& nBad_225	&& nBad_226	// bad bars in MC vs data v2.0 04-21-2021 WITHOUT theta cut 
				//&& nBad_424	&& nBad_426	&& nBad_521	&& nBad_236	
				//&& nBad_336	&& nBad_431	&& nBad_532	&& nBad_241
				//&& nBad_242	&& nBad_244	&& nBad_346	&& nBad_445
				//&& nBad_541	&& nBad_542	&& nBad_255	&& nBad_256
				//&& nBad_257	&& nBad_351;


	// Conditions for a final accepted event neutron in signal region
	TCut nToF	= Form("nHits[nleadindex]->getTof() > %f && nHits[nleadindex]->getTof() < %f"	,NCUT_Tof_min, NCUT_Tof_max);
	TCut nPn	= Form("tag[nleadindex]->getMomentumN().Mag() > %f && tag[nleadindex]->getMomentumN().Mag() < %f", NCUT_Pn_min, NCUT_Pn_max );
	TCut nPnNaN	= Form("tag[nleadindex]->getMomentumN().Mag() == tag[nleadindex]->getMomentumN().Mag()");
	TCut nWp	= Form("tag[nleadindex]->getWp() > %f && tag[nleadindex]->getWp() < %f",	NCUT_Wp_min,	NCUT_Wp_max);
	TCut nAs	= Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",	NCUT_As_min,	NCUT_As_max);
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
		double N_mixed		= inTree->GetEntries();
		
		bacnorm.SetXYZ( N_mixed , 0 , 0 );

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
	if( MC_DATA_OPT == 1 ){ // if this is data file
		hToF_bac->Write();
		hToF_far->Write();
	}
	outTree->Write();
	bacnorm.Write("bacnorm");	
	outFile->Close();
	return 0;
}
