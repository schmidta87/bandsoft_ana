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
		cerr << "\t\t[<MC,DATA> = <0, 1> \n";
		cerr << "\t\t[inputFile] = ____.root ____.root ____.root ...\n\n";
		return -1;
	}
	int MC_DATA_OPT = atoi(argv[2]);

	// Conditions for a final accepted event neutron - signal or background
	TCut nGood	= Form("goodneutron == %i",					NCUT_goodneutron);
	TCut nLeadIdx	= Form("nleadindex != %i",					NCUT_leadindex);
	TCut nStatus	= Form("nHits[nleadindex]->getStatus() == %i",			NCUT_status);
	TCut nEdep;
	if( MC_DATA_OPT == 0 )
		nEdep	= Form("nHits[nleadindex]->getPmtLadc() > %f",			NCUT_Edep*SimAdcToMeVee);
	else if( MC_DATA_OPT == 1)
		nEdep	= Form("nHits[nleadindex]->getEdep() > %f",			NCUT_Edep*DataAdcToMeVee);
		// kill any bad bars:
	TCut nBad_122	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==2 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_132	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_133	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==3)");
	TCut nBad_135	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==5)");
	TCut nBad_136	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==6)");
	TCut nBad_142	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_146	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==6)");
	TCut nBad_152	= Form("!(nHits[nleadindex]->getLayer()==1 && nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_251	= Form("!(nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getComponent()==1)");
	TCut nBad_252	= Form("!(nHits[nleadindex]->getLayer()==2 && nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_333	= Form("!(nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==3)");
	TCut nBad_342	= Form("!(nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_346	= Form("!(nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==6)");
	TCut nBad_352	= Form("!(nHits[nleadindex]->getLayer()==3 && nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_411	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1)");
	TCut nBad_435	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==5)");
	TCut nBad_442	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_443	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==3)");
	TCut nBad_444	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==4)");
	TCut nBad_445	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==5)");
	TCut nBad_452	= Form("!(nHits[nleadindex]->getLayer()==4 && nHits[nleadindex]->getSector()==5 && nHits[nleadindex]->getComponent()==2)");
	TCut nBad_511	= Form("!(nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1)");
	TCut nBad_531	= Form("!(nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==1)");
	TCut nBad_535	= Form("!(nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getSector()==3 && nHits[nleadindex]->getComponent()==5)");
	TCut nBad_541	= Form("!(nHits[nleadindex]->getLayer()==5 && nHits[nleadindex]->getSector()==4 && nHits[nleadindex]->getComponent()==1)");
	TCut nHole_1 	= Form("!(nHits[nleadindex]->getSector()==2 && (nHits[nleadindex]->getComponent()==4 || nHits[nleadindex]->getComponent()==5 || nHits[nleadindex]->getComponent()==6 || nHits[nleadindex]->getComponent()==7) && (nHits[nleadindex]->getX()>90 || nHits[nleadindex]->getX()<-110))");
	TCut nHole_2	= Form("!(nHits[nleadindex]->getSector()==3 && (nHits[nleadindex]->getComponent()==1 || nHits[nleadindex]->getComponent()==2 ) && (nHits[nleadindex]->getX()>80 || nHits[nleadindex]->getX() < 45) )");
	TCut neutron 	= nGood && nLeadIdx && nStatus && nEdep && nHole_1 && nHole_2 
				&& nBad_122 && nBad_132 && nBad_133 && nBad_135 && nBad_136 && nBad_142 && nBad_146
				&& nBad_152 && nBad_251 && nBad_252 && nBad_333 && nBad_342 && nBad_346 && nBad_352 
				&& nBad_411 && nBad_435 && nBad_442 && nBad_443 && nBad_444 && nBad_445 && nBad_452 
				&& nBad_511 && nBad_531 && nBad_535 && nBad_541;

	// Final TCut:
	TString cut = Form("%s",neutron.GetTitle());

	// Load input files
	TChain* inTree = new TChain("neutrons");
	for( int i = 3 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree->Add(argv[i]);
	}

	// Create an output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");
	inTree->LoadTree(0);
	TTree * outTree = inTree->CloneTree(0);


	// get background normalization level for this neutron PID
	TH1D * hToF_bac = new TH1D("hToF_bac","hToF_bac",1000,-25,75);
	inTree->Draw("nHits[nleadindex]->getTofFadc() / (nHits[nleadindex]->getDL().Mag()/100.) >> hToF_bac",cut);
	TVector3 bacnorm;
	if( MC_DATA_OPT == 1 ){
		TFitResultPtr fit = (TFitResultPtr)hToF_bac->Fit("pol0","QESR","",-20,0);
		double norm_per_bin = fit->Parameter(0);
			// Given our momentum max and min, solve for bins in ToF/m
		double beta_min = 1./sqrt(1.+ pow(mN/NCUT_Pn_min,2));
		double beta_max = 1./sqrt(1.+ pow(mN/NCUT_Pn_max,2));
			// max beta = min ToF and vice versa
		double TofpM_max = 1./(cAir*beta_min)*100;
		double TofpM_min = 1./(cAir*beta_max)*100;
		int TofpM_min_bin = hToF_bac->FindBin( TofpM_min );
		int TofpM_max_bin = hToF_bac->FindBin( TofpM_max );
		int nBins = (TofpM_max_bin - TofpM_min_bin); 	
		double background_counts = norm_per_bin * nBins;
	
		bacnorm.SetXYZ(background_counts,0,0);
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


	// Write and close
	outFile->cd();
	outTree->Write();
	bacnorm.Write("bacnorm");	
	hToF_bac->Write();
	outFile->Close();
	return 0;
}
