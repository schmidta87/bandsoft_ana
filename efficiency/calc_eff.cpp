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
#include "genpart.h"

// For processing data

using namespace std;
bool pointsToBand(double theta,double phi,double z_m);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [inputFiles] \n\n";
		cerr << "\t\t[inputFile] = ____.root ____.root ____.root ...\n\n";
		return -1;
	}

	// Conditions for a final accepted event neutron - signal or background
	TCut nGood	= Form("goodneutron == %i",					NCUT_goodneutron);
	TCut nLeadIdx	= Form("nleadindex != %i",					NCUT_leadindex);
	TCut nStatus	= Form("nHits[nleadindex]->getStatus() == %i",			NCUT_status);
	TCut nEdep	= Form("nHits[nleadindex]->getPmtLadc() > %f",			NCUT_Edep);
	TCut nToF	= Form("nHits[nleadindex]->getTofFadc() > %f",			NCUT_Tofabove0);
	//TCut nPn	= Form("tag[nleadindex]->getMomentumN().Mag() > %f && tag[nleadindex]->getMomentumN().Mag() < %f", NCUT_Pn_min, NCUT_Pn_max );
	TCut nPn	= Form("mcParts[1]->getMomentum() > %f && mcParts[1]->getMomentum() < %f",	NCUT_Pn_min,	NCUT_Pn_max);
	TCut nPnNaN	= Form("tag[nleadindex]->getMomentumN().Mag() == tag[nleadindex]->getMomentumN().Mag()");
	TCut good_neutrons = nGood && nLeadIdx && nStatus && nEdep && nToF && nPn && nPnNaN;

	// Final TCut:
	TString cut = Form("%s",good_neutrons.GetTitle());

	// Load input files
	TChain* inTree 		= new TChain("tagged");
	for( int i = 1 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree->Add(argv[i]);
	}
	// Setup the input branches I'll use
	TClonesArray * mcParts 	= new TClonesArray("genpart");
	inTree->SetBranchAddress("mcParts"	,&mcParts	);

	// Create the good event list from our cuts defined above
	inTree->Draw(">>goodEvents",cut,"",100000);
	TEventList * goodEvents = (TEventList*) gDirectory->Get("goodEvents");
	int nEvents = goodEvents->GetN();


	// Loop over all the good events and write the output tree
	int PTB = 0;
	for( int ev = 0 ; ev < 100000 ; ev++ ){
		mcParts->Clear();
		if( ev % 100000 == 0 ) cout << "\ton event " << ev << "\n";
		
		inTree->GetEntry(ev);

		genpart * gen_neutron = (genpart*)mcParts->At(1);
		double thetaN 	= gen_neutron->getTheta();
		double phiN	= gen_neutron->getPhi();
		bool ptb = pointsToBand(thetaN,phiN,0.);
		if( ptb && gen_neutron->getMomentum() > NCUT_Pn_min && gen_neutron->getMomentum() < NCUT_Pn_max ) PTB++;
	}
	cout << "generated events that point-to-band: " << PTB <<"\n";

	cout << "events that are reconstructed and pass neutron PID: " << nEvents << "\n";

	cout << "efficiency: " << (double)nEvents/PTB << "\n";



	return 0;
}



bool pointsToBand(double theta,double phi,double z_m){
	//double z = z_m*100; // from m to cm
	double z = z_m;

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = fabs(zDown/cos(theta));
	double xDown = rho*sin(theta)*cos(phi);
	double yDown = rho*sin(theta)*sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
	  )
		return 1;

	return 0;
}
