#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TClonesArray.h"

#include "constants.h"
#include "genpart.h"
#include "bandhit.h"


using namespace std;

bool pointsToBand(double theta,double phi,double z_m);
int main( int argc, char ** argv ){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [SimFile]\n";
		return -1;
	}

	TFile * inSimFile = new TFile(argv[1]);
	TTree * inSimTree = (TTree*) inSimFile->Get("tagged");
	const int nEvents = inSimTree->GetEntries();

	TClonesArray * mcParts 	= new TClonesArray("genpart");
	TClonesArray * nHits 	= new TClonesArray("bandhit");
	int nMult = 0;
	inSimTree->SetBranchAddress("mcParts"	,&mcParts	);
	inSimTree->SetBranchAddress("nHits"	,&nHits		);
	inSimTree->SetBranchAddress("nMult"	,&nMult		);



	int PTB = 0;
	int INTR = 0;
	int processed = 0;
	int CNT0 = 0;
	int CNT1 = 0;
	int CNT2 = 0;
	double thres = 5;
	double time_thres = 0.3;
	int recoFlag = 1;
	for( int event = 0 ; event < nEvents ; event++ ){
		if( event > 100000 ) break;
		mcParts->Clear();
		nHits->Clear();
		nMult = 0;
		inSimTree->GetEntry(event);
		processed++;

		genpart * gen_neutron = (genpart*)mcParts->At(1);
		double thetaN 	= gen_neutron->getTheta();
		double phiN	= gen_neutron->getPhi();
		bool ptb = pointsToBand(thetaN,phiN,0.);

		if(ptb) PTB++;
		
		if( nMult == 1 ){
			bandhit * this_neutron = (bandhit*)nHits->At(0);
			int status 	= this_neutron->getStatus();
			double Edep 	= this_neutron->getEdep()/1E4;
			double ToF	= this_neutron->getTofFadc();
			if( ToF != 0 && status == 0 && Edep > 5) INTR++;
		}

		// Loop over all hits in BAND and get the earliest hit in time
		// but filter out hits that have low Edep
		bool accept = false;
		bool thresPass = false;
		bool vetoHit = false;
		double fastestTime = 1E5;
		int fastestTimeIdx = -1;
		for( int thishit = 0; thishit < nMult ; thishit++){
			bandhit * neutron = (bandhit *) nHits->At(thishit);
			if( neutron->getPmtLadc()/1E4 < 2 ) continue; // if the Edep in this PMT is less than 2MeVee, do not count it
			if( neutron->getTof() < fastestTime ){
				fastestTime = neutron->getTof();
				fastestTimeIdx = thishit;
			}
		}

		// Once we have the fastest hit, ask how many other hits there
		// are in the same event and flag how far in time they are,
		// but keep this earliest hit.
		// Accept the event based on if there is another hit > 300ps
		if( fastestTimeIdx != -1 ){
			accept = true;
			for( int thishit = 0; thishit < nMult ; thishit++){

				bandhit * neutron = (bandhit *) nHits->At(thishit);
				// if we have a veto bar hit and if the Edep is > 0.25MeVee, then count it as a true veto hit
				if( neutron->getLayer() == 6 && neutron->getPmtLadc()/1E4 > 0.25 ) vetoHit = true; 
				// if this bar does not have > 2MeVee, do not count the hit
				if( neutron->getPmtLadc()/1E4 < 2 ) continue;

				// if this bar does not have > X MeVee, do not count the hit (X = threshold set by user)
				if( neutron->getPmtLadc()/1E4 > thres ) thresPass = true;

				double tdiff = neutron->getTof() - fastestTime;
				if( tdiff == 0 ){
					continue;
				}
				bool adjLayer = false;
				bool adjComp = false;
				bool adjSpace = false;
				bandhit * lead = (bandhit *) nHits->At(fastestTimeIdx);

				int layerDiff = lead->getLayer() - neutron->getLayer();
				if( fabs(layerDiff) <= 1 ) adjLayer = true;
				if( lead->getSector() == neutron->getSector() ){
				// if hits are in the same sector, then just ask for comp diff
					int compDiff = lead->getComponent() - neutron->getComponent();
					if( fabs(compDiff) <= 1 ) adjComp = true;
				}
				else{
				// have to implement some fancy component check due to sec comp differences
					int yDiff = lead->getY() - neutron->getY();
					if( fabs(yDiff) <= 10 ) adjComp = true;
				}
				if( adjComp && adjLayer ) adjSpace = true;

				if( tdiff < time_thres && tdiff != 0 ) accept = false;
				if( recoFlag == 1 ) // Use improved recostruction:
					if( tdiff < time_thres && tdiff != 0 && adjSpace ) accept = true;
			}

			if( thresPass == true ) CNT0++;
			if( vetoHit == false && thresPass == true ) CNT1++;
			if( vetoHit == false && thresPass == true && accept == true) CNT2++;
		}
		

	}
	cout << (double)PTB/processed << "\n";
	cout << (double)INTR/PTB << "\n";
	cout << (double)CNT0/PTB << "\n";
	cout << (double)CNT1/PTB << "\n";
	cout << (double)CNT2/PTB << "\n";

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
