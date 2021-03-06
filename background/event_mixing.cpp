#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <vector>

#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVectorT.h"
#include "TRandom3.h"
#include "constants.h"
#include "kinematic_cuts.h"
#include "TClonesArray.h"


using namespace std;

//clashit (const clashit &input) {x = p2.x; y = p2.y; }


int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./event_mixing [outputRootfile] [inputElectronSkim] [inputNeutronSkim]\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("tagged","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	bool goodneutron = false;
	int nleadindex = -1;
	double weight		= 0;
	//	MC info:
	int genMult		= 0;
	TClonesArray * mcParts = new TClonesArray("genpart");
	TClonesArray &saveMC = *mcParts;
	// 	Neutron info:
	int nMult		= 0;
	TClonesArray * nHits = new TClonesArray("bandhit");
	TClonesArray &saveHit = *nHits;
	//	Electron info:
	clashit eHit;
	//	Tagged info:
	TClonesArray * tags = new TClonesArray("taghit");
	TClonesArray &saveTags = *tags;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("weight"	,&weight		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//	Tagged branches:
	outTree->Branch("tag"		,&tags			);
	//	MC branches:
	outTree->Branch("genMult"	,&genMult		);
	outTree->Branch("mcParts"	,&mcParts		);

	// Setup rand for picking particles
	TRandom3 * myRand = new TRandom3(0);


	// Define 4D binned ToF plots and bin edges and kin cuts from kinematic_cuts.h file
	const double max_CosTheta_nq 	= cos(NCUT_THETANQ_min);	// bit non-intuitive but min angle = 150 which is a "max" costhetanq of -0.8
	const double min_CosTheta_nq 	= cos(NCUT_THETANQ_max);	// and then max angle = 180 which is a "min" costhetanq of -1

	// Define background and signal edges from from kinematic_cuts.h file
	const double bkgrd_min 		= NCUT_BACK_TofpM_min;
	const double bkgrd_max 		= NCUT_BACK_TofpM_max;
	const double signal_min 	= NCUT_TofpM_min;
	const double signal_max 	= NCUT_TofpM_max;

	// Loop over the neutron and electron skim files given
	TFile * inFile_e = new TFile(argv[2]);
	TFile * inFile_n = new TFile(argv[3]);
	TTree * inTree_e = (TTree*)inFile_e->Get("electrons");
	TTree * inTree_n = (TTree*)inFile_n->Get("neutrons");

	clashit * input_eHit 		= new clashit;
	TClonesArray* input_nHit 	= new TClonesArray("bandhit");
	int input_nleadindex		= -1;
	double input_Ebeam		= 0;

	inTree_e->SetBranchAddress("eHit",		&input_eHit);
	inTree_e->SetBranchAddress("Ebeam",		&input_Ebeam);
	inTree_n->SetBranchAddress("nHits",		&input_nHit);
	inTree_n->SetBranchAddress("nleadindex",	&input_nleadindex);



	// Read entire tree into memory due to issue with jumping around in tree as we read it..
	vector<clashit> electron_list;
	vector<bandhit> neutron_list;
	vector<double>  Ebeam_list;
	cout << "Saving neutrons to mem...\n";
	for( int neutron = 0 ; neutron < inTree_n->GetEntries() ; neutron++ ){
		input_nHit->Clear();
		input_nleadindex = -1;

		inTree_n->GetEntry(neutron);

		// This is AFTER final neutron, so all lead neutrons are good events
		// so just get the lead neutron and save it
		bandhit * lead_n = (bandhit*) input_nHit->At(input_nleadindex);
		bandhit copy_n;
		copy_n = *lead_n;
		neutron_list.push_back( copy_n );
	}

	cout << "Saving electrons to mem...\n";
	for( int electron = 0 ; electron < inTree_e->GetEntries() ; electron++ ){
		input_eHit->Clear();
		input_Ebeam = 0;

		inTree_e->GetEntry(electron);

		// This is AFTER final inclusive, so all electrons are good events
		// so just get it and save it
		clashit copy_e;
		copy_e = *input_eHit;
		electron_list.push_back( copy_e );
		Ebeam_list.push_back( input_Ebeam );
	}

	
	// Loop over all neutron events
	for( int neutron = 0 ; neutron < neutron_list.size() ; neutron++ ){
		if( neutron % 1000 == 0 ) cout << "working on neutron " << neutron << "\n";

		int nPairs = 0; // counts how many pairs were made for this neutron
		int nFails = 0; // counts how many times a pair was tried for
		int maxPairs = 10;
		int maxFails = 10;
		while( nPairs < maxPairs && nFails < maxFails ){

			int electron = myRand->Rndm() * inTree_e->GetEntries();

			// Read in the stored variables
			double this_Ebeam 	= Ebeam_list	[electron];
			double p_e		= electron_list	[electron].getMomentum();
			double theta_e		= electron_list	[electron].getTheta();
			double phi_e		= electron_list	[electron].getPhi();
			double vtz 		= electron_list	[electron].getVtz();
			double theta_n		= neutron_list	[neutron].getDL().Theta();
			double phi_n		= neutron_list	[neutron].getDL().Phi();
			double dL		= neutron_list	[neutron].getDL().Mag();

			// Re-drawn in ToF per meter within our signal region for neutrons:
			double TofpM = myRand->Rndm() * (NCUT_TofpM_max - NCUT_TofpM_min) + NCUT_TofpM_min;
			double beta = (1./TofpM) * (1./cAir) * (100./1);
			double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
			double ToF = TofpM * (dL/100.);

			// Create vectors for calculating angles
			TVector3 beamVec(0,0,this_Ebeam);
			TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
			TVector3 qVec;	qVec = beamVec - eVec;
			TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);

			double q 	= electron_list[electron].getQ();
			double theta_q  = electron_list[electron].getThetaQ();
			double phi_q 	= electron_list[electron].getPhiQ();
			double nu 	= electron_list[electron].getOmega();
			double Q2 	= electron_list[electron].getQ2();
			double xB	= electron_list[electron].getXb();
			double W2	= electron_list[electron].getW2();
			double E_n 	= sqrt( mN*mN + p_n*p_n );

			// Calculate the nq angles
			TVector3 norm_scatter = qVec.Cross( beamVec );
			norm_scatter 	= norm_scatter.Unit();
			TVector3 norm_reaction = qVec.Cross( nVec );
			norm_reaction 	= norm_reaction.Unit();
			double phi_nq 	= norm_scatter.Angle( norm_reaction );
			double theta_nq = nVec.Angle( qVec );
			double CosTheta_nq = cos(theta_nq);
			TVector3 direction = norm_scatter.Cross(norm_reaction);
			if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
			}
			else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
				phi_nq *= (-1);
			}
			double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
			double Wp = sqrt(W_primeSq);
			double Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
			double As = (E_n - p_n*CosTheta_nq)/mN;
			double Xp2 = Q2/(W_primeSq - mN*mN + Q2);

			TVector3 Pt;
			TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
			Pt = nVec - pN_par_q;

			// Only cut not implemented in the 
			// final skim is the cosThetaNQ because neutron had no electron
			if( CosTheta_nq > min_CosTheta_nq 	&& 
			    CosTheta_nq < max_CosTheta_nq 	){
				// Valid pair
				nFails = 0;
				nPairs++;
				
				// Clear ehit branches
				eHit.Clear();
				Ebeam 		= 0;
				gated_charge 	= 0;
				livetime 	= 0;
				starttime 	= 0;
				current 	= 0;
				weight 		= 1;
				// Clear nhit branches
				nHits->Clear();
				nMult 		= 1;
				goodneutron 	= NCUT_goodneutron;
				nleadindex 	= 0;
				// Clear tag branches
				tags->Clear();

				// Save electron info
				Ebeam = this_Ebeam;
				eHit = electron_list[electron];

				// Save neutron info with new ToF
				neutron_list[neutron].setTof(ToF);
				neutron_list[neutron].setTofFadc(ToF);
				new(saveHit[0]) bandhit;
				saveHit[0] = &neutron_list[neutron];

				// Save tag info with new vectors
				taghit new_tag;
				new_tag.setMomentumE	(eVec 		);
				new_tag.setMomentumN	(nVec		);
				new_tag.setMomentumQ	(qVec		);
				new_tag.setMomentumB	(beamVec	);
				new_tag.setPhiNQ	(phi_nq		);
				new_tag.setThetaNQ	(theta_nq	);
				new_tag.setWp		(Wp		);
				new_tag.setXp		(Xp		);
				new_tag.setAs		(As		);
				new_tag.setPt		(Pt		);
				new_tag.setXp2		(Xp2		);
				new(saveTags[0]) taghit;
				saveTags[0] = &new_tag;


				outTree->Fill();


			}
			else{ nFails++; }


		} // end loop over finding electrons

	}
	
	// Write output file:
	outFile->cd();
	outTree->Write();
	outFile->Close();


	return 0;
}
