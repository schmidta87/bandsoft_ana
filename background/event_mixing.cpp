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
#include "TClonesArray.h"


using namespace std;

int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./event_mixing [outputRootfile] [inputElectronSkim] [inputNeutronSkim]\n";
		return -1;
	}

        // Create output tree
        TFile * outFile = new TFile(argv[1],"RECREATE");
        TTree * outTree = new TTree("tagged","BAND Neutrons and CLAS Electrons");
        //      Event info:
	int Runno               = 1;
	double Ebeam            = 0;
	double gated_charge     = 1;
	double livetime         = 1;
	double starttime        = 1;
	double current          = 1;
        //      Neutron info:
	int output_nMult = 1;
	bandhit output_nHit;
	//      Electron info:
        clashit output_eHit;
        //      Tagged info:
        taghit  output_tag;
        //      Event branches:
	outTree->Branch("Runno"         ,&Runno                 );
        outTree->Branch("Ebeam"         ,&Ebeam                 );
        outTree->Branch("gated_charge"  ,&gated_charge          );
        outTree->Branch("livetime"      ,&livetime              );
        outTree->Branch("starttime"     ,&starttime             );
        outTree->Branch("current"       ,&current               );
        //      Neutron branches:
	outTree->Branch("nMult"		,&output_nMult			);
        outTree->Branch("nHits"          ,&output_nHit                   );
        //      Electron branches:
        outTree->Branch("eHit"          ,&output_eHit                  );
        //      Tagged branches:
        outTree->Branch("tag"           ,&output_tag                    );

	TRandom3 * myRand = new TRandom3(0);


	// Define 4D binned ToF plots and bin edges and kin cuts
	double const min_Q2 		= 1.5;
	double const max_Q2 		= 11.;
	double const min_CosTheta_nq 	= -1.1;
	double const max_CosTheta_nq 	= -0.5;

/*
	int nBins_Q2		= 1;
	int nBins_CosTheta_nq 	= 1;
	TH1D *** hToF_4D	= new TH1D**[nBins_Q2];
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		hToF_4D[bin_Q2]	= new TH1D*[nBins_CosTheta_nq];
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			hToF_4D[bin_Q2][bin_CosTheta_nq] = new TH1D(Form("hToF_%i_%i",bin_Q2,bin_CosTheta_nq),"",1600,-100,300);
		}
	}
*/

	// Define background and signal edges
	const double bkgrd_min = -50;
	const double bkgrd_max = 0;
	const double signal_min = 13;
	const double signal_max = 63;

	// Loop over the neutron and electron skim files given
	TFile * inFile_e = new TFile(argv[2]);
	TFile * inFile_n = new TFile(argv[3]);
	TTree * inTree_e = (TTree*)inFile_e->Get("electrons");
//	TTree * inTree_n = (TTree*)inFile_n->Get("neutrons");
	TTree * inTree_n = (TTree*)inFile_n->Get("tagged");


	clashit* input_eHit = new clashit;
	TClonesArray* input_nHit = new TClonesArray("bandhit");
//	bandhit* input_nHit = new bandhit[maxNeutrons];
        int nMult;
	double input_ebeam;

	inTree_e->SetBranchAddress("Ebeam", &input_ebeam);
	inTree_e->SetBranchAddress("eHit",	&input_eHit);
	inTree_n->SetBranchAddress("nHits",	&input_nHit);
	inTree_n->SetBranchAddress("nMult",	&nMult);

	vector<double> n_theta;
	vector<double> n_phi;
	vector<double> n_dL;
	vector<double> n_Edep;
	vector<int> n_status;
	vector<double> e_theta;
	vector<double> e_phi;
	vector<double> e_mom;
	vector<double> e_beam;
	vector<double> e_vtz;
	// Read entire tree into memory due to issue with jumping around in tree as we read it..
	cout << "Saving neutrons to mem...\n";
	for( int neutron = 0 ; neutron < inTree_n->GetEntries() ; neutron++ ){

		double input_theta_n = 0;
		double input_phi_n = 0;
		double input_Edep_n = 0;
		double input_dL = 0;
		int input_status = 0;

		inTree_n->GetEntry(neutron);

		if (nMult != 1) continue;

		bandhit* this_nHit = (bandhit*)input_nHit->At(0);

		input_theta_n = this_nHit->getDL().Theta();
		input_phi_n = this_nHit->getDL().Phi();
		input_Edep_n = this_nHit->getEdep();
		input_dL = this_nHit->getDL().Mag();
		input_status = this_nHit->getStatus();

/*
		input_theta_n = input_nHit[0].getDL().Theta();
		input_phi_n = input_nHit[0].getDL().Phi();
		input_Edep_n = input_nHit[0].getEdep();
		input_dL = input_nHit[0].getDL().Mag();
		input_status = input_nHit[0].getStatus();
*/
		n_theta.push_back(input_theta_n);
		n_phi.push_back(input_phi_n);
		n_Edep.push_back(input_Edep_n);
		n_status.push_back(input_status);
		n_dL.push_back(input_dL);
	}
	cout << "Saving electrons to mem...\n";
	for( int electron = 0 ; electron < inTree_e->GetEntries() ; electron++ ){
		double input_p_e = 0;
		double input_theta_e = 0;
		double input_phi_e = 0;
		double input_vtz_e = 0;

		inTree_e->GetEntry(electron);

		input_p_e = input_eHit->getMomentum();
		input_theta_e = input_eHit->getTheta();
		input_phi_e = input_eHit->getPhi();
		input_vtz_e = input_eHit->getVtz();

		//if( electron > inTree_e->GetEntries()/100 ) break;

		e_theta.push_back(input_theta_e);
		e_phi.push_back(input_phi_e);
		e_mom.push_back(input_p_e);
		e_vtz.push_back(input_vtz_e);
		e_beam.push_back(input_ebeam);
	}

	// Loop over all neutron events
	for( int neutron = 0 ; neutron < n_status.size() ; neutron++ ){
		if( neutron % 1000 == 0 ) cout << "working on neutron " << neutron << "\n";

		int nPairs = 0; // counts how many pairs were made for this neutron
		int nFails = 0; // counts how many times a pair was tried for
		int maxPairs = 10;
		int maxFails = 10;
		while( nPairs < maxPairs && nFails < maxFails ){

			int electron = myRand->Rndm() * inTree_e->GetEntries();

			// Read in the stored variables
			double fixed_Ebeam 	= e_beam	[electron];
			double p_e		= e_mom		[electron];
			double theta_e		= e_theta	[electron];
			double phi_e		= e_phi		[electron];
			double vtz 		= e_vtz		[electron];
			double theta_n		= n_theta	[neutron];
			double phi_n		= n_phi		[neutron];
			double dL		= n_dL		[neutron];
			double Edep 		= n_Edep	[neutron];
			int status 		= n_status	[neutron];

			// Redraw my ToF and calculate all quantities to save
			double ToF = myRand->Rndm()*( signal_max - signal_min ) + signal_min ;
			double beta = dL / (ToF*cAir);
			double p_n = mN / sqrt( 1./pow(beta,2) - 1. );

			// Create vectors for calculating angles
			TVector3 beamVec(0,0,fixed_Ebeam);
			TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
			TVector3 qVec;	qVec = beamVec - eVec;
			TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);

			double E_e 	= sqrt( mE*mE + p_e*p_e );
			double q 	= qVec.Mag();
			double theta_q  = qVec.Theta();
			double phi_q 	= qVec.Phi();
			double nu 	= fixed_Ebeam - E_e;
			double Q2 	= q*q - nu*nu;
			double xB	= Q2 / (2.*mP*nu);
			double W2	= mP*mP - Q2 + 2*nu*mP;
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

			// Ask if passes cuts, and if so, count valid pair and fill histogram/tree
			if( Q2 > min_Q2 && Q2 < max_Q2 && CosTheta_nq > min_CosTheta_nq && CosTheta_nq < max_CosTheta_nq ){
				// Valid pair
				nFails = 0;
				nPairs++;
				Ebeam 	= 		fixed_Ebeam 	;

				output_eHit.setMomentum(p_e);
				output_eHit.setTheta(theta_e);
				output_eHit.setPhi(phi_e);
				output_eHit.setW2(W2);
				output_eHit.setQ2(Q2);
				output_eHit.setXb(xB);
				output_eHit.setOmega(nu);
				output_eHit.setQ(q);
				output_eHit.setThetaQ(theta_q);
				output_eHit.setPhiQ(phi_q);

				output_eHit.setPID(11);
				output_eHit.setCharge(-1);
				output_eHit.setEpcal(E_e);
				output_eHit.setVtz(vtz);

				output_eHit.setEoP(0.2);
				output_eHit.setU(10);
				output_eHit.setV(10);
				output_eHit.setW(10);

				output_nHit.setTof(ToF);
				output_nHit.setTofFadc(ToF);
				output_nHit.setX(dL * sin(theta_n) * cos(phi_n));
				output_nHit.setY(dL * sin(theta_n) * sin(phi_n));
				output_nHit.setZ(dL * cos(theta_n));
				output_nHit.setEdep(Edep);
				output_nHit.setStatus(status);

				output_tag.setWp(Wp);
				output_tag.setXp(Xp);
				output_tag.setAs(As);
				output_tag.setThetaNQ(theta_nq);
				output_tag.setPhiNQ(phi_nq);
				output_tag.setPt(Pt);
				output_tag.setXp2(Xp2);

				output_tag.setMomentumE(eVec);
				output_tag.setMomentumN(nVec);
				output_tag.setMomentumQ(qVec);
				output_tag.setMomentumB(beamVec);

				outTree->Fill();


			}
			else{ nFails++; }


		} // end loop over finding electrons

	}

	outFile->cd();
/*
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			hToF_4D[bin_Q2][bin_CosTheta_nq] -> Write();
		}
	}
*/
	outTree->Write();
	outFile->Close();



	return 0;
}
