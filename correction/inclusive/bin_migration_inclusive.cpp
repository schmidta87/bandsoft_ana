#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"

#include "genpart.h"
#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;
void calculate_AsPt( double & As, double & Pt , const double Ebeam, genpart * const electron, genpart * const neutron );
void fillHist( double Q2, double Xb, TH1D* hist_lowQ2 , TH1D* hist_highQ2 );

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File]\n";
		return -1;
	}

	TFile inFile(argv[1]);
	TTree * inTree = (TTree*) inFile.Get("electrons");

	TClonesArray* gen_particles 	= new TClonesArray("genpart");
	clashit * rec_electron 		= new clashit;
	int 	genMult 	= 0;
	double	Ebeam		= 0;
	inTree->SetBranchAddress("mcParts"		, &gen_particles		);
	inTree->SetBranchAddress("genMult"		, &genMult			);
	inTree->SetBranchAddress("eHit"			, &rec_electron			);
	inTree->SetBranchAddress("Ebeam"		, &Ebeam			);

	TFile * outFile = new TFile("test_binmigration.root","RECREATE");
	TH1D * gen_xb_lowQ2  = new TH1D("gen_xb_lowQ2", "gen_xb_lowQ2", bins_Xb,Xb_min,Xb_max);
	TH1D * rec_xb_lowQ2  = new TH1D("rec_xb_lowQ2", "rec_xb_lowQ2", bins_Xb,Xb_min,Xb_max);
	TH1D * gen_xb_highQ2 = new TH1D("gen_xb_highQ2","gen_xb_highQ2",bins_Xb,Xb_min,Xb_max);
	TH1D * rec_xb_highQ2 = new TH1D("rec_xb_highQ2","rec_xb_highQ2",bins_Xb,Xb_min,Xb_max);

	for( int event = 0 ; event < inTree->GetEntries() ; ++event ){
		gen_particles	->Clear();
		rec_electron	->Clear();
		genMult 	= 0;
		Ebeam 		= 0;
		inTree->GetEntry(event);

		// Get the reconstructed variables
		double rec_Q2		= rec_electron->getQ2();
		double rec_Xb		= rec_electron->getXb();

		// Form As, ThetaNQ, W, X with generated:
		genpart * gen_electron 	= (genpart*) gen_particles->At(0);
		double	gen_Q2		= gen_electron->getQ2();
		double 	gen_Xb		= gen_electron->getXb();

		// sort the generated values:
		fillHist( gen_Q2, gen_Xb, gen_xb_lowQ2 , gen_xb_highQ2 );

		// sort the reconstructed values:
		fillHist( rec_Q2, rec_Xb, rec_xb_lowQ2 , rec_xb_highQ2 );


	}

	// Print the background subtracted distributions to a pdf file:
	TCanvas * c 	= new TCanvas("c",	"",800,600);
	c ->Divide(1,2);
	outFile->cd();

	gen_xb_lowQ2->Write();
	rec_xb_lowQ2->Write();

	rec_xb_lowQ2->Divide( gen_xb_lowQ2 );
	c->cd(1);

	rec_xb_lowQ2->SetLineWidth(3);
	rec_xb_lowQ2->SetMarkerColor(4);
	rec_xb_lowQ2->SetMaximum(2);
	rec_xb_lowQ2->SetMinimum(0);
	rec_xb_lowQ2->Draw("*P");
	TLine * hline = new TLine(Xb_min,1,Xb_max,1);
	hline->SetLineWidth(2);
	hline->SetLineColor(1);
	hline->SetLineStyle(2);
	hline->Draw("SAME");


	gen_xb_highQ2->Write();
	rec_xb_highQ2->Write();

	rec_xb_highQ2->Divide( gen_xb_highQ2 );
	c->cd(2);

	rec_xb_highQ2->SetLineWidth(3);
	rec_xb_highQ2->SetMarkerColor(4);
	rec_xb_highQ2->SetMaximum(2);
	rec_xb_highQ2->SetMinimum(0);
	rec_xb_highQ2->Draw("*P");
	hline->Draw("SAME");

	c->Print("binmigration_inclusive.pdf");

	outFile->Close();

	inFile.Close();

	return 1;
}


void calculate_AsPt( double & As, double & Pt , const double Ebeam, genpart * const electron, genpart * const neutron ){

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; 
	eVec.SetMagThetaPhi( electron->getMomentum(), electron->getTheta(), electron->getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;
	TVector3	nVec; 
	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	nVec.Unit();


	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_nq 		= norm_scatter.Angle( norm_reaction );
	double theta_nq 	= nVec.Angle( qVec );
	double gen_CosThetaNQ 	= cos(theta_nq);

	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_nq *= (-1);
	}

	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	double p_n = neutron->getMomentum();

	double E_n 	= sqrt( mN*mN + p_n*p_n );
	double W_primeSq = mD*mD - electron->getQ2() + mN*mN + 2.*mD*(electron->getOmega()-E_n) - 2.*electron->getOmega()*E_n + 2.*electron->getQ()*p_n*cos(theta_nq);
	double Wp = sqrt(W_primeSq);
	double Xp = electron->getQ2()/(2.*( electron->getOmega()*(mD-E_n) + p_n*electron->getQ()*gen_CosThetaNQ));
	double Xp2 = electron->getQ2()/(W_primeSq - mN*mN + electron->getQ2());

	TVector3 Pt_vec;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;

	As = (E_n - p_n*gen_CosThetaNQ)/mN;
	Pt_vec = nVec - pN_par_q;
	Pt = Pt_vec.Mag();
}

void fillHist( double Q2, double Xb, TH1D* hist_lowQ2 , TH1D* hist_highQ2 ){

	if( Q2 >= Q2Bins[0] && Q2 < Q2Bins[1] ){
		hist_lowQ2->Fill( Xb );
		return;
	}
	else if(Q2 >= Q2Bins[1] && Q2 < Q2Bins[2] ){
		hist_highQ2->Fill( Xb );
		return;
	}
	else return;
}
