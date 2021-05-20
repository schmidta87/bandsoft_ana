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
void fillHist( double Q2, double As, double Pt, double Xb, TH1D** hists_lowQ2 , TH1D** hists_highQ2 );

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File]\n";
		return -1;
	}

	TFile inFile(argv[1]);
	TTree * inTree = (TTree*) inFile.Get("tagged");

	TClonesArray* gen_particles 	= new TClonesArray("genpart");
	clashit * rec_electron 		= new clashit;
	TClonesArray* rec_tagged 	= new TClonesArray("taghit");
	TClonesArray* rec_neutrons 	= new TClonesArray("bandhit");
	int 	genMult 	= 0;
	int	nleadindex 	= 0;
	double	Ebeam		= 0;
	inTree->SetBranchAddress("mcParts"		, &gen_particles		);
	inTree->SetBranchAddress("genMult"		, &genMult			);
	inTree->SetBranchAddress("eHit"			, &rec_electron			);
	inTree->SetBranchAddress("tag"			, &rec_tagged			);
	inTree->SetBranchAddress("nHits"		, &rec_neutrons			);
	inTree->SetBranchAddress("nleadindex"		, &nleadindex			);
	inTree->SetBranchAddress("Ebeam"		, &Ebeam			);

	TFile * outFile = new TFile("test_binmigration.root","RECREATE");
	TH1D ** gen_xb_lowQ2 	= new TH1D*[lowQ2_bins];
	TH1D ** gen_xb_highQ2 	= new TH1D*[highQ2_bins];
	TH1D ** rec_xb_lowQ2 	= new TH1D*[lowQ2_bins];
	TH1D ** rec_xb_highQ2 	= new TH1D*[highQ2_bins];
	for( int i = 0 ; i < lowQ2_bins ; ++i){
		gen_xb_lowQ2[i] = new TH1D(Form("gen_xb_lowQ2_%i",i),Form("gen_xb_lowQ2_%i",i),bins_Xb,Xb_min,Xb_max);
		rec_xb_lowQ2[i] = new TH1D(Form("rec_xb_lowQ2_%i",i),Form("rec_xb_lowQ2_%i",i),bins_Xb,Xb_min,Xb_max);
	}
	for( int i = 0 ; i < highQ2_bins; ++i){
		gen_xb_highQ2[i] = new TH1D(Form("gen_xb_highQ2_%i",i),Form("gen_xb_highQ2_%i",i),bins_Xb,Xb_min,Xb_max);
		rec_xb_highQ2[i] = new TH1D(Form("rec_xb_highQ2_%i",i),Form("rec_xb_highQ2_%i",i),bins_Xb,Xb_min,Xb_max);
	}

	for( int event = 0 ; event < inTree->GetEntries() ; ++event ){
		gen_particles	->Clear();
		rec_electron	->Clear();
		rec_neutrons	->Clear();
		rec_tagged	->Clear();
		genMult 	= 0;
		nleadindex	= 0;
		Ebeam 		= 0;
		//if(event > 1000) break;
		inTree->GetEntry(event);

		// Get the reconstructed variables
		taghit * rec_tag	= (taghit*)  rec_tagged->At(nleadindex);
		bandhit * rec_neutron	= (bandhit*) rec_neutrons->At(nleadindex);
		double rec_As		= rec_tag->getAs();
		double rec_Pt 		= rec_tag->getPt().Mag();
		double rec_Q2		= rec_electron->getQ2();
		double rec_Xb		= rec_electron->getXb();

		// Form As, ThetaNQ, W, X with generated:
		genpart * gen_electron 	= (genpart*) gen_particles->At(0);
		genpart * gen_neutron 	= (genpart*) gen_particles->At(1);
		double	gen_Q2		= gen_electron->getQ2();
		double 	gen_Xb		= gen_electron->getXb();

		double gen_As, gen_Pt;
		calculate_AsPt( gen_As, gen_Pt, Ebeam, gen_electron, gen_neutron );

		// sort the generated values:
		fillHist( gen_Q2, gen_As, gen_Pt, gen_Xb, gen_xb_lowQ2 , gen_xb_highQ2 );

		// sort the reconstructed values:
		fillHist( rec_Q2, rec_As, rec_Pt, rec_Xb, rec_xb_lowQ2 , rec_xb_highQ2 );


	}

	// Print the background subtracted distributions to a pdf file:
	TCanvas * c_lowQ2 	= new TCanvas("c_lowQ2",	"",800,600);
	TCanvas * c_highQ2 	= new TCanvas("c_highQ2",	"",800,600);
	c_lowQ2 ->Divide(5,5);
	c_highQ2->Divide(5,5);
	outFile->cd();

	for( int i = 0 ; i < lowQ2_bins ; ++i){
		gen_xb_lowQ2[i]->Write();
		rec_xb_lowQ2[i]->Write();

		if( gen_xb_lowQ2[i]->Integral() < 1 || rec_xb_lowQ2[i]->Integral() < 1 ) continue;
		rec_xb_lowQ2[i]->Divide( gen_xb_lowQ2[i] );
		c_lowQ2->cd(i+1);

		rec_xb_lowQ2[i]->SetLineWidth(3);
		rec_xb_lowQ2[i]->SetMarkerColor(4);
		rec_xb_lowQ2[i]->SetMaximum(2);
		rec_xb_lowQ2[i]->SetMinimum(0);
		rec_xb_lowQ2[i]->Draw("*P");
		TLine * hline = new TLine(Xb_min,1,Xb_max,1);
		hline->SetLineWidth(2);
		hline->SetLineColor(1);
		hline->SetLineStyle(2);
		hline->Draw("SAME");
	}
	c_lowQ2->Print("binmigration_lowQ2.pdf");

	for( int i = 0 ; i < highQ2_bins ; ++i){
		gen_xb_highQ2[i]->Write();
		rec_xb_highQ2[i]->Write();

		if( gen_xb_highQ2[i]->Integral() < 1 || rec_xb_highQ2[i]->Integral() < 1 ) continue;
		rec_xb_highQ2[i]->Divide( gen_xb_highQ2[i] );
		c_highQ2->cd(i+1);

		rec_xb_highQ2[i]->SetLineWidth(3);
		rec_xb_highQ2[i]->SetMarkerColor(4);
		rec_xb_highQ2[i]->SetMaximum(2);
		rec_xb_highQ2[i]->SetMinimum(0);
		rec_xb_highQ2[i]->Draw("*P");
		TLine * hline = new TLine(Xb_min,1,Xb_max,1);
		hline->SetLineWidth(2);
		hline->SetLineColor(1);
		hline->SetLineStyle(2);
		hline->Draw("SAME");
	}
	c_highQ2->Print("binmigration_highQ2.pdf");

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

void fillHist( double Q2, double As, double Pt, double Xb, TH1D** hists_lowQ2 , TH1D** hists_highQ2 ){

	if( Q2 >= Q2Bins[0] && Q2 < Q2Bins[1] ){
		for( int this_bin = 0 ; this_bin < lowQ2_bins ; ++this_bin ){
			double As_low  = lowQ2_AsPtBins[this_bin][0];
			double As_high = lowQ2_AsPtBins[this_bin][1];
			double Pt_low  = lowQ2_AsPtBins[this_bin][2];
			double Pt_high = lowQ2_AsPtBins[this_bin][3];
			
			if( As < As_low 	) continue;
			if( As >= As_high	) continue;
			if( Pt < Pt_low		) continue;
			if( Pt >= Pt_high	) continue;

			hists_lowQ2[this_bin]->Fill( Xb );
			return;
		}
	}
	else if(Q2 >= Q2Bins[1] && Q2 < Q2Bins[2] ){
		for( int this_bin = 0 ; this_bin < highQ2_bins ; ++this_bin ){
			double As_low  = highQ2_AsPtBins[this_bin][0];
			double As_high = highQ2_AsPtBins[this_bin][1];
			double Pt_low  = highQ2_AsPtBins[this_bin][2];
			double Pt_high = highQ2_AsPtBins[this_bin][3];
			
			if( As < As_low 	) continue;
			if( As >= As_high	) continue;
			if( Pt < Pt_low		) continue;
			if( Pt >= Pt_high	) continue;

			hists_highQ2[this_bin]->Fill( Xb );
			return;
		}

	}
	else return;
}
