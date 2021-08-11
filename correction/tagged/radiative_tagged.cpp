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
void fillHist( double Q2, double Pt, double Xb, double As, double weight, TH1D**** hist , TH1D**** counts );
void setError( double * err , TH1D * ravg , TH1D * cnts );

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File]\n";
		return -1;
	}

	TFile inFile(argv[1]);
	TTree * inTree = (TTree*) inFile.Get("tagged");

	TClonesArray* gen_particles 	= new TClonesArray("genpart");
	int 	genMult 	= 0;
	double	Ebeam		= 0;
	double	rad_weight	= 0;
	inTree->SetBranchAddress("mcParts"		, &gen_particles		);
	inTree->SetBranchAddress("genMult"		, &genMult			);
	inTree->SetBranchAddress("Ebeam"		, &Ebeam			);
	inTree->SetBranchAddress("weight"		, &rad_weight			);

	TFile * outFile = new TFile("radcorrection_tagged.root","RECREATE");

	TH1D **** h4_gen_as 	= new TH1D***[bins_Q2];
	TH1D **** h4_cnt_as	= new TH1D***[bins_Q2];

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_gen_as[i]	= new TH1D**[bins_Pt];
		h4_cnt_as[i]	= new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_gen_as[i][j] 	= new TH1D*[bins_Xb];
			h4_cnt_as[i][j] 	= new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){ // bins in Xb
				h4_gen_as[i][j][k] = new TH1D(Form("radcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),Form("radcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),bins_As,As_min,As_max);  // histogram in As
				h4_cnt_as[i][j][k] = new TH1D(Form("cnt_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),Form("cnt_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),bins_As,As_min,As_max);  // histogram in As
			}
		}
	}

	for( int event = 0 ; event < inTree->GetEntries() ; ++event ){
		gen_particles	->Clear();
		genMult 	= 0;
		Ebeam 		= 0;
		rad_weight	= 0;
		//if(event > 1000) break;
		inTree->GetEntry(event);

		// Form As, ThetaNQ, W, X with generated:
		genpart * gen_electron 	= (genpart*) gen_particles->At(0);
		genpart * gen_neutron 	= (genpart*) gen_particles->At(1);
		double	gen_Q2		= gen_electron->getQ2();
		double 	gen_Xb		= gen_electron->getXb();

		double gen_As, gen_Pt;
		calculate_AsPt( gen_As, gen_Pt, Ebeam, gen_electron, gen_neutron );

		// sort the generated values:
		fillHist( gen_Q2, gen_Pt, gen_Xb, gen_As, rad_weight, h4_gen_as , h4_cnt_as );

	}

	// Print the background subtracted distributions to a pdf file:
	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		c_Q2[i]->Divide( bins_Pt, bins_Xb );
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bins

				if( h4_gen_as[i][j][k]->Integral() < 1 || h4_cnt_as[i][j][k]->Integral() < 1 ) continue;
				double * errors = new double[bins_As];
				setError( errors , h4_gen_as[i][j][k] , h4_cnt_as[i][j][k] );
				h4_gen_as[i][j][k]->Divide( h4_cnt_as[i][j][k] );

				for( int bin = 1 ; bin < h4_gen_as[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					h4_gen_as[i][j][k]->SetBinError( bin , errors[bin-1] );
				}

				c_Q2[i]->cd( (k*2)+1 + j );

				h4_gen_as[i][j][k]->SetLineWidth(3);
				h4_gen_as[i][j][k]->SetMarkerColor(4);
				h4_gen_as[i][j][k]->SetMaximum(2);
				h4_gen_as[i][j][k]->SetMinimum(0);
				h4_gen_as[i][j][k]->Draw("*P");
				TLine * hline = new TLine(As_min,1,As_max,1);
				hline->SetLineWidth(2);
				hline->SetLineColor(1);
				hline->SetLineStyle(2);
				hline->Draw("SAME");

				h4_gen_as[i][j][k]->Write();
				//h4_cnt_as[i][j][k]->Write();
				delete[] errors;

			}
		}
		c_Q2[i]->Print(Form("radiative_tagged_Q2_%i.pdf",i));
	}
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

void fillHist( double Q2, double Pt, double Xb, double As, double weight, TH1D**** hist , TH1D**** counts ){
	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Need to figure out which bin of Q2, Pt, Xb (then fill as As):
	int this_bin_q2 = -1;
	int this_bin_pt = -1;
	int this_bin_xb = -1;
	
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}

	if( Pt < Pt_min			) return;
	if( Pt > Pt_max			) return;
	if( Xb < Xb_min[this_bin_q2]			) return;
	if( Xb > Xb_min[this_bin_q2]+Xb_step*6		) return;

	this_bin_pt = (int) ((Pt - Pt_min)/Pt_step);
	this_bin_xb = (int) ((Xb - Xb_min[this_bin_q2])/Xb_step);

	// Safety clauses
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }
	if( this_bin_pt == -1 ){ cerr << "how\n"; return; }
	if( this_bin_xb == -1 ){ cerr << "how\n"; return; }
	if( this_bin_q2 > bins_Q2-1 ){ cerr << "how\n"; return; }
	if( this_bin_pt > bins_Pt-1 ){ cerr << "how\n"; return; }
	if( this_bin_xb > bins_Xb-1 ){ cerr << "how\n"; return; }

	// For radiative weight, we cummulatively add the content in the same bin:
	int fill_bin = hist[this_bin_q2][this_bin_pt][this_bin_xb]->GetXaxis()->FindBin( As );
	double hist_curr_content = hist[this_bin_q2][this_bin_pt][this_bin_xb]->GetBinContent( fill_bin );
	double counts_curr_content = counts[this_bin_q2][this_bin_pt][this_bin_xb]->GetBinContent( fill_bin );
	
	hist[this_bin_q2][this_bin_pt][this_bin_xb]->SetBinContent( fill_bin, hist_curr_content + weight );
	counts[this_bin_q2][this_bin_pt][this_bin_xb]->SetBinContent( fill_bin , counts_curr_content + 1 );

}

void setError( double * err , TH1D * ravg , TH1D * cnts ){
	
	for( int bin = 1 ; bin < ravg->GetXaxis()->GetNbins(); ++bin ){

		double R = ravg->GetBinContent(bin);
		double N = cnts->GetBinContent(bin);

		double e = sqrt( R/pow(N,3./2) );
		cout << R << " " << N << " " << e << "\n";
		if( e!=e ) e = 0;
		err[bin-1] = e;
	}


	return;
}
