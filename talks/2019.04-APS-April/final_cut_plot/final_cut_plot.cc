///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////		final_cut_plot.cxx
//////		Nov 2018, Brian Clark 
//////		Preparing to make the slanted cut data
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TF1.h"
#include "TExec.h"
#include "TColor.h"
#include "TPaletteAxis.h"

using namespace std;

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(11);

	if(argc!=3){
		printf("Invalid usage! Usage is ./SlantedCut <slope> <cut_val> \n");
		return -1;
	}

	TH2D *PeakCorr_vs_SNR[2];
	PeakCorr_vs_SNR[0]=new TH2D("","V",50,0,1,75,0,30);
	PeakCorr_vs_SNR[1]=new TH2D("","H",50,0,1,75,0,30);

	TH2D *PeakCorr_vs_SNR_sim[2];
	PeakCorr_vs_SNR_sim[0]=new TH2D("","V",50,0,1,75,0,30);
	PeakCorr_vs_SNR_sim[1]=new TH2D("","H",50,0,1,75,0,30);

	TH1D *Counts_vs_SNR[2];
	Counts_vs_SNR[0]=new TH1D("","V",200,0,20);
	Counts_vs_SNR[1]=new TH1D("","H",200,0,20);

	//load the data

	TFile *inputFile_data = TFile::Open("data_vals_for_cuts_A2_2013.root");
	if(!inputFile_data){
		cout<<"Can't open vals for cuts file data!"<<endl;
		return -1;
	}
	TTree *tree_data = (TTree*) inputFile_data->Get("CutTree");
	double corr_val_data[2];
	double snr_val_data[2];
	int pass_ana_data[2];
	tree_data->SetBranchAddress("corr_val", &corr_val_data);
	tree_data->SetBranchAddress("snr_val", &snr_val_data);
	tree_data->SetBranchAddress("pass_ana", &pass_ana_data);

	//load the sim

	TFile *inputFile_sim = TFile::Open("sim_vals_for_cuts.root");
	if(!inputFile_sim){
		cout<<"Can't open vals for cuts file sim!"<<endl;
		return -1;
	}
	TTree *tree_sim = (TTree*) inputFile_sim->Get("CutTree");
	double corr_val_sim[2];
	double snr_val_sim[2];
	int pass_ana_sim[2];
	double weight_out;
	double energy_out;
	double VPeak_out[16];
	tree_sim->SetBranchAddress("corr_val", &corr_val_sim);
	tree_sim->SetBranchAddress("snr_val", &snr_val_sim);
	tree_sim->SetBranchAddress("pass_ana", &pass_ana_sim);
	tree_sim->SetBranchAddress("weight", &weight_out);
	tree_sim->SetBranchAddress("energy", &energy_out);
	tree_sim->SetBranchAddress("pass_ana", &pass_ana_sim);
	tree_sim->SetBranchAddress("VPeak", &VPeak_out);

	vector <double> corr_vals_V;
	vector <double> snr_vals_V;
	vector <bool> cut_already_V;

	vector <double> corr_vals_H;
	vector <double> snr_vals_H;
	vector <bool> cut_already_H;

	//first, loop data
	int numEntries = tree_data->GetEntries();
	for(int entry=0; entry<numEntries; entry++){

		tree_data->GetEvent(entry);

		for(int pol=0; pol<2; pol++){
			if(pass_ana_data[pol]){
				PeakCorr_vs_SNR[pol]->Fill(corr_val_data[0],snr_val_data[0]);
				if(pol==0){
					corr_vals_V.push_back(corr_val_data[pol]);
					snr_vals_V.push_back(snr_val_data[pol]);
					cut_already_V.push_back(false);
				}
				else if(pol==1){
					corr_vals_H.push_back(corr_val_data[pol]);
					snr_vals_H.push_back(snr_val_data[pol]);
					cut_already_H.push_back(false);
				}
			}
		}
	}

	//first, loop sim
	int numEntries_sim = tree_data->GetEntries();
	cout<<"num entrie sim "<<numEntries_sim<<endl;
	for(int entry=0; entry<numEntries_sim; entry++){
		tree_sim->GetEvent(entry);
		for(int pol=0; pol<2; pol++){
			if(pass_ana_sim[pol]){
				double snr = snr_val_sim[0];
				if(snr>20) snr=19.75;
				PeakCorr_vs_SNR_sim[pol]->Fill(corr_val_sim[0],snr,weight_out);
			}
		}
	}

	double y_int = 7.;
	double slope= double(atof(argv[1]));
	printf("Slope this time is %.2f \n",slope);

	//now, slide the line up from y-int of zero, counting cuts as we go
	for(double test_y_int=0.;test_y_int<15; test_y_int+=0.1){
		for(int count=0; count<snr_vals_V.size(); count++){
			if(!cut_already_V[count]){
				double this_y_val = (slope * corr_vals_V[count] ) + test_y_int;
				if(snr_vals_V[count]<this_y_val){
					cut_already_V[count]=true;
					Counts_vs_SNR[0]->Fill(test_y_int);
				}
			}
		}
	}
	int max_bin = Counts_vs_SNR[0]->GetMaximumBin();
	double start_of_fit = Counts_vs_SNR[0]->GetBinCenter(max_bin+2);
	int last_filled_bin = Counts_vs_SNR[0]->FindLastBinAbove(0.,1);
	double last_bin  = Counts_vs_SNR[0]->GetBinCenter(last_filled_bin);
	printf("Max bin is %d , two bins away bin center is %.2f , last filled bin center is %.2f \n", max_bin, start_of_fit,last_bin);

	//now to do exponential fit
	char equation[150];
	sprintf(equation,"(exp([0] - [1]*x))");
	TF1 *fit = new TF1("ExpFit",equation,start_of_fit,last_bin);
	fit->SetParameter(1,3.58);
	Counts_vs_SNR[0]->Fit("ExpFit","R");
	printf("Chi-Square/NDF %.2f / %.2f \n",fit->GetChisquare(),double(fit->GetNDF()));

	gStyle->SetOptStat(0);
	TCanvas *c_exp = new TCanvas("","",1100,850);
		Counts_vs_SNR[0]->Draw("");
		gPad->SetLogy();
		Counts_vs_SNR[0]->GetYaxis()->SetTitle("Events");
		Counts_vs_SNR[0]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		Counts_vs_SNR[0]->GetXaxis()->SetRangeUser(4,10);
		Counts_vs_SNR[0]->SetTitle("");
	c_exp->SaveAs("exp_dist.pdf");

	double cut = double(atof(argv[2]));
	double amplitude = exp(fit->GetParameter(0))*50.;
	double index= fit->GetParameter(1);
	double back_estimate = amplitude/index * exp(-index * cut);
	printf("Background estimate: %.2f \n", back_estimate);

	vector <double> x_vals_for_line;
	vector <double> y_vals_for_line;
	for(double x=0; x<1; x+=0.01){
		double y_val = (slope * x ) + cut;
		x_vals_for_line.push_back(x);
		y_vals_for_line.push_back(y_val);
	}
	TGraph *cut_line = new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]);

	char save_title[300];

	TCanvas *c3 = new TCanvas("aa","aa",3*1000,3*670);
	PeakCorr_vs_SNR_sim[0]->Draw("colz");
		TExec *ex1 = new TExec("ex1","gStyle->SetPalette(kDeepSea); ");
		ex1->Draw();
		PeakCorr_vs_SNR_sim[0]->Draw("colz same");
	PeakCorr_vs_SNR_sim[0]->GetZaxis()->SetRangeUser(1.,200.);
	PeakCorr_vs_SNR_sim[0]->GetYaxis()->SetRangeUser(2.,20.);
	PeakCorr_vs_SNR_sim[0]->GetXaxis()->SetRangeUser(0.,0.8);
	gPad->SetLogz();
	gPad->Modified();
	gPad->Update();

	TPaletteAxis *palette1 = (TPaletteAxis*)PeakCorr_vs_SNR_sim[0]->GetListOfFunctions()->FindObject("palette");
	palette1->SetX1NDC(0.66);
	palette1->SetX2NDC(0.71);
	palette1->SetY1NDC(0.1);
	palette1->SetY2NDC(0.95);

	// gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.35);
	gPad->SetTopMargin(0.05);

	PeakCorr_vs_SNR[0]->Draw("same colz");
		TExec *ex2 = new TExec("ex2","gStyle->SetPalette(kSolar);"); // TColor::InvertPalette();
		ex2->Draw();
		PeakCorr_vs_SNR[0]->Draw("same colz");
	gPad->SetLogz();
	gPad->Modified();
	gPad->Update();

	TPaletteAxis *palette2 = (TPaletteAxis*)PeakCorr_vs_SNR[0]->GetListOfFunctions()->FindObject("palette");
	palette2->SetX1NDC(0.86);
	palette2->SetX2NDC(0.91);
	palette2->SetY1NDC(0.1);
	palette2->SetY2NDC(0.95);

	cut_line->Draw("same");
	cut_line->SetLineColor(kRed);
	cut_line->SetLineWidth(4);

	PeakCorr_vs_SNR_sim[0]->SetTitle("");
	PeakCorr_vs_SNR_sim[0]->GetXaxis()->SetTitle("Peak Correlation Value");
	PeakCorr_vs_SNR_sim[0]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");

	PeakCorr_vs_SNR_sim[0]->GetZaxis()->SetTitle("Simulated Neutrinos (Kotera Max)");
	PeakCorr_vs_SNR[0]->GetZaxis()->SetTitle("A2 Measured Background (2013 10\%)");

	PeakCorr_vs_SNR_sim[0]->GetYaxis()->SetTitleSize(0.05);
	PeakCorr_vs_SNR_sim[0]->GetXaxis()->SetTitleSize(0.05);
	PeakCorr_vs_SNR_sim[0]->GetYaxis()->SetLabelSize(0.045);
	PeakCorr_vs_SNR_sim[0]->GetXaxis()->SetLabelSize(0.045);

	PeakCorr_vs_SNR_sim[0]->GetZaxis()->SetTitleSize(0.045);
	PeakCorr_vs_SNR[0]->GetZaxis()->SetTitleSize(0.045);
	PeakCorr_vs_SNR_sim[0]->GetZaxis()->SetLabelSize(0.04);
	PeakCorr_vs_SNR[0]->GetZaxis()->SetLabelSize(0.04);
	
	gStyle->SetOptStat(0);
	sprintf(save_title, "corr_vs_snr_data_and_noise.pdf",year_now, month_now, day_now);
	c3->SaveAs(save_title);

}
