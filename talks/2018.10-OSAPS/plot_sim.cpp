///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////		Making some nice plots for OSPAS talk; this is for *simulation*
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
#include "TMath.h"

//AraRoot includes
//#include "PlottingFns.h"
//#include "RecoFns.h"

using namespace std;

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(0);
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <joined filename>"<<endl;;
		return 0;
	}
	int station = atoi(argv[1]);

	//just to have the cut parameters up front and easy to find
	int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

	TH2D *peakcorrv_vs_snr_final[2];
	peakcorrv_vs_snr_final[0]=new TH2D("","V",20,0,20,50,0,0.5);
	peakcorrv_vs_snr_final[1]=new TH2D("","H",20,0,20,50,0,0.5);

	TH2D *snr_vs_wrms[2];
	snr_vs_wrms[0]=new TH2D("","V",90,-5,4,400,0,40);
	snr_vs_wrms[1]=new TH2D("","H",90,-5,4,400,0,40);
	
	int num_total=0;
	int num_WFRMS[2]={0};

	for(int file_num=2; file_num<argc; file_num++){
		
		cout << "Run " << file_num << " :: " << argv[file_num] << endl;
		
		//first, load in the data file
		//this shoud be a "joined" file
		//meaning it should contain "filter" trees and "reco" trees
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		
		//next, we need to load the filter tree
		ss.str("");
		ss << "OutputTree_filter";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		double thirdVPeakOverRMS[3]; //the third highest vpeak over RMS
		double rms_pol_thresh_face[2][15][12];
		bool isCalPulser;
		bool isSoftTrigger;
		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face", &rms_pol_thresh_face);
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);

		//next, load the reco tree
		TTree *inputTree_reco[35];
		double peakCorr[35][2];
		int peakTheta[35][2];
		int peakPhi[35][2];
		int recoBinSelect = 19; //300 m map
		for(int i=0; i<35; i++){
			if(i==recoBinSelect){
				ss.str("");
				ss << "OutputTree_recoRadius_" << i;
				inputTree_reco[i] = (TTree*) inputFile->Get(ss.str().c_str());
				if(!inputTree_reco[i]) {
					std::cout << "Can't find OutputTree: " << i << "\n";
					return -1;
				}
				inputTree_reco[i]->SetBranchAddress("peakCorr_single", &peakCorr[i]);
				inputTree_reco[i]->SetBranchAddress("peakTheta_single", &peakTheta[i]);
				inputTree_reco[i]->SetBranchAddress("peakPhi_single", &peakPhi[i]);
			}
		}

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries%20000;
		if(starEvery==0) starEvery++;
		
		//now to loop over events
		for(int event=0; event<inputTree_filter->GetEntries(); event++){
		//for(int event=0; event<25; event++){

			if(event%starEvery==0) {
				std::cout << "	On event "<<event<<endl;
			}

			inputTree_filter->GetEvent(event);

			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			num_total++;
			
			for (int i = 0; i < 35; i++){
				if (i == recoBinSelect){
					inputTree_reco[i]->GetEntry(event);
				}
			}

			//figure out which reconstruction map (vpol or hpol) is best
			//in the present analysis, this only matters for the 300m bin

			double bestCorr[] = {0., 0., 0.};
			int bestCorrRadiusBin[3];
			int bestPol = 2;
			int bestTheta[3];
			int bestPhi[3];

			for(int pol=0; pol<2; pol++){
				for(int i=0; i<35; i++){
					if(i==recoBinSelect){
						if(peakCorr[i][pol] > bestCorr[pol]){
							bestCorr[pol]=peakCorr[i][pol];
							bestCorrRadiusBin[pol]=i;
							bestTheta[pol]=peakTheta[i][pol];
							bestPhi[pol]=peakPhi[i][pol];
						}
						if(peakCorr[i][pol] > bestCorr[2]){
							bestCorr[2]=peakCorr[i][pol];
							bestCorrRadiusBin[2]=i;
							bestTheta[2]=peakTheta[i][pol];
							bestPhi[2]=peakPhi[i][pol];
							bestPol=pol;
						}
					}//300m bin check
				}//loop over reco bins
			}//loop over polarizations

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			if(SNRs[0]>29.) SNRs[0]=29.;
			if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double>  rms_faces_V;
			rms_faces_V.resize(12);
			vector <double> rms_faces_H;
			rms_faces_H.resize(12);

			//now, we must loop over the faces
			for(int i=0; i<12; i++){
				rms_faces_V[i] = rms_pol_thresh_face[0][thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
				rms_faces_H[i] = rms_pol_thresh_face[1][thresholdBin_pol[1]][i];
			}

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];

			if(log(bestFaceRMS[0])/log(10) >= wavefrontRMScut[0]) failWavefrontRMS[0]=true;
			if(log(bestFaceRMS[1])/log(10) >= wavefrontRMScut[1]) failWavefrontRMS[1]=true;

			if(failWavefrontRMS[0]) num_WFRMS[0]++;
			if(failWavefrontRMS[1]) num_WFRMS[1]++;

			for(int pol=0; pol<2; pol++){
				peakcorrv_vs_snr_final[pol]->Fill(SNRs[pol],bestCorr[pol]);
				snr_vs_wrms[pol]->Fill(TMath::Log10(bestFaceRMS[pol]),SNRs[pol]);
				
				// //cut WRMS
				// if(!failWavefrontRMS[pol]){
				// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Fill(SNRs[pol],bestCorr[pol]);	
				// }//loop over polarizations
			
			}//loop over events
		}
		inputFile->Close();
		delete inputFile;
	}

	printf("Num total is %d \n", num_total);

	gStyle->SetOptStat(11);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	TCanvas *cfinal = new TCanvas("","",1100,850);
	peakcorrv_vs_snr_final[0]->Draw("colz");
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetTitle("Correlation Value");
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetTitle("Signal-to-Noise Ratio");
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetTitle("Number of Events");	
	peakcorrv_vs_snr_final[0]->SetTitle("");
	
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetTitleFont(43);
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetTitleSize(40);
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetTitleOffset(1.0);

	peakcorrv_vs_snr_final[0]->GetXaxis()->SetLabelFont(43);
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetLabelSize(35);
	peakcorrv_vs_snr_final[0]->GetXaxis()->SetLabelOffset(0.005);
	
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetTitleFont(43);
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetTitleSize(40);
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetTitleOffset(1.3);

	peakcorrv_vs_snr_final[0]->GetYaxis()->SetLabelFont(43);
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetLabelSize(35);
	peakcorrv_vs_snr_final[0]->GetYaxis()->SetLabelOffset(0.015);
	
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetTitleFont(43);
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetTitleSize(40);
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetTitleOffset(1.0);

	peakcorrv_vs_snr_final[0]->GetZaxis()->SetLabelFont(43);
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetLabelSize(35);
	peakcorrv_vs_snr_final[0]->GetZaxis()->SetLabelOffset(0.005);


	gPad->SetLogz();
	TVirtualPad* p1 = cfinal->cd();
	p1->SetRightMargin(0.15);
	p1->SetTopMargin(0.05);	
	p1->SetLeftMargin(0.13);	

	cfinal->SaveAs("./results/SIM_peak_corr_vs_snr_final.png");
	delete cfinal;


	TCanvas *cfinal2 = new TCanvas("","",1100,850);
	snr_vs_wrms[0]->Draw("colz");
	snr_vs_wrms[0]->GetXaxis()->SetTitle("log_{10}(Wavefront RMS)");
	snr_vs_wrms[0]->GetYaxis()->SetTitle("Signal-to-Noise Ratio");
	snr_vs_wrms[0]->GetZaxis()->SetTitle("Number of Events");	
	snr_vs_wrms[0]->SetTitle("");
	
	snr_vs_wrms[0]->GetXaxis()->SetTitleFont(43);
	snr_vs_wrms[0]->GetXaxis()->SetTitleSize(40);
	snr_vs_wrms[0]->GetXaxis()->SetTitleOffset(1.0);

	snr_vs_wrms[0]->GetXaxis()->SetLabelFont(43);
	snr_vs_wrms[0]->GetXaxis()->SetLabelSize(35);
	snr_vs_wrms[0]->GetXaxis()->SetLabelOffset(0.005);
	
	snr_vs_wrms[0]->GetYaxis()->SetTitleFont(43);
	snr_vs_wrms[0]->GetYaxis()->SetTitleSize(40);
	snr_vs_wrms[0]->GetYaxis()->SetTitleOffset(1.3);

	snr_vs_wrms[0]->GetYaxis()->SetLabelFont(43);
	snr_vs_wrms[0]->GetYaxis()->SetLabelSize(35);
	snr_vs_wrms[0]->GetYaxis()->SetLabelOffset(0.015);
	
	snr_vs_wrms[0]->GetZaxis()->SetTitleFont(43);
	snr_vs_wrms[0]->GetZaxis()->SetTitleSize(40);
	snr_vs_wrms[0]->GetZaxis()->SetTitleOffset(1.0);

	snr_vs_wrms[0]->GetZaxis()->SetLabelFont(43);
	snr_vs_wrms[0]->GetZaxis()->SetLabelSize(35);
	snr_vs_wrms[0]->GetZaxis()->SetLabelOffset(0.005);

	gPad->SetLogz();
	TVirtualPad* p2 = cfinal2->cd();
	p2->SetRightMargin(0.15);
	p2->SetTopMargin(0.05);
	p2->SetLeftMargin(0.13);
	p2->SetBottomMargin(0.12);

	cfinal2->SaveAs("./results/SIM_snr_vs_wrms_final.png");
	delete cfinal2;


}