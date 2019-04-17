///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////		2018.07.13_NewCW_ApplyGeoFilter_InitialDistributions.cxx
//////		June 2018, Brian Clark 
//////		Making Initial distributions, filtering with geometric filter
//////		the frequencies identified by phase variance going forwards and backwards
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
	
	TH2D *PeakCorr_vs_SNR_all[2];
	PeakCorr_vs_SNR_all[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_all[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal[2];
	PeakCorr_vs_SNR_cutCal[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft[2];
	PeakCorr_vs_SNR_cutCal_cutSoft[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *peakcorrv_vs_snr_final[2];
	peakcorrv_vs_snr_final[0]=new TH2D("","V",20,0,20,50,0,0.5);
	peakcorrv_vs_snr_final[1]=new TH2D("","H",20,0,20,50,0,0.5);

	TH2D *snr_vs_wrms[2];
	snr_vs_wrms[0]=new TH2D("","V",90,-5,4,400,0,40);
	snr_vs_wrms[1]=new TH2D("","H",90,-5,4,400,0,40);
	
	int num_total=0;
	int num_cal=0;
	int num_soft=0;
	int num_short=0;
	int num_CW=0;
	int num_newbox=0;
	int num_surf=0;
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
		int waveformLength[16];
		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face", &rms_pol_thresh_face);
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);

		//next, load the reco tree
		TTree *inputTree_reco[35];
		double peakCorr[35][2];
		int peakTheta[35][2];
		int peakPhi[35][2];
		int recoBinSelect = 19; //300 m map
		int recoBinCalpulser = 6; //41 m map
		for(int i=0; i<35; i++){
			if(i==recoBinSelect||i==recoBinCalpulser){
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

		//next, load the CW tree
		ss.str("");
		ss << "OutputTree_CW";
		TTree *inputTree_CW = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_CW){
			cout<<"Can't open CW tree"<<endl;
			return -1;
		}
		vector<vector<double> > *badFreqs = 0;
		vector<vector<double> > *badSigmas = 0;
		inputTree_CW->SetBranchAddress("badFreqs",&badFreqs);
		inputTree_CW->SetBranchAddress("badSigmas",&badSigmas);

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		char summary_file_name[400];
		sprintf(summary_file_name,"/data/user/brianclark/A23Diffuse/NewCWAllEv/A2/NewCW_run_%d.root",runNum);
		TFile *NewCWFile = TFile::Open(summary_file_name);
		if(!NewCWFile) {
			std::cerr << "Can't open new CW file\n";
			return -1;
		}
		TTree* NewCWTree = (TTree*) NewCWFile->Get("NewCWTree");   
		if(!NewCWTree) {
			std::cerr << "Can't find NewCWTree\n";
			return -1;
		}
		vector<vector<double> > *badFreqs_fwd =0;
		vector<vector<double> > *badFreqs_back=0;
		vector<vector<double> > *badSigmas_fwd=0;
		vector<vector<double> > *badSigmas_back=0;

		NewCWTree->SetBranchAddress("badFreqs_fwd",&badFreqs_fwd);
		NewCWTree->SetBranchAddress("badSigmas_fwd",&badSigmas_fwd);
		NewCWTree->SetBranchAddress("badFeqs_back",&badFreqs_back);
		NewCWTree->SetBranchAddress("badSigmas_back",&badSigmas_back);

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

			bool isShort=false;
			bool isSurf=false;
			bool isCP5=false;
			bool isCP6=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			num_total++;

			for(int i=0;i<16;i++){ if(waveformLength[i]<550) isShort=true; }

			
			for (int i = 0; i < 35; i++){
				if (i == recoBinSelect || i == recoBinCalpulser){
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


			for(int pol=0; pol<2; pol++){
				if(bestTheta[pol] > 37) isSurf=true;
			}

			//figure out which reconstruction map (vpol or hpol) is best
			//for the 41m bin
			double bestCorr_pulser[] = {0., 0., 0.};
			int bestCorrRadiusBin_pulser[3];
			int bestPol_pulser = 2;
			int bestTheta_pulser[3];
			int bestPhi_pulser[3];

			for(int pol=0; pol<2; pol++){
				for(int i=0; i<35; i++){
					if (i == recoBinCalpulser){
						if (peakCorr[i][pol] > bestCorr_pulser[pol]){
							bestCorr_pulser[pol] = peakCorr[i][pol];
							bestCorrRadiusBin_pulser[pol] = i;
							bestTheta_pulser[pol] = peakTheta[i][pol];
							bestPhi_pulser[pol] = peakPhi[i][pol];
							//cout<<"best phi pulser "<<bestPhi_pulser[pol]<<endl;
							//cout<<"best peak corr "<<bestCorr_pulser[pol]<<endl;
						}
						if (peakCorr[i][pol] > bestCorr_pulser[2]){
							bestCorr_pulser[2] = peakCorr[i][pol];
							bestCorrRadiusBin_pulser[2] = i;
							bestTheta_pulser[2] = peakTheta[i][pol];
							bestPhi_pulser[2] = peakPhi[i][pol];
							bestPol_pulser = pol;
						}
					}//cal pulser (41m) bin check
				}//loop over reco bins
			}//loop over polarizations
			
			//draw a box around the cal pulser
			for (int pol = 0; pol < 2; pol++){
				if (bestPhi_pulser[pol] > -30 && bestPhi_pulser[pol] < -20 && bestTheta_pulser[pol] > -25 && bestTheta_pulser[pol] < -10){
					isCP5=true;
				}
				//if (bestPhi_pulser[pol] > 60 && bestPhi_pulser[pol] < 70 && bestTheta_pulser[pol] > 10 && bestTheta_pulser[pol] < 25){
				if (bestPhi_pulser[pol] > 60 && bestPhi_pulser[pol] < 70 && bestTheta_pulser[pol] > 0 && bestTheta_pulser[pol] < 15){
					isCP6=true;
				}
			}

			//deal w/ CW cut
			//inputTree_CW->GetEntry(event);
			NewCWTree->GetEntry(event);
			bool isCutonCW_fwd[2]; isCutonCW_fwd[0]=false; isCutonCW_fwd[1]=false;
			bool isCutonCW_back[2]; isCutonCW_back[0]=false; isCutonCW_back[1]=false;
			
			double threshCW = 1.0;
			vector<double> badFreqList_fwd;
			vector<double> badSigmaList_fwd;
			//first, loop over 
			for(int pol=0; pol<badFreqs_fwd->size(); pol++){
				badFreqList_fwd=badFreqs_fwd->at(pol);
				badSigmaList_fwd=badSigmas_fwd->at(pol);
				for(int iCW=0; iCW<badFreqList_fwd.size(); iCW++){
					if(
						badSigmaList_fwd[iCW] > threshCW 
						&& 
						abs(300. - badFreqList_fwd[iCW]) > 2.
						&&
						abs(500. - badFreqList_fwd[iCW]) > 2.
						&&
						abs(125. - badFreqList_fwd[iCW]) > 2.
					){
						isCutonCW_fwd[pol] = true;
					}
				}
			}
			vector<double> badFreqList_back;
			vector<double> badSigmaList_back;
			for(int pol=0; pol<badFreqs_back->size(); pol++){
				badFreqList_back=badFreqs_back->at(pol);
				badSigmaList_back=badSigmas_back->at(pol);
				for(int iCW=0; iCW<badFreqList_back.size(); iCW++){
					if(
						badSigmaList_back[iCW] > threshCW 
						&& 
						abs(300. - badFreqList_back[iCW]) > 2.
						&&
						abs(500. - badFreqList_back[iCW]) > 2.
						&&
						abs(125. - badFreqList_back[iCW]) > 2.
					){
						isCutonCW_back[pol] = true;
					}
				}
			}

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

			if(isCalPulser) num_cal++;
			if(isSoftTrigger) num_soft++;
			if(isShort) num_short++;
			if(failWavefrontRMS[0]) num_WFRMS[0]++;
			if(failWavefrontRMS[1]) num_WFRMS[1]++;
			if(isCP5 || isCP6 ) num_newbox++;
			if(isSurf) num_surf++;

			for(int pol=0; pol<2; pol++){
				PeakCorr_vs_SNR_all[pol]->Fill(SNRs[pol],bestCorr[pol]);

					snr_vs_wrms[pol]->Fill(TMath::Log10(bestFaceRMS[pol]),SNRs[pol]);

				// if(!isSoftTrigger 
				// 	&& !isShort 
				// 	&& !isCutonCW_fwd[pol] 
				// 	&& !isCutonCW_back[pol]
				// 	&& !isSurf)
				// {
				// } 

				
				//cut cal pulsers
				if(!isCalPulser){
					PeakCorr_vs_SNR_cutCal[pol]->Fill(SNRs[pol],bestCorr[pol]);
					
					//cut software triggers
					if(!isSoftTrigger){
						PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Fill(SNRs[pol],bestCorr[pol]);
						
						//cut short
						if(!isShort){
							PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Fill(SNRs[pol],bestCorr[pol]);
							
							//cut WRMS
							if(!failWavefrontRMS[pol]){
								PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Fill(SNRs[pol],bestCorr[pol]);
								
								//cut cal box
								if(!isCP5 && !isCP6){
									PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Fill(SNRs[pol],bestCorr[pol]);

									//cut events which have CW
									if(!isCutonCW_fwd[pol] && !isCutonCW_back[pol]){
										PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[pol]->Fill(SNRs[pol],bestCorr[pol]);
										
										//cut events which are surface
										if(!isSurf){
											PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[pol]->Fill(SNRs[pol],bestCorr[pol]);
											peakcorrv_vs_snr_final[pol]->Fill(SNRs[pol],bestCorr[pol]);
											//save these out to a tree for us to optimize over later
										}//not surface
									}//not cut on CW
								}//not cal pulser
							}//now wavefront RMS
						}//not short
					}//not software
				}//not cal pulser
			}//loop over polarizations
		}//loop over events
		inputFile->Close();
		NewCWFile->Close();
		delete inputFile;
	}

	gStyle->SetOptStat(11);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	//save out SNR vs WavefrontRMS plot
	char graph_title[2][300];
	char title[300];

	int cal=0;
	int soft=0;
	int Short=0;
	int wrms=0;
	int box = 0;
	int surf=0;
	int cw=0;

	//save out the Corr vs SNR plot for all 
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c2 = new TCanvas("","",2.1*850,850);
	c2->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
		PeakCorr_vs_SNR_all[pol]->Draw("colz");
		PeakCorr_vs_SNR_all[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_all[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_all[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c2->SaveAs(title);
	delete c2;
	delete PeakCorr_vs_SNR_all[0]; delete PeakCorr_vs_SNR_all[1];

	//turn on cal
	cal=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c3 = new TCanvas("","",2.1*850,850);
	c3->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c3->cd(pol+1);
		PeakCorr_vs_SNR_cutCal[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c3->SaveAs(title);
	delete c3;
	delete PeakCorr_vs_SNR_cutCal[0]; delete PeakCorr_vs_SNR_cutCal[1];

	//turn on cal, soft
	soft=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c4 = new TCanvas("","",2.1*850,850);
	c4->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c4->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c4->SaveAs(title);
	delete c4;
	delete PeakCorr_vs_SNR_cutCal_cutSoft[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft[1];

	//turn on cal, soft, short
	Short=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c5 = new TCanvas("","",2.1*850,850);
	c5->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c5->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c5->SaveAs(title);
	delete c5;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1];

	//turn on cal, soft, short, wmrs
	wrms=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c6 = new TCanvas("","",2.1*850,850);
	c6->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c6->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c6->SaveAs(title);
	delete c6;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1];

	//turn on cal, soft, short, wmrs, box
	box=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c7 = new TCanvas("","",2.1*850,850);
	c7->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c7->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c7->SaveAs(title);
	delete c7;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1];

	//turn on cal, soft, short, wmrs, box, surf
	cw=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, CW %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, CW %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c9 = new TCanvas("","",2.1*850,850);
	c9->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c9->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c9->SaveAs(title);
	delete c9;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW[1];

	//turn on cal, soft, short, wmrs, box, surf
	surf=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, CW %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, CW %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c13 = new TCanvas("","",2.1*850,850);
	c13->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c13->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "./results/%d.%d.%d_A%d_2013_%dEvents_Correlation_vs_SNR_cal%d_soft%d_short%d_wrms%d_newbox%d_CW%d_surf%d.png",year_now, month_now, day_now,station,num_total,cal,soft,Short,wrms,box,cw,surf);
	c13->SaveAs(title);
	delete c13;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutCW_cutSurf[1];

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

	cfinal->SaveAs("./results/peak_corr_vs_snr_final.png");
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

	cfinal2->SaveAs("./results/snr_vs_wrms_final.png");
	delete cfinal2;


}