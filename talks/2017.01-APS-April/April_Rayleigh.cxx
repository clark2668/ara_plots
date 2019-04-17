////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  analysis.cxx 
////
////    Dec 2016,  clark.2668@osu.edu
////	Pass one of the solar flare analysis
////
////    This code will implement the following cuts
////    Correlation peak map value
////    Second highest vpeak/ rms value
////    Area cut/ reconstruction quality cut
////    Saturation cut
////    Slanted cut
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <memory>
#include <boost/scoped_ptr.hpp>
#include <vector>
//#include <array>
#include <cmath>
#include <algorithm>
//#include <complex.h>
#include <fftw3.h>
#include <ctime>
#include "time.h" // for time convert

//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "RawAraStationEvent.h"
//#include "AraEventCorrelator.h"
#include "AraAntennaInfo.h"
#include "AraGeomTool.h"
#include "FFTtools.h"
#include "AraStationInfo.h"
#include "AraAntennaInfo.h"
#include "/home/clark.2668/workspace/TrunkSolarFlares/source_AraRoot_Trunk/include/RayTraceCorrelator.h"
//#include "RayTraceCorrelator.h"

//AraSim Includes
#ifndef ___CINT___ //shielded from CINT
#include "Settings.h"
//#include "Event.h"
//#include "Detector.h"
//#include "Report.h"
//#include "Vector.h"
#endif

//ROOT Includes
#include "TTreeIndex.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TText.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TString.h"
#include "TVector3.h"
#include "TLegend.h"

using namespace std;
const double DEG2RAD=0.017453292;   // radians/degree (convert to radians by multiplying degrees by rad/deg)
const double RAD2DEG=57.2957795;    // degree/rad (convert to degrees by multiplying radians by deg/rad)


RayTraceCorrelatorType::RayTraceCorrelatorType_t raytrace = RayTraceCorrelatorType::kRayTrace;
RayTraceCorrelatorType::RayTraceCorrelatorType_t planewave = RayTraceCorrelatorType::kPlaneWave;

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

//some helper functions
void fillAntennaPositions( int antenna_number, Double_t (&antenna_position)[3], Double_t &antenna_rho, Double_t &antenna_phi); 
Double_t calcDeltaT(Double_t ant1[3],Double_t rho1, Double_t phi1, Double_t ant2[3],Double_t rho2, Double_t phi2, Double_t phiWave, Double_t thetaWave);
void GetCorrMapPeakAngle_1deg( TH2D *theCorrMap_input, int &peaktheta, int &peakphi, double &peakvalue);
void GetCorrMapPeakAngle_1deg_SolarAnalysis( TH2D *theCorrMap_input, int &peaktheta, int &peakphi, double &peakvalue, int eventType) ;
void CheckCorrMapRange( int &theta, int &phi );
int GetMaxVector( vector <int> &array );
int GetMinVector( vector <int> &array );
TH2D *GetPeakBoundaryPlot ( TH2D *theCorrMap_input, int peaktheta, int peakphi, double peakfrac, int &max_theta, int &min_theta, int &max_phi, int &min_phi, int bin_x, int bin_y );
TH2D *GetPeakAreaPlot ( TH2D *theCorrMap, TH2D *theCorrBoundary, int peaktheta, int peakphi, double peakfrac, int max_theta, int min_theta, int max_phi, int min_phi, double &area, int bin_x, int bin_y );
void computeFFT(TGraph *inputWaveform, TGraph* (& SpectrumPhase)[3]);
TH2D *makeSpectrogram(TGraph *waveform, double ideal_dT_ns, int spectrogram_freq_bin_width);
TGraph *makeSpectrum_mVPerRootHz(TGraph *grWave);
TGraph *makeCoherentSum(UsefulIcrrStationEvent *realIcrrEvPtr, int antenna_1_number, int num_antennas_used, int source_phi, int source_theta, Double_t antenna_positions[8][3], Double_t antenna_pol[8], Double_t antenna_rho[8], Double_t antenna_phi[8]);
TGraph *normalizeWaveform(TGraph *inputGraph);


//AraGeomTool *geomTool;
//RawIcrrStationEvent *rawIcrrEvPtr=0;
//RawAtriStationEvent *rawAtriEvPtr=0;
//RawAraStationEvent *rawEvPtr=0;
//UsefulIcrrStationEvent *realIcrrEvPtr=0;
//UsefulAtriStationEvent *realAtriEvPtr=0;

/* Channel mappings for the testbed
Channel 0: H Pol
Channel 1: H Pol
Channel 2: V Pol
Channel 3: V Pol
Channel 4: V Pol
Channel 5: H Pol
Channel 6: V Pol
Channel 7: H Pol
Channel 8: V Pol
Channel 9: H Pol
Channel 10: V Pol
Channel 11: H Pol
Channel 12: H Pol
Channel 13: H Pol
Channel 14: Surface
Channel 15: Surface
*/

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
    
	std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input file 1> <input file 2> ....<input file n>" << std::endl;
		return -1;
	}

	bool MakeAverageSpectra=true;
	bool MakeRayleigh=true;
	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;
	
	//stealing this setup for the Ray Trace Correlator from Carl
	Settings *settings1 = new Settings(); //make a settings file
	string setupfile = "setup.txt";
	settings1 -> ReadFile(setupfile);
	settings1->NOFZ=1;
	int select_station = 0; //make a testbed map
	double select_radius = 3000.; //make a map with a 30m punitive source direction
	double select_ang_bin = 1.; //1 degree angular bins
	int test_mode =4;
	RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(select_station, select_radius, settings1, select_ang_bin, test_mode);
	Detector *detector = 0;
	AraGeomTool *geomTool = new AraGeomTool();
	const int num_antennas_used=8;
	
	int spectralBins[16]={}; //to hold the spectral bins of interest
	TH1D *rayleighPlotsFeb15V[16]={};
 	TH1D *rayleighPlotsFeb15H[16]={};
	TH1D *rayleighPlotsFeb11V[16]={};
	TH1D *rayleighPlotsFeb11H[16]={};
	for(int i=0; i<16; i++){
		spectralBins[i]=((double) (i+6))*15; //50 MHz bins
		if(i==1 || i==2 || i==3 || i==4){
			rayleighPlotsFeb15V[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb15H[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb11V[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb11H[i]=new TH1D("","",100,0,20);
		}
		else{
			rayleighPlotsFeb15V[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb15H[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb11V[i]=new TH1D("","",100,0,20);
			rayleighPlotsFeb11H[i]=new TH1D("","",100,0,20);
		}	
	}
	double spectralValues[16]={}; //to hold the spectral values of interest
	
	
	//first, I need to deal with the Feb 15 data, which is already event selected on
		TChain chain("eventTree"); //this for the events for the exterior loop
		TChain CutsChain("TreeCutValues"); //get the correlation maps
		for(int file=1	; file<argc; file++){ //start at 2 when there is a pedestal file
			TString fileKey(argv[file]); //a file key
			chain.Add(fileKey); //add files to the chain
			CutsChain.Add(fileKey);
		}
     
		double CorrMapPeak_value = 0.;
		int PeakTheta_value=0;
		int PeakPhi_value=0;
		double VpeakRMS_value = 0.;
		double AreaTotal_value =0.;
		double AreaContour_value =0.;
		TH2D *VMap_30m = 0;
     
		double CorrMapPeakCut_value = 0.13;
		double VpeakRMSCut_value = 4.;
		double AreaMinCut_value =1.;
		double AreaMaxCut_value = 50.;
		double SlantedCut_yint_value= 8.8;
		double SlantedCut_slope_value = -14.;
		double PeakFracCut_value = 0.85;
		double RatioCut_value = 1.5;

		RawIcrrStationEvent *rawIcrrEvPtr=0; //make the raw event pointer
		chain.SetBranchAddress("event",&rawIcrrEvPtr); //set the branch address
		CutsChain.SetBranchAddress("CorrMapPeak",&CorrMapPeak_value); //set the branch address
		CutsChain.SetBranchAddress("PeakTheta",&PeakTheta_value); //set the branch address
		CutsChain.SetBranchAddress("PeakPhi",&PeakPhi_value); //set the branch address
		CutsChain.SetBranchAddress("VpeakRMS",&VpeakRMS_value); //set the branch address
		CutsChain.SetBranchAddress("AreaTotal",&AreaTotal_value); //set the branch address
		CutsChain.SetBranchAddress("AreaContour",&AreaContour_value); //set the branch address
		CutsChain.SetBranchAddress("VMap_30m",&VMap_30m); //set the branch address
     	
		vector<double> peakPhi;
		vector<double> peakTheta;
		vector<double> reportedpeakPhi;
		vector<double> reportedpeakTheta;
		vector<double> reportedpeakPhi_tshift;
		vector<double> reportedpeakTheta_tshift;
		vector<double> eventTime;
		vector<double> eventTime_tshift;
		
		Int_t numEntries = chain.GetEntries(); //get the number of entries
		
		//get ready to average over everything
		chain.GetEvent(0); //just pick something for now...
		CutsChain.GetEvent(0); //just pick something for now...
		
		UsefulIcrrStationEvent *firstrealIcrrEvPtr = new UsefulIcrrStationEvent(rawIcrrEvPtr, AraCalType::kLatestCalib); //make a real event pointer for the Atri station event with a run dependent calibration
		//make coherent sums
		TGraph *firstcoherentSumRayTraceV = theCorrelator->makeCoherentSum(firstrealIcrrEvPtr,2,PeakPhi_value,PeakTheta_value);
		TGraph *firstcoherentSumRayTraceH = theCorrelator->makeCoherentSum(firstrealIcrrEvPtr,0,PeakPhi_value,PeakTheta_value);
		TGraph *RayleighPrepVFeb15 = makeSpectrum_mVPerRootHz(firstcoherentSumRayTraceV);	
		TGraph *RayleighPrepHFeb15 = makeSpectrum_mVPerRootHz(firstcoherentSumRayTraceH);		
		//just get a baseline starting spectra to add to
		Double_t *AverageSpectraVFeb15 = RayleighPrepVFeb15->GetY(); //prepare to average the spectrum
		Double_t *AverageSpectraHFeb15 = RayleighPrepHFeb15->GetY(); //prepare to average the spectrum
		//know what the frequency values we needt to loop over are
		Double_t *FreqVals = RayleighPrepVFeb15->GetX();  //get the spectral values
		//how many frequency bins to loop over; debugging only
		Int_t numSampsPrep = RayleighPrepVFeb15 -> GetN();
		
		TH2D *TwoCuts = new TH2D("","",100,0,1,100,0,20);
		
		int num_found=0;
		int num_found_15=0;
		
		for(int event=0; event<numEntries; event++){
			chain.GetEvent(event);
			CutsChain.GetEvent(event);
			int unixtime_1;
			unixtime_1 = rawIcrrEvPtr->head.unixTime; // get  event unixTime
			time_t test_time_1 = unixtime_1;
			tm *test_time_tm_1 = gmtime( &test_time_1 );
			int hour = test_time_tm_1->tm_hour;
			if(hour!=2) {continue;} //impose a time constraint
			UsefulIcrrStationEvent *realIcrrEvPtr = new UsefulIcrrStationEvent(rawIcrrEvPtr, AraCalType::kLatestCalib); //make a real event pointer for the Atri station event with a run dependent calibration
			double area_ratio = AreaTotal_value/AreaContour_value;
			double slantedcut = (SlantedCut_slope_value * CorrMapPeak_value) + SlantedCut_yint_value;
			
			TwoCuts->Fill(CorrMapPeak_value,VpeakRMS_value);
			
			//cout<<"on event "<<event<<endl;
			
			if( //2==2 //make a whole bunch of cuts
				(CorrMapPeak_value >= CorrMapPeakCut_value)
				&&
				(VpeakRMS_value > VpeakRMSCut_value)
				&&
				(AreaContour_value >= AreaMinCut_value)
				&&
				(AreaContour_value <= AreaMaxCut_value)
				&&
				(VpeakRMS_value >= slantedcut)
				&&
				(area_ratio<RatioCut_value)
				&&
				(realIcrrEvPtr->isCalPulserEvent() != true) 
				&& 
				(rawIcrrEvPtr -> trig.trigType !=68)
			){ //if the event passes all cuts, do something with it

				//cout<<"Phi, Theta  "<<PeakPhi_value<<" , "<<PeakTheta_value<<endl;
				if(PeakPhi_value==180) continue;
				num_found++;
				num_found_15++;
				TGraph *coherentSumRayTraceV = theCorrelator->makeCoherentSum(realIcrrEvPtr,2,PeakPhi_value,PeakTheta_value);
				TGraph *coherentSumRayTraceH = theCorrelator->makeCoherentSum(realIcrrEvPtr,0,PeakPhi_value,PeakTheta_value);
				TGraph *forRayleighV = makeSpectrum_mVPerRootHz(coherentSumRayTraceV);
				TGraph *forRayleighH = makeSpectrum_mVPerRootHz(coherentSumRayTraceH);
				Int_t numSamps = forRayleighV -> GetN(); //get number of frequency bins to loop over
				Double_t *yValsV = forRayleighV -> GetY();
				Double_t *yValsH = forRayleighH -> GetY();
				Double_t *xVals = forRayleighV -> GetX();
				if((numSamps!=numSampsPrep)){cout<<"Something is wrong between the base and new spectrum"<<endl;} //just check this, though it should never happen the way I've designed the analysis
				for(int i=0; i<numSamps; i++){ //for averaged spectra
					AverageSpectraVFeb15[i]+=yValsV[i];
					AverageSpectraHFeb15[i]+=yValsH[i];
				}
				for(int bin=0; bin<16; bin++){ //iterate over all the bins //for Rayleigh
					double spectral_value_this_run_V = yValsV[spectralBins[bin]]; //get the spectra value for this V event
					double spectral_value_this_run_H = yValsH[spectralBins[bin]]; //get the spectra value for this V event
					spectralValues[bin]=xVals[spectralBins[bin]];
					rayleighPlotsFeb15V[bin]->Fill(spectral_value_this_run_V); //add it to the Rayleigh plot
					rayleighPlotsFeb15H[bin]->Fill(spectral_value_this_run_H); //add it to the Rayleigh plot
				}
				char thisevent[150];
				/*
				TCanvas *thisOne = new TCanvas("","",1100,850);
				thisOne->Divide(1,2);
				thisOne->cd(1);
				forRayleighV->Draw();
				thisOne->cd(2);
				forRayleighH->Draw();
				sprintf(thisevent,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/April_SingleSpectra_Event%d.pdf",event);
				delete thisOne;
				*/
				delete forRayleighV;
				delete forRayleighH;
				delete coherentSumRayTraceV;
				delete coherentSumRayTraceH;
			}
			delete realIcrrEvPtr;
		}
		
		TCanvas *twoCuts = new TCanvas("two","two",1100,850);
		TwoCuts->Draw("colz");
		TwoCuts->GetXaxis()->SetTitle("Max Peak Correlation");
		TwoCuts->GetYaxis()->SetTitle("Signal Strength (arb units)");
		TwoCuts->SetTitle("");
		twoCuts->SaveAs("/home/clark.2668/workspace/TrunkSolarFlares/results/today/April_RotatedCut.pdf");
		delete twoCuts;
		delete TwoCuts;
		
	//okay, now to deal with February 11, the background window
	
	/*
	//need to find when the feb 11 events start
	int event_find =0;
	for(Int_t find=0; find<numEntries; find++){
		//if( (realIcrrEvPtr->isCalPulserEvent() != true) && (rawIcrrEvPtr->trig.trigType==68)) { //Plot only the RF events                                                                                                                                                   
		int unixtime_find;
		unixtime_find = rawIcrrEvPtr->head.unixTime; // get  event unixTime
		time_t test_time_find = unixtime_find;
		tm *test_time_tm_find = gmtime( &test_time_find );
		int day_find = test_time_tm_find->tm_mday;
		if(day_find ==11){event_find = find; break;}
	}//find the first one
	sumsChain.GetEvent(event_find);
	Int_t testnumpoints2 = thecoherentSum->GetN();
	TGraph *RayleighPrep2 = makeSpectrum_mVPerRootHz(thecoherentSum);
	Double_t *AverageSpectra_Feb11 = RayleighPrep2 ->GetY();	
	int num_found_11=0;
	*/
	
	cout<<"Found "<<num_found_15<<" events from Feb 15"<<endl;

	char save_title[150];
	char graph_title[150];	
	if(MakeAverageSpectra==true){
		for(int i=0; i<numSampsPrep; i++){
			AverageSpectraVFeb15[i]/=(double) num_found_15; //make this the averaged spectra specifically
			AverageSpectraHFeb15[i]/=(double) num_found_15;
		}
		TGraph *AveragedSpectraVFeb15 = new TGraph(numSampsPrep,FreqVals,AverageSpectraVFeb15);
		TGraph *AveragedSpectraHFeb15 = new TGraph(numSampsPrep,FreqVals,AverageSpectraHFeb15);
		TCanvas *AveragedSpectraCanvas = new TCanvas("","",2*1100,2*850);
			sprintf(graph_title,"Averaged Spectra, 2-3AM",num_found);
			AveragedSpectraVFeb15->Draw();
			AveragedSpectraVFeb15->SetLineWidth(3);
			AveragedSpectraHFeb15->Draw("same");
			AveragedSpectraHFeb15->SetLineColor(kBlue);
			AveragedSpectraHFeb15->SetLineWidth(3);
			AveragedSpectraVFeb15->SetTitle("");
			AveragedSpectraVFeb15->GetXaxis()->SetRangeUser(0.,1200.);
			AveragedSpectraVFeb15->GetYaxis()->SetTitle("Spectral Amplitude (V/MHz)");
			AveragedSpectraVFeb15->GetYaxis()->SetTitleOffset(1.5);
			AveragedSpectraVFeb15->GetXaxis()->SetTitle("Frequency (MHz)");
		sprintf(save_title,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_April_AveragedSpectra.pdf",year_now,month_now,day_now);
		AveragedSpectraCanvas->SaveAs(save_title);
		delete AveragedSpectraCanvas;
		delete AveragedSpectraVFeb15;
		delete AveragedSpectraHFeb15;
	}
		
	if(MakeRayleigh==true){
		TCanvas *RayleighCanvas = new TCanvas("rayleigh","rayleigh",1100,850);
		RayleighCanvas->Divide(4,4);
		char plots_title[150];
		for(int i=0; i<16; i++){
			RayleighCanvas->cd(i+1);
			rayleighPlotsFeb15V[i]->Draw();
			sprintf(plots_title,"distribution of spectral coefficients for coherent sums, %f MHz",spectralValues[i]);
			rayleighPlotsFeb15V[i]->SetTitle(plots_title);
			rayleighPlotsFeb15V[i]->GetXaxis()->SetTitle("Spectral Amplitude");
			rayleighPlotsFeb15V[i]->GetYaxis()->SetTitle("Counts");
			rayleighPlotsFeb15V[i]->SetLineColor(kBlack);
			rayleighPlotsFeb15H[i]->Draw("same");
			rayleighPlotsFeb15H[i]->SetLineColor(kRed);
		}
		sprintf(save_title,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_April_RayleighPlots.pdf(",year_now,month_now,day_now);
		RayleighCanvas->Print(save_title,"pdf");
		delete RayleighCanvas; //delete this
		gStyle->SetOptStat(kFALSE);//clear everything out
		
		
		vector <double> sigmasV;
		vector <double> sigmasH;
		TF1 *RayleighFitFeb15V=0;
 		TF1 *RayleighFitFeb15H=0;
		TLegend *legend2=0;
		TCanvas *NormalizedRayleighCanvas = new TCanvas("rayleighfit","rayleighfit",4*1100,4*850);
		NormalizedRayleighCanvas->Divide(4,4);
		for(int i=1; i<16; i++){ //for the moment, don't do all 16 frequencies I pulled out
			legend2 = new TLegend(0.4,0.8,0.9,0.4); 
			NormalizedRayleighCanvas->cd(i);
			char infofeb15V[150];
			char infofeb15H[150];
			rayleighPlotsFeb15V[i]->Sumw2();
			rayleighPlotsFeb15H[i]->Sumw2();
			rayleighPlotsFeb15V[i]->Scale(1/rayleighPlotsFeb15V[i]->Integral("width"));
			rayleighPlotsFeb15H[i]->Scale(1/rayleighPlotsFeb15H[i]->Integral("width"));
			rayleighPlotsFeb15V[i]->Draw();
			rayleighPlotsFeb15H[i]->Draw("same");
			rayleighPlotsFeb15H[i]->SetLineColor(kRed);
			sprintf(plots_title,"normalized distribution of spectral coefficients, %f MHz",spectralValues[i]);
			rayleighPlotsFeb15V[i]->SetTitle(plots_title);
			rayleighPlotsFeb15V[i]->GetXaxis()->SetTitle("Magnitude Spectral Coefficient");
			rayleighPlotsFeb15V[i]->GetYaxis()->SetTitle("Counts");
			rayleighPlotsFeb15V[i]->SetLineColor(kBlack);
			RayleighFitFeb15V = new TF1("RayleighFitFeb15V","(x/([0]*[0]))*exp(((-0.5)*x*x)/([0]*[0]))",-10,200); //create a function to fit with
				RayleighFitFeb15V->SetParameter(0,rayleighPlotsFeb15V[i]->GetRMS()); //this should be a good starting place for sigma
				RayleighFitFeb15V->SetLineColor(kBlack);
				RayleighFitFeb15V->SetLineWidth(2);
				rayleighPlotsFeb15V[i]->Fit("RayleighFitFeb15V");
				double sigma15V = RayleighFitFeb15V->GetParameter(0); //get the sigma
				double chisquare15V = RayleighFitFeb15V->GetChisquare();
				int ndf15V = RayleighFitFeb15V->GetNDF();
				sprintf(infofeb15V,"2/15 V, chi^2/NDF = %d/%d",(int) chisquare15V,ndf15V);
				sigmasV.push_back(sigma15V);
			 	legend2->AddEntry(RayleighFitFeb15V,infofeb15V,"lp");
			RayleighFitFeb15H = new TF1("RayleighFitFeb15H","(x/([1]*[1]))*exp(((-0.5)*x*x)/([1]*[1]))",-10,200); //create a function to fit with
				RayleighFitFeb15H->SetParameter(1,rayleighPlotsFeb15H[i]->GetRMS()); //this should be a good starting place for sigma
				RayleighFitFeb15H->SetLineColor(kRed);
				RayleighFitFeb15H->SetLineWidth(2);
				rayleighPlotsFeb15H[i]->Fit("RayleighFitFeb15H");
				double sigma15H = RayleighFitFeb15H->GetParameter(1); //get the sigma
				double chisquare15H = RayleighFitFeb15H->GetChisquare();
				int ndf15H = RayleighFitFeb15H->GetNDF();
				sprintf(infofeb15H,"2/15 H, chi^2/NDF = %d/%d",(int) chisquare15H,ndf15H);
			 	legend2->AddEntry(RayleighFitFeb15H,infofeb15H,"lp");
			 	sigmasH.push_back(sigma15H);
			 legend2->SetBorderSize(0);  //no border for legend
			 legend2->SetFillColor(0);  //fill color is white                        
			 legend2->SetTextSize(0.04);
			 legend2->DrawClone("same");  
			 delete RayleighFitFeb15V;
			 delete RayleighFitFeb15H;
			 delete legend2;
		}
		
		for(int i=0; i<16; i++){
			cout<<"V sigma "<<sigmasV[i]<<endl;
			cout<<"H sigma "<<sigmasH[i]<<endl;
		}
		sprintf(save_title,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_April_RayleighPlots.pdf)",year_now,month_now,day_now);
		NormalizedRayleighCanvas->Print(save_title,"pdf");
		delete NormalizedRayleighCanvas;
	}
	
	for(int i=0; i<16;i++){delete rayleighPlotsFeb15H[i]; delete rayleighPlotsFeb15V[i];} //delete these	     
}//close the main program


//a function to rescale a waveform so that the previously maximum value is one
TGraph *normalizeWaveform(TGraph *inputGraph){

		Int_t numSamps = inputGraph -> GetN(); //get the number of points
		Double_t *xVals = inputGraph -> GetX(); //get the x values of this specific waveform
		Double_t *yVals = inputGraph -> GetY(); //get the y values of this specific waveform
		double max_value = TMath::MaxElement(numSamps,yVals);
		double min_value = TMath::MinElement(numSamps,yVals);		
		double normalize_value=-1050; //something wierd so that I can identify a problem if it exists

		if(abs(min_value)>abs(max_value)){ normalize_value = abs(min_value);} //either the min value is abs_value larger, and in which case we want to noramlize by it
		else{normalize_value = max_value;} //otherwise the max is larger, or they are equal, and in which case, why not just use the max_value

		Double_t yVals_norm[numSamps];
		for (int iter=0; iter<numSamps; iter++){
			yVals_norm[iter]=yVals[iter]/normalize_value; //divide by the maximum value
		}

		TGraph *outputGraph = new TGraph(numSamps,xVals,yVals_norm); //make the new graph with the normalized y points
		return outputGraph; //return the final product
} //end function to make the max vale of a waveform 1

//a function to make a coherent sum from a given event
//Needs to take the
//antenna_1_number is the antenna against which all of the others are being shifted and summed to
//num_antennas_used is the number of antennas used in this coherent sum
//source_theta
//source_phi
//the array of antenna_positions, be careful here, the 8 in antenna_positions[8][3] depends on how many antennas you are using, and is not fixed to 8 in the main body of the code
//the array of antenna_pol
TGraph *makeCoherentSum(UsefulIcrrStationEvent *realIcrrEvPtr, int antenna_1_number, int num_antennas_used, int source_phi, int source_theta, Double_t antenna_positions[8][3], Double_t antenna_pol[8], Double_t antenna_rho[8], Double_t antenna_phi[8]){

	Double_t antenna_locations_first_station[3]; // a holder variable for the location of a station in this iteration; pick one antenna to shift everything else with respect to
	Double_t antenna_locations_second_station[3]; //a holder variable for the location of a second staion

	double interpolation_step=0.1; //this is not needed unless we are interpolating
	//Now, lots of things have to happen in a loop over the REST of the antennas
	/*
	vector <TGraph*> Waveforms;
	vector <TGraph*> Waveforms_Interpolated;
	vector <TGraph*> Waveforms_Shifted;
	vector <TGraph*> Waveforms_Padded;
	vector <Double_t> time_delays;
	*/
	vector <TGraph*> Waveforms_Cropped;
						
	//the flow here is: get the waveform, interpolated it, perform its time shift, pad it out, and then cut it off (get, interpolate, shift, pad, crop)
	//added the interpolated step so that we treat all of the antennas on even footing, otherwise the code depends on the time sampling of the waveforms to which all others are being shifted

	/* //old code, allows you to save the interpolation, cropping, etc as you go
	for(int antenna_iterator=0; antenna_iterator<num_antennas_used; antenna_iterator++){ //begin the loop over all of the other antennas							
		if(antenna_iterator==antenna_1_number){ //because we're using pushback, this antenna should get shoved in the right place in the line
			for(int coordinate_index=0; coordinate_index<3; coordinate_index++){ //loop over all coordinates (x,y,z) for the first station
				antenna_locations_first_station[coordinate_index]=antenna_positions[antenna_1_number][coordinate_index]; //fill in the x,y,z, values for the first station
			}
			
			double time_delay = 0;
			time_delays.push_back(0.); //because this antenna should experience no time delay
			Waveforms.push_back(realIcrrEvPtr->getGraphFromRFChan(antenna_1_number)); //the waveform from antenna_1_number gets it's place
			Waveforms_Interpolated.push_back(FFTtools::getInterpolatedGraph(Waveforms[antenna_1_number],interpolation_step));
			Waveforms_Shifted.push_back(FFTtools::translateGraph(Waveforms_Interpolated[antenna_1_number],0)); //no time delay; doing this explicitly to avoid pointing this at something else that might get deleted
			Waveforms_Padded.push_back(FFTtools::padWaveToLength(Waveforms_Shifted[antenna_1_number], Waveforms_Shifted[antenna_1_number]->GetN()+6000)); //zero pad by a lot
			Waveforms_Cropped.push_back(FFTtools::cropWave(Waveforms_Padded[antenna_1_number],-300.,300.)); //crop the waveform down from -250 to 250
		} //finish if statement over the reference antenna
				
		for(int coordinate_index=0; coordinate_index<3; coordinate_index++){ //loop over all coordinates (x,y,z) for the second antenna
			antenna_locations_second_station[coordinate_index]=antenna_positions[antenna_iterator][coordinate_index]; //fill in the x,y,z values for the second antenna
		}
		time_delays.push_back(calcDeltaT(antenna_locations_first_station, antenna_rho[antenna_1_number], antenna_phi[antenna_1_number], antenna_locations_second_station, antenna_rho[antenna_iterator],antenna_phi[antenna_iterator], source_phi, source_theta));//get the appropriate time delays
		Waveforms.push_back(realIcrrEvPtr->getGraphFromRFChan(antenna_iterator)); //extract the waveform from the second antenna
		Waveforms_Interpolated.push_back(FFTtools::getInterpolatedGraph(Waveforms[antenna_iterator],interpolation_step));
		Waveforms_Shifted.push_back(FFTtools::translateGraph(Waveforms_Interpolated[antenna_iterator],time_delays[antenna_iterator])); //time shift them
		Waveforms_Padded.push_back(FFTtools::padWaveToLength(Waveforms_Shifted[antenna_iterator], Waveforms_Shifted[antenna_iterator]->GetN()+6000)); //pad them
		Waveforms_Cropped.push_back(FFTtools::cropWave(Waveforms_Padded[antenna_iterator],-300.,300.)); //crop them
					
	}//end loop over all of the other antennas
	*/
	
	for(int antenna_iterator=0; antenna_iterator<num_antennas_used; antenna_iterator++){ //begin the loop over all of the other antennas							
		if(antenna_iterator==antenna_1_number){ //because we're using pushback, this antenna should get shoved in the right place in the line
			for(int coordinate_index=0; coordinate_index<3; coordinate_index++){ //loop over all coordinates (x,y,z) for the first station
				antenna_locations_first_station[coordinate_index]=antenna_positions[antenna_1_number][coordinate_index]; //fill in the x,y,z, values for the first station
			}	
			double time_delay = 0.;
			TGraph *Waveform = realIcrrEvPtr->getGraphFromRFChan(antenna_1_number);
			TGraph *Waveform_Interpolated = FFTtools::getInterpolatedGraph(Waveform,interpolation_step);
			delete Waveform;
			TGraph *Waveform_Shifted = FFTtools::translateGraph(Waveform_Interpolated,time_delay);
			delete Waveform_Interpolated;
			TGraph *Waveform_Padded = FFTtools::padWaveToLength(Waveform_Shifted, Waveform_Shifted->GetN()+6000);
			delete Waveform_Shifted;
			Waveforms_Cropped.push_back(FFTtools::cropWave(Waveform_Padded,-300.,300.));
			delete Waveform_Padded;
		} //finish if statement over the reference antenna
		else{		
			for(int coordinate_index=0; coordinate_index<3; coordinate_index++){ //loop over all coordinates (x,y,z) for the second antenna
				antenna_locations_second_station[coordinate_index]=antenna_positions[antenna_iterator][coordinate_index]; //fill in the x,y,z values for the second antenna
			}
			double time_delay_2 = calcDeltaT(antenna_locations_first_station, antenna_rho[antenna_1_number], antenna_phi[antenna_1_number], antenna_locations_second_station, antenna_rho[antenna_iterator],antenna_phi[antenna_iterator], source_phi, source_theta);//get the appropriate time delays
			TGraph *Waveform_2=(realIcrrEvPtr->getGraphFromRFChan(antenna_iterator)); //extract the waveform from the second antenna
			TGraph *Waveform_Interpolated_2=(FFTtools::getInterpolatedGraph(Waveform_2,interpolation_step));
			delete Waveform_2;
			TGraph *Waveform_Shifted_2=(FFTtools::translateGraph(Waveform_Interpolated_2,time_delay_2)); //time shift them
			delete Waveform_Interpolated_2;
			TGraph *Waveform_Padded_2=(FFTtools::padWaveToLength(Waveform_Shifted_2, Waveform_Shifted_2->GetN()+6000)); //pad them
			delete Waveform_Shifted_2;
			Waveforms_Cropped.push_back(FFTtools::cropWave(Waveform_Padded_2,-300.,300.)); //crop them
			delete Waveform_Padded_2;
		}
					
	}//end loop over all of the other antennas
	
	//to be clear, the whole point of the pad and crop exercise was to generate waveforms that all had the same mutual start time of -250 and end time of 250
	//TGraph *summed_waveforms = Waveforms_Cropped[antenna_1_number]->Clone(); //clone this graph, we want to add everything else to it
						
	Double_t *summed_x_values = Waveforms_Cropped[antenna_1_number]->GetX(); //get the x values of the summed object
	Int_t num_summed_points = Waveforms_Cropped[antenna_1_number]->GetN(); //get the number of data points to be in the final summed object
	vector<Double_t> summed_y_values;
	vector<Double_t> summed_y_values_vpol;
	for(Int_t x_iter=0; x_iter<num_summed_points; x_iter++){
		double new_y_value=0.;
		double new_y_value_vpol=0.;
		int  num_contributing_waveforms=0.; //this is to count the number of waveforms with non-zero contribution
		//at this x-value, we need to go into each waveform, pull down it's y value at this x-value, and add it to the total y-value
		
		for(int waveform_iterator = 0; waveform_iterator<num_antennas_used; waveform_iterator++){ //loop over all the waveforms, which I really mean antennas
			if(antenna_pol[waveform_iterator]==0){ //it's a v-pol //zero is Vpol, 1 is Hpol...
				new_y_value=Waveforms_Cropped[waveform_iterator]->Eval(summed_x_values[x_iter]);
				if(abs(new_y_value)>0.0){ //just want to know if there is something there, positive or negative are both fine, predicated on the assumption that we never measure zero
					 new_y_value_vpol+=new_y_value;
					 num_contributing_waveforms++;
				}
			}
			else{continue;} //it's an hpol  don't do anything
		}
		if(num_contributing_waveforms==0) num_contributing_waveforms=1; //if this number is not bigger than zero, then none of the waveforms contributed anything, and we need to set this equal to zero to avoid an undefined quantity
		summed_y_values_vpol.push_back(new_y_value_vpol/sqrt((double) num_contributing_waveforms)); //divide by the sqrt of the number of waveforms
	} //end loop over the waveforms to make the coherent sum

	//need to clear the variables before starting the next iteration
	/*
	for(int i=0; i<num_antennas_used; i++){
		delete Waveforms[i];
		delete Waveforms_Interpolated[i];
		delete Waveforms_Shifted[i];
		delete Waveforms_Padded[i];
		delete Waveforms_Cropped[i];
	}
	*/

	//cout<<"                        Just before actually making the new waveform"<<endl;
	TGraph *summed_waveforms_vpol= new TGraph(num_summed_points, summed_x_values, &summed_y_values_vpol[0]); //actually make the coherent sum
	//cout<<"                        The new waveform has been made"<<endl;
	for(int i =0; i<num_antennas_used; i++){
		delete Waveforms_Cropped[i];
	}
	return summed_waveforms_vpol; //return this final product


}//end coherent sum function


//this snippet of code is used for making Rayleigh plots, and I'm basically taking it straight from Eugene
TGraph *makeSpectrum_mVPerRootHz(TGraph *grWave){

	double *oldY = grWave->GetY();
	double *oldX = grWave->GetX();
	double deltaT=(oldX[1]-oldX[0])*1.e-9;
	int length=grWave->GetN();
	FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

	int newLength=(length/2)+1;

	double *newY = new double [newLength];
	double *newX = new double [newLength];

	//    double fMax = 1/(2*deltaT);  // In Hz
	double deltaF=1/(deltaT*length); //Hz
	deltaF*=1e-6; //MHz
	//    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t"
	//<< length << std::endl;

	double tempF=0;
	for(int i=0;i<newLength;i++) {
		newY[i] = FFTtools::getAbs(theFFT[i])*1.e-3;

		newX[i]=tempF;
		// newY[i]=power;
		tempF+=deltaF;
	}

	TGraph *grPower = new TGraph(newLength,newX,newY);
	delete [] theFFT;
	delete [] newY;
	delete [] newX;
	return grPower;

}


TH2D *makeSpectrogram(TGraph *waveform, double ideal_dT_ns, int spectrogram_freq_bin_width){
	//borrowing the following snipet of code from Eugen "AnalyzeAraData_tmp6_copy.cxx"
     	//ideal_dT_ns = 0.1; // in unit ns //I think this should be the same as the interpolation size
	//int spectrogram_freq_bin_width = 25; // 25 MHz bin //work with a 25 MHz bin //replaces "SP_time_Freq_width"
	//*waveform is a pointer to the waveform we want to make the spectrogram of
	
	TGraph * interpolated_waveform = FFTtools::getInterpolatedGraph(waveform, ideal_dT_ns); //get the interpolated version of the summed_waveforms (which should already have this spacing, but just to be sure...)
	int num_spectrogram_freq_bins= (int)(1000. / (double)spectrogram_freq_bin_width); // number of freq bins in a GHz //replaces "SP_time_Freq_Bin"
	//cout<<"The number of spectrogram frequency bins is"<<num_spectrogram_freq_bins<<endl;
	int num_spectrogram_time_bins = (int)( 1. / ( ideal_dT_ns * 1.e-9 * (double)spectrogram_freq_bin_width * 1.e6 ) ); //get the number of time bins  //replaces "SP_time_Bin"
	//cout<<"Number of spectrogram time bins "<<num_spectrogram_time_bins<<endl;
	int num_waveform_points = interpolated_waveform -> GetN(); //get the number of points in the summed waveform
	int num_waveform_pieces = (int) (num_waveform_points / num_spectrogram_time_bins); //get the number of pieces into which the waveform needs to be chopped, also the number of time slices in the spectrogram //this replaces "N_sub") //ie, the number of time bins in the spectrogram
	//cout<<"Number of summed waveform pieces "<<num_summed_waveform_pieces<<endl;
	char title_test[150];
	sprintf(title_test, "histSP_Time");
	TH2D *Spectrogram = new TH2D (title_test, "",num_waveform_pieces,0,num_waveform_pieces,num_spectrogram_freq_bins,0,1000); //the spectrogram//this replaces "histSP_Time_EachWF"
	
		
	for( int t_bin=0; t_bin < num_waveform_pieces; t_bin++){ //loop over all of the pieces of the summed waveform
		//cout<<"I'm at time bin "<<t_bin<<endl;	
		TGraph *waveform_piece = new TGraph(); //the sub-piece of the summed waveform //replaces "chan_sub"
		double x_tmp,  y_tmp; //temporary holders for the x and y values of this sub-piece of the waveform
			for( int bin=0; bin<num_spectrogram_time_bins; bin++){ //loop over the number of spectrogram time bins
				interpolated_waveform -> GetPoint(t_bin*num_spectrogram_time_bins + bin, x_tmp, y_tmp); //get the data out of the waveform
				waveform_piece -> SetPoint(bin, x_tmp, y_tmp); //put it inside the "piece" of the waveform we are interested in in this time slice
			}
														
			TGraph *waveform_piece_spectrum = FFTtools :: makePowerSpectrumMilliVoltsNanoSeconds(waveform_piece); //make a spectrum from this piece //the spectgrum of the sub-piece of the summed waveform //replaces "chan_sub_sp"
							
			int getn_sp = waveform_piece_spectrum -> GetN(); //get the number of points in this piece
			double *gety_sp = waveform_piece_spectrum -> GetY(); //get the y values of this piece of the waveform
			double *getx_sp = waveform_piece_spectrum -> GetX(); //get the x values of this piece of the waveform
							
			for (int bin=0; bin<getn_sp; bin++){
				Spectrogram -> Fill(t_bin, getx_sp[bin], gety_sp[bin]); //fill the spectrogram for this time slice of the waveform, and for all of the frequency bins
				/*
				cout<<""<<endl;
				cout<<"Time bin "<<t_bin<<endl;
				cout<<"Spectrogram bin "<<bin<<endl;
				cout<<"X bin content "<<getx_sp[bin]<<endl;
				cout<<"Y bin content "<<gety_sp[bin]<<endl;
				*/
			}
			
			delete waveform_piece;
			delete waveform_piece_spectrum;
	} //end loop over waveform pieces
     
     return Spectrogram;
    
}//end makeSpectrogram

void computeFFT(TGraph *inputWaveform, TGraph* (& SpectrumPhase)[3]){
		

		if(!inputWaveform) cout<<"Something is wrong!"<<endl; //get out of the function if something is wrong with the input waveform
		
		//We have to go through several steps to successfully take the FFT
		//We must
		//1) Interpolate the graph to make sure we have even spacing between points
		//2) Pad the waveform with zeros so that it has a power-of-two number of zeros
		//3) Do the FFT
		//4) Compute the magnitude
		//5) Make the Graph and pass it back out

		Double_t interpolationStep=0.5; //this is the interpolation distance for the waveform
		Int_t maxSamps = 1024; //the number of samples to be in the final waveform; here I'm using 2^10 for safety

	//////////////////////////////////////////
	// 1) interpolate the waveform
	//////////////////////////////////////////
		
		TGraph *interpolatedWaveform = FFTtools::getInterpolatedGraph(inputWaveform,interpolationStep); //get an interpolated graph
		
	//////////////////////////////////////////
	// 2) pad the waveform
	//////////////////////////////////////////
		
		Int_t interpNumSamps = interpolatedWaveform -> GetN(); //get the number of data points in the interpolated waveform
		Double_t *interpXVals = interpolatedWaveform -> GetX();  //get the x values of the interpolated waveforms
		Double_t *interpYVals = interpolatedWaveform -> GetY();  //get the y values of the interpolated waveforms
		
		Double_t paddedXVals[maxSamps]; //a array of size "maxSamps" to hold the new X-values
		Double_t paddedYVals[maxSamps]; //a array of size "maxSamps" to hold the new Y-valuesf
		
		for(int i=0; i<maxSamps; i++){ //loop over all of the data points in the final set
			if(i<interpNumSamps){ //if the iterator is less than the number of samples in the waveform
				paddedXVals[i]=interpXVals[i]; //then the new vector can hold the old values
				paddedYVals[i]=interpYVals[i]; //then the new vecotr can hold the old values
			} //end if loop
			else{
				paddedXVals[i]=paddedXVals[i-1]+interpolationStep; //otherwise, add another data point
				paddedYVals[i]=0; //make its y-value zero, which adds a bunch of zeros to the length of the waveform
			} //end else loop
			
		} //end padding for loop
		

		
	//////////////////////////////////////////
	// 3) now to do the FFT
	//////////////////////////////////////////
		
		
		//if(paddedYVals.size()!=paddedXVals.size()){//stop to make sure the X and Y arrays are of the same length
			//std::cerr << "The X and Y FFT arrays are of differing length! Stop!"<<std::endl; //something is wrong, the arrays for the FFT are of differing length
			//return -1; //and we should get out!
		//}
		
		//int length = paddedYVals.size(); //get the length of what we're feeding to the FFTW
		int length = maxSamps; //I think this should just be maxSamps? the prior syntax was a hold over from someone else's function
		double deltaT = paddedXVals[1]-paddedXVals[0]; //get the size of the timestep
		double deltaF = 1/(deltaT * length); //get the size of the frequency bin
		deltaF*=1e3; //convert this from GHz to MHz   
		deltaT*=1e3; //convert this from ns to us (nanoseconds to microseconds)
		//cout<<"The value of deltaF is "<<deltaF<<endl;

		fftw_complex *FFTWoutput = new fftw_complex [(length/2)+1]; //create an output array of size (N/2)+1
		double *FFTWinput = new double [length]; //create an input array for the FFT
		//fftw_complex *FFTWinput, *FFTWoutput;
		//fftw_complex *FFTWoutput
		//FFTWinput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
		//fftw_complex *FFTWoutput;
		//FFTWoutput = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * ((length/2)+1));
		fftw_plan thePlan = fftw_plan_dft_r2c_1d(length, FFTWinput, FFTWoutput, FFTW_MEASURE); //create a plan, with "newInput" as the input, "theOutput" as the output, and "FFTW_MEASURE" which tells FFTW to find an optimized plan
		
		if(!thePlan) { //check the plan before proceeding
				std::cerr<<"Something is wrong with the plan! Stop!"<<std::endl; //get out if something is wrong
		}
		

		for(int i=0; i<length; i++){
			FFTWinput[i]=paddedYVals[i];
		}
		
		for(int i=0; i<(length/2)+1; i++){//instantiate every element of the complex output with zeros
			FFTWoutput[i][0]=0.0; //instantiate the real part of the complex number to zero
			FFTWoutput[i][1]=0.0; //instantiate the imaginary part of the complex number to zero
		}
		
		fftw_execute(thePlan); //execute the plan, so actually do the FFT
		fftw_destroy_plan(thePlan); //destroy the plan, we don't need it anymore either
		


		
	//////////////////////////////////////////
	// 4) extract the magnitude
	//////////////////////////////////////////
		
		
		
		double *magnitude = new double [(length/2)+1]; //the y values, or magnitudes
		double *frequencies = new double [(length/2)+1]; //the x values, or frequencies
		double *phase = new double [(length/2)+1];
	   	double *PowerSpectralDensity = new double [(length/2)+1]; //y values for the power spectral density
		double tempF =0.; //the first frequency bin
		
		for(int j=0; j<(length/2)+1; j++){
				double realPart = FFTWoutput[j][0]; //get the real part
				double imagPart = FFTWoutput[j][1]; //get the imaginary part
				//cout<<"The real part is "<<realPart<< " and the imaginary part is "<<imagPart<<endl;
				//cout<<"The magnitude is "<<TMath::Sqrt((realPart * realPart)+(imagPart*imagPart))<<endl;
				//cout<<FFTWoutput[j][0]<<endl;
				magnitude[j] = TMath::Sqrt((realPart * realPart)+(imagPart*imagPart)); //fill in the magnitude
				phase[j]= TMath::ATan2(imagPart,realPart);
				PowerSpectralDensity[j] = magnitude[j]*magnitude[j]; //square the magnitude to get the PSD
				frequencies[j]=tempF; //fill in the frequency bin
				//cout<<"I just made entry (Spectrum,Frequency)= ("<<magnitude[j]<<" , "<<tempF<<")"<<endl;
				tempF+=deltaF; //increment up the frequency bin
				//cout<<"I raised to frequency bin "<<tempF<<endl;
		}
		fftw_free(FFTWoutput);
		
		//lastly, just need to multiply by deltaT to finish the conversion from discrete coefficients to continuous estimate
		for(int j=0; j<(length/2)+1; j++){
			magnitude[j]*=deltaT; //multiply by deltaT
			PowerSpectralDensity[j]*=(deltaT/( (double) length)); //multiply by deltaT and divide by the length of the number of points (see numerical recipes pg 653), and I think this should be N and not (N/2)+1... 
		}
	//////////////////////////////////////////
	// 5) make the graph
	//////////////////////////////////////////
		
		TGraph *theSpectrum = new TGraph((length/2)+1,frequencies,magnitude);
		TGraph *thePhase = new TGraph((length/2)+1, frequencies, phase);
		TGraph *thePSD = new TGraph((length/2)+1, frequencies, PowerSpectralDensity);
		SpectrumPhase[0]=theSpectrum;
		SpectrumPhase[1]=thePhase;
		SpectrumPhase[2]=thePSD;

	
	//////////////////////////////////////////
	// clean up: need to delete some stuff
	//////////////////////////////////////////

		//delete [] FFTWoutput;
		delete [] FFTWinput;
		delete [] magnitude;
		delete [] frequencies;
		delete [] PowerSpectralDensity;
		
		
	//////////////////////////////////////////
	// return the final product
	//////////////////////////////////////////
		
	//return theFFT; //return what we spent this time computing
		
	
}

//this function will return the time delay between between ant1 and ant2, given their locations (rho1,phi1) and (rho2,phi2), and the 
//direction of the incoming plane wave, (phiWave,thetaWave)
Double_t calcDeltaT(Double_t ant1[3],Double_t rho1, Double_t phi1, Double_t ant2[3],Double_t rho2, Double_t phi2, Double_t phiWave, Double_t thetaWave){
	phiWave*=TMath::DegToRad(); //convert phi in degrees to radians
	thetaWave*=TMath::DegToRad(); //convert theta in degress to radians
     
    Double_t nTopOfIce = 1.4;
	Double_t d1=TMath::Cos(thetaWave)*(ant1[2]*TMath::Tan(thetaWave)+rho1*TMath::Cos(phi1-phiWave)); //compute the distance from the source to antenna 1
	Double_t d2=TMath::Cos(thetaWave)*(ant2[2]*TMath::Tan(thetaWave)+rho2*TMath::Cos(phi2-phiWave)); //compute the distance from the source to antenna 2
	Double_t t1t2=(d2-d1)*AraGeomTool::nTopOfIce/TMath::C(); //compute the time it would take
	t1t2*=1e9; //convert it to nanoseconds
	return t1t2; //return it
}

//need a function to find the position of the antenna, and the rho, and the phi
//it should take an integer (which antenna), and return a vector of positions, a double rho and the double phi for this antenna
void fillAntennaPositions( int antenna_number, Double_t (& antenna_position)[3], Double_t &antenna_rho, Double_t &antenna_phi){ 
		
	AraGeomTool *araGeom=AraGeomTool::Instance(); //need an instance of the AraGeomTool

	antenna_position[0]=araGeom->getStationInfo(0)->getAntennaInfo(antenna_number)->getLocationXYZ()[0]; //get the x-position of the antenna
	antenna_position[1]=araGeom->getStationInfo(0)->getAntennaInfo(antenna_number)->getLocationXYZ()[1]; //get the y-position of the antenna
	antenna_position[2]=araGeom->getStationInfo(0)->getAntennaInfo(antenna_number)->getLocationXYZ()[2]; //get the z-position of the antenna
	
	antenna_rho=TMath::Sqrt(antenna_position[0]*antenna_position[0]+antenna_position[1]*antenna_position[1]); //compute the rho of the antenna (in cylindrical coordinates), rho = sqrt(x^2 + y^2)
	antenna_phi=TMath::ATan2(antenna_position[1],antenna_position[0]); //compute the phi of the antenna (in cylindrical coordinates), arctan(y/x)
}

// get the correlation map peak location (in 1 degree bins)
void GetCorrMapPeakAngle_1deg( TH2D *theCorrMap_input, int &peaktheta, int &peakphi, double &peakvalue) {
	
	double corr_peak = -1.;
	//for (int theta=1; theta<=180; theta++) { //loop over theta
	for (int theta=110; theta<=180; theta++) { //right now, only look at peaks above theta = 20
		for (int phi=1; phi<=360; phi++) { //loop over phi
			if(   (phi>30) && (phi<330)){ 
				continue; //if you're outside of where the sun should be, don't do anything
			}
			else{ //if you're where the sun should be, then proceed
				if (theCorrMap_input->GetBinContent(phi,theta) > corr_peak) {
					corr_peak = theCorrMap_input->GetBinContent(phi,theta); //get the bin content
					peaktheta = theta - 90; //map from bin number to angle in degrees
					peakphi = phi - 180; //map from bin number to angle in degrees
					peakvalue = corr_peak; //return the peak value
				}
			}
		}
	}
}

// get the correlation map peak location (in 1 degree bins) for the Solar Flare analysis, which needs special conditioning on allowed peaks
void GetCorrMapPeakAngle_1deg_SolarAnalysis( TH2D *theCorrMap_input, int &peaktheta, int &peakphi, double &peakvalue, int eventType) {

	double corr_peak = -1.;
	
	if(eventType==0){ //a flare event
		for (int theta=110; theta<=180; theta++) { //right now, only look at peaks above theta = 20
			for (int phi=1; phi<=360; phi++) { //loop over phi
				if(   (phi>30) && (phi<330) ){
				continue; //if you're outside of where the sun should be, don't do anything
				}
				else{ //if you're where the sun should be, then proceed
					if (theCorrMap_input->GetBinContent(phi,theta) > corr_peak) {
						//cout<<"                    I did get inside the if statement"<<endl;
						corr_peak = theCorrMap_input->GetBinContent(phi,theta); //get the bin content
						//cout<<"                    I got thorugh getting the corr peak"<<endl;
						peaktheta = theta - 90; //map from bin number to angle in degrees
						peakphi = phi - 180; //map from bin number to angle in degrees
						peakvalue = corr_peak; //return the peak value
					}
				}
			}
		}
	}//end flare event
	else if(eventType==1){ //a flare background; same code as for a flare
		for (int theta=110; theta<=180; theta++) { //right now, only look at peaks above theta = 20
			for (int phi=1; phi<=360; phi++) { //loop over phi
				if(   (phi>30) && (phi<330) ){ //temporarily saying if >500 so that this always gets satisfied; this is done in the particular case for the cal pulsers
				continue; //if you're outside of where the sun should be, don't do anything
				}
				else{ //if you're where the sun should be, then proceed
					if (theCorrMap_input->GetBinContent(phi,theta) > corr_peak) {
						corr_peak = theCorrMap_input->GetBinContent(phi,theta); //get the bin content
						peaktheta = theta - 90; //map from bin number to angle in degrees
						peakphi = phi - 180; //map from bin number to angle in degrees
						peakvalue = corr_peak; //return the peak value
					}
				}
			}
		}
		//cout<<"Finished the background event"<<endl;
	}//end flare background
	
	else if(eventType==2){ //a a cal pulser
		for (int theta=1; theta<=180; theta++) { //loop over theta
			for (int phi=1; phi<=360; phi++) { //loop over phi
				if (theCorrMap_input->GetBinContent(phi,theta) > corr_peak) {
					corr_peak = theCorrMap_input->GetBinContent(phi,theta); //get the bin content
					peaktheta = theta - 90; //map from bin number to angle in degrees
					peakphi = phi - 180; //map from bin number to angle in degrees
					peakvalue = corr_peak; //return the peak value
				}
			}
		}
		//cout<<"Finished the Cal Pulser check"<<endl;
	}//end cal pulser loop
}

//a function for dealing with the boundary cases of where you have reached the edge of the histogram (theta = 91 or -91, and phi = 181 or -181) and need to correclty re-wrap the iterator
void CheckCorrMapRange( int &theta, int &phi ) {

	if(theta > 90) theta -= 181; //if you've incremented over the 90 high mark, wrap all the way around to -90
	else if ( theta < -90 ) theta += 181; //if you've incremented below -91, wrap all the way around to 90

	if (phi > 180) theta -= 361; //if you've incremented over the 180 high mark, wrap all the way back around to -180
	else if ( phi < -180 ) theta += 361; //if you've incremented below -180, wrap all the way around to 180
}


//find the maximum value in an array
int GetMaxVector( vector <int> &array ) {
	int max = -1000;
	for (int bin=0; bin<(int)array.size(); bin++) {
		if (array[bin] > max) max = array[bin];
	}
	return max;

}

//find the minimum value in an array
int GetMinVector( vector <int> &array ) {
	int min = 1000;
	for (int bin=0; bin<(int)array.size(); bin++) {
		if (array[bin] < min) min = array[bin];
	}
	return min;
}


//a function that will return the boundary of a peak
TH2D *GetPeakBoundaryPlot ( TH2D *theCorrMap_input, int peaktheta, int peakphi, double peakfrac, int &max_theta, int &min_theta, int &max_phi, int &min_phi, int bin_x, int bin_y ) {

	//////////////////////////////////////////
	// define variables needed for this routine, including the histogram to be returned at the end
	//////////////////////////////////////////
	
				double peakbinvalue = theCorrMap_input->GetBinContent(peakphi+180,peaktheta+90); //the value of the peak passed to the function
				double peakrange = peakbinvalue * peakfrac; //a fraction of the peak value that defines the "width" of the peak you're looking for
				double binvalue_tmp; //a variable to hold the temporary bin value during scanning over candidates
				TH2D *histPeakBoundary = new TH2D("","",bin_x,-180,180,bin_y,-90,90); //the PeakBoundary histogram that will be returned at the end of the routine
				vector <int> thetas; //a vector of theta candidates
				vector <int> phis; // a vector of phi candidates

	//////////////////////////////////////////
	// find the rightmost theta boundary point
	//////////////////////////////////////////
    
				binvalue_tmp = peakbinvalue; //start at the peak value
				int theta = peaktheta; //start at the peak theta (in degrees)
				int phi = peakphi; //start at the peak phi (in degrees)
				
				while ( binvalue_tmp >= peakrange ) { //sweep to the right and look for the edge of the peak--ie, keep moving over until the bin contents are smaller than the "width" you wanted
					theta += 1; //increase the theta value
					CheckCorrMapRange( theta, phi ); //handle the boundary case
					binvalue_tmp = theCorrMap_input->GetBinContent(phi+180,theta+90); //get the peak value
				}

				theta -= 1; //the peak value for which this loop will stop will be LESS than peak range (ie, it will be outside your "width"), so step back one theta degree to get back to the actual edge
				CheckCorrMapRange( theta, phi ); //again, account for the edge case if you have to
				histPeakBoundary->SetBinContent( phi+180, theta+90, 1 ); //fill this entry in the boundary map
				thetas.push_back(theta); //record this theta value
				phis.push_back(phi); //record this phi value
				
				int org_phi = phi; //the phi of the original boundary
				int org_theta = theta; //the theta of the original boundary

    
    
	//////////////////////////////////////////
	// define variables that allow us to directionally creep across the histogram bins
	//////////////////////////////////////////
    
				// call bitset for direction define
				int pre_dir = 2; //down 0, left 1, up 2, right 3
				int backward = (pre_dir+2)%4;
				int dir=0; //just instantiate this value to be used later to indicate direction
				
				int search_cnt=0; // a counter; if we can't find the direction after 4 times, we are done
				int theta_tmp=0; //temporary theta value
				int phi_tmp=0; //temporary phi value
				binvalue_tmp = 0; // reset the bin value to 0 to continue the boundary search

	
	//////////////////////////////////////////
	// test the first cells immediately around the peak value
	//////////////////////////////////////////
		
				/*
				everytime we change the direction (to get the rotation right, we have to do %4 to get correct direction)
			
				left rotate = add 1 (example: imagine you are going down, or 0. Add 1, gets you to 1, and 1%4 = 1, which means you are now going left. This corresponds to a left turn from "south" to "east")
				right rotate = add 3 (example: you are going left, or 1. Add 3, gets you to 4, and  4%4=0, which means you are now going down. This coresponds to a right turn from "east" to "south"
			
			
				(is this wrong?) the very first rotate should be right, while all other rotates should be left
				*/
			
				while ( binvalue_tmp < peakrange && search_cnt<5) {

					//dir = (backward+3)%4; // right rotate
					dir = (pre_dir+1)%4; // left rotate starting from upward
					
					if (dir==0) { //search downward
						theta_tmp = theta-1; //move down in theta
						phi_tmp = phi; //hold phi constant
						CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
						binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
					}
					else if (dir==1) { //search leftward
						theta_tmp = theta; //hold theta constant
						phi_tmp = phi-1; //move phi to the left
						CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
						binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
					}
					else if (dir==2) { //search upward
						theta_tmp = theta+1; //move up in theta
						phi_tmp = phi; //hold phi constant
						CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
						binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
					}
					else if (dir==3) { //search rightward
						theta_tmp = theta; //hold theta constant
						phi_tmp = phi+1; //move phi to the right
						CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
						binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
					}

					//backward+=3; // for the next round
					pre_dir+=1; // turn left for the next round
					search_cnt++; //increase the counter
				}
				
				// shift that direction
				if (dir==0) theta-=1; // downward
				else if (dir==1) phi-=1; // leftward
				else if (dir==2) theta+=1; // upward
				else if (dir==3) phi+=1; // rightward

				CheckCorrMapRange( theta, phi ); //handle the boundary cases


	
	
	//////////////////////////////////////////
	// store what we learned from the first step
	//////////////////////////////////////////
  
				if ( search_cnt > 4 || (theta==org_theta && phi==org_phi) ) { //  check if we exceeded the search count or if we reached the original location; if so, call it a day and leave the loop
					max_theta = GetMaxVector( thetas ); //find the maximum theta
					min_theta = GetMinVector( thetas ); //find the minimum theta
					max_phi = GetMaxVector( phis ); //find the maximum phi
					min_phi = GetMinVector( phis ); //find the minimum phi
					thetas.clear(); //clear the theta array
					phis.clear(); //clear the phi array
					return histPeakBoundary; //return the boundary map
				}
				else { //otherwise, store the information you found; the first step locates the adjacent bin which fits within the width
					histPeakBoundary->SetBinContent( phi+180, theta+90, 1 ); //fill in this entry in the histogram
					thetas.push_back(theta); //store the theta value
					phis.push_back(phi); //store the phi value
					pre_dir = dir; //store pre_dir
					backward = (pre_dir+2)%4; //store backward
				}

	//////////////////////////////////////////
	// now we did the first step, now we do all other steps in left handed search until it meets the very begining bin
	//////////////////////////////////////////

				int cnt_loop = 0; //a loop counter variable

				// search until we find the original location
				while ( theta!=org_theta || phi!=org_phi ) { // don't look at the original points; also, we don't need to worry about single peak point case as it should have selected in the first run

					search_cnt = 0; //reset the counter to 0
					binvalue_tmp = 0; // reset to 0

					while ( binvalue_tmp < peakrange && search_cnt < 5) {
							
						dir = (backward+1)%4; // left rotate
						if (dir==0) { //search downward
							theta_tmp = theta-1; //move down in theta	
							phi_tmp = phi; //hold phi constant
							CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
							binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value	
						}
						else if (dir==1) { //search leftward
							theta_tmp = theta; //hold theta constant
							phi_tmp = phi-1; //move phi to the left
							CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
							binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
						}
						else if (dir==2) { //search upward
							theta_tmp = theta+1; //move up in theta
							phi_tmp = phi; //hold phi constant
							CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
							binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
						}
						else if (dir==3) { //search rightward
							theta_tmp = theta; //hold theta constant
							phi_tmp = phi+1; //move phi to the right
							CheckCorrMapRange( theta_tmp, phi_tmp ); //handle the boundary cases
							binvalue_tmp = theCorrMap_input->GetBinContent(phi_tmp+180,theta_tmp+90); //get the temporary value
						}

					
						backward+=1; // for the next round
						search_cnt++; //increase the counter
					}

					// shift that direction
					if (dir==0) theta-=1; // downward
					else if (dir==1) phi-=1; // leftward
					else if (dir==2) theta+=1; // upward
					else if (dir==3) phi+=1; // rightward

					CheckCorrMapRange( theta, phi ); //handle the boundary cases

					// check if we exceeded search count or we reached original location
					if ( search_cnt > 4 || (theta==org_theta && phi==org_phi) ) { // we are done
						break; // get out from this while loop to avoid infinite loop
					}
					else {
						// fill in histogram
						histPeakBoundary->SetBinContent( phi+180, theta+90, 1 ); //fill in that boundary point
						thetas.push_back(theta); //store the theta
						phis.push_back(phi); //store the phi
						
						pre_dir = dir; //store pre_dir
						backward = (pre_dir+2)%4; //store backward
					}

					cnt_loop++; //increment up the loop counter

				} // while loop that's running until you find the original location
			
	
	
	
	//////////////////////////////////////////
	// write out the final values of interest
	//////////////////////////////////////////
			
				max_theta = GetMaxVector( thetas ); //find the maximum theta
				min_theta = GetMinVector( thetas ); //find the minimum theta
				max_phi = GetMaxVector( phis ); //find the maximum phi
				min_phi = GetMinVector( phis ); //find the minimum phi
				thetas.clear(); //clear the theta array
				phis.clear(); //clear the phi array
				return histPeakBoundary; //return the histogram containing the Peak Boundary



}

// get peak range area plot and number of bins
TH2D *GetPeakAreaPlot ( TH2D *theCorrMap, TH2D *theCorrBoundary, int peaktheta, int peakphi, double peakfrac, int max_theta, int min_theta, int max_phi, int min_phi, double &area, int bin_x, int bin_y ) {

	double peakbinvalue = theCorrMap->GetBinContent(peakphi+180,peaktheta+90); //the value of the peak passed to the function
	double peakrange = peakbinvalue * peakfrac;  //a fraction of the peak value that defines the "width" of the peak you're looking for
	double binvalue_tmp; //a variable to hold the temporary bin value during scanning over candidates
	
	TH2D *histPeakArea = new TH2D("","",bin_x,-180,180,bin_y,-90,90); // a histogram to hold the plot of the peak Area

	area = 0.; //a variable to hold the computed area

    // next bin values
	int next_theta=0;
	int next_phi=0;

    int in_area=0; // flag for if the found spot is inside the peak in peak area
    // ok, starting from very low theta, phi bins, sweep the angle range and get area

	for (int theta=min_theta; theta<=max_theta; theta++) {
		in_area = 0; //reset the value
		for (int phi=min_phi; phi<=max_phi; phi++) {
			// if we meet boundary
			if ( theCorrBoundary->GetBinContent(phi+180,theta+90) > 0 ) { //see if you've found a boundary
				histPeakArea->SetBinContent(phi+180,theta+90,1); //fill in that boundary spot
				area += cos( (double)theta * DEG2RAD ); //add it's value to the area amount

				// calculate next bin location, and set it's "in_area" flag
				if ( phi+1 <= max_phi ) { //try and move up, if that spot is still allowed within the boundary we've identified
					if ( theCorrMap->GetBinContent(phi+1+180,theta+90) >= peakrange ) in_area = 1; //if you're in the area, set the flag to yes (1)
					else if ( theCorrMap->GetBinContent(phi+1+180,theta+90) < peakrange ) in_area = 0; //if you're not in the area, set the flag to no (0)
                }
            }
			else { //if you're not at the boundary you must be inside!
				if ( in_area == 1 ) {
					histPeakArea->SetBinContent(phi+180,theta+90,1); //fill in the area
					area += cos( (double)theta * DEG2RAD ); //add it to the area count
				}
            }
        } // loop over phi
    } // loop over theta

    return histPeakArea;

}

