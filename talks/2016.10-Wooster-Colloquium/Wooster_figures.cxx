////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  analysis.cxx 
////  Making coherent sum waveforms of the Feb 15 flare
////
////    September 2016,  clark.2668@osu.edu
////    Figures for Wooster colloquium talk
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
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
#include "AraEventCorrelator.h"
#include "AraAntennaInfo.h"
#include "AraGeomTool.h"
#include "FFTtools.h"
#include "AraStationInfo.h"
#include "AraAntennaInfo.h"

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

using namespace std;
const double DEG2RAD=0.017453292;   // radians/degree (convert to radians by multiplying degrees by rad/deg)
const double RAD2DEG=57.2957795;    // degree/rad (convert to degrees by multiplying radians by deg/rad)

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
TGraph *getCorrelationGraph_timeDomain_SqrtN(TGraph *gr1, TGraph *gr2);
//TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2);
TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]);
TGraph *getCorrelationGraph_timeDomain_powerNormalization(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]);

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
	
	gStyle -> SetOptStat(111111111);
	bool PrintMode=false; //if false, don't print stuff out, if true, print stuff out
	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;
	
	int num_found =0; //a variable to hold the number of events we've found in the first file
	int num_found_max=2; //a variable to hold the maximum number of waveforms you want to allow in the first input file loop
	
	AraEventCorrelator *theCorrelator = AraEventCorrelator::Instance(); //make a correlator	
	AraGeomTool *geomTool = new AraGeomTool();

	//TString dataPath("/home/brianclark/SolarFlares/"); //the location of the files
	TString dataPath("/home/clark.2668/workspace/SolarFlares/"); //the location of the files
	
	TChain chain("eventTree"); //this for the events for the exterior loop
	for(int file=1; file<argc; file++){
		TString fileKey(argv[file]); //a file key
		chain.Add(dataPath+fileKey); //add files to the chain
	}

	RawIcrrStationEvent *rawIcrrEvPtr=0; //make the raw event pointer
	chain.SetBranchAddress("event",&rawIcrrEvPtr); //set the branch address
	Int_t numEntries = chain.GetEntries(); //get the number of entries
	
	for(Int_t event=0;event<numEntries;event++) { //loop over all events
	//for(Int_t event=0;event<10;event++) { //loop over only this number of events
		if(num_found>num_found_max) {continue;} //bounce out if you have found more than this many events
		chain.GetEvent(event);  //This line gets the RawIcrr or RawAtri Event
		UsefulIcrrStationEvent *realIcrrEvPtr = new UsefulIcrrStationEvent(rawIcrrEvPtr, AraCalType::kLatestCalib); //make a real event pointer for the Icrr station		
		
		if( realIcrrEvPtr->isCalPulserEvent() != true) { //plot only the non-Calpulser, non RF events
			if(rawIcrrEvPtr -> trig.trigType == 68){
				delete realIcrrEvPtr; //delete this
				continue; //bounce out of doing any other calculation
			} //if this is a RF trigger, move on //trig.trigType!=68 asks for non-RF/ software trigger stuff (bit 2 high = 4) + (bit 6 high = 64) = 68
			                                                                                                                                             
			int unixtime_1;
			unixtime_1 = rawIcrrEvPtr->head.unixTime; // get  event unixTime
			time_t test_time_1 = unixtime_1;
			tm *test_time_tm_1 = gmtime( &test_time_1 );
			int year = test_time_tm_1->tm_year+1900;
			int month = test_time_tm_1->tm_mon+1;
			int day = test_time_tm_1->tm_mday;
			int hour = test_time_tm_1->tm_hour;
			int min = test_time_tm_1->tm_min;
			int sec = test_time_tm_1->tm_sec;
			
			num_found++; //increment up the number we have found that pass our cuts
			char histName[150];
			char saveName[150];
			
			TGraph *Waveform = realIcrrEvPtr->getGraphFromRFChan(2);
			TGraph *Wavedform_Box = FFTtools::getBoxCar(Waveform,5);
			
			TCanvas *save = new TCanvas("","",850,1100);
			save -> Divide(2,1);
			save->cd(1);
			Waveform -> Draw();
			save->cd(2);
			Waveform_Box -> Draw();
			sprintf(saveName,"/home/clark.26682/workspace/SolarFlares/results/today/2016.09_Wooster_%d",event);
			save -> SaveAs(saveName);
			
			
				//if(hour==2){continue;} //use this if you want to exclude the 2 o clock hour

				/*
				TH2D *V_map = theCorrelator->getInterferometricMapHilbert(realIcrrEvPtr,AraAntPol::kVertical,AraCorrelatorType::kPlaneWave); //Get the V-pol correlation map
				int V_corr_map_peak_theta=0; //a variable to hold the correlation map peak theta value
				int V_corr_map_peak_phi=0; //a variable to hold the correlation map peak phi value
				double V_corr_map_peak_value=0; //a variable to hold the value of the correlation map peak value
				//GetCorrMapPeakAngle_1deg(V_map,V_corr_map_peak_theta,V_corr_map_peak_phi,V_corr_map_peak_value);
				GetCorrMapPeakAngle_1deg_SolarAnalysis(V_map,V_corr_map_peak_theta,V_corr_map_peak_phi,V_corr_map_peak_value,eventType);
				delete V_map;
				*/
				
				/*
				TGraph *correlation_with_self = FFTtools::getCorrelationGraph(coherentSum_normalized,coherentSum_normalized);
				Int_t numSamps_self = correlation_with_self -> GetN();
				Double_t *yVals_self = correlation_with_self -> GetY();
				double selfMax = 0.;
				double max_value_self = TMath::MaxElement(numSamps_self,yVals_self);
				double min_value_self = TMath::MinElement(numSamps_self,yVals_self);
				if(abs(min_value_self)>abs(max_value_self)){ selfMax = abs(min_value_self);} //either the min value is abs_value larger, and in which case it's the biggest correlation, albiet negative
				else{selfMax = max_value_self;} //otherwise the max is larger, or they are equal, and in which case, why not just use the max_value //this should not happen very often
				*/   
				     
				//now, to initiate another loop over files
				cout<<"Working on loop 1 event "<<event<<endl;

		}//close cal pulser loop		
		delete realIcrrEvPtr; //delete this
	}//close the outer loop over events
	
}//close the main program


//a function to do a time domain correlation of two tgraphs
// the graphs must have the same interpolation time step, and start at the same time window
//this function will normalize the cross correlation coefficient according to the amount of power in the overlapping waveforms
//this is done to follow A. Romero Wolf, "An Interferometric Analysis Method for Radio Impulses from Ulra-high Energy Particle Showers", arXiv 1304.5663v1
//particularly equation 10
TGraph *getCorrelationGraph_timeDomain_powerNormalization(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]){

	int length1= gr1->GetN(); //get the number of points in graph 1
	int length2= gr2->GetN(); //get the number of points in graph 2

	if(length1!=length2){cout<<"The Waveforms are not of the same length! This many not work..."<<endl;}

	double x1,y1=0.; //starting values for the first data point
	double x2,y2=0.; //starting values for the second data point
	
	int padLength = 3*length1; // x3 to allow for one whole waveform to go before, middle, and after during the correlation
	
	TGraph *gr1Padded = FFTtools::padWaveToLength(gr1,padLength); //pad them to the same length for ease
	TGraph *gr2Padded = FFTtools::padWaveToLength(gr2,padLength); //pad them to the same length for ease

	gr1Padded->GetPoint(1,x2,y2); //get this first point
	gr1Padded->GetPoint(0,x1,y1); //get this second point
	double deltaT = x2-x1; //find the time step for the waveform
	double testx,testy; //for checking if the two first time points are the same
	gr2Padded->GetPoint(0,testx,testy); //just get the first time point of the other array
	double waveOffset = x1-testx; //get the time offset of wave 2 from wave 1
	if(abs(waveOffset)>0.0001){cout<<"The two graphs don't have the same starting point, this could be a problem..."<<endl;}
	int maxDelay = length1; //the maximum time delay should be to delay one waveform by the entire lenght of the other
	
	Double_t *oldX1 = gr1Padded->GetX();
	Double_t *oldY1 = gr1Padded->GetY();
	
	Double_t *oldX2 = gr2Padded->GetX();
	Double_t *oldY2 = gr2Padded->GetY();

	vector <double> time_lag;
	vector <double> correlationValues;
	vector <double> numContributing;

	//cout<<"About to start here"<<endl;
	for(int delay = -maxDelay; delay<maxDelay; delay++){
	//for(int delay = -maxDelay; delay<-5998; delay++){
	    double sum=0.;
	    double integral1=0.;
	    double integral2=0.;
		int num_contributing_points = 0.;
		double denominator=0.;
		for(int i=length1; i<2*length1; i++){//integrate over just the center region of waveform 1
		//for(int i=length1; i<6002; i++){//integrate over just the center region of waveform 1
			int j=i+delay; //but, integrate over all possible delayed regions of waveform 2 //for that given i step, you want to look at i+delay in the other waveform
			double yVal1 = oldY1[i]; //one of the old voltages
			double yVal2 = oldY2[j]; //the other old voltage
			double product = yVal1*yVal2/**deltaT*/; //the product of the two with deltaT to make the time integral, added 8/5, removed 8/7
			if((abs(yVal1)>0)&&(abs(yVal2)>0)){ //only if both are non-zero should you do anything with it
				num_contributing_points++; //increment up the number of contributing points
				sum+=product; //increase the sum
				//denominator+=((yVal1*yVal1)+(yVal2*yVal2))/2; //8/7, trying in the denominator (x^2 + y^2)/ 2
				integral1+=(yVal1*yVal1/**deltaT*/); //add to the integral for this step //8/7 deltaT removed to avoid units conflict
				integral2+=(yVal2*yVal2/**deltaT*/); //add to the integral for this step //8/7, deltaT removed to avoid units conflict
			} //increase the number of counted overlap, the thing added to the sum, and the integrals for the denominator
			//integral1 += FFTtools::integrateVoltageSquared(gr1,i,1+length1); //integrate the waveform from i to length1 away from i now
			//integral2 += FFTtools::integrateVoltageSquared(gr2,j,j+length1); //integrate the second (shifted) waveform from i to length1 away from i now
		}
		time_lag.push_back(-delay*deltaT); //save this time lag value for plotting
		if(num_contributing_points==0){
			correlationValues.push_back(0.); 
			numContributing.push_back(0.);
		} //if there is no content, input some zeros
		else{
			//integral1= FFTtools::integrateVoltageSquared(gr1);
			//integral2=FFTtools::integrateVoltageSquared(gr2); 
			//integral1=FFTtools::sumVoltageSquared(gr1,0,(3*length1)-1); //when the deltaT is absent above, we need to divide by sum of squared voltages, not the square*time which would be power; this is done to keep the units the same
			//integral2=FFTtools::sumVoltageSquared(gr2,0,(3*length1)-1);
			correlationValues.push_back(sum/TMath::Sqrt(integral1)/TMath::Sqrt(integral2));  //according to equation 10 in arXiv 1304.5663v1
			//correlationValues.push_back(sum/denominator);
			numContributing.push_back((double) num_contributing_points);
		} //otherwise, input these values
	}
	
	if(time_lag.size()!=correlationValues.size()){cout<<"Mismatch in size of time lag and correlation values, something is wrong!"<<endl;}
	
	TGraph *output = new TGraph(time_lag.size(),&time_lag[0],&correlationValues[0]);
	final[0]=new TGraph(time_lag.size(),&time_lag[0],&numContributing[0]);

	delete gr1Padded; //delete this
	delete gr2Padded; //delete this
	return output; //return the function
}


//a function to do a time domain correlation of two tgraphs
// the graphs must have the same interpolation time step, and start at the same time window

TGraph *getCorrelationGraph_timeDomain_SqrtN(TGraph *gr1, TGraph *gr2){

	int length1= gr1->GetN(); //get the number of points in graph 1
	int length2= gr2->GetN(); //get the number of points in graph 2

	if(length1!=length2){cout<<"The Waveforms are not of the same length! This many not work..."<<endl;}

	double x1,y1=0.; //starting values for the first data point
	double x2,y2=0.; //starting values for the second data point
	
	int padLength = 3*length1; // x3 to allow for one whole waveform to go before, middle, and after during the correlation
	
	TGraph *gr1Padded = FFTtools::padWaveToLength(gr1,padLength); //pad them to the same length for ease
	TGraph *gr2Padded = FFTtools::padWaveToLength(gr2,padLength); //pad them to the same length for ease

	gr1Padded->GetPoint(1,x2,y2); //get this first point
	gr1Padded->GetPoint(0,x1,y1); //get this second point
	double deltaT = x2-x1; //find the time step for the waveform
	gr2Padded->GetPoint(0,x2,y2);
	double waveOffset = x1-x2; //get the time offset of wave 2 from wave 1
	if(abs(waveOffset)>0.0001){cout<<"The two graphs don't have the same starting point, this could be a problem..."<<endl;}
	//double maxDelay = deltaT * length1;
	int maxDelay = length1;
	
	Double_t *oldX1 = gr1Padded->GetX();
	Double_t *oldY1 = gr1Padded->GetY();
	
	Double_t *oldX2 = gr2Padded->GetX();
	Double_t *oldY2 = gr2Padded->GetY();
	
	//Double_t *time_lag[padLength]={};
	//Double_t *correlationValues[padLength]={};

	vector <double> time_lag;
	vector <double> correlationValues;
	
	//cout<<"     MaxDelay is "<<maxDelay<<endl;
	//cout<<"     Length1 is "<<length1<<endl;

	for(int delay = -maxDelay; delay<maxDelay; delay++){
	//for(int delay = -maxDelay; delay<-5998; delay++){
		//cout<<"I'm working on delay "<<delay<<endl;
	    double sum=0.;
		int num_contributing_points = 0;
		for(int i=length1; i<2*length1; i++){//integrate over just the center region of waveform 1
		//for(int i=length1; i<6002; i++){//integrate over just the center region of waveform 1
			//cout<<"I'm working on i= "<<i<<endl;
			int j=i+delay; //but, integrate over all possible delayed regions of waveform 2
			//cout<<"     I'm working on j= "<<j<<endl;
			double yVal1 = oldY1[i];
			//cout<<"          First Y value is "<<yVal1<<endl;
			double yVal2 = oldY2[j];
			//cout<<"          Second Y value is "<<yVal2<<endl;
			double product = yVal1*yVal2*deltaT; //deltaT added 8/5
			//cout<<"          The product of the two is "<<product<<endl;
			if((abs(yVal1)>0)&&(abs(yVal2)>0)){ num_contributing_points++;} //increase the number of counted overlap
			//cout<<"The number of contributing points now is "<<num_contributing_points<<endl;
			sum+=product;
		}
		//cout<<"For time delay "<<delay<<" I'm working with time lag "<<(-maxDelay*deltaT)<<" and sum value "<<sum/((double) num_contributing_points)<<endl;
		time_lag.push_back(-delay*deltaT);
		if(num_contributing_points==0){correlationValues.push_back(0.) /*num_contributing_points=1*/;}
		else{correlationValues.push_back(sum/sqrt((double) num_contributing_points));}		
		//correlationValues.push_back(sum/((double) num_contributing_points));
	}
	
	if(time_lag.size()!=correlationValues.size()){cout<<"Mismatch in size of time lag and correlation values, something is wrong!"<<endl;}
	
	TGraph *output = new TGraph(time_lag.size(),&time_lag[0],&correlationValues[0]);

	delete gr1Padded;
	delete gr2Padded;

	//delete output;
	return output;
}



//a function to do a time domain correlation of two tgraphs
// the graphs must have the same interpolation time step, and start at the same time window

TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]){

	int length1= gr1->GetN(); //get the number of points in graph 1
	int length2= gr2->GetN(); //get the number of points in graph 2

	if(length1!=length2){cout<<"The Waveforms are not of the same length! This many not work..."<<endl;}

	double x1,y1=0.; //starting values for the first data point
	double x2,y2=0.; //starting values for the second data point
	
	int padLength = 3*length1; // x3 to allow for one whole waveform to go before, middle, and after during the correlation
	
	TGraph *gr1Padded = FFTtools::padWaveToLength(gr1,padLength); //pad them to the same length for ease
	TGraph *gr2Padded = FFTtools::padWaveToLength(gr2,padLength); //pad them to the same length for ease

	gr1Padded->GetPoint(1,x2,y2); //get this first point
	gr1Padded->GetPoint(0,x1,y1); //get this second point
	double deltaT = x2-x1; //find the time step for the waveform
	gr2Padded->GetPoint(0,x2,y2);
	double waveOffset = x1-x2; //get the time offset of wave 2 from wave 1
    if(abs(waveOffset)>0.0001){cout<<"The two graphs don't have the same starting point, this could be a problem..."<<endl;}
	//double maxDelay = deltaT * length1;
	int maxDelay = length1;
	
	Double_t *oldX1 = gr1Padded->GetX();
	Double_t *oldY1 = gr1Padded->GetY();
	
	Double_t *oldX2 = gr2Padded->GetX();
	Double_t *oldY2 = gr2Padded->GetY();
	
	//Double_t *time_lag[padLength]={};
	//Double_t *correlationValues[padLength]={};

	vector <double> time_lag;
	vector <double> correlationValues;
	vector <double> numContributing;
	
	//cout<<"     MaxDelay is "<<maxDelay<<endl;
	//cout<<"     Length1 is "<<length1<<endl;

	for(int delay = -maxDelay; delay<maxDelay; delay++){
	//for(int delay = -maxDelay; delay<-5998; delay++){
		//cout<<"I'm working on delay "<<delay<<endl;
	    double sum=0.;
		int num_contributing_points = 0;
		for(int i=length1; i<2*length1; i++){//integrate over just the center region of waveform 1
		//for(int i=length1; i<6002; i++){//integrate over just the center region of waveform 1
			//cout<<"I'm working on i= "<<i<<endl;
			int j=i+delay; //but, integrate over all possible delayed regions of waveform 2
			//cout<<"     I'm working on j= "<<j<<endl;
			double yVal1 = oldY1[i];
			//cout<<"          First Y value is "<<yVal1<<endl;
			double yVal2 = oldY2[j];
			//cout<<"          Second Y value is "<<yVal2<<endl;
			double product = yVal1*yVal2/**deltaT*/; //deltaT added 8/5, not sure if wanted on 8/7
			//cout<<"          The product of the two is "<<product<<endl;
			if((abs(yVal1)>0)&&(abs(yVal2)>0)){ num_contributing_points++;} //increase the number of counted overlap
			//cout<<"The number of contributing points now is "<<num_contributing_points<<endl;
			sum+=product;
		}
		//cout<<"For time delay "<<delay<<" I'm working with time lag "<<(-maxDelay*deltaT)<<" and sum value "<<sum/((double) num_contributing_points)<<endl;
		time_lag.push_back(-delay*deltaT);
		if(num_contributing_points==0){correlationValues.push_back(0.) /*num_contributing_points=1*/; numContributing.push_back(0.);}
        //else{correlationValues.push_back(sum/((double) num_contributing_points));}
		else{correlationValues.push_back(sum/((double) num_contributing_points)); numContributing.push_back((double) num_contributing_points);}                                                                                                                            	}
	
	if(time_lag.size()!=correlationValues.size()){cout<<"Mismatch in size of time lag and correlation values, something is wrong!"<<endl;}
	
	TGraph *output = new TGraph(time_lag.size(),&time_lag[0],&correlationValues[0]);
	final[0]=new TGraph(time_lag.size(),&time_lag[0],&numContributing[0]);

	delete gr1Padded;
	delete gr2Padded;

	//delete output;
	return output;
}


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

