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
TGraph *getCorrelationGraph_timeDomain_SqrtN(TGraph *gr1, TGraph *gr2);
//TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2);
TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]);
TGraph *getCorrelationGraph_timeDomain_powerNormalization(TGraph *gr1, TGraph *gr2,TGraph* (& final)[1]);
TH2D *ShiftedMap ( TH2D *inputMap, int ThetaShift, int PhiShift);
void CheckRange_byBin( int &theta, int &phi );
double TransformPredictedAzimuth(double inputAzimuth, double station_long, bool ReRange =true);
double TransformMapPeakToGlobalFrame(double inputAzimuth);
bool SaturationCut(UsefulIcrrStationEvent *realIcrrEvPtr); // a function to check and see if the event has any saturated waveforms
double GetTotalArea ( TH2D *theCorrMap, int peaktheta, int peakphi, double peakfrac);

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
	
	gStyle->SetOptStat(0);                                                                                                                                                                    
	
	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	int num_found =0; //a variable to hold the number of events we've found in the first file
	int num_found_max=2; //a variable to hold the maximum number of waveforms you want to allow in the first input file loop
		
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
     
	//set up the station geometry, and establish the longitude and latitude
	int station = 0;
	TVector3 stationVector = geomTool->getStationVector(station); //get the vector of easting, northings, and upping for station 0 (testbed) in ARRAY coordinates
	double station_long = AraGeomTool::getLongitudeFromArrayCoords(stationVector[1], stationVector[0], 2011); //get the longitude
	double station_lat = AraGeomTool::getGeographicLatitudeFromArrayCoords(stationVector[1], stationVector[0], 2011); //get the latitude
     
	Int_t numEntries = chain.GetEntries(); //get the number of entries
     
	for(Int_t event=0;event<numEntries;event++) { //loop over all events
     
		if(event%50==0) cout<<"On event "<<event<<endl;
          //if(num_found>num_found_max) break;
     
		chain.GetEvent(event);  //This line gets the RawIcrr or RawAtri Event
		CutsChain.GetEvent(event);  //This line gets the RawIcrr or RawAtri Event
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
          
		if(hour!=2) {continue;} //impose a time constraint
		
		UsefulIcrrStationEvent *realIcrrEvPtr = new UsefulIcrrStationEvent(rawIcrrEvPtr, AraCalType::kLatestCalib); //make a real event pointer for the Atri station event with a run dependent calibration

          
          //now, to compute cut values and make cuts
          
			double area_ratio = AreaTotal_value/AreaContour_value;
			bool saturatedcut = SaturationCut(realIcrrEvPtr); //check to see if this event saturated
			double slantedcut = (SlantedCut_slope_value * CorrMapPeak_value) + SlantedCut_yint_value;
          
			if( //2==2
				(CorrMapPeak_value >= CorrMapPeakCut_value)
				&&
				(VpeakRMS_value > VpeakRMSCut_value)
				&&
				(AreaContour_value >= AreaMinCut_value)
				&&
				(AreaContour_value <= AreaMaxCut_value)
				&&
				(!saturatedcut)
				&&
				(VpeakRMS_value >= slantedcut)
				&&
				(area_ratio<RatioCut_value)
			){ //if the event passes all cuts, do something with it
          
				num_found++;
               
				/*
				//store the event answer
				double peak_phi_corrected = TransformMapPeakToGlobalFrame(PeakPhi_value); //V_corr_map_peak_phi;
				eventTime.push_back((double) min + ((double)sec)/60. + ((double)(hour-2))*60.); //store the time in minutes
				peakTheta.push_back((double) PeakTheta_value);
				peakPhi.push_back((double) peak_phi_corrected );
               
				//now, get the predicted answer for the current time
				double this_event_phi = AraGeomTool::getSunAzimuthLongLat(station_long,station_lat,unixtime_1);
				double this_event_phi_corrected = TransformPredictedAzimuth(this_event_phi,station_long);
				double this_event_theta = AraGeomTool::getSunZenithLongLat(station_long, station_lat, unixtime_1);
				reportedpeakTheta.push_back(this_event_theta);
				reportedpeakPhi.push_back(this_event_phi_corrected);

				double this_event_phi_tshift = AraGeomTool::getSunAzimuthLongLat(station_long,station_lat,unixtime_tshift);
				double this_event_phi_corrected_tshift = TransformPredictedAzimuth(this_event_phi,station_long);
				double this_event_theta_tshift = AraGeomTool::getSunZenithLongLat(station_long, station_lat, unixtime_tshift);
				//eventTime_tshift.push_back((double) min_tshift + ((double)sec_tshift)/60. + ((double)(hour_tshift-2))*60.); //store the time in minutes
				reportedpeakTheta_tshift.push_back(this_event_theta_tshift);
				reportedpeakPhi_tshift.push_back(this_event_phi_corrected_tshift);
				*/
               
               
				if(event%10==0){
					char vtitle[150];
					char htitle[150];
					if(min<10){
						sprintf(vtitle,"VPol Interferometric Map, %d:0%d GMT",hour,min);
						sprintf(htitle,"HPol Interferometric Map, %d:0%d GMT",hour,min);
					}
					else{
						sprintf(vtitle,"VPol Interferometric Map, %d:%d GMT",hour,min);
						sprintf(htitle,"HPol Interferometric Map, %d:%d GMT",hour,min);
					}
				
					double TBx[3] = {-0.5986468972,0.8010128835,-0.000399302}; //this is from the ARA coordinates document, and contains the easting, northing, and upping components of the TB origin relative to the array (and therefore global) origin
					double angle_from_0 = TMath::ATan2(TBx[1],TBx[0]) * TMath::RadToDeg(); //I can compute this as the angle from zero because ATan2 is intelligent, and knows I mean the angle from zero!
				
					TH2D *ShiftedMapV= ShiftedMap(VMap_30m,0,(int) angle_from_0);

					TH2D * HMap_30m = theCorrelator -> getInterferometricMapHilbert_RT_Icrr(settings1, detector, realIcrrEvPtr, Hpol, false, 0, 2); //get an interferometric map with just BH channels
					TH2D *ShiftedMapH= ShiftedMap(HMap_30m,0,(int) angle_from_0);

				
					char save_title[150];
				
					TCanvas *VCanvas = new TCanvas("","",4*1600,4*900);
						ShiftedMapV->Draw("colz");
						ShiftedMapV->GetXaxis()->SetTitle("Phi (deg)");
						ShiftedMapV->GetYaxis()->SetTitle("Theta (deg)");
						ShiftedMapV->GetYaxis()->SetTitleOffset(0.9);
						ShiftedMapV->GetXaxis()->SetTitleSize(0.05);
						ShiftedMapV->GetXaxis()->SetLabelSize(0.05);
						ShiftedMapV->GetYaxis()->SetTitleSize(0.05);
						ShiftedMapV->GetYaxis()->SetLabelSize(0.05);
						ShiftedMapV->GetZaxis()->SetLabelSize(0.05);
						ShiftedMapV->SetTitle(vtitle);
						gPad->SetRightMargin(0.12);
						gPad->SetBottomMargin(0.12);
						sprintf(save_title,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/AprilMaps/%d.%d.%d_April_VSunMaps_Event%d_%d.%d.jpg",year_now,month_now,day_now,event,hour,min);
						VCanvas->SaveAs(save_title);
					delete VCanvas;
					delete ShiftedMapV;
					
					TCanvas *HCanvas = new TCanvas("h","h",4*1600,4*900);
						ShiftedMapH->Draw("colz");
						ShiftedMapH->GetXaxis()->SetTitle("Phi (deg)");
						ShiftedMapH->GetYaxis()->SetTitle("Theta (deg)");
						ShiftedMapH->GetYaxis()->SetTitleOffset(0.9);
						ShiftedMapH->GetXaxis()->SetTitleSize(0.05);
						ShiftedMapH->GetXaxis()->SetLabelSize(0.05);
						ShiftedMapH->GetYaxis()->SetTitleSize(0.05);
						ShiftedMapH->GetYaxis()->SetLabelSize(0.05);
						ShiftedMapH->GetZaxis()->SetLabelSize(0.05);
						ShiftedMapH->SetTitle(htitle);
						gPad->SetRightMargin(0.12);
						gPad->SetBottomMargin(0.12);
						sprintf(save_title,"/home/clark.2668/workspace/TrunkSolarFlares/results/today/AprilMaps/%d.%d.%d_April_HSunMaps_Event%d_%d.%d.jpg",year_now,month_now,day_now,event,hour,min);
						HCanvas->SaveAs(save_title);
					delete HCanvas;
					delete ShiftedMapH;
				
				}
          
			}
          
		delete realIcrrEvPtr;
          
	}

	//delete theCorrelator;
	//delete settings1;
	
}//close the main program



bool SaturationCut(UsefulIcrrStationEvent *realIcrrEvPtr){
	bool saturated = false;
	for(int antenna_iterator=0; antenna_iterator<16; antenna_iterator++){ //begin the loop over all of the other antennas
		if(saturated == true) break;
		TGraph *Waveform=realIcrrEvPtr->getGraphFromRFChan(antenna_iterator); //extract the waveform
		int numPoints= Waveform->GetN(); //get the number of points in the graph
		Double_t *yVals = Waveform->GetY();
		for(int i =0; i<numPoints; i++){
			if( abs(yVals[i]) > 995.){
				saturated = true;
				delete Waveform;
				break;
			}
		}
		delete Waveform;
	}
	return saturated;
}

//takes two inputs, the azimuth angle returned by getSunAzimuthLongLat (once that function has been corrected, and no longer adds the longitude before returning,
//Also takes the station longitude, and a boolean. True for returning a number between -180 and 180, and false for 0 to 360, default to True
//it returns the location of the sun, in the CCW coordinate system were phi = 0 = East
double TransformPredictedAzimuth(double inputAzimuth, double station_long, bool ReRange){
	//double Jonathan_azimuth = inputAzimuth + station_long; //this is what Jonathan used to have in the Trunk version of the code before I fixed it, and ported the longitude correction over here so I could do it correctly. //just here for posterity

	//To prepare the azimuth to be used in analysis, we need to do 3 things
	//Thing 1: Rotate the answer out of the testbed frame, and into the antarctic frame
	//Thing 2 : Reverse the clock-wise-ness (handedness?) of the system from N = 0, E = 90, S = 180, W = 180 (a CW system) to E=0, N=90, W=180, S=270 (a CCW system)
	//Thing 3: Can possibly convert this from the 0->360 to the -180->180 range, but that's optional
	
	////////
	//Thing 1:	
	///////
	//Okay, so the function getSunAzimuthLongLat returns the phi of the sun, relative to a phi=0 that is defined to be the -74 meridian
	//To be clear, this means that phi=0 for this function is the equivalent of phi=-74 for someone standing perfectly at south pole
	//So, we need to correct for that, but do it carefully. The function always returns a number between 0 and 360. 
	//Special cases are:
	// (1) if we add a negative longitude, and this brings our angle below 0, then we need to correct this to be some number that is less than positive 360
		//ie, -10 -> 350
	// (2) if we add a positive longitude, and this brings our angle above 360, then we need to correct this to be some number that is bigger than positive 0
		//ie, 380 -> 20
	//This tranformed coordiante system will have N = 0, E = 90, S = 180, W = 180
	if((inputAzimuth + station_long)<0.) {inputAzimuth+=360.; inputAzimuth +=station_long;}
	else if((inputAzimuth + station_long)>360.){ inputAzimuth-=360.; inputAzimuth-=station_long;}
	else {inputAzimuth +=station_long;}
		
	////////
	//Thing 2:	
	///////
	//Now, we need to reverse the clock-wise-ness of this sytem. Currently, it goes N = 0, E = 90, S = 180, W = 180, which is the geography convention, and is also a CW system
	//To agree with a physics convention, it actually needs to go E=0, N=90, W=180, S=270, which is a CCW system, and is the same convention we calculat the interferometric maps with
	//Let me demonstrate four cases. Call "Frame 1" = N = 0, E = 90, S = 180, W = 180 and "Frame 2" = E=0, N=90, W=180, S=270
	//Case 1: Theta_F1 = 30, then Theta_F2 = 60, can be obtained by doing 90 - (Theta_F2-0) = 90- Theta_F2. This means that if Theta_F1<90, then we should do 90-Theta_F2
	//Case 2: Theta_F1 = 290, then Theta_F2 = 160, can be obtaine dby doing 180- (Theta_F2-270) = 450 - Theta_F2. This means that if Theta_F1 < 360 and Theta_F1 > 270, we should do 180-(Theta_F2 - 270) = 450 - Theta_F2
	//Case 3: Theta_F1 = 210, then Theta_F2 = 240, can be obtaine dby doing 270- (Theta_F2-180) = 450 - Theta_F2. This means that if Theta_F1 < 270 and Theta_F1 > 180, we should do 270-(Theta_F2 - 180) = 450 - Theta_F2
	//Case : Theta_F1 = 120, then Theta_F2 = 330, can be obtaine dby doing 360- (Theta_F2-90) = 450 - Theta_F2. This means that if Theta_F1 < 180 and Theta_F1 > 90, we should do 360-(Theta_F2 - 90) = 450 - Theta_F2

	//So, you could just say that if Theta_F1<90 then do Theta_F2 = 90 - ThetaF1, else 450 - ThetaF1
	//Or, you could equivalently do (http://math.stackexchange.com/questions/926226/conversion-from-azimuth-to-counterclockwise-angle):
		//double new_azimuth = 450 - Theta_F1;
		//if(Theta_F1> 360) Theta_F1-=360;
		//else return Theta_F1;
	//Both should work just fine, but I'll go for the latter for clarity, and becasue someone else suggested it to
		
	double new_azimuth = 450. - inputAzimuth;
	if(new_azimuth > 360.) new_azimuth -=360;
		
		
	//At the end of this whole do-hicky, what's happened is that we've corrected for the longitude offset in the zenith angle, and also fixed the clock-wise-ness of the frame
	//now, Azimuth is measured going CCW from a phi = 0 = E
		
	////////
	//Thing 3:
	///////
	//Now, we need to domain correct it if desired, so switch from a frame where we have the boundary as [0,360] to a system where the domain is [-180,180]
	//This functionally means that if the azimuth is <180, leave it alone, but if it's bigger, do -360 + azimuth
		
	double new_domain_azimuth = new_azimuth;
	if(ReRange){ //this means correct the range from 0->360 to -180->180
		if(new_azimuth > 180.) new_domain_azimuth = -360. + new_azimuth;		
	}
	
	return new_domain_azimuth;//this will only have been changed if ReRange as true
}

double TransformMapPeakToGlobalFrame(double inputAzimuth){
     
	double TBx[3] = {-0.5986468972,0.8010128835,-0.000399302}; //this is from the ARA coordinates document, and contains the easting, northing, and upping components of the TB origin relative to the array (and therefore global) origin
	double angle_from_0 = TMath::ATan2(TBx[1],TBx[0]) * TMath::RadToDeg(); //I can compute this as the angle from zero because ATan2 is intelligent, and knows I mean the angle from zero!
	double outputAzimuth = inputAzimuth + angle_from_0; //correct this
	if(outputAzimuth > 180.) outputAzimuth -=360.; //this accounts for domain correcting the answer (ie, spinning from -180 -> 180, to 0->360
	return outputAzimuth; //return the answer
	
}




//A function that will take a map (inputMap) and shift the bin contents around by a specified theta and phi
TH2D *ShiftedMap ( TH2D *inputMap, int ThetaShift, int PhiShift) {
	TH2D *outputMap = new TH2D("","",360,-180,180,180,-90,90); //make an out put map
	for (int theta=1; theta<=180; theta++) { //loop over theta
		for (int phi=1; phi<=360; phi++) { //loop over phi
			int temp_phi = phi + PhiShift; //make these temporary variables
			int temp_theta= theta + ThetaShift; //make these temporary variables
			CheckRange_byBin(temp_theta,temp_phi);
			double content = inputMap -> GetBinContent(phi,theta);
			outputMap -> SetBinContent(temp_phi, temp_theta, content);
		}
	}
	
	return outputMap;
}

//a function for dealing with the boundary cases of where you have reached the edge of the histogram and need to correclty re-wrap the iterator
//adapted for wrappign by bins (not degrees) from CheckCorrMapRange from Eugene
void CheckRange_byBin( int &theta, int &phi ) {
	
	if(theta > 180) theta -= 180; //if you've incremented over the 90 high mark, wrap all the way around to -90
	else if ( theta < 0 ) theta += 180; //if you've incremented below -91, wrap all the way around to 90
	
	if (phi > 360) phi -= 360; //if you've incremented over the 180 high mark, wrap all the way back around to -180
	else if ( phi < 0 ) phi += 360; //if you've incremented below -180, wrap all the way around to 180
}


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
	    double integral1=0;
	    double integral2=0;
		int num_contributing_points = 0;
		for(int i=length1; i<2*length1; i++){//integrate over just the center region of waveform 1
		//for(int i=length1; i<6002; i++){//integrate over just the center region of waveform 1
			int j=i+delay; //but, integrate over all possible delayed regions of waveform 2 //for that given i step, you want to look at i+delay in the other waveform
			double yVal1 = oldY1[i]; //one of the old voltages
			double yVal2 = oldY2[j]; //the other old voltage
			double product = yVal1*yVal2*deltaT; //the product of the two with deltaT to make the time integral
			if((abs(yVal1)>0)&&(abs(yVal2)>0)){ //only if both are non-zero should you do anything with it
				num_contributing_points++; //increment up the number of contributing points
				sum+=product; //increase the sum
				//integral1+=(yVal1*yVal1*deltaT); //add to the integral for this step
				//integral2+=(yVal2*yVal2*deltaT); //add to the integral for this step
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
			integral1= FFTtools::integrateVoltageSquared(gr1);
			integral2=FFTtools::integrateVoltageSquared(gr2);
			correlationValues.push_back(sum/TMath::Sqrt(integral1)/TMath::Sqrt(integral2));  //according to equation 10 in arXiv 1304.5663v1
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
			double product = yVal1*yVal2*deltaT; //deltaT added 8/5
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

	//double interpolation_step=0.1; //this is not needed unless we are interpolating
	//Now, lots of things have to happen in a loop over the REST of the antennas
	/*
	vector <TGraph*> Waveforms;
	vector <TGraph*> Waveforms_Interpolated;
	vector <TGraph*> Waveforms_Shifted;
	vector <TGraph*> Waveforms_Padded;
	vector <Double_t> time_delays;
	*/
	vector <TGraph*> Waveforms_Cropped;
	double interpolation_step  = 0.5; //1/2 nanosecond interpolation step
						
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
}

// get the correlation map peak location (in 1 degree bins) for the Solar Flare analysis, which needs special conditioning on allowed peaks
void GetCorrMapPeakAngle_1deg_SolarAnalysis( TH2D *theCorrMap_input, int &peaktheta, int &peakphi, double &peakvalue, int eventType) {

	double corr_peak = -1.;
	
	if(eventType==0){ //a flare event
		for (int theta=90; theta<=180; theta++) { //right now, only look at peaks above theta = 20 //actually, just above the horizon for now
			for (int phi=0; phi<=360; phi++) { //loop over phi
				if(   (phi>30) && (phi<330) ){ continue;} //if you're outside of where the sun should be, don't do anything
				//if(2==1);
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

// get peak range area plot and number of bins
double GetTotalArea ( TH2D *theCorrMap, int peaktheta, int peakphi, double peakfrac) {

	double area=0;
	double peakbinvalue = theCorrMap->GetBinContent(peakphi+180,peaktheta+90); //the value of the peak passed to the function
	double peakrange = peakbinvalue * peakfrac;  //a fraction of the peak value that defines the "width" of the peak you're looking for

	for (int theta=0; theta<=180; theta++) {
		for (int phi=0; phi<=360; phi++) {
			if( theCorrMap->GetBinContent(phi+180,theta+90) > peakrange)  area += cos( (double)theta * DEG2RAD );
		}
	}  
    return area;
}
