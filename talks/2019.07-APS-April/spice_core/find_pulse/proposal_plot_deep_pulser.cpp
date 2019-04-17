////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////		proposal_plot_deep_puler.cxx
////		Plot deep pulser event
////
////		Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <fstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraEventCalibrator.h"
#include "FFTtools.h"

//ROOT Includes
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<6){
		cout<<"Usage: "<<argv[0]<<" <station> <ant1> <event> <ped file> <run file 1> "<<endl;
		return 0;
	}
	gStyle->SetOptStat(1111111);

	int stationId = atoi(argv[1]);
	int ant1 = atoi(argv[2]);
	int event = atoi(argv[3]);
  
	char pedFileName[200];
	sprintf(pedFileName, "%s", argv[4]);
  
	printf("------------------------------------------------------------------------\n");
	printf("%s \n", argv[0]);
	printf("ant1: %d \n",ant1);
	printf("event: %d \n",event);
	printf("pedFileName %s\n", pedFileName);
	for(int i=5; i<argc; i++) printf("runFileName %s\n", argv[i]);
	printf("------------------------------------------------------------------------\n");

	AraEventCalibrator *calib = AraEventCalibrator::Instance(); //make a calibrator
	calib->setAtriPedFile(pedFileName,stationId); //set the pedestal file
	
	TFile *fp = TFile::Open(argv[5]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree = (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}	
		
	RawAtriStationEvent *rawAtriEvPtr = 0; //empty pointer
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	
	int numEntries=eventTree->GetEntries();
	
	printf("stationId %i, numEntries %i\n", stationId,numEntries);

	eventTree->GetEntry(event);
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib); //full calib

	TGraph *wvfm = realAtriEvPtr->getGraphFromRFChan(ant1);
	delete realAtriEvPtr;

	fp->cd();
	fp->Close();

	//now make pretty pdf version
	gStyle->SetOptStat(111111111);
	TCanvas *canvas = new TCanvas("","",3*500,3*500);
		wvfm->Draw("AL");
		wvfm->GetXaxis()->SetTitle("Time (ns)");
		wvfm->GetYaxis()->SetTitle("Voltage (mV)");
		// wvfm->GetXaxis()->SetRangeUser(50.,500.);
		wvfm->GetYaxis()->SetRangeUser(-1900.,1900.);
		wvfm->SetLineWidth(4);
		wvfm->SetTitle("");
		wvfm->GetYaxis()->SetTitleSize(0.05);
		wvfm->GetXaxis()->SetTitleSize(0.05);
		wvfm->GetYaxis()->SetLabelSize(0.04);
		wvfm->GetXaxis()->SetLabelSize(0.04);
		wvfm->GetYaxis()->SetTitleOffset(1.4);

	gPad->SetRightMargin(0.035);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.15);

	char title[200];
	sprintf(title,"deep_pulser_ant%d_event%d.png",ant1,event);
	canvas->SaveAs(title);
	delete canvas;

	ofstream myfile;
  	myfile.open ("waveform.txt");
  	myfile << "time, volts "<<endl;
  	for(int i=0; i<wvfm->GetN(); i++){
  		myfile<<wvfm->GetX()[i]<<" , "<<wvfm->GetY()[i]<<endl;
  	}
  	myfile.close();

	delete wvfm;
	
}
