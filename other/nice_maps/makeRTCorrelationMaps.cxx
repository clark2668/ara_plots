#include <iostream>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

// ARA Includes
#include "AraGeomTool.h"
#include "RayTraceCorrelator.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
RawAtriStationEvent *rawAtriEvPtr;

//
// ./makeRTCorrelationmaps 2 /data/wipac/ARA/2013/filtered/burnSample1in10/ARA02/root/run1942/event1942.root

void
set_plot_style()
{
     const Int_t NRGBs = 5;
     const Int_t NCont = 255;

     Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
     Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
     Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
     Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
     TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
     gStyle->SetNumberContours(NCont);
}

int main(int argc, char **argv)
{
    double interpV = 0.4;
    double interpH = 0.625;

    if(argc<3) {
        std::cout << "Usage\n" << argv[0] << " <station> <input file>\n";
        std::cout << "e.g.\n" << argv[0] << " 2 http://www.hep.ucl.ac.uk/uhen/ara/monitor/root/run1841/event1841.root\n";
        return 0;
    }

    int station = atoi(argv[1]);

    /////////////////////////////////////////////////
    /////////////////////////////////////////////////
    //// Initialize the correlator
    /////////////////////////////////////////////////
    /////////////////////////////////////////////////

    // setup the paths to our ray tracing tables
    double radius = 300.;
    double angular_size = 1.;
    int iceModel = 0;
    char dirPath[500];
    char refPath[500];
    std::string topDir = "/cvmfs/ara.opensciencegrid.org/data/raytrace_tables/";
    sprintf(dirPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_0.root",
        topDir.c_str(), station, iceModel, radius, angular_size
    );
    sprintf(refPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_1.root",
        topDir.c_str(), station, iceModel, radius, angular_size
    );

    int numAntennas = 16;
    // initialize a correlator
    RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(station, numAntennas,
        radius, angular_size, dirPath, refPath
    );

    // and tell it to load up the arrival times tables
    theCorrelator->LoadTables();

    // How you set up the pairs is up to you!
    // There are a few helper functions;
    // for example, here we can load all of the VPol pairs.

    AraGeomTool *geomTool = AraGeomTool::Instance();
    // std::vector<int> excludedChannels = {1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::vector<int> excludedChannels = {15};
    std::map< int, std::vector<int> > pairs = theCorrelator->SetupPairs(station, geomTool, AraAntPol::kVertical, excludedChannels);
    std::cout<<"Number of pairs "<<pairs.size()<<std::endl;
    for(int i=0; i<pairs.size(); i++){
        printf("Pair %d: %d, %d\n",i,pairs.find(i)->second[0], pairs.find(i)->second[1]);
    }


    printf("------------------\n");


    /////////////////////////////////////////////////
    /////////////////////////////////////////////////
    //// Actually use it on some data
    /////////////////////////////////////////////////
    /////////////////////////////////////////////////

    TFile *fp = TFile::Open(argv[2]);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    eventTree->SetBranchAddress("event", &rawAtriEvPtr);
    Long64_t numEntries=eventTree->GetEntries();

    char ped_file_name[400];
    sprintf(ped_file_name, "/data/user/mkim/OMF_filter/ARA03/ped_full/ped_full_values_A3_R1814.dat");
    // sprintf(ped_file_name, "/data/user/brianclark/ARA/ara5_analysis/peds/2013/A3/reped_run_001814.dat.gz");
    AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
    calibrator->setAtriPedFile(ped_file_name, 3);


    numEntries=10;
    for(Long64_t event=0;event<numEntries;event++) {
        eventTree->GetEntry(event);

        bool isCalpulser = rawAtriEvPtr->isCalpulserEvent();
        if(!isCalpulser) continue;

        std::cout<<"Looking at event number "<<event<<std::endl;

        UsefulAraStationEvent * realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);

        std::map<int, TGraph*> interpolatedWaveforms;
        for(int i=0; i<16; i++){
            TGraph *gr = realAtriEvPtr->getGraphFromRFChan(i);
            TGraph *grInt = FFTtools::getInterpolatedGraph(gr,i<8?interpV:interpH);
            interpolatedWaveforms[i] = grInt;
            delete gr;
        }
        std::vector<TGraph*> corrFunctions = theCorrelator->GetCorrFunctions(pairs, interpolatedWaveforms); // apply Hilbert envelope is default

        // get the map
        TH2D *dirMap = theCorrelator->GetInterferometricMap(pairs, corrFunctions, 0); // direct solution

        gStyle->SetOptStat(0);
        TCanvas *c = new TCanvas("","", 1200, 950);
        dirMap->Draw("z aitoff"); // aitoff projection
        dirMap->GetXaxis()->SetTitle("Phi [deg]");
        dirMap->GetYaxis()->SetTitle("Theta [deg]");

        dirMap->GetYaxis()->SetTitleSize(0.05);
        dirMap->GetYaxis()->SetLabelSize(0.03);
        dirMap->GetYaxis()->SetTitleOffset(0.6);
        
        dirMap->GetXaxis()->SetTitleSize(0.05);
        dirMap->GetXaxis()->SetLabelSize(0.03);
        dirMap->GetXaxis()->SetTitleOffset(0.6);

        gStyle->SetPalette(kRainBow);
        // gStyle->SetNumberContours(100);
        // set_plot_style();

        // dirMap->GetZaxis()->SetTitleSize(0.05);
        // dirMap->GetZaxis()->SetLabelSize(0.03);	
        // dirMap->GetZaxis()->SetTitleOffset(0.95);
        
        dirMap->GetXaxis()->CenterTitle();
        dirMap->GetYaxis()->CenterTitle();
        // dirMap->GetZaxis()->CenterTitle();

        gPad->SetRightMargin(0.15);
        char title[500];
        sprintf(title,"maps_ev%d.png", event);
        c->SaveAs(title);
        delete c;

        // cleanup
        delete dirMap;
        for(int i=0; i<16; i++){
            delete interpolatedWaveforms[i];
        }
        for(int i=0; i<corrFunctions.size(); i++){
            delete corrFunctions[i];
        }
        delete realAtriEvPtr;

    }
}
