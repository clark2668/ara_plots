/////////////
///// do fit
/////////////

// C++ Includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>

// ROOT includes
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"

using namespace std;

int main(int argc, char **argv)
{




	double mydepths[] = {3,8,12,18,22,28,35,45,55,65,75,85,95,105,115,125,135,145,400,500,600,700,800};
	double myindexs[] = {1.35,1.38,1.45,1.52,1.46,1.54,1.50,1.57,1.53,1.56,1.56,1.63,1.68,1.73,1.72,1.72,1.75,1.77,1.78,1.78,1.78,1.78,1.78};

	vector<double> depths;
	vector<double> indexs;

	for(int i=0; i<23; i++){
		depths.push_back(mydepths[i]);
		indexs.push_back(myindexs[i]);
	}

	// vector<double> depths(mydepths, mydepths+sizeof(mydepths)/sizeof(double));
	// vector<double> indexs(myindexs, myindexs+sizeof(myindexs)/sizeof(double));

	TGraph *g = new TGraph(depths.size(), &depths[0], &indexs[0]);

	char equation[150];
	// sprintf(equation,"[0] + ([1])*exp([2]*-x)");
	sprintf(equation,"[0] - ([0]-[1])*exp([2]*-x)");
	TF1 *fit = new TF1("expFit",equation,0,800);
	fit->SetParameter(0,1.7);
	fit->SetParameter(1,1.3);
	fit->SetParameter(2,0.014);
	g->Fit("expFit","MER");

	TCanvas *c = new TCanvas("","",1100,850);
	
	char graph_title[200];
	sprintf(graph_title,"n_{d}=%.3f #\pm %.3f, n_{s}=%.3f #\pm %.3f, n_{c}=%.4f #\pm %.4f",
		fit->GetParameter(0),fit->GetParError(0),
		fit->GetParameter(1),fit->GetParError(1),
		fit->GetParameter(2),fit->GetParError(2)
	);
	g->SetTitle(graph_title);
	g->Draw("AP");
	g->SetMarkerStyle(kFullCircle);	
	g->GetYaxis()->SetRangeUser(1.35,1.8);

	g->GetXaxis()->SetTitle("Depth [m]");
	g->GetYaxis()->SetTitle("Index of Refraction");


	c->SaveAs("index_vs_depth.png");

	// gStyle->SetLineWidth(4); //set some style parameters
	// TFile *stripe1file = new TFile("hstripe1snr.root"); //import the file containing the histogram to be fit
	// TH1D *hstripe1 = (TH1D*)stripe1file->Get("stripe1snr"); //get the histogram to be fit
	// TCanvas c1("c1","c1",1000,800); //make a canvas
	// hstripe1->Draw(""); //draw it
	// c1.SetLogy(); //set a log axis

	// //need to declare an equation
	// //I want to fit for two paramters, in the equation these are [0] and [1]
	// //so, you'll need to re-write the equation to whatever you're trying to fit for
	// //but ROOT wants the variables to find to be given as [0], [1], [2], etc.
	
	// char equation1[150]; //declare a container for the equation
	// sprintf(equation1,"([0]*(x^[1]))"); //declare the equation
	// TF1 *fit1 = new TF1("PowerFit",equation1,20,50); //create a function to fit with, with the range being 20 to 50

	// //now, we need to set the initial parameters of the fit
	// //fit->SetParameter(0,H->GetRMS()); //this should be a good starting place for a standard deviation like variable
	// //fit->SetParameter(1,H->GetMaximum()); //this should be a good starting place for amplitude like variable
	// fit1->SetParameter(0,60000.); //for our example, we will manually choose this
	// fit1->SetParameter(1,-3.);

	// hstripe1->Fit("PowerFit","R"); //actually do the fit;
	// fit1->Draw("same"); //draw the fit

	// //now, we want to print out some parameters to see how good the fit was
	// cout << "par0 " << fit1->GetParameter(0) << " par1 " << fit1->GetParameter(1) << endl;
	// cout<<"chisquare "<<fit1->GetChisquare()<<endl;
	// cout<<"Ndf "<<fit1->GetNDF()<<endl;
	// cout<<"reduced chisquare "<<double(fit1->GetChisquare())/double(fit1->GetNDF())<<endl;
	// cout<<"   "<<endl;
	// c1.SaveAs("stripe1snrfit.png");

}