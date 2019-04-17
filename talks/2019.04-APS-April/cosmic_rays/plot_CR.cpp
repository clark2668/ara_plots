#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char **argv)
{
	TFile *f = TFile::Open("uzair_CR.root", "READ");
	TGraph *graphs_in[16];
	// std::vector<TGraph*> graphs;
	for(int i=0; i<16; i++){
		char this_graph_title[20];
		sprintf(this_graph_title,"gr%d",i);
		// TGraph *g;
		f->GetObject(this_graph_title,graphs_in[i]);
		// graphs.push_back((TGraph*)g->Clone());
	}

	vector<TGraph*> graphs;
	for(int i=0; i<16; i++){
		vector <double> thisX, thisY;
		for(int samp=0; samp<graphs_in[i]->GetN(); samp++){
			thisX.push_back(graphs_in[i]->GetX()[samp]);
			thisY.push_back(graphs_in[i]->GetY()[samp]);
		}
		graphs.push_back(new TGraph(graphs_in[i]->GetN(), &thisX[0], &thisY[0]));
	}


	gStyle->SetOptStat(0);
	// vector<TH2D*> dummy;
	// for(int i=0; i<16; i++){
	// 	dummy.push_back(new TH2D("","",2,-50,250,2,-100,100));

	// }

	// TCanvas *c = new TCanvas("","",2*1500,2*500);
	// c->Divide(4,2);
	// for(int i=8; i<16; i++){
	// 	c->cd(i-8+1);
	// 	gPad->SetTopMargin(0.06);
	// 	gPad->SetRightMargin(0.06);
	// 	gPad->SetLeftMargin(0.20);
	// 	gPad->SetBottomMargin(0.20);
	// 	dummy[i]->Draw("");
	// 	graphs[i]->Draw("Lsame");
	// 	dummy[i]->GetXaxis()->SetTitle("Time (ns)");
	// 	dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
	// 	dummy[i]->GetXaxis()->SetLabelSize(0.10);
	// 	dummy[i]->GetXaxis()->SetTitleSize(0.10);
	// 	dummy[i]->GetYaxis()->SetLabelSize(0.10);
	// 	dummy[i]->GetYaxis()->SetTitleSize(0.10);
	// 	dummy[i]->GetYaxis()->SetLabelOffset(0.03);
	// 	dummy[i]->GetYaxis()->SetNdivisions(4,0,0,false);
	// 	dummy[i]->GetXaxis()->SetNdivisions(3,0,0,false);
	// }
	// c->SaveAs("cr.pdf");

	vector<TGraph*> graphs_new;
	graphs_new.push_back(graphs[8]);
	graphs_new.push_back(graphs[9]);
	graphs_new.push_back(graphs[0]);
	graphs_new.push_back(graphs[1]);
	graphs_new.push_back(graphs[12]);
	graphs_new.push_back(graphs[13]);
	graphs_new.push_back(graphs[4]);
	graphs_new.push_back(graphs[5]);

	vector<TH2D*> dummy;
	for(int i=0; i<8; i++){
		if(i==2 || i==3 || i==6 || i==7){
			dummy.push_back(new TH2D("","",2,-50,250,2,-5,5));
		}
		else 
			dummy.push_back(new TH2D("","",2,-50,250,2,-100,100));
	}

	TCanvas *c = new TCanvas("","",3*500,4*500);
	c->Divide(2,4);
	for(int i=0; i<8; i++){
		c->cd(i+1);
		gPad->SetTopMargin(0.06);
		gPad->SetRightMargin(0.06);
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.20);
		dummy[i]->Draw("");
		graphs_new[i]->Draw("Lsame");
		dummy[i]->GetXaxis()->SetTitle("Time (ns)");
		dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
		dummy[i]->GetXaxis()->SetLabelSize(0.10);
		dummy[i]->GetXaxis()->SetTitleSize(0.10);
		dummy[i]->GetYaxis()->SetLabelSize(0.10);
		dummy[i]->GetYaxis()->SetTitleSize(0.10);
		dummy[i]->GetYaxis()->SetLabelOffset(0.03);
		dummy[i]->GetYaxis()->SetTitleOffset(1.);
		dummy[i]->GetYaxis()->SetNdivisions(4,0,0,false);
		dummy[i]->GetXaxis()->SetNdivisions(3,0,0,false);
		if(i==2 || i==3 || i==6 || i==7)
			graphs_new[i]->SetLineColor(kRed+1);
		else
			graphs_new[i]->SetLineColor(kBlue+1);
		dummy[i]->GetYaxis()->CenterTitle();
		// if(i!=6 || i!=7)
		// 	dummy[i]->GetXaxis()->SetLabelSize(0);

	}
	c->SaveAs("cr.pdf");
}