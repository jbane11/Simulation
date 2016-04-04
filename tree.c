#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;
void tree() {
	ifstream in;
	in.open("data/ESroot.txt");
	Float_t lambda,qsquared,theta,Xb,efinal,Ge,P,Diffcross,deltaErest,deltaElab;
	Int_t nlines = 0;
	remove("/lustre/expphy/volatile/halla/triton/Bane/Rootfiles/basic.root");
	TFile *f = new TFile("/lustre/expphy/volatile/halla/triton/Bane/Rootfiles/basic.root","RECREATE");
// Histograms for different angles.
/*
	TH1F *XbW = new TH1F("Xb_Weighted","Counts of Xb weighted with cross section",100,0,2.5);
TH1F *XbW1 = new TH1F("Xb_Weighted1","Counts of Xb weighted with cross section at Forward Angle",100,0,2.5);
	TH1F *XbW2 = new TH1F("Xb_Weighted2","Counts of Xb weighted with cross section at Backward Angle",100,0,2.5);
	TH1F *XbW3 = new TH1F("Xb_Weighted3","Counts of Xb weighted with cross section",1000,0,2.5);
// Coloring the different histograms.
	XbW1-> SetLineColor(1);
	XbW2-> SetLineColor(2);
	XbW3-> SetLineColor(14);

//Unweighted histograms
	TH1F *xb = new TH1F("Xb","Counts of Xb",100,0,2.5);
	TH1F *xb1 = new TH1F("Xb1","Counts of Xb p=0.25",100,0,2.5);
	TH1F *xb2 = new TH1F("Xb2","Counts of Xb p=0.5",100,0,2.5);

//Q^2 histogram
	TH1F *QQ = new TH1F("Qsquared", "Counts in Q^2", 100, 0, 15);
	TH1F *QQW= new TH1F("Qsquared_weighted", "Counts in Q^2 weighted with Sigma", 100, 0, 15);	

//Histograms in Xb for different values of Q^2, weighted
	TH1F *XbWa = new TH1F("Xb_Weighted_a","Counts of Xb weighted with cross section with Q^2 central value of 11.5",100,0,2.5);
	TH1F *XbWb = new TH1F("Xb_Weighted_b","Counts of Xb weighted with cross section with Q^2 central value of 12.5",100,0,2.5);
	TH1F *XbWc = new TH1F("Xb_Weighted_c","Counts of Xb weighted with cross section with Q^2 central value of 13.5",100,0,2.5);
	TH1F *XbWd = new TH1F("Xb_Weighted_d","Counts of Xb weighted with with Q^2 central value of 14.5",100,0,2.5);
// Coloring the different histograms.
	XbWa-> SetLineColor(1);
	XbWb-> SetLineColor(2);
	XbWc-> SetLineColor(14);
	XbWd-> SetLineColor(4);
//Histograms in Xb for different values of Q^2, unweighted
	TH1F *Xba = new TH1F("Xb_a","Counts of Xb with Q^2 central value of 5",100,0,2.5);
	TH1F *Xbb = new TH1F("Xb_b","Counts of Xb with Q^2 central value of 7",100,0,2.5);
	TH1F *Xbc = new TH1F("Xb_c","Counts of Xb with Q^2 central value of 9",100,0,2.5);
	TH1F *Xbd = new TH1F("Xb_d","Counts of Xb with Q^2 central value of 11",100,0,2.5);
*/
//Histograms for deltaE 
	TH1F *deltaEr = new TH1F("deltaEr","Counts of deltaE for the rest frame",100,0,10);
	TH1F *deltaEl = new TH1F("deltaEl","Counts of deltaE for the lab frame",100,0,10);



// Coloring the different histograms.
//	Xba-> SetLineColor(1);
//	Xbb-> SetLineColor(2);
//	Xbc-> SetLineColor(14);
//	Xbd-> SetLineColor(4);


	TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","lambda:qsquared:theta:xb:efinal:p:deltaErest:deltaElab");




int j=0;
int i=0;
while (1) {
in >>lambda>>qsquared>>theta>>Xb>>efinal>>Ge>>P>>Diffcross>>deltaElab>>deltaErest;
if (!in.good())break;
//QQ->Fill(qsquared);
//QQW->Fill(qsquared,Diffcross);

//Cut in Q^2


i++;
//if(theta >= 0 && theta <= 360){ XbW->Fill(Xb,Diffcross);}
//if(theta >0 && theta < 90){ XbW1->Fill(Xb,Diffcross);}//Green
//if(theta >-90 && theta < -45){ XbW1->Fill(Xb,Diffcross);}
//if(theta >90 && theta < 270){ XbW2->Fill(Xb,Diffcross);}//Red
//if(theta >270 && theta < 360){ XbW1->Fill(Xb,Diffcross);}



//xb->Fill(Xb);
//XbW->Fill(Xb, Diffcross);


//Q^2 central value of 2.
//if(qsquared >= 11 && qsquared <=12){ XbWa->Fill(Xb,Diffcross); Xba->Fill(Xb);}
//Q^2 central value of 4.
//if(qsquared >= 12 && qsquared <=13){ XbWb->Fill(Xb,Diffcross); Xbb->Fill(Xb);}
//Q^2 central value of 5.
//if(qsquared >= 13 && qsquared <=14){ XbWc->Fill(Xb,Diffcross); Xbc->Fill(Xb);}
//Q^2 central value of 6.
//if(qsquared >= 14 && qsquared <=15){ XbWd->Fill(Xb,Diffcross); Xbd->Fill(Xb);}


deltaEl -> Fill(deltaElab);
deltaEr -> Fill(deltaErest);

	ntuple->Fill(lambda,qsquared,theta,Xb,efinal,Ge,P,Diffcross,deltaErest,deltaElab);
nlines++;
}
deltaEr->Draw();
TCanvas *c = new TCanvas;
deltaEl->Draw();
//xb->Draw();
 // TCanvas *aaa = new TCanvas;
/*QQW->Draw();
  TCanvas *aaaa = new TCanvas;
QQ->Draw();*/
 // TCanvas *c = new TCanvas;
//XbW->Draw();

printf(" found %d points\n",nlines);

// Now I begin making the 2*2 canvas where the different angle cuts  get displayed together.
//	   TCanvas *c = new TCanvas;
//		c->Divide (2,2);

/*	c->cd(1);
	TPad *F = new TPad;
	F->SetLogy();
	F->Draw();
	c->Update();
	F->cd();
	XbW->Draw();

	c->cd(2);
	TPad *Try = new TPad;
	Try->SetLogy();
	Try->Draw(); 
	c->Update();
	Try->cd();
	XbW->Draw();
	XbW1->Draw("SAME");
	XbW2->Draw("SAME");

	c->cd(3);
	TPad *H = new TPad;
	H->SetLogy();
	H->Draw();
	c->Update();
	H->cd();
	XbW1->Draw();
	
	c->cd(4);
	TPad *G = new TPad;
	G->SetLogy();
	G->Draw();
	c->Update();
	G->cd();
	XbW2->Draw();
*/
//
/*
	TImage *img =TImage::Create();
	img->FromPad(b);
	img->WriteImage("images/Xb_Qcut_2.png");
	delete img;*/
/*
TCanvas *a= new TCanvas;
Xba->Draw();

TCanvas *b= new TCanvas;
Xbb->Draw();

TCanvas *d= new TCanvas;
Xbc->Draw();

TCanvas *e= new TCanvas;
Xbd->Draw();

TCanvas *ab= new TCanvas;
XbWa->Draw();

TCanvas *bq= new TCanvas;
XbWb->Draw();

TCanvas *da= new TCanvas;
XbWc->Draw();

TCanvas *ea= new TCanvas;
XbWd->Draw();
*/


/*TCanvas *quad= new TCanvas("ALL Q^2","Xb Counts for 4 different cuts in Q^2");
Xba->Draw();
Xbb->Draw("SAME");
Xbc->Draw("SAME");
Xbd->Draw("SAME");
*/

	

	printf("%d \n",i);

	
in.close();
f->Write();
}
