#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;
void Ftwo() {
	ifstream in;
	in.open("data/ISrest.txt");
	double x,QQ,F_two,W_two,mew,xb,new;
	Int_t nlines = 0;
	remove("ISrest.root");
	TFile *f = new TFile("/lustre/expphy/volatile/halla/triton/Bane/Rootfiles/ISrest.root","RECREATE");
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",0,0,700,700);
	TH2F *F_2x = new TH2F("F_2vXb", "2D histo with F2 v Xb",100,0,1.5,100,0,.4);
	TH2F *F_2xb = new TH2F("F_2vXbb", "2D histo with F2 v Xb",100,0,1.5,100,0,.4);

	F_2xb-> SetMarkerColor(3);
	TH2F *F_2Q = new TH2F("F_2vQ", "2D histo with F2 v Q^2",100,0,11,100,0,.4);
	TH2F *F_2v = new TH2F("F_2vnu_restframe", "2D histo with F2 v nu",100,0,7,100,0,.4);
	TH2F *F_2vb = new TH2F("F_2vnu_labframe", "2D histo with F2 v nu",100,0,7,100,0,.4);
F_2vb-> SetMarkerColor(3);
	TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","Xb:qsquared:F_two:W_two:mew:xb:new");
while (1) {
in >>x >> QQ >> F_two >> W_two >>mew>>xb>>new;
if (!in.good())break;

F_2x->Fill(x,F_two);
F_2xb->Fill(xb,F_two);
F_2Q->Fill(QQ,F_two);
F_2v->Fill(mew,F_two);
F_2vb->Fill(new,F_two);

ntuple->Fill(x,QQ,F_two,W_two,mew,xb,new);

nlines++;}
	

F_2xb->Draw();
F_2x->Draw("same");
  // TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",700,0,500,500);
//F_2Q->Draw();
   //TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",1200,0,500,500);
//F_2v->Draw();
  // TCanvas *c4 = new TCanvas("c4","A Simple Graph Example",1200,600,500,500);	
//F_2vb->Draw();

printf("%d \n",nlines);


	
in.close();
f->Write();
}
