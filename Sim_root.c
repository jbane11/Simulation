#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;
void Sim_root() {
		time_t start = time(0) ;
	cout << "\n" << "\n";
int run =0;
cout << "Please input the run number of the disturbution used." <<endl;
cin >> run;

  char* data_dir;
  data_dir = getenv ("OUTPUT_DIR");
  char* root_dir;
  root_dir = getenv ("OUT_DIR");


	ifstream in;
	in.open(Form("%s/ISroot_run_%d.txt",data_dir,run));
	Float_t lambda,qsquared,theta,Xb,efinal,F_two,P,Diffcross,deltaErest,deltaElab,DeltaE,Qsquared_Rf,xb_RF,invarM_RF,invarM_LF;
	Int_t nlines = 0;

remove(Form("%s/IS_%d.root",root_dir,run));

	TFile *f = new TFile(Form("%s/IS_%d.root",root_dir,run),"RECREATE");

TH1F *XbW1_1 = new TH1F("Xb_Weighted_1",Form("Counts of Xb weighted with cross section focousing on EMC region with dist %d",run),100,0,1);


	TH1F *XbW = new TH1F("Xb_Weighted","Counts of Xb weighted with cross section",100,0.1,1.0);
	TH1F *XbW1 = new TH1F("Xb_Weighted1","Counts of Xb weighted with cross section",100,0.1,1.0);
	
	TH1F *XbW1_2 = new TH1F("Xb_Weighted_2","Counts of Xb weighted with cross section focousing on EMC region",100,0,1);
	TH1F *XbW2 = new TH1F("Xb_Weighted2","Counts of Xb weighted with cross section at Backward Angle",100,0,1);
	TH1F *XbW3 = new TH1F("Xb_Weighted3","Counts of Xb weighted with cross section",1000,0,2.5);

// Coloring the different histograms.
	XbW1-> SetLineColor(1);
	XbW2-> SetLineColor(2);
	XbW3-> SetLineColor(14);

//Unweighted histograms
	TH1F *xb = new TH1F("Xb","Counts of Xb",100,0.1,1.0);
	TH1F *xb_1 = new TH1F("Xb_1", Form("Counts of Xb focousing on EMC region: dist. %d",run),100,0,1);


	TH1F *xb_2 = new TH1F("Xb_2", "Counts of Xb focousing on EMC region: dist 2",100,0,1);
	TH1F *xb1 = new TH1F("Xb1","Counts of Xb p=0.25",100,0,2.5);
	TH1F *xb2 = new TH1F("Xb2","Counts of Xb p=0.5",100,0,2.5);

//Q^2 histogram
	TH1F *QQ = new TH1F("Qsquared", "Counts in Q^2", 100, 0, 20);
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

//Histograms for deltaE 
	TH1F *deltaEr = new TH1F("deltaEr","Counts of deltaE for the rest frame",100,0,10);
	TH1F *deltaEl = new TH1F("deltaEl","Counts of deltaE for the lab frame",100,0,10);
//deltaE vrs xb
	TH2F *dEXr = new TH2F("deltaEvxbr","Counts of deltaE vrs Xb for the rest frame",100,0,2.5,100,0,10);
	TH2F *dEXl = new TH2F("deltaEvxbl","Counts of deltaE vrs Xb for the lab frame",100,0,2.5,100,0,10);
//	dEXl-> SetMarkerColor(3);;
//historgram for momentum disturbution
	TH1F *Pmom = new TH1F("Pmom","Counts for momentum disturbution",100,0,1.0);

	TH2F *Histo2 = new TH2F("Histo2" ,"cross " ,1000,0,2.5,1000, 0,0.02);
	TH1F *Xbftwo = new TH1F("Xbftwo" ,"Xb counts weighted with Ftwo",100,0,2.5);
// Coloring the different histograms.
//	Xba-> SetLineColor(1);
//	Xbb-> SetLineColor(2);
//	Xbc-> SetLineColor(14);
//	Xbd-> SetLineColor(4);

//histograms for invarient Mass
	TH1F *Wmass  = new TH1F("Wmass" ,"invariant Mass " ,50,0,5);
	TH1F *WmassW = new TH1F("WmassW" ,"invariant Mass weighted with Cross section" ,50,0,5);
	TH2F *WDiff  = new TH2F("WDiff" ,"invariant Mass vs Cross section",50,0,5,50,0.0001,0.02);
//Invariant Mass v Xb
	TH2F *W_xb  = new TH2F("W_xb" ,"invariant Mass vs Xb " ,1000,0,1.2,100,0,5);
//Q^2 V Xb
	TH2F *QQ_xb  = new TH2F("QQ_xb" ,"Q^2 vs Xb " ,100,0,1.2,100,0,20.5);	
	TH2F *QQ_W  = new TH2F("QQ_W" ,"Q^2 vs W " ,100,0,5,100,0,20.5);	

	TH2F *xb_bin[6];
	for(int i=0;i<7;i++){
	xb_bin[i] = new TH2F(Form("QQ_W%d",i) ,"Q^2 vs W " ,100,0,5,100,0,20.5);
	}
xb_bin[0]->SetMarkerColor(1);
xb_bin[1]->SetMarkerColor(kGreen+3);
xb_bin[2]->SetMarkerColor(kBlue+1); //EMC effect
xb_bin[3]->SetMarkerColor(kMagenta+1);
xb_bin[4]->SetMarkerColor(kRed+1);

	TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","lambda:qsquared:theta:xb:efinal:Proton_mom:Wmass");



int k=0;
int j=0;
int i=0;
while (1) { //Xb-> lab frame, xb-> rest fram

 in >>lambda>>qsquared>>theta>>Xb>>efinal>>F_two>>P>>Diffcross>>deltaElab>>deltaErest>>Qsquared_Rf>>xb_RF>>invarM_RF>>invarM_LF;
double jjk = i *100.0/100.0;

//cout <<lambda<< "  "<<qsquared<< "  "<<theta<<"  "<<Xb<< "  "<<efinal<< "  "<< F_two<<"  "<< P <<"  " <<Diffcross<< " " <<deltaElab<< "  "<<deltaErest<< "  " <<Qsquared_Rf<< " "<<xb_RF<< " "<<invarM_RF<<"  " <<invarM_LF<< endl;


//if(jjk/10000 == jjk/10000.0) {cout << jjk << " "<<in.good()<<endl;}


//if(invarM_LF < 2){ continue;} // Only lets DIS events stay.
if (!in.good()) { break;} 
//QQ->Fill(qsquared);
//QQW->Fill(qsquared,Diffcross);

//Cut in Q^2

//if(qsquared < 3.5 || qsquared > 15.5){continue;}
i++;
//if(theta >= 0 && theta <= 360){ XbW->Fill(Xb,Diffcross);}
//if(theta >0 && theta < 90){ XbW1->Fill(Xb,Diffcross);}//Green
//if(theta >-90 && theta < -45){ XbW1->Fill(Xb,Diffcross);}
//if(theta >90 && theta < 270){ XbW2->Fill(Xb,Diffcross);}//Red
//if(theta >270 && theta < 360){ XbW1->Fill(Xb,Diffcross);}

int xbin; //trash catcher

	if(Xb < 0.06){xbin = 0;}
		else if(Xb < 0.3){xbin=1;}
			else if(Xb < 0.71){xbin=2;}
				else if(Xb < 1){xbin=3;}
					else{ xbin=4;}
	

xb_bin[xbin]->Fill(invarM_LF, qsquared);

if(qsquared > 2 && qsquared < 14 && invarM_LF > 2.5 && invarM_LF < 4){
 XbW1_1->Fill(Xb, Diffcross); 
 xb_1->Fill(Xb);  }


//1D
Pmom->Fill(P);
xb->Fill(Xb);
XbW->Fill(Xb, Diffcross);
QQ->Fill(qsquared);
QQW->Fill(qsquared,Diffcross);
Wmass->Fill(invarM_RF);
WmassW->Fill(invarM_RF,Diffcross);

//2D
dEXr->Fill(xb_RF,deltaErest);
dEXl->Fill(Xb,deltaElab);
W_xb->Fill(Xb,invarM_LF);
WDiff->Fill(invarM_LF,Diffcross);
QQ_xb->Fill(Xb, qsquared);
QQ_W->Fill(invarM_LF, qsquared);

//3D


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


	ntuple->Fill(lambda,qsquared,theta,Xb,efinal,P,invarM_LF);
nlines++;

if(i/250000 == i/250000.0){
	if(i >= 1000000){cout << "Event count  " << i/1000000.0<< " "<<"mil" << "\n";}  
	else cout << "Event count  " << i << "\n";	}


}
/*TCanvas *aawa = new TCanvas("asas","A Simple Graph Example",0,0,1200,900);
	QQ_W->Draw();
	xb_bin[0]->Draw("same");
	xb_bin[1]->Draw("same");
	xb_bin[2]->Draw("same");
	xb_bin[3]->Draw("same");
	xb_bin[4]->Draw("same");*/

//TCanvas *aa2wa = new TCanvas("as2as","A Simple Graph Example",0,0,700,700);
//	QQ_W->Draw();

//TCanvas *aaww = new TCanvas("d","A Simple Graph Example",0,0,700,700);
//XbW1_1->Draw();
//aaww->cd();
//aaww->SetLogy();
//aaww->update;

//TCanvas *a = new TCanvas("a","A Simple Graph Example",0,0,700,700);
//WDiff->Draw();
//TCanvas *c = new TCanvas("c","A Simple Graph Example",700,0,700,700);
//dEXl->Draw();

//TCanvas *aaa = new TCanvas("aaa","A Simple Graph Example",0,0,700,700);
//xb->Draw();
//TCanvas *aaa = new TCanvas;
//QQW->Draw();
// TCanvas *aaaa = new TCanvas;
// QQ->Draw();
//TCanvas *d = new TCanvas("b","A Simple Graph Example",0,0,700,700);
//XbW->Draw();

//TCanvas *aaaww = new TCanvas("ad","A Simple Graph Example",0,0,700,700);
//dEXr->Draw();


printf(" found %d points\n",nlines);
	
in.close();
f->Write();

	time_t finish = time(0) ;
 cout << "This program took " << finish-start << " seconds to run"<<endl;
}


//////////////////////////////////JUNK///////////////////////////////////////


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

