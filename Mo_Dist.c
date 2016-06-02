#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;
void modist_2() {
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
	in.open(Form("%s/pmom_%d.txt",data_dir,run));

	Float_t P,random;
	
	Int_t nlines = 0;
	

	TFile *f = new TFile(Form("%s/Pmom%d.root",root_dir,run),"recreate");
	
	TH1F *Pmom  = new TH1F(Form("Pmom_%d",run) ,Form("Momentum disturbution for run %d",run),100,0,1.0);
	TH1F *Pmom1 = new TH1F("Pmom2","Counts for momentum disturbution",100,0,0.8);
	TH1F *Pmom2 = new TH1F("Pmom3","Counts for momentum disturbution",100,0,4.5);		
	TH1F *Pmom3 = new TH1F("Pmom4","Counts for momentum disturbution",100,0,4.5);	
	TH1F *Pmom4 = new TH1F("Pmom5","Counts for momentum disturbution",100,0,4.5);
//		Pmom1-> SetLineColor(2);
//		Pmom2-> SetLineColor(3);
//		Pmom3-> SetLineColor(4);
//		Pmom4-> SetLineColor(7);

int j=0;
int i=0;
while (1) {
in >> P>>random;
if (!in.good()) break; 

if(1==1){Pmom -> Fill(P);}
//	else if (random <= 40){Pmom1 -> Fill(P);}
//	else if (random <= 60){Pmom2 -> Fill(P);}
//	else if (random <= 80){Pmom3 -> Fill(P);}
//	else 				   {Pmom4 -> Fill(P);}

i++;
nlines++;

if(i/500000 == i/500000.0){
	if(i >= 1000000){cout << "Event count  " << i/1000000.0<< " "<<"mil" << "\n";}  
	else cout << "Event count  " << i << "\n";	}


}

TCanvas *aaww = new TCanvas("d","A Simple Graph Example",0,0,1000,700);
aaww->cd();
aaww->SetLogy();
Pmom->Draw();
//TH1F *h_sum = Pmom->Add(*Pmom1, 1);
/*
Pmom->Add(Pmom1);
Pmom->Add(Pmom2);
Pmom->Add(Pmom3);
Pmom->Add(Pmom4);
*/

Double_t norm = Pmom->GetEntries();
 norm = norm + Pmom1->GetEntries();
 norm = norm + Pmom2->GetEntries();
 norm = norm + Pmom3->GetEntries();
 norm = norm + Pmom4->GetEntries();
Pmom->Scale(1/norm);





/*
/Pmom->Draw();
Pmom1->Draw("same");
Pmom2->Draw("same");
Pmom3->Draw("same");
Pmom4->Draw("same");
*/
	printf("%d \n",i);

f->Write();	
in.close();


//cout<< "intergral   " << Pmom->Integral()<<endl;

	time_t finish = time(0) ;
 cout << "This program took " << finish-start << " seconds to run"<<endl;
}

