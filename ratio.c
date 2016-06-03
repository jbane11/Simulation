#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;
void ratio() {
		time_t start = time(0) ;


	cout << "\n" << "\n";
int run1 =0;
int run2 =0;
cout << "Please input the run numbers of the disturbution you like to look at the ratios of (1st/2nd)." <<endl;
cin >> run1 ;
cin >> run2 ;
  char* data_dir;
  data_dir = getenv ("OUTPUT_DIR");
  char* root_dir;
  root_dir = getenv ("OUT_DIR");

 TFile *file_1 = new TFile(Form("%s/IS_%d.root",root_dir,run1));
 TFile *file_2 = new TFile(Form("%s/IS_%d.root",root_dir,run2));

TH1F *XbW1 = file_1->Get("Xb_Weighted_1");
XbW1->SetTitle(Form("Weighted Counts for runs %d & %d",run1,run2));


//TCanvas *b= new TCanvas("b","B",0,500,500,500);

TH1F *XbW_2  = file_2->Get("Xb_Weighted_1");
XbW_2->SetLineColor(2);


TCanvas *c= new TCanvas("c","C",0,0,900,700);

TH1F * rXbW = new TH1F("Ratio_XbW","The ratio of weighted counts in X",100,0,1);

rXbW =  (TH1F*)XbW1->Clone();
rXbW->SetTitle(Form("Ratio of Weighted counts in Xb for runs %d/%d",run1,run2));
rXbW->Divide(XbW_2);
rXbW->Draw();

TH1F *Xb1 = file_1->Get("Xb_1");
TH1F *Xb2 = file_2->Get("Xb_1");
TH1F * rXb = new TH1F("Ratio_Xb","The ratio of weighted counts in X",100,0,1);

TCanvas *c2= new TCanvas("c2","C2",1100,0,900,700);
rXb = (TH1F*)Xb1->Clone();

rXb->Divide(Xb2);
rXb->SetTitle(Form("Ratio of counts in Xb for runs %d/%d",run1,run2));
rXb->Draw();

//TCanvas *D= new TCanvas("D","D",1000,0,500,500);
//Xb1->Draw();

//TCanvas *F= new TCanvas("F","F",1000,700,500,500);
//Xb2->Draw();
TH1F *XbWa = file_1->Get("Xb_Weighted");
TH1F *XbWb = file_2->Get("Xb_Weighted");
TH1F * rXbWuncut = new TH1F("Ratio_XbW_uncut","The ratio of weighted counts in X",100,0,1);
rXbWuncut = (TH1F*)XbWa->Clone();

rXbWuncut->Divide(XbWb);
//rXbWuncut->Draw();
TCanvas *a= new TCanvas("a","A",600,0,500,500);
XbW1->Draw();
XbW_2 ->Draw("same");
}
