#include <iomanip>

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
const int steps =20000;
const int normalizer=1000;
using namespace std;

void av18(){
  char* data_dir;
  data_dir = getenv ("OUTPUT_DIR");

cout <<endl<<endl;

ofstream function;
function.open("function.tex");

double 	min = 0;
double 	max = 5;
//int   steps = 10000;
double stepsize = (max-min)/steps;
double B[steps];
double C[steps];
double D[steps];
double S_Rule[steps];
double CDF[steps];
double org[steps];
double org1[steps];
B[0]=1.0;  
C[0]=0.004;  /// coiefeiant for second half orginally 0.0107
D[0]=1.0;       /// coiefeiciant for first half
//B[0]=C[0]+D[0];
double F=0.5;

double x[steps];
double sum= 0;
double intersection[2];
double sum_check; 
int run;
cout << "Please input the run number" <<endl;
cin >> run;
cout << "Please input the maxium number for the second exp. The orginal value is 0.0107." <<endl;
cin >> C[0];
cout << "Please input the slope modifier, between 0 and 1. The orginal value is 1." <<endl;
cin >> F;
cout <<endl<<endl;
TGraph *av[normalizer];
 // TCanvas *c_av = new TCanvas("c_av", "Normalization",0,0,1000,900);;
//c_av->cd();
//c_av->SetLogy();
//c_av->Update();

int quit =0;
while(1){

	double sumrule =0;
	double sum= 0;
	int j =0;

	for(int i=1;i<=steps;i++){
	 	x[i] = i*stepsize;
		org[i] = 0.0107*exp(-1.63*x[i]);
		org1[i] = 1*exp(-5.52*x[i]);
		C[i] = C[0]*exp(F*-1.63*x[i]);		
		D[i] = D[0]*exp(-5.52*x[i]);
		B[i] = (C[0]*exp(F*-1.63*x[i])+  D[0]*exp(-5.52*x[i]));
		function << x[i] << " " << B[i]<<endl;
		sum= sum + stepsize*((B[i-1]+B[i])/2.0);
		CDF[i] =sum;
//		cout <<C[i] << " "<< D[i] <<" "<< (C[i] + D[i])/D[i] << " " <<stepsize<< endl;
		if(abs((C[i] + D[i])/D[i] -2) <= stepsize){j++; intersection[0] = x[i]; intersection[1]=D[i];}
		

		sumrule = sumrule + stepsize*(x[i]*x[i]*B[i] +x[i-1]*x[i-1]*B[i-1])/2;

		S_Rule[i]=sumrule;
}

double diff= sumrule - 0.01677;



av[quit] = new TGraph(steps,x,D);
TGraph *avc = new TGraph(steps,x,C);
avc->SetLineColor(2);
//if(quit == 0){avc->Draw();}
//if(quit/10 == quit/10.0){

int color = (quit/normalizer)*19 - 9; 
 
av[quit]->SetLineColor(kViolet+color);
//av[quit]->Draw("same");
//av[quit]->SetMaximum(2);}
//cout << sumrule <<"   "<< D[0] <<endl;

//		cout<<setprecision(4) <<"intergral of P(k) = "<<sum<< ":  intergral of k^2 * P(k) = "<< sumrule<<  ": "<< endl << "intersections for k = " << intersection[0] << ", intersection for P(k)  "<< intersection[1]<<endl<<"  number of tries to find the intersection   "<< j<< "  D[0] = " <<D[0] <<endl<<endl<<endl;

if(quit == normalizer){break;}
if(abs(sumrule - 0.01677) <= 0.0005){cout<< "We had to edit B(1) "<< quit<<" times to normalize." <<endl;break;}
if(sum = sum_check){cout<< "not working" <<endl ;}

if(D[0] <= C[0] + diff){cout << "This function can not be normalized with the orginal by decreasing the y intercept, we have tried "<< quit<< " times" << endl<< endl<< endl; cout << "Warning: This is not normalized to 0.01677. The current value of the int. is "<< sumrule <<endl; break;}
 sum_check = sum;
D[0] = D[0]- diff;
quit++;

}

//cout << sum <<"   "<< D[0] <<endl;
cout<<setprecision(4) <<"intergral of P(k) = "<<sum<< ":  intergral of k^2 * P(k) = "<< sumrule<<  ": "<< endl << "intersections for k = " << intersection[0] << ", intersection for P(k)  "<< intersection[1]<<endl<<endl<<endl<<endl;
TGraph *av18c = new TGraph(steps,x,C);
TGraph *av18d = new TGraph(steps,x,D);
TGraph *av18 = new TGraph(steps,x,B);
TGraph *av18_1 = new TGraph(steps,B,x);
TGraph *av18_11 = new TGraph(steps,B,x);
TGraph *av18_12 = new TGraph(steps,B,x);
TGraph *av18c_1 = new TGraph(steps,C,x);
TGraph *av18d_1 = new TGraph(steps,D,x);
TGraph *orig 	= new TGraph(steps,x,org);
TGraph *orig1 	= new TGraph(steps,x,org1);
orig->SetLineStyle(2);
orig1->SetLineStyle(2);
av18c->SetLineColor(9);
av18c_1->SetLineColor(9);
av18d->SetLineColor(2);
av18d_1->SetLineColor(2);



  TCanvas *c = new TCanvas("c", "orignal function",0,0,700,700);
av18->Draw();
av18c->Draw("same");
av18d->Draw("same");
orig->Draw("same");
orig1->Draw("same");
c->cd();
c->SetLogy();
c->Update();

/*
 TCanvas *CDFsum = new TCanvas("s", "Intregal of P(k)",0,0,700,700);
TGraph *CDF_P = new TGraph(steps,CDF,x);
CDF_P->Draw();
CDFsum->cd();
CDFsum->SetLogx();
CDFsum->Update();
*/




/* TCanvas *c_1 = new TCanvas("c_1", "orignal function projected onto the x",0,0,700,700);
av18_1->Draw();
av18c_1->Draw("same");
av18d_1->Draw("same");
c_1->SetLogx();
c_1->Update();
*/
TF1 *frac = new TF1("frac" ,"([0] + [1]*x)/(1 + [2]*x + [3]*x^2)",0,5);
TF1 *frac1 = new TF1("frac1" ,"([0] + [1]*x)/(1 + [2]*x + [3]*x^2)",0,5);
TF1 *frac2 = new TF1("frac2" ,"([0] + [1]*x)/(1 + [2]*x + [3]*x^2) + ([4] + [5]*x)/(1 + [6]*x + [7]*x^2)",0,5);
TF1 *log   = new TF1("log" ,"log(x/[0])/[1]",0,5);

 TCanvas *Fit = new TCanvas("Fit", "Fit",700,0,700,700);
av18_11 ->Fit("frac","R","",0,intersection[1]);
av18_12 ->Fit("frac1","R+","",intersection[1],1/0.187721);
av18_11->Draw();
av18_12->Draw("same");	
Fit->SetLogx();
Fit->Update();


function.close();

Double_t param[7];
Double_t params[4];

frac->GetParameters(&param[0]);
frac1->GetParameters(&params[0]);

for(int i=0;i<=3;i++){cout <<param[i]<< ",  ";}
for(int i=0;i<=3;i++){cout <<params[i]<< ",  ";}
cout<< "intersection for P(k)  "<< intersection[1];
cout <<endl<<endl;



	ofstream dist;
	dist.open(Form("%s/av18_data_%d.txt",data_dir,run));
for(int i=0;i<=3;i++){dist <<param[i]<< "  ";}
for(int i=0;i<=3;i++){dist <<params[i]<< "  ";}
dist << intersection[1]<< " " << D[0];

dist.close();



}








