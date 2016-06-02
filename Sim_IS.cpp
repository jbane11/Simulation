#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
//#include <Tstring.h> 
#include <time.h>
#include "nr3.h"
#include "ran.h"
//#include "/u/home/jbane/headers/trapz-simp.h"
using namespace std;



const 	int number_of_electrons =10; 
//Function protype.
double distribution(double A ,double Fit[8],double intersect, double max);
double rotationZ(double A[4], double angle,int i);
void printvectors(double A[4], double B[4],int f);
double Ftwo(double Beam, double Eprime, double theta);
double Boostx(double A[4], double gamma, double beta, int i);
double Random(double A);
double RANDOM(double A);

int main(void){

	time_t start = time(0) ;

	cout << "\n" << "\n";
int run =0;
  char* data_dir;
  data_dir = getenv ("OUTPUT_DIR");
 



cout << "Please input the run number of the momentum disturbution you would like to use." <<endl;
cin >> run;

char input_file[100];

int m = sprintf(input_file, "%s/av18_data_%d.txt",data_dir,run);

	cout << "\n" << "\n";
	ifstream av18;
	av18.open(input_file);


char output_file[100];

int n = sprintf(output_file, "%s/ISroot_run_%d.txt",data_dir,run);
char pmom_file[100];
int nm = sprintf(pmom_file, "%s/pmom_%d.txt",data_dir,run);
	ofstream root;
	root.open(output_file);
//	ofstream rest;
//	rest.open("/lustre/expphy/volatile/halla/triton/Bane/data/ISrest.txt");
	ofstream Pmom;
	Pmom.open(pmom_file);
// Constants


	double Fit[8];
	double intersect;
	double random_max;
    double mp = 0.989 ;            // Mass of Proton in GeV/c^2
	double me = .0005 ;            // Mass of Electron in GeC/c^2
	double c  = 3.0 * pow(10,8);   // Speed of light in m/s
	double pi = 3.14159265359;     // Pi
	double rad = pi/180;           // Convert from degrees to radians
    double e[4];                   //Electron vector array
    double p[4];                   //Proton vector array
    double e_prime[4];             //ELectron vector array in prime frame
    double p_prime[4];             //Proton vector array in prime frame
    double e_prime2[4];            //ELectron vector array in prime frame
    double p_prime2[4];            //Proton vector array in prime frame
	double e_final[4],p_final[4];	//Final vectors for the electron and proton
    double elec_Beam_momentum = 10;//Beam momentum in GeV
    double Prot_momentum = 0.250     ;  //Proton momentum in GeV   
//	int number_of_electrons =3000000;  
	double min_angle = 0;// By defaults scattered angle covers 360 degrees
	double angle_range = 360;  
	char check;	
	int N = 0;
	int M = 0;
	int print = 0; //If you want to check the 4 vectors during the rotations set this to 10;
	if(print == 10){cout << "Electron 4 vector" << "\t"<<  "Proton 4 vector " <<endl;}
	

//input for momentum disturbution
	av18 >> Fit[0] >> Fit[1] >> Fit[2] >> Fit[3] >> Fit[4] >> Fit[5] >> Fit[6] >> Fit[7] >> intersect >>random_max;


cout << Fit[0] << " " << intersect <<endl;


	int jj = 0;
	int i;
	cout << "\n" << "\n";
// Loop for the number of electrons
	for(i=1; i <= number_of_electrons; i++){

// Random proton momentum from a disturbution	
	Prot_momentum = distribution(i, Fit,intersect,random_max)*0.1973;
	double random_Pmom = Random(i);

//					if(random_Pmom <= 0.0012){	jj++; 
//				cout << jj << " " << i << " " << random_Pmom << " " << Prot_momentum<<endl;}


// Random angle in degrees used for Lambda the incoming proton angle
	srand(i);
	double random_lam = rand();
			random_lam = random_lam/RAND_MAX;
	
	double Lambda = (random_lam)*360*rad;
//	if(i <= 100) {cout<<Prot_momentum << " "<<Lambda/rad<< endl ;}
// incomeing electron and proton 4 vectors.
    e[0] = elec_Beam_momentum;      //electron energy
    e[1] = elec_Beam_momentum;      //Electron Beam Momentum
    e[2] = 0 ;
    e[3] = 0 ;
    p[1] = Prot_momentum*cos(Lambda);   //Proton x momentum
    p[2] = Prot_momentum*sin(Lambda);   //Proton y momentum
    p[3] = 0;
    p[0] = sqrt( mp*mp+Prot_momentum*Prot_momentum); //Proton energy
if(print == 10){cout << setw(25)<< "Electron 4 vector" <<"\t" << "\t"<<  "Proton 4 vector " <<endl;}
printvectors(e, p, print);

// Rotating by Lambda,  clockwise!! This sets Py to 0.
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e , -Lambda, k);
    p_prime[k] = rotationZ(p , -Lambda, k); }
	printvectors(e_prime, p_prime, print);

// Lorentz factors!
	double gamma = p[0]/mp;
	double beta  = p_prime[1]/p[0];

//Lorentz boost in the x direction.
for(int k = 0; k < 4; k++){ 
    e_prime2[k] = Boostx(e_prime , gamma, beta, k);
    p_prime2[k] = Boostx(p_prime , gamma, beta, k); }
	printvectors(e_prime2, p_prime2, print);
// Calculate new angle Delta, angle between e and the x axis. 
    double Delta = atan(e_prime2[2]/e_prime2[1]);

// Rotate the electron to the x axis. 
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , -Delta, k);
    p_prime[k] = rotationZ(p_prime2 , -Delta, k); }
	printvectors(e_prime, p_prime, print);

	double eEri = e_prime[0]; //Electron energy in rest frame before

// Due to limitations of atan, forceing the e[1](x momentum) to be postive. This makes the beam always head to the right, in my internal drawling. This will also let me know if we need to flip back.
	int flip =0;
	if(e_prime[1] < 0){flip = 10; e_prime[1] = abs(e_prime[1]);}  

	printvectors(e_prime, p_prime, print);

//WE ARE NOW IN THE REST FRAME OF THE TARGET IN THE ORENTATION OF A FIXED TARGET SCATTERING
//Randomly select scattered angle
	double scatter_theta = Random(i*(i+2))*360*rad;

// Scatter!! Scatter the electron to a random angle, using conservation of momentum to calulate all varibles.
    e_prime2[0] = -1*(e_prime[0]*mp/(cos(scatter_theta)*e_prime[0] - e_prime[0]-mp));

// Inelastic scattering! Some energy will be lost!
	double Eloss = e_prime2[0]*Random(i*((3+i)/4)*3);

//Recalculate energy for the e`
	double ISE = e_prime2[0]-Eloss;
	double nuw = e_prime[0] - ISE;
	e_prime2[0] =  ISE ;
    e_prime2[1] =  e_prime2[0]*cos(scatter_theta);
    e_prime2[2] =  e_prime2[0]*sin(scatter_theta);
    e_prime2[3] =  0.0 ;   

//Calculate x,Q,F_2. Call a function, sending
double F_two = Ftwo(e_prime[0], e_prime2[0],scatter_theta);
	double QQ = 4*e_prime[0]*e_prime2[0]*sin(scatter_theta/2)*sin(scatter_theta/2);
	double x = QQ/(2*0.989*((e_prime[0]-e_prime2[0])));


//Using results for e, calculate p's info. 
double scattered_p_theta;
double scattered_p_momentum;

    scattered_p_momentum = sqrt(e_prime[0]*e_prime[0] + e_prime2[0]*e_prime2[0] + 2*mp*e_prime[0] - 2*e_prime[0]*e_prime2[0] - 2*mp*e_prime2[0]);
    p_prime2[0] = sqrt(scattered_p_momentum*scattered_p_momentum +mp*mp );
   scattered_p_theta = asin(-e_prime2[0]*sin(scatter_theta)/scattered_p_momentum);
    p_prime2[1] = scattered_p_momentum*cos(scattered_p_theta);
    p_prime2[2] = scattered_p_momentum*sin(scattered_p_theta);
    p_prime2[3] = 0.0;  

printvectors(e_prime2, p_prime2, print);

double eErf = e_prime2[0]; // Electron energy after scattering in rest frame!

// We begin transforming back to the lab frame!
// If we needed to flip before, we unflip.
if(flip == 10){for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , 180*rad, k);
    p_prime[k] = rotationZ(p_prime2 , 180*rad, k); }   
    for(int k = 0; k < 4; k++) {
        e_prime2[k]=e_prime[k];
        p_prime2[k]=p_prime[k];}}

// Rotated back by delta.
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , Delta, k);
    p_prime[k] = rotationZ(p_prime2 , Delta, k); }

printvectors(e_prime, p_prime, print);

//Boost back!   
for(int k = 0; k < 4; k++){ 
    e_prime2[k] = Boostx(e_prime , gamma, -beta, k);
    p_prime2[k] = Boostx(p_prime , gamma, -beta, k); }

printvectors(e_prime2, p_prime2, print);

//Rotate back by Lambda
for(int k = 0; k < 4; k++){ 
    e_final[k] = rotationZ(e_prime2 , Lambda, k);
    p_final[k] = rotationZ(p_prime2 , Lambda, k); }

	double theta = atan(e_final[2]/e_final[1]); // Final lab frame scattered angle
// To determine the right quadrent. 
if(e_final[1] < 0 && e_final[2] < 0){theta = pi + theta;}	
if(e_final[1] < 0 && e_final[2] > 0){theta =  (pi+theta);}
if(e_final[1] > 0 && e_final[2] < 0){theta = 2*pi + theta;}
   
printvectors(e_final, p_final,  print);

	double deltaE_lab  = elec_Beam_momentum -e_final[0];
	double deltaE_rest = eEri-eErf;
	double Qsquared = 4 * elec_Beam_momentum * e_final[0] * sin(theta/2) *sin(theta/2);
	double xb	= Qsquared/(2*mp*((elec_Beam_momentum-e_final[0])));
	double invar_mass = sqrt(mp*mp - Qsquared + 2*mp*deltaE_lab);
	double invar_mass_RF = sqrt(mp*mp - Qsquared + 2*mp*deltaE_rest);

//##################################################

//Inelastic Cross Section
	double alpha = (1.0/137.0);
	double A = 4*alpha*alpha*e_final[0]*e_final[0]*pow(cos(theta/2.0),2)/pow(Qsquared,2);
	double B = F_two/mp;
	double C = 2 * (F_two/xb) *pow(tan(theta/2),2);
	double diffcross = A*(B+C);

// If the scattered angle is really close to the beam, I get NaN fo P(0).
	double L = p_final[0];

if(L != L){ N++; continue;   } //Forces the code to not record any NaN
if(Qsquared < 0.6){M++; continue;} //Does not record any unwanted q^2 values.
if(invar_mass != invar_mass){M++; continue;} 

// Output files!!!
root << Lambda/rad << " "<< Qsquared<<" "<<  theta/rad<<" "<< xb <<" "<<e_final[0]<<" "<<F_two<< " "<< Prot_momentum<< " "<<diffcross<<" "<<deltaE_lab<<" " <<deltaE_rest<<" "<< QQ << " "<<x<< " "<<invar_mass_RF<< " "<< invar_mass << "\n";

Pmom << Prot_momentum << " "<< random_Pmom <<endl;

// Counter for long runs.
	if(i/250000 == i/250000.0){  
		if(i >= 1000000){cout << "Event count  " << i/1000000.0<< " "<<"mil" << "\n";}  
		else cout << "Event count  " << i << "\n";	}

}// End of the beam loop! ////////////////////////////////////////

 root.close();
 rest.close();
 Pmom.close();

//Output statment for general information
 time_t finish = time(0) ;
 cout << i << " " << "electrons were sent." <<"\n";
 cout << M+N << " events have been thrown out due to q^2 Value or unphysical results." << "\n";

if(finish -start >= 60){
 	cout << "This program took " << floor((finish-start)/60) << " minutes and ";
	double secs = (finish-start)/60.0 - floor((finish-start)/60);
	cout << secs*60 << " seconds to run" <<endl; }
else{ cout << "This program took " << finish-start << " seconds to run"<<endl;}

cout <<i-(M+N)<< "  Electrons were scattered."<<endl<<endl<<endl;

}//End of main program
  
////////////////////////////////////////////////////////////////////////////////////////
////Functions!!!

// Rotation about the Z axis! Rotate CCW from the postive X axis.
double rotationZ(double A[], double angle,int i) { 
    double Prime[4] = {0.0};
    double rot[4][4] = {0.0};
    rot[1][1] = rot[2][2] = cos(angle);
    rot[1][2] = - sin(angle);
    rot[2][1] = sin(angle);
    rot[0][0] = rot[3][3] =1.0;
        for (int j = 0; j < 4; j++) {
          Prime[i] = Prime[i]+ A[j]*rot[i][j]  ; }
return Prime[i];
}

// Function to print out two four vectors E|x|y|z
void printvectors(double A[], double B[],int f ){
if(f==10){
for(int i=0; i < 4; i++){ cout <<setprecision(4)<< setw(9)<< A[i] << "|"  ;}
cout << "          ";
for(int j=0; j < 4; j++){ cout <<setprecision(4)<< setw(9)<< B[j] << "|"  ;}
cout << "\n";
}}

// Boost in the x axis!
double Boostx(double A[], double gamma, double beta, int i){
    double Prime[4] = {0.0};
    double B[4][4] = {0.0};
    B[0][1] = B[1][0] = -beta*gamma;
    B[0][0] = B[1][1] =gamma;
    B[2][2]=B[3][3] = 1.0;
        for (int j = 0; j < 4; j++) {
          Prime[i] = Prime[i]+ A[j]*B[i][j]  ; }
return Prime[i];
}

// random number generator: Thinking about making this a seperate called function in order to make the seed random: WIP
double Random(double A){
	Ran yourran(A);  
	double B = yourran.doub();
return B; }

//Currently unused
double RANDOM(double A){
	Ran yourran(A);
	double C = yourran.doub();
	double B = 0.4 - 0.05*C; 
return B;}

//Calculate Ftwo-Need Beam energy,scattered electron energy and angle
double Ftwo(double Beam, double Eprime, double theta){
	double QQ = 4*Beam*Eprime*sin(theta/2)*sin(theta/2);
	double xb = QQ/(2*0.989*((Beam-Eprime)));  
	double A = 1.22*exp(3.2*xb);
	double Ftwo_thr = 0;
	double lambda[2] = {0};
	double C[12] = {5.226,  1.225*pow(10,05),  4.837*pow(10,4),  3.415*pow(10,7),  1.904,  405.5,  527.6,  4756, };
	double beta =1;
	for(int i=0; i < 4; i++){lambda[0] = lambda[0] + C[i+8]*pow(xb,i); }
	if(QQ < A){lambda[1] = C[5] + C[6]*xb + C[7]*xb*xb;}
		else{lambda[1] = 0;} 
	for(int i=1; i < 6; i++){Ftwo_thr = Ftwo_thr + C[i-1]*pow((1-xb),i+2);}
	double Ftwo = beta*Ftwo_thr*(1+lambda[0]*log(QQ/A) + lambda[1]*pow(log(QQ/A),2));
return Ftwo; }

//Produces a momemtum disturbution to roughly match av18. WIP
double distribution(double A,double Fit[],double intersect, double max){
/* This function brings;
 A, a seed for the random number generator. It is currently the electron number.
Fit[8]: an eight part array that contains the coefficients from the av18 fit functions.
intersect: The intersection between the two exp. functions from av18.
max: The max number for our random number generator from av18.
*/
	// Generate a random number
	srand(A);
	double x = rand();
			x = (x/RAND_MAX)*max; 

	double B;
	double xlow  =  intersect - intersect*0.1;
	double xhigh =  intersect + intersect*0.1 ;
	double frac;
	double cor;

	 if(x <= xlow  ){B= (Fit[0] + Fit[1]*x)/(1+ Fit[2]*x + Fit[3]*x*x);}
		else if(x <= xhigh){frac = (x -xlow)/(xhigh-xlow); 
				B = (Fit[0] + Fit[1]*x)/(1+ Fit[2]*x + Fit[3]*x*x)*frac + (1-frac)*(Fit[0] + Fit[1]*x)/(1+ Fit[2]*x + Fit[3]*x*x);}
			else{B =(Fit[4] + Fit[5]*x)/(1+ Fit[6]*x + Fit[7]*x*x);}
		

	double p = B;
return p; }


//////// Old methonds


//	B = log(x/13.2412)/-6.761 + log(x/0.013428)/-1.6886; // equation extracted from avx 18? for first function
	//if(B > 1.36  ){ B = log(x/0.013428)/-1.6886 ;} // equation extracted from avx 18? for second function

//, B[1] =0.05:::    5.202,  2.475*pow(10,4),  9866,  1.554*pow(10,6) ,0.9987,  -4.38,  22.36,  11.3 , intersect at 0.03
//*pow(10,4) , B[1]=0.107   5.226,  1.225*pow(10,5),  4.837*pow(10,4),  3.415*pow(10,7),  1.904,  405.5,  527.6,  4756,  intersect at 0.0015

//	double a[8] = { 5.202,  2.475*pow(10,4),  9866,  1.554*pow(10,6) ,0.9987,  -4.38,  22.36,  11.3 };
//		 if(x <= 0.03){B= (a[0] + a[1]*x)/(1+ a[2]*x + a[3]*x*x);}
//		else{B =(a[4] + a[5]*x)/(1+ a[6]*x + a[7]*x*x);}}


/*
	double lambda = 13.2412; //constants for first function
	double beta   = -6.761;
	double gamma  = 0.013428; //constants for second function
	double delta  = -1.6886; */

/*else if(count == 0){	
	double b[8] = { 5.226,  1.225*pow(10,5),  4.837*pow(10,4),  3.415*pow(10,7),  1.904,  405.5,  527.6,  4756};
	 	if(x <= 0.0015){B= (b[0] + b[1]*x)/(1+ b[2]*x + b[3]*x*x);}
		else{B =(b[4] + b[5]*x)/(1+ b[6]*x + b[7]*x*x);}}

else{
	if(A <= half){
	double a[8] = {75.85,  7.604*pow(10,6),  1.293*pow(10,6),  3.097*pow(10,9), 6274, 5.493*pow(10,6),7.169*pow(10,6),2.098*pow(10,8)  };
			 if(x <=  0.001618){B= (a[0] + a[1]*x)/(1+ a[2]*x + a[3]*x*x);}
			else{B =(a[4] + a[5]*x)/(1+ a[6]*x + a[7]*x*x);}}
	else{
			double b[8] = { 5.226,  1.225*pow(10,5),  4.837*pow(10,4),  3.415*pow(10,7),  1.904,  405.5,  527.6,  4756};
	 		if(x <= 0.0015){B= (b[0] + b[1]*x)/(1+ b[2]*x + b[3]*x*x);}
			else{B =(b[4] + b[5]*x)/(1+ b[6]*x + b[7]*x*x);}}
		}

	double p = B;
return p; }




*/
