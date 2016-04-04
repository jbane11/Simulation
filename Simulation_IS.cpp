#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include "nr3.h"
#include "ran.h"
using namespace std;

//Function protype.
double rotationZ(double A[4], double angle,int i);
void printvectors(double A[4], double B[4],int f);
double Ftwo(double Beam, double Eprime, double theta);
double Boostx(double A[4], double gamma, double beta, int i);
double Random(double A);
double RANDOM(double A);
int main(void){

	ofstream root;
	root.open("data/ISroot.txt");
	ofstream rest;
	rest.open("data/ISrest.txt");
// Constants

    double mp = 0.989 ;             // Mass of Proton in GeV/c^2
	double me = .0005 ;             // Mass of Electron in GeC/c^2
	double c  = 3.0 * pow(10,8);    // Speed of light in m/s
	double pi = 3.14159265359;      // Pi
	double rad = pi/180;            // Convert from degrees to radians
    double e[4];                   //Electron vector array
    double p[4];                   //Proton vector array
    double e_prime[4];             //ELectron vector array in prime frame
    double p_prime[4];             //Proton vector array in prime frame
    double e_prime2[4];             //ELectron vector array in prime frame
    double p_prime2[4];             //Proton vector array in prime frame
    double elec_Beam_momentum = 5;  //Beam momentum in GeV
    double Prot_momentum = 0.50     ;  //Proton momentum in GeV   
	int number_of_electrons = 100000;  
	double min_angle = 0;// By defaults scattered angle covers 360 degrees
	double angle_range = 360;  
	char check;	
	int N = 0;
	int M = 0;
	int print = 0; //If you want to check the 4 vectors during the rotations set this to 10;
// Loop for the number of electrons
	for(int i=1; i <= number_of_electrons; i++){

// Random angle in degrees used for Lambda
	double Lambda = Random(i)*360*rad;

// incomeing electron and proton 4 vectors.
    e[0] = elec_Beam_momentum;      //electron energy
    e[1] = elec_Beam_momentum;      //Electron Beam Momentum
    e[2] = 0 ;
    e[3] = 0 ;
    p[1] = Prot_momentum*cos(Lambda);   //Proton x momentum
    p[2] = Prot_momentum*sin(Lambda);   //Proton y momentum
    p[3] = 0;
    p[0] = sqrt( mp*mp+Prot_momentum*Prot_momentum); //Proton energy
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
// Due to limitations of atan, forceing the e[1] to be postive. This makes the beam always head to the right, in my internal drawling. This will also let me know if we need to flip back.
int flip =0;
if(e_prime[1] < 0){flip = 10; e_prime[1] = abs(e_prime[1]);}  

printvectors(e_prime, p_prime, print);
//WE ARE NOW IN THE REST FRAME OF THE TARGET IN THE ORENTATION OF A FIXED TARGET SCATTERING
//Randomly select scattered angle
	double scatter_theta = Random(i*(i+2))*360*rad;

// Scatter!! Scatter the electron to a random angle, using conservation of momentum to calulate all varibles.

    e_prime2[0] = -1*(e_prime[0]*mp/(cos(scatter_theta)*e_prime[0] - e_prime[0]-mp));

// Inelastic scattering! Some energy will be lost!
	double DeltaE = e_prime2[0]*Random(i*((3+i)/4)*3);
//Recalculate energy for the e`
	double ISE = e_prime2[0]-DeltaE;
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

// We begin transforming back to the lab frame!
// If we needed to flip before, we unflip.
if(flip == 10){for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , 180*rad, k);
    p_prime[k] = rotationZ(p_prime2 , 180*rad, k); }   
    for(int k = 0; k < 4; k++) {
        e_prime2[k]=e_prime[k];
        p_prime2[k]=p_prime[k];}
	}

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

double e_final[4],p_final[4];

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
//cout << theta/rad ;

	double Qsquared = 4 * elec_Beam_momentum * e_final[0] * sin(theta/2) *sin(theta/2);
	double xb	= Qsquared/(2*mp*((elec_Beam_momentum-e_final[0])));

//##################################################

//Cross section
double	Mottdiff =  (1.0/137.0)*(1.0/137.0)*pow(cos(theta/2.0),2)/(4.0*elec_Beam_momentum*e_final[0]*pow(sin(theta/2.0),4));
double Tau = Qsquared/(4*mp*mp);
double Ge = pow((1+Qsquared)/.71,-2);
double Gm = 2.79*Ge;
double frec = 1 + (2*elec_Beam_momentum/mp)*pow(sin(theta/2),2);
double diffcross = (Mottdiff/frec)*((Ge*Ge+Tau*Gm*Gm)/(1+Tau) + 2*Tau*Gm*Gm*pow(tan(theta/2),2));

// If the scattered angle is really close to the beam, I get NaN fo P(0).
double L = p_final[0];

if(L != L){ N++; continue;   } //Forces the code to not record any NaN
if(Qsquared < 0.6){M++; continue;} //Does not record any unwanted q^2 values.

root <<setprecision(3)<< setw(5)<< Lambda/rad<< " "<< setw(7)<< Qsquared<< " "<< setw(7)<< theta/rad<< " " <<  setw(8)<<  xb << " " << setw(7)<<  e_final[0]<< " " <<  setw(9)<< F_two<< " "<<  setw(9)<< DeltaE<< " " <<  setw(7)<< diffcross<<"\n";
rest <<setprecision(3)<< setw(7)<<x << " "<<setprecision(3)<< setw(7)<< QQ << " "<<setprecision(3)<< setw(8)<< F_two<< " "<<setprecision(3)<< setw(8) << F_two/nuw << " "<<setprecision(3)<< setw(7)<< nuw << " "<< xb <<" "<< elec_Beam_momentum-e_final[0]<<"\n" ;

if(i/500000 == i/500000.0){cout << "Event count  " << i << "\n";}  
//cout <<setprecision(3)<< setw(5)<< x <<"\t" << xb << "\t" <<nuw << "\t" <<elec_Beam_momentum-e_final[0] << "\n";
}// End of the beam loop!



root.close();
rest.close();
}
  
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

double Random(double A){
	Ran yourran(A);  
	double B = yourran.doub();
return B; }

double RANDOM(double A){
	Ran yourran(A);
	double C = yourran.doub();
	double B = 0.4 - 0.05*C; 
	
return B;}

double Ftwo(double Beam, double Eprime, double theta){
	double QQ = 4*Beam*Eprime*sin(theta/2)*sin(theta/2);
	double xb = QQ/(2*0.989*((Beam-Eprime)));  
	double A = 1.22*exp(3.2*xb);
	double Ftwo_thr = 0;
	double lambda[2] = {0};
	double C[12] = {1.417,-0.108,1.486,-5.979,3.524,-0.011,-0.619,1.385,0.270,-2.179, 4.722, -4.363};
	double beta =1;
	for(int i=0; i < 4; i++){lambda[0] = lambda[0] + C[i+8]*pow(xb,i); }
	if(QQ < A){lambda[1] = C[5] + C[6]*xb + C[7]*xb*xb;}
		else{lambda[1] = 0;} 
	for(int i=1; i < 6; i++){Ftwo_thr = Ftwo_thr + C[i-1]*pow((1-xb),i+2);}
	double Ftwo = beta*Ftwo_thr*(1+lambda[0]*log(QQ/A) + lambda[1]*pow(log(QQ/A),2));
return Ftwo; }




