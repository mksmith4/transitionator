#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream> 
#include <stdlib.h>
using namespace std;

// Coefficients
// doubles have 15 digits
double cE1=1.5906248778818252500e15;
double cE2=1.2255069133707089424e9;
double cE3=570.94451255476974438;
double cE4=0.00016970979189532430951;
double cE5=3.4579563256864936680e-11;

double cM1=1.7588286310724052734e13;
double cM2=1.3551005499699665233e7;
double cM3=6.3132016190533493116;
double cM4=1.8765608730848246391e-6;
double cM5=3.8236247120153468948e-13;

double wE1=0.064457751952217604008;
double wE2=0.059404264199163746285;
double wE3=0.059404264199163746285;
double wE4=0.062847286858870379800;
double wE5=0.069289133761904597719;

double wM1=1.7904931097838225895;
double wM2=1.6501184499767707070;
double wM3=1.6501184499767707070;
double wM4=1.7457579683019552874;
double wM5=1.9246981600529053935;

//Units
double Eu, Tu, half_or_mean, delt;

//Variables
double Energy, cE, cM, wE, wM;
double Egamma, Eerr, Lifeerr, l, Life, A, a, branch, frac;	
string Eunits, Tunits, EM, h_or_m, mixtype, TEXtable;

//Arrays
double mass[1000], energy[1000], enerr[1000], br[1000], halflife[1000], terr[1000], alph[1000], ang[1000], fract[1000], delta[1000];
string mult[1000], UNITStime[1000], UNITSenergy[1000]; //, WeissE[1000], WeissM[1000];

//File for latex output
ofstream outputFile;

//Reduced Transition Probability
double B (double E, double T, double BR, double alpha, double coeff, double lam, double fr)
{
  double b;
  if (lam == 1 ) {coeff == cE1;}
  b=fr*BR*log(2)/(T*(1+alpha)*coeff*pow(E, 2*lam+1));
  return b;
}

double errB (double E, double E_err, double T, double T_err, double BR, double alpha, double coeff, double lam, double fr)
{
  double errb;
  double db_E, db_T, db_a, db_BR, db_fr;
  db_fr = BR*log(2)/(T*(1+alpha)*coeff*pow(E, 2*lam+1));
  db_BR = fr*log(2)/(T*(1+alpha)*coeff*pow(E, 2*lam+1));
  db_T = fr*BR*log(2) / (pow(T, 2)*(1 + alpha)*coeff*pow(E, 2*lam+1));
  db_a = fr*BR*log(2) / (T*pow((1+alpha), 2)*coeff*pow(E, 2*lam+1));
  db_E = fr*BR*log(2)*(-2*lam-1) / ( T*(1 + alpha)*coeff*pow(E, 2+2*lam));
  
  errb= sqrt(pow(db_E*E_err, 2) + pow(db_T*T_err, 2));
  return errb;
}


//Weisskopf Single Particle transition probabilities
double BWeissE(double mass, double lam, double Wcoeff)
{
  double bwe;
  bwe = Wcoeff * pow(mass, 2*lam/3);
  return bwe;
}

double BWeissM(double mass, double lam, double Wcoeff)
{
  double bwm;
  bwm = Wcoeff * pow(mass, (2*lam-2)/3);
  return bwm;
}


//Main Program - two ways to run - with an input file, or interactively
int main(int argc, char* argv[])
{
  cout << "--------------------------------------------------------------------------\n";
  cout << "-   A C++ calculator to find Gamma-ray transition probabilities          -\n";
  cout << "-   [B(E/M l) values] in Weisskopf units.                                -\n";
  cout << "- If you want to calculate the transition probability                    -\n";
  cout << "- for more than one state, you can use an input file.                    -\n";
  cout << "- Use: > ./TransProb.exe -help                                           -\n";
  cout << "-   to read the full instructions / options                              -\n";
  cout << "- Author: Mallory Smith, msmith40@nd.edu, July 2014                      -\n";

if (argc == 2 && argc != 3 && argc != 1) 
{ 
        if (strcmp(argv[1], "-help") == 0);
	{
        cout << "---------------------------- Full Instructions----------------------------\n";        
        cout << "- This is compiler specific - to compile from source:                    -\n";
        cout << "- g++ -o TransProb.exe TransProb.cpp                                     -\n";
	cout << "-                                                                        -\n";
	cout << "- The input file should have the following format:                       -\n";
	cout << "A   Gam_En Err Units B.R. T1/2 Err Units Alpha(ICC) E/M? AngMom Fraction -\n";
        cout << "156 100    1   keV   1    20   0.2 ps    0.082      E    2      0.10     -\n";
	cout << "- Spaces don't matter. (Spreadsheet compatible)                          -\n";
	cout << "-                                                                        -\n";
	cout << "- There are different options to evaluate:                               -\n";
        cout << "- ./TransProb.exe                 -- run interactively                   -\n";
        cout << "- ./TransProb.exe -f file.dat                                            -\n";
	cout << "- -- evaluate 'file'assuming T = HALF-life; output results to terminal   -\n";
        cout << "- ./TransProb.exe -m file.dat                                            -\n";
	cout << "- -- evaluate 'file'assuming T = MEAN-life; output results to terminal   -\n";
        cout << "- ./TransProb.exe -f file.dat -latex                                     -\n";
	cout << "-     -latex flag = output results to latex-ready table, and to terminal -\n";
	cout << "-                                                                        -\n";
        cout << "- this can calculate up to E5/M5                                         -\n";
        cout << "-     Energy must be in units, from eV to MeV and lifetimes              -\n";
	cout << "-     Lifetimes must be from days to femtoseconds                        -\n";
        cout << "-     Branching Ratio must be a fractional percent!!                     -\n";
        cout << "- Careful! Branching ratio = I_gamma (1+alpha) / [SUM I_gamma*(1+alpha)] -\n";
	cout << "-                                                                        -\n";
        cout << "- Also note - in the case of a mixed transition, ie, M1+E2, the mixing   -\n";
        cout << "- ratio, 'delta', must be converted to fraction E2 and M1.               -\n";
        cout << "- This is delta^2 / (1 + delta^2) for the L+1 transition, and    -\n";
        cout << "- 1 / (1 + delta^2) for the L transition (holds for M1+E2 and higher     -\n";
        cout << "- order mixed transitions, E3+M2, etc.                                   -\n";
        cout << "- Thus, the E2 and M1 contributions must be listed as fractions of the   -\n";
	cout << "- total intensity                                                        -\n";
        cout << "- One last thing - all errors here are assumed to be symmetric!          -\n";
	} 
} 
 cout << "--------------------------------------------------------------------------\n";
if ( argc == 3 || argc == 4) 
 {
    cout << "Opening" << " "  << argv[2] << endl;
	if( argc == 4) {    

	if (strcmp(argv[3], "-latex") == 0) {   TEXtable == "yes"; 
						cout << "Output saved to file results.tex" << endl;	
						ofstream outputFile;
						outputFile.open("results.tex");
						outputFile << "\\documentclass[preprint,aps]{revtex4-1}\n \\usepackage[margin=0.4in]{geometry}\n \\begin{document}\n \\begin{table}[h]\n \\setlength{\\tabcolsep}{4pt}\n \\begin{tabular}{| c | c | c | c | c | c | c |}\n \\hline E$_{\\gamma}$ & Branching Ratio & Lifetime & $\\alpha$ & Multipolarity & Fraction & B(EM) [W.u.] \\\\ \\hline \\hline" << endl;
					     };
	if (strcmp(argv[3], "-latex") == 1 || argv[3] == NULL) { TEXtable == "no";}
			}
//**************************************************************************************************************
//This part is for the interactive calculator... The main bit for the interactive part is below.
double half_or_mean;
string ans;
    
        if (strcmp(argv[1], "-f") == 0) {half_or_mean=1; cout << "HALF LIFE given" << endl;}
        if (strcmp(argv[1], "-m") == 0) {half_or_mean=0.693; cout << "mean-life given" << endl;}
//**************************************************************************************************************

   string ifile;
   ifile = argv[2];

//opening the data file, name provided in command line...
     ifstream Chicago(ifile.c_str());
	     if(Chicago. fail())
      		{
		cerr << "could not open input file \n" << ifile << endl;
		return 0;
      		}

     	     for(int j = 0; j < 100; ++j)
		{
		Chicago >> mass[j] >> energy[j] >> enerr[j] >> UNITSenergy[j] >> br[j] >> halflife[j] >> terr[j] >> UNITStime[j] >> alph[j]  >> mult[j] >> ang[j] >> fract[j];
		}

double g[1000], h[1000], gerr[1000], herr[1000], WeissE[1000], WeissM[1000];


// Starting the calculations... L=1, L=2, etc each calculated separately, using the same B(E/M) formula as defined above. The coefficients are assigned first, then the units are corrected, then the corrected numbers are passed to the B(E/M) formula.

        printf ("%s %s %3s %11s %6s %6s %3s \n", "Mass", "Energy (keV)", "ICC", "Lifetime ", "Mult.", "Fraction","B(E/Ml) [w.U.]");

            for(int m = 0; m < 87; ++m)
               {
		if ( ang[m] == 1 ) 
		{	
//cM, wM are the coefficients for the B(E/M) formula and the Weisskopf estimate, respectively.
	        cM = cM1;	
		wM = wM1;
		cE = cE1;
		wE = wE1;	
			if( UNITSenergy[m] == "MeV") { Eu=1; }	
			if( UNITSenergy[m] == "keV") { Eu=1e-3;}		
			if( UNITSenergy[m] == "eV") { Eu=1e-6;}
			if( UNITStime[m] == "days" || UNITStime[m] == "day" || UNITStime[m] == "d" ) { Tu = half_or_mean*3600*24;}
			if( UNITStime[m] == "hours" || UNITStime[m] == "hour" || UNITStime[m] == "hrs" || UNITStime[m] == "hr" || UNITStime[m] == "h"  ) { Tu = half_or_mean*3600;}
			if( UNITStime[m] == "min" || UNITStime[m] == "m" ) {Tu = half_or_mean*60;}
			if( UNITStime[m] == "sec" || UNITStime[m] == "s" ) {Tu = half_or_mean*1;}
			if( UNITStime[m] == "ms" ) {Tu = half_or_mean*1e-3;}
			if( UNITStime[m] == "us" ) {Tu = half_or_mean*1e-6;}
			if( UNITStime[m] == "ns" ) {Tu = half_or_mean*1e-9;}
			if( UNITStime[m] == "ps" ) {Tu = half_or_mean*1e-12;}
			if( UNITStime[m] == "fs" ) {Tu = half_or_mean*1e-15;}
			
	//	if( delta[m] != 0 ) { frac = 1/(1 + pow(delta[m], 2));} else {frac=1;}

		g[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cE, 1, fract[m]) / BWeissE(mass[m], ang[m], wE);
		gerr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cE, 1, 1) / BWeissE(mass[m], ang[m], wE);
		WeissE[m] = BWeissE(mass[m], ang[m], wE);
		h[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cM, 1, fract[m]) / BWeissM(mass[m], ang[m], wM);
		herr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cM, 1, 1) / BWeissM(mass[m], ang[m], wM);
		WeissM[m] = BWeissM(mass[m], ang[m], wM);
		}

		if ( ang[m] == 2 )
		{	
	        cE = cE2;
		cM = cM2;	
		wE = wE2;
		wM = wM2;
			if( UNITSenergy[m] == "MeV") { Eu=1; }	
			if( UNITSenergy[m] == "keV") { Eu=1e-3;}		
			if( UNITSenergy[m] == "eV") { Eu=1e-6;}
			if( UNITStime[m] == "days" || UNITStime[m] == "day" || UNITStime[m] == "d" ) { Tu = half_or_mean*3600*24;}
			if( UNITStime[m] == "hours" || UNITStime[m] == "hour" || UNITStime[m] == "hrs" || UNITStime[m] == "hr" || UNITStime[m] == "h"  ) { Tu = half_or_mean*3600;}
			if( UNITStime[m] == "min" || UNITStime[m] == "m" ) {Tu = half_or_mean*60;}
			if( UNITStime[m] == "sec" || UNITStime[m] == "s" ) {Tu = half_or_mean*1;}
			if( UNITStime[m] == "ms" ) {Tu = half_or_mean*1e-3;}
			if( UNITStime[m] == "us" ) {Tu = half_or_mean*1e-6;}
			if( UNITStime[m] == "ns" ) {Tu = half_or_mean*1e-9;}
			if( UNITStime[m] == "ps" ) {Tu = half_or_mean*1e-12;}
			if( UNITStime[m] == "fs" ) {Tu = half_or_mean*1e-15;}

	//	if( delta[m] != 0 ) { frac = pow(delta[m], 2)/(1 + pow(delta[m], 2));} else {frac=1;}

		g[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cE, 2, fract[m]) / BWeissE(mass[m], ang[m], wE);
		gerr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu, terr[m]*Tu, br[m], alph[m], cE, 2, 1) / BWeissE(mass[m], ang[m], wE);
		WeissE[m] = BWeissE(mass[m], ang[m], wE);
		h[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cM, 2, fract[m]) / BWeissM(mass[m], ang[m], wM);
		herr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cM, 2, 1) / BWeissM(mass[m], ang[m], wM);
		WeissM[m] = BWeissM(mass[m], ang[m], wM);
		}

		if ( ang[m] == 3 )
		{	
        	cE = cE3;
		cM = cM3;
		wE = wE3;
		wM = wM3;
			if( UNITSenergy[m] == "MeV") { Eu=1; }	
			if( UNITSenergy[m] == "keV") { Eu=1e-3;}		
			if( UNITSenergy[m] == "eV") { Eu=1e-6;}
			if( UNITStime[m] == "days" || UNITStime[m] == "day" || UNITStime[m] == "d" ) { Tu = half_or_mean*3600*24;}
			if( UNITStime[m] == "hours" || UNITStime[m] == "hour" || UNITStime[m] == "hrs" || UNITStime[m] == "hr" || UNITStime[m] == "h"  ) { Tu = half_or_mean*3600;}
			if( UNITStime[m] == "min" || UNITStime[m] == "m" ) {Tu = half_or_mean*60;}
			if( UNITStime[m] == "sec" || UNITStime[m] == "s" ) {Tu = half_or_mean*1;}
			if( UNITStime[m] == "ms" ) {Tu = half_or_mean*1e-3;}
			if( UNITStime[m] == "us" ) {Tu = half_or_mean*1e-6;}
			if( UNITStime[m] == "ns" ) {Tu = half_or_mean*1e-9;}
			if( UNITStime[m] == "ps" ) {Tu = half_or_mean*1e-12;}
			if( UNITStime[m] == "fs" ) {Tu = half_or_mean*1e-15;}
		g[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cE, 3, fract[m]) / BWeissE(mass[m], ang[m], wE);
		gerr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cE, 3, 1) / BWeissE(mass[m], ang[m], wE);
		WeissE[m] = BWeissE(mass[m], ang[m], wE);
		h[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cM, 3, fract[m]) / BWeissM(mass[m], ang[m], wM);
		herr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cM, 3, 1) / BWeissM(mass[m], ang[m], wM);
		WeissM[m] = BWeissM(mass[m], ang[m], wM);
		}

		if ( ang[m] == 4 )
		{	
    	        cE = cE4;
		cM = cM4;
		wE = wE4;
		wM = wM4;
			if( UNITSenergy[m] == "MeV") { Eu=1; }	
			if( UNITSenergy[m] == "keV") { Eu=1e-3;}		
			if( UNITSenergy[m] == "eV") { Eu=1e-6;}
			if( UNITStime[m] == "days" || UNITStime[m] == "day" || UNITStime[m] == "d" ) { Tu = half_or_mean*3600*24;}
			if( UNITStime[m] == "hours" || UNITStime[m] == "hour" || UNITStime[m] == "hrs" || UNITStime[m] == "hr" || UNITStime[m] == "h"  ) { Tu = half_or_mean*3600;}
			if( UNITStime[m] == "min" || UNITStime[m] == "m" ) {Tu = half_or_mean*60;}
			if( UNITStime[m] == "sec" || UNITStime[m] == "s" ) {Tu = half_or_mean*1;}
			if( UNITStime[m] == "ms" ) {Tu = half_or_mean*1e-3;}
			if( UNITStime[m] == "us" ) {Tu = half_or_mean*1e-6;}
			if( UNITStime[m] == "ns" ) {Tu = half_or_mean*1e-9;}
			if( UNITStime[m] == "ps" ) {Tu = half_or_mean*1e-12;}
			if( UNITStime[m] == "fs" ) {Tu = half_or_mean*1e-15;}
		g[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cE, 4, fract[m]) / BWeissE(mass[m], ang[m], wE);
		gerr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cE, 4, 1) / BWeissE(mass[m], ang[m], wE);
		WeissE[m] = BWeissE(mass[m], ang[m], wE);
		WeissM[m] = BWeissM(mass[m], ang[m], wM);
		h[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cM, 4, fract[m]) / BWeissM(mass[m], ang[m], wM);
		herr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cM, 4, 1) / BWeissM(mass[m], ang[m], wM);
		}

		if ( ang[m] == 5 )
		{	
   	   	cE = cE5;
		cM = cM5;
		wE = wE5;
		wM = wM5;
			if( UNITSenergy[m] == "MeV") { Eu=1; }	
			if( UNITSenergy[m] == "keV") { Eu=1e-3;}		
			if( UNITSenergy[m] == "eV") { Eu=1e-6;}
			if( UNITStime[m] == "days" || UNITStime[m] == "day" || UNITStime[m] == "d" ) { Tu = half_or_mean*3600*24;}
			if( UNITStime[m] == "hours" || UNITStime[m] == "hour" || UNITStime[m] == "hrs" || UNITStime[m] == "hr" || UNITStime[m] == "h"  ) { Tu = half_or_mean*3600;}
			if( UNITStime[m] == "min" || UNITStime[m] == "m" ) {Tu = half_or_mean*60;}
			if( UNITStime[m] == "sec" || UNITStime[m] == "s" ) {Tu = half_or_mean*1;}
			if( UNITStime[m] == "ms" ) {Tu = half_or_mean*1e-3;}
			if( UNITStime[m] == "us" ) {Tu = half_or_mean*1e-6;}
			if( UNITStime[m] == "ns" ) {Tu = half_or_mean*1e-9;}
			if( UNITStime[m] == "ps" ) {Tu = half_or_mean*1e-12;}
			if( UNITStime[m] == "fs" ) {Tu = half_or_mean*1e-15;}
		g[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cE, 5, fract[m]) / BWeissE(mass[m], ang[m], wE);
		gerr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cE, 5, 1) / BWeissE(mass[m], ang[m], wE);
		h[m] = B(energy[m]*Eu, halflife[m]*Tu*half_or_mean, br[m], alph[m], cM, 5, fract[m]) / BWeissM(mass[m], ang[m], wM);
		herr[m] = errB(energy[m]*Eu, enerr[m]*Eu, halflife[m]*Tu*half_or_mean, terr[m]*Tu*half_or_mean, br[m], alph[m], cM, 5, 1) / BWeissM(mass[m], ang[m], wM);
		WeissE[m] = BWeissE(mass[m], ang[m], wE);
		WeissM[m] = BWeissM(mass[m], ang[m], wM);
		}


        if( mult[m] == "E" )
	{
	if (TEXtable == "no"){
	cout << mass[m] << "  " << energy[m] << " "  << enerr[m] << " " << br[m] << " " << halflife[m] << " " << terr[m] << " " << UNITStime[m] << " " << alph[m] << " " << mult[m] << ang[m] << " " << fract[m] << " "  << g[m] << " " << gerr[m] << " " << WeissE[m] << endl;
	}
        if (TEXtable == "yes") {
//	outputFile << energy[m] << " $\\pm$ "  << enerr[m] << " " << UNITSenergy[m] << " & " << br[m] << " & " << halflife[m] << "$\\pm$ " << terr[m] << " " << UNITStime[m] << " & " << alph[m] << " & " << mult[m] << ang[m] << " & " << fract[m] << " & "  << g[m] << " $\\pm$ " << gerr[m] << "\\\\ \\hline" << endl;
		}
 	}


        if( mult[m] == "M" )
	{
	if (argv[3] == NULL){
	cout << mass[m] << "  " << energy[m] << " "  << enerr[m] << " " << br[m] << " " << halflife[m] << " " << terr[m] << " " << UNITStime[m] << " " <<  alph[m] << " " << mult[m] << ang[m] << " " << fract[m] << " "  << h[m] << " " << herr[m] << " " << WeissM[m] << endl;
}
        if (TEXtable == "no") {
//	outputFile << energy[m] << " $\\pm$ "  << enerr[m] << " " << UNITSenergy[m] << " & " << br[m] << " & " << halflife[m] << " $\\pm$ " << terr[m] << " " << UNITStime[m] << " & " << alph[m] << " & " << mult[m] << ang[m] << " & " << fract[m] << " & "  << g[m] << " $\\pm$ " << gerr[m] << "\\\\ \\hline" << endl;
		}
	}
 }
outputFile << "\\end{tabular}\n \\caption{}\n \\end{table}\n \\end{document}\n";
 }

//**************************************************************************************************************
// Interactive Mode - needs more work, currently can not compute B(E/M) values.

  else { 
	if ( argc == 1 && argc != 2 && argc !=3 )
	{
		cout << "Running in interactive mode " << endl;
	
// Inputs *********************************************        
	cout << "Do you want to use Mean lifetime (mean) or Half-life (half) ?: " << endl;
	cin >> h_or_m;	
	cout << "Mass Number, A " << endl;
	cin >> A;
 	cout << "Transition energy + Err + Units (ie, keV):" << endl;
        cin >> Egamma >> Eerr >> Eunits;
	cout << "Lifetime + Err + Units (ie, ps):" << endl;
	cin >> Life >> Lifeerr >> Tunits;
        cout << "Electric (E) or Magnetic (M) " << endl;
        cin >> EM;
	cout << "Angular momentum  " << endl;
	cin >> l;
        cout << "Internal Conversion (alpha) " << endl;
        cin >> a;
	cout << "Is there mixing? If not, enter zero\n delta = ";
	cin >> delt;
	if(delt != 0){cout << "What type of mixing? Enter E2+M1, or E3+M2\n If it's some other kind of mixing, this program will not calculate it\n but can be easily modified to handle other types of mixing." << endl;
	cin >> mixtype;}
	cout << "Branching Ratio " << endl;
	cin >> branch;
// ***************************************************        
	}

	//if( delt != 0 ) { frac = pow(delt, 2)/(1 + pow(delt, 2));} else {frac=1;}

	if( Eunits == "MeV") { Eu=1; }	
	if( Eunits == "keV") { Eu=1e-3;}			
	if( Eunits == "eV") { Eu=1e-6;}
	if( Tunits == "days" || Tunits == "day" || Tunits == "d" ) { Tu = 3600*24;}
	if( Tunits == "hours" || Tunits == "hour" || Tunits == "hrs" || Tunits == "hr" || Tunits == "h"  ) { Tu = 3600;}
	if( Tunits == "min" || Tunits == "m" ) {Tu = 60;}
	if( Tunits == "sec" || Tunits == "s" ) {Tu = 1;}
	if( Tunits == "ms" ) {Tu = 1e-3;}
	if( Tunits == "us" ) {Tu = 1e-6;}
	if( Tunits == "ns" ) {Tu = 1e-9;}
	if( Tunits == "ps" ) {Tu = 1e-12;}
	if( Tunits == "fs" ) {Tu = 1e-15;}

 	if ( l == 1 )
 	{	
        cE = cE1;
	cM = cM1;
	wE = wE1;
	wM = wM1;
	if(mixtype=="E2+M1"){frac = 1/(1+pow(delt,2));} else{frac=1;}
 	}

	if ( l == 2 )
	{	
        cE = cE2;
	cM = cM2;
	wE = wE2;
	wM = wM2;
	if(mixtype=="E2+M1"){frac = pow(delt, 2)/(1+pow(delt,2));} else{frac=1;}
	if(mixtype=="E3+M2"){frac = 1/(1+pow(delt,2));} else{frac=1;}
	}

	if ( l == 3 )
	{	
        cE = cE3;
	cM = cM3;
	wE = wE3;
	wM = wM3;
	if(mixtype=="E3+M2"){frac = pow(delt, 2)/(1+pow(delt,2));} else{frac=1;}
	}

	if ( l == 4 )
	{	
        cE = cE4;
	cM = cM4;
	wE = wE4;
	wM = wM4;
	}

	if ( l == 5 )
	{	
        cE = cE5;
	cM = cM5;
	wE = wE5;
	wM = wM5;
	}

    if ( h_or_m == "half" )
	{
	half_or_mean=1;
	}
    if ( h_or_m == "mean" )
	{
	half_or_mean=0.693;
	}


double z, x, fraction=1;
        if( EM == "E" )
	{
//	z = B(energy[m]*Eu, halflife[m]*Tu, br[m], alph[m], cE, 2, frac) / BWeissE(mass[m], ang[m], wE);
	z = B(Egamma*Eu, Life*half_or_mean*Tu, branch, a, cE, l, frac) / BWeissE(A, l, wE);
	cout << Energy << " " << Life << " " << a << " E" << l << " " << z << endl;
	}

        if( EM == "M" )
	{
//	x = B(Egamma, Life, branch, a, cM, l) / BWeissM(A, l, wM);
	x = B(Egamma*Eu, Life*half_or_mean*Tu, branch, a, cE, l, frac) / BWeissE(A, l, wE);
	cout << "M" << l << " (w. U.) " << x << endl;
	}

 }

//**************************************************************************************************************
    return 0;
}
