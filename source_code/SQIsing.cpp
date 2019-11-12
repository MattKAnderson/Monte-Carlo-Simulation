/* Program that uses the MC functions to perform a Monte Carlo simulation of
 * the square lattice Ising model. */

#include "MC_functions.h"
#include "outputFunctions.h"
#include <iostream>
#include <fstream>
using namespace std;

double calcg(double M4, double M2);
vector <int> LB, LF;

int main()
{
	int L, NWarmup, NSweep, Nmeas, Tpoints;
	double Averages[6], error[4], T, Tmin, Tmax, Tstep;
	vector <double> table;

	cout << "Please input the lattice size: ";
	cin >> L;
	cout << "Please input the minimum temperature in units of J/k: ";
	cin >> Tmin;
	cout << "Please input the maximum temperature in units of J/k: ";
	cin >> Tmax;
	cout << "Please input the desired number of temperature points: ";
	cin >> Tpoints;
	cout << "Please input the number of Warm - up MC sweeps: ";
	cin >> NWarmup;
	cout << "Please input the number of MC sweeps between measurements: ";
	cin >> NSweep;
	cout << "Please input the number of measurements to take: ";
	cin >> Nmeas;
	
	//creating, looped boundary, Spin Configuration and Measurement vectors	
	LB.resize(L + 1);
	LF.resize(L + 1);
	vector < vector <int> > SpinConf(L, vector <int> (L));
	vector < vector <long> > Measurements(6, vector <long> (Nmeas));
	
	//vector for average energy(0,1), X(2,3), Cv/k(4,5) and their error bar and Temperature(6)
	vector < vector <double> > CalculatedValues(8, vector <double> (Tpoints));

	//initializing looped boundary vectors
	for (int i = 1; i < L; i++)
	{
		LB[i] = i - 1;
	}
	LB[0] = L - 1;
	for (int i = 0; i < (L - 1); i++)
	{
		LF[i] = i + 1;
	}
	LF[L - 1] = 0;
	
	//determining temperature step size
	Tstep = (Tmax - Tmin) / (Tpoints - 1);

	//running simulation for k Tpoints
	for(int k = 0; k < Tpoints; k++)
	{
		//determining value of T for iteration of simulation
		T = Tmin + Tstep * k;		
		CalculatedValues[6][k] = T;
		//initializing Spin conformation
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < L; j++)
			{	
				if ( getRandom() < 0.5)
				{
					SpinConf[i][j] = -1;
					//cout << "Position (" << i << ", " << j << ") is " << SpinConf[i][j] << endl;
				}
				else
				{
					SpinConf[i][j] = 1;
					//cout << "Position (" << i << ", " << j << ") is " << SpinConf[i][j] << endl;
				}
			}
		}
		//generating exponent(delta E) table
		Exponent_table_calculator(L, T, table);		
	
		//Warm up sweeps
		Measurements[0][0] = MCSweeps(SpinConf, table, L,T, NWarmup); 
		//cout << "finished warm - up fine\n";
		//Running simulation and getting measurements
		for (int i = 0; i < Nmeas; i++)
		{
			Measurements[0][i] = MCSweeps(SpinConf, table, L,T, NSweep);
			make_measurement(SpinConf, Measurements, i, L);
		}
		calculate_averages(Averages, Nmeas, Measurements);
		calculate_error_bar(Averages, error, Nmeas);
		CalculatedValues[0][k] = Averages[0];
		CalculatedValues[1][k] = error[0];
		CalculatedValues[2][k] = calcX(Averages[4], T, L);
		CalculatedValues[3][k] = Xerror(error[3], T, L);
		CalculatedValues[4][k] = CvCalc(Averages[0], Averages[1], T, L);
		CalculatedValues[5][k] = CvError(error[0], error[1], Averages[0], T, L);
		CalculatedValues[6][k] = T;  	
		CalculatedValues[7][k] = calcg(Averages[5],Averages[4]);
		cout << k << "  M4 = " << Averages[5] << endl;
		cout << k << "  M2 = " << Averages[4] << endl;	
	}
	Printout(Averages, error, L);
	file_output(Measurements[3], Nmeas);
	calc_values_file_output(CalculatedValues, Tpoints);
}

double calcg(double M4, double M2)
{
	double g = 1.0 / 2.0 * (3.0 - M4 / (M2*M2));
	return (g);
}
