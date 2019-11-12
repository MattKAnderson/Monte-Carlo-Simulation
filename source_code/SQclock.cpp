/* main source code file for monte carlo q=3 clock model. written by matt 
 * anderson in december 2018, as part of the bonus assignment for phys 4g03 at 
 * mcmaster university */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <omp.h>
#include "MCsim.h"
#include "dataFunctions.h"
using namespace std;

typedef vector<double> Row;
typedef vector<Row> Column;

void To_file(Row &Observable_Avgs, Row &error_bars, double T, int L);
void prep_file();

int main()
{
    //simulation variable declarations
    int L, Nmeas, Nstep, NWarmup, Tpoints, nthreads;
    double T, Tmin, Tmax, Tstep;
    Row Observable_Avgs(7);   //E, E^2, E^4, M^2, M^4, X, Cv
    Row error_bars(5);        //E, E^2, M^2, X, Cv
    vec Magnetization;
    Column Observables;   
    Observables.resize(5);

    //getting variable values from terminal prompts
    cout << "Enter a lattice size: ";
    cin >> L;
    cout << "Enter the max temperature: ";
    cin >> Tmax;
    cout << "Enter the min temperature: ";
    cin >> Tmin;
    cout << "Enter the number of temperature points: ";
    cin >> Tpoints;
    cout << "Enter the number of warm-up MC sweeps: ";
    cin >> NWarmup;
    cout << "Enter the number of MC sweeps between measurements: ";
    cin >> Nstep;
    cout << "Enter the number of measurements to take: ";
    cin >> Nmeas;
    cout << "Enter the number of threads to use for calculation: ";
    cin >> nthreads;
    
    //resizing the rows of Observables for the number of measurements to be taken
    for (int i=0; i<5; i++)
    {
      Observables[i].resize(Nmeas);
    }

    //determining Tstep
    Tstep = (Tmax - Tmin) / (Tpoints - 1);
    
    //clearing all contents in file to print output to
    prep_file();
    //iterating simulation for Tpoints number of temperatures
    for (int i=0; i<Tpoints; i++)
    {
      T = Tmin + Tstep * i;

      //initializing the Monte Carlo simulation & doing warm-up sweeps
      MCsim sim1(L, T);
      sim1.MCSweeps(NWarmup);
     
      //iterating simulation and taking measurements
      for (int i=0; i<Nmeas; i++)
      {
        sim1.MCSweeps(Nstep);
        Observables[0][i] = sim1.Energy();
        Magnetization = sim1.Magnetization();
        Calc_Observables(Observables, Magnetization, i);
      }
      //processing recorded observations
      calc_averages(Observables, Observable_Avgs);
      error_bar_calc(Observable_Avgs, error_bars, Nmeas);
      X_calc(Observable_Avgs, error_bars, T, L);
      Cv_calc(Observable_Avgs, error_bars, T, L); 

      //sending data to file
      To_file(Observable_Avgs, error_bars, T, L);
    }
    
    cout << "Simulaton complete. Data output to file: 2Dclock.dat" << endl;  
    //sending results to terminal
    /*
    cout << "Results are:" << endl;
    cout << "<E> =\t" << Observable_Avgs[0] << endl;
    cout << "<E^2> =\t" << Observable_Avgs[1] << endl;
    cout << "<E^4> =\t" << Observable_Avgs[2] << endl;
    cout << "<M^2> =\t" << Observable_Avgs[3] << endl;
    cout << "<M^4> =\t" << Observable_Avgs[4] << endl;
    cout << "X =\t" << Observable_Avgs[5] << endl;
    cout << "Cv =\t" << Observable_Avgs[6] << endl; */
}

void To_file(Row &Observable_Avgs, Row &error_bars, double T, int L)
{
    ofstream output;
    output.open("2Dclock.dat", ios_base::app);
    output << setprecision(5) << Observable_Avgs[0] << "\t\t" << error_bars[0] 
           << "\t\t" << Observable_Avgs[1] << "\t\t" << error_bars[1] << "\t\t" 
           << Observable_Avgs[2] << "\t\t" << Observable_Avgs[3] << "\t\t" << error_bars[2]
           << Observable_Avgs[4] << "\t\t" << Observable_Avgs[5] << "\t\t" << error_bars[3]
           << Observable_Avgs[6] << "\t\t" << error_bars[4] << "\t\t" << T << "\t\t"
           << L << endl;

            /*"\t\t" << Observable_Avgs[7]
           << "\t\t" << error_bars[5] << "\t\t" << T << "\t\t" << L << endl;*/
}

void prep_file()
{
    ofstream output;
    output.open("2Dclock.dat");
    output << "2D, q=3 clock model data. Columns are <E>, error <E>, <E^2>, error <E^2>, <E^4>, <M^2>, error <M^2>, X, Cv, T, L" << endl;
    output.close();
}
