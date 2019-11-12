/* main source code file for monte carlo q=3 clock model. Parallelization
 * added to speed up code run time. written by Matt Anderson in december 2018,
 * as part of the bonus assignment for phys 4g03 at McMaster University */

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
void G_calc(Row &Observable_Avgs);

int main()
{
    //simulation variable declarations
    int L, j, Nmeas, Nstep, NWarmup, Tpoints, nthreads;
    double Tmin, Tmax, Tstep;
    Row T;
    Column Observable_Avgs;   //E, E^2, E^4, M^2, M^4, X, Cv
    Column error_bars;        //E, E^2, M^2, X, Cv
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

    //each row is a different measurement, resizing the rows to have the necessary number of components
    T.resize(Tpoints);
    Observable_Avgs.resize(Tpoints);
    error_bars.resize(Tpoints);
    for (int i=0; i<Tpoints; i++)
    {
      Observable_Avgs[i].resize(8);
      error_bars[i].resize(5);
    }

    //determining Tstep
    Tstep = (Tmax - Tmin) / (Tpoints - 1);
    
    //clearing all contents in file to print output to
    prep_file();

    //setting number of threads and starting the parallel region
    omp_set_num_threads(nthreads);
    #pragma omp parallel private(Magnetization) firstprivate(Observables) shared(Observable_Avgs, error_bars, T)
    {
      //iterating simulation for Tpoints number of temperatures with multiple threads
      #pragma omp for nowait private(j) 
      for (int i=0; i<Tpoints; i++)
      {
        T[i] = Tmin + Tstep * i;
        //initializing the Monte Carlo simulation & doing warm-up sweeps
        MCsim sim(L, T[i]);
        sim.MCSweeps(NWarmup);
        //cout << "Simulation object successfully created" << endl;
        //iterating simulation and taking measurements
        for (j=0; j<Nmeas; j++)
        {
          sim.MCSweeps(Nstep);
          //cout << "did an MCSweep" << endl;
          Observables[0][j] = sim.Energy();
          //cout << "saved to observables" << endl;
          Magnetization = sim.Magnetization();
          Calc_Observables(Observables, Magnetization, j);
        }
        //cout << "Finished the simulation" <<endl;
        //processing recorded observations
        calc_averages(Observables, Observable_Avgs[i]);
        error_bar_calc(Observable_Avgs[i], error_bars[i], Nmeas);
        X_calc(Observable_Avgs[i], error_bars[i], T[i], L);
        Cv_calc(Observable_Avgs[i], error_bars[i], T[i], L);
        G_calc(Observable_Avgs[i]); 
      }
    }

    //sending calculated data to file
    for (int i=0; i<Tpoints; i++)
    {
      To_file(Observable_Avgs[i], error_bars[i], T[i], L); 
    }
    cout << "Simulaton complete. Data output to file: 2Dclock.dat" << endl;  
    //sending results to terminal
    
    cout << "Results are:" << endl;
    cout << "For T =   "  << T[Tpoints / 2] << endl;
    cout << "<E> =     " << Observable_Avgs[Tpoints / 2][0] << "\t  error = " << error_bars[Tpoints / 2][0] << endl;
    cout << "<E^2> =   " << Observable_Avgs[Tpoints / 2][1] << "\t  error = " << error_bars[Tpoints / 2][1] << endl;
    cout << "<E^4> =   " << Observable_Avgs[Tpoints / 2][2] << endl;
    cout << "<M^2> =   " << Observable_Avgs[Tpoints / 2][3] << "\t  error = " << error_bars[Tpoints / 2][2] << endl;
    cout << "<M^4> =   " << Observable_Avgs[Tpoints / 2][4] << endl;
    cout << "X =       " << Observable_Avgs[Tpoints / 2][5] << "\t  error = " << error_bars[Tpoints / 2][3] << endl;
    cout << "Cv =      " << Observable_Avgs[Tpoints / 2][6] << "\t  error = " << error_bars[Tpoints / 2][4] << endl;
    cout << "G(T,L) =  " << Observable_Avgs[Tpoints / 2][7] << endl; 
}

void To_file(Row &Observable_Avgs, Row &error_bars, double T, int L)
{
    ofstream output;
    output.open("2Dclock.dat", ios_base::app);
    output << fixed << setprecision(5) << Observable_Avgs[0] << "  " << error_bars[0] 
           << "  " << Observable_Avgs[1] << "  " << error_bars[1] << "  " 
           << Observable_Avgs[2] << "  " << Observable_Avgs[3] << "  " << error_bars[2]
           << "  " << Observable_Avgs[4] << "  " << Observable_Avgs[5] << "  " 
           << error_bars[3] << "  " << Observable_Avgs[6] << "  " << error_bars[4] 
           << "  " << Observable_Avgs[7] << "  " << T << "  " << L << endl;

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

void G_calc(Row &Observable_Avgs)
{
    Observable_Avgs[7] = Observable_Avgs[6] / (Observable_Avgs[5] * Observable_Avgs[5]);
    Observable_Avgs[7] = 3.0 - Observable_Avgs[7];
    Observable_Avgs[7] /= 2.0;
}
