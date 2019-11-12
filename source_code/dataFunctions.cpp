/* Source file for functions that do operations on the data produced from the
 * SQclock MC simulation program. Written by Matt Anderson Dec 2018  */

#include <iostream>
#include "MCsim.h"
#include "dataFunctions.h"

void Calc_Observables(Column &Observables, vec Mag, int i)
{
    Observables[1][i] = Observables[0][i] * Observables[0][i];
    Observables[2][i] = Observables[1][i] * Observables[1][i];
    Observables[3][i] = sqrt(Mag.x * Mag.x + Mag.y * Mag.y);
    //std::cout << "value of M^2: " << Observables[3][i] << "   for measurement: " << i << std::endl;
    Observables[4][i] = Observables[3][i] * Observables[3][i];
}

void calc_averages(Column &Observables, Row &Observable_Avgs)
{
    for (int i=0; i<5; i++)
    {
      Observable_Avgs[i] = 0.0;
      for (int j=0; j<Observables[i].size(); j++)
      {
        Observable_Avgs[i] += Observables[i][j];
      }
      Observable_Avgs[i] /= Observables[i].size();
    }
}

void error_bar_calc(Row &Observable_Avgs, Row &error_bars, int Nmeas)
{
    //error bar on average energy:
    error_bars[0] = sqrt((Observable_Avgs[1] - Observable_Avgs[0] * Observable_Avgs[0]) / (Nmeas - 1));

    //error bar on average energy squared:
    error_bars[1] = sqrt((Observable_Avgs[2] - Observable_Avgs[1] * Observable_Avgs[1]) / (Nmeas - 1));

    //error bar on average Magnetization squared:
    error_bars[2] = sqrt((Observable_Avgs[4] - Observable_Avgs[3] * Observable_Avgs[3]) / (Nmeas - 1));
}

void X_calc(Row &Avgs, Row &errorbars, double T, int L)
{
    Avgs[5] = Avgs[3] / (T * L * L);    
    errorbars[3] = errorbars[2] / (T * L * L);
}

void Cv_calc(Row &Avgs, Row &errorbars, double T, int L)
{
    double error = (1 / (L * L * T * T)) * (1 / (L * L * T * T));
    double relError1 = errorbars[0] / Avgs[0];
    double relError2 = errorbars[1] / Avgs[1];
    
    relError1 *= relError1;
    relError2 *= relError2;
    error = 4*Avgs[1]*Avgs[1]*error*relError1 + error*relError2;
    errorbars[4] = sqrt(error);
    
    Avgs[6] = (Avgs[1] - (Avgs[0] * Avgs[0])) / (L * L * T * T);
}


