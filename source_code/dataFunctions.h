/* Header file for the data functions used in the SQclock program
 * Written by Matt Anderson Dec 2018  */

#include "MCsim.h"
void Calc_Observables(Column &Observables, vec Mag, int i);
void calc_averages(Column &Observables, Row &Observable_Avgs);
void error_bar_calc(Row &Obervable_Avgs, Row &error_bars, int Nmeas);
void X_calc(Row &Avgs, Row &errorbars, double T, int L);
void Cv_calc(Row &Avgs, Row &errorbars, double T, int L);
