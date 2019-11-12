/* Header file for the MCsim class. An MCsim object is a single 2D, q=3 clock
 * model Monte Carlo simulation */

#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#define _USE_MATH_DEFINES

typedef std::vector<double> Row;
typedef std::vector<int> iRow;
typedef std::vector<Row> Column;
typedef std::vector<iRow> iColumn; 
typedef std::chrono::high_resolution_clock myclock;

struct vec
{
    double x, y;
};
class MCsim 
{
    private:
      //variables
      iColumn SpinConf;        //2D vector grid to hold the spin conformation
      iRow F, B;              //forward and back vectors to shift indices by 1 with looped boundaries
      double delta_E, E_tot;  //Energy change between previous and current and total energy
      double T;               //temperature of the system
      double ExpTable[5];     //table to hold precalculated exp(delta_E) values for delta_E > 0
      double SineTable[3];    //table to hold possible sine(theta) values
      double CosTable[3];     //table to hold possible cosine(theta) values
      int L;                  //lattice size of the system
      int SweepSize;          //number of iterations in one MCS
      int changepos[2], old_s;//coordinates of position "being changed", and it's previous value
      std::mt19937 rng;       //MT random number generator (initialized in constructor)
      std::uniform_real_distribution<double> dis{0.0, 1.0};   //distribution for random number generator
      std::uniform_int_distribution<int> gen0toL;
      
      //functions
      void initialize_conf();      //called by constructor to initialize a conformation & F,B vectors
      void initialize_rand();      //called by constructor to initialize random num generator
      void initialize_energy();    //called by constructor to detemine initial energy 
      void precalculate_table();   //calculate the values for ExpTable, sineTable, cosTable
      
      int random_state();         //returns at random 1, 2, or 3. called by initialize_conf & change_conf
      double energy_contribution(int s1, int s2); //returns energy contribution of a single nearest neighbour pair 

      void change_conf();          //randomly selects a position to assign a new random state, s
      void calc_delta_E();         //calculate the associated change in energy
      void check_pass();           //check to see if new conformation passes
      void revert();               //revert to old conformation on failure
      void advance();              //update E_tot in the case of a pass

    public:
      //variables
      
      //functions
      MCsim(int L, double T);    //constructor
      void MCSweeps(int Nmcs);   //performs Ncms MC sweeps on the Spin Conformation
      vec Magnetization();       //calculate and return magnetization of Spin Conformation
      double Energy();           //calculates and return Energy of Spin Conformation
      
};
