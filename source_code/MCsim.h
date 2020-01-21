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
typedef std::vector<int> intRow;
typedef std::vector<Row> Column;
typedef std::vector<intRow> intMat; 
typedef std::chrono::high_resolution_clock myclock;

struct vec
{
    double x, y;
};

class MCsim 
{
    private:
      intMat SpinConf;   //2D grid for spin conformation
      intRow F, B;       //forward and back vectors to facilitate looped boundaries
      double delta_E;    //energy difference between old and proposed states
	  double E_tot;      //total energy of system
	  int changepos[2];  //coordinates of position to change 
	  int old_s;         //old state of changepos, used if reversion is needed
      
	  //simulation parameters
	  double T;          //temperature of the system
      int L;             //lattice size of the system
      int SweepSize;     //number of iterations in one MCS

	  //tables to hold re-used computations
      double ExpTable[5];     
      double SineTable[3];    
      double CosTable[3];     
       
	  //rng seeded, and int distribution initialized in constructor
	  std::mt19937 rng;      
	  std::uniform_real_distribution<double> dis{0.0, 1.0};
      std::uniform_int_distribution<int> gen0toL;
      
      //functions to initialize simulation 
      void initialize_conf();  	  
      void initialize_rand();   
      void initialize_energy();  
      void precalculate_table();  
      
	  double energy_contribution(int s1, int s2); 
	  int random_state(); //randomly select a state for a node in the grid      

	  //utility functions
      void change_conf();          //assign a random state to a random position
      void calc_delta_E();         //calculate the associated change in energy
      void check_pass();           //check to see if new conformation passes
      void revert();               //revert to old conformation on failure
      void advance();              //update E_tot in the case of a pass

    public:
      MCsim(int L, double T);    
      void MCSweeps(int Nmcs);   //performs Ncms MC sweeps on the Spin Conformation
      vec Magnetization();       //returns magnetization of Spin Conformation
      double Energy();           //returns Energy of Spin Conformation
      
};
