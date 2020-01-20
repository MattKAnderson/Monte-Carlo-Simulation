/* Source files for the MCsim class functions. An MCsim object is a single
 * 2D, q=3 clock model Monte Carlo simulation */

#include <iostream>
#include "MCsim.h"

MCsim::MCsim(int Length, double Temp)
{
    L = Length;
    T = Temp;
    std::uniform_int_distribution<int> temp(0, (L-1));
    gen0toL.param(temp.param());
    SweepSize = L * L;
    initialize_rand();
    initialize_conf();
    initialize_energy();
    precalculate_table(); 
}

void MCsim::MCSweeps(int Nmcs)  
{
    for(int i=0; i<Nmcs; i++)
    {
      for(int j=0; j<SweepSize; j++)
      {
        change_conf();
        calc_delta_E();
        check_pass();
      }
    }
}
      
vec MCsim::Magnetization()
{
    vec Mag;
    Mag.x = 0.0;
    Mag.y = 0.0;
    for(int i=0; i<L; i++)
    {
      for(int j=0; j<L; j++)
      { 
        Mag.x += CosTable[SpinConf[i][j] - 1];
        Mag.y += SineTable[SpinConf[i][j] - 1];
      }
    }
    return Mag;
}

double MCsim::Energy()
{
    return(E_tot);
}

void MCsim::initialize_conf()
{
    SpinConf.resize(L);
    F.resize(L);
    B.resize(L);
    for(int i=0; i<L; i++)
    {
      SpinConf[i].resize(L);
      for(int j=0; j<L; j++)
      {
        SpinConf[i][j] = random_state(); 
      }
    }
    for(int i=0; i<(L-1); i++)  
    {
      F[i] = i + 1;
      B[L-i-1] = L-i-2;
    }
    F[L-1] = 0;                 
    B[0] = L - 1;
}

void MCsim::initialize_rand()
{
    myclock::time_point begin = myclock::now();
    myclock::duration dur = myclock::now() - begin;
    unsigned rngseed = dur.count();
    rng.seed(rngseed);
}

void MCsim::initialize_energy()
{
    E_tot = 0;
    for(int i=0; i<L; i++)
    {
      for(int j=0; j<L; j++)
      {
        E_tot += energy_contribution(SpinConf[i][j], SpinConf[F[i]][j]);
        E_tot += energy_contribution(SpinConf[i][j], SpinConf[i][F[j]]);
      }
    } 
}

void MCsim::precalculate_table()
{
    for (int i=0; i<5; i++)
    {
      ExpTable[i] = exp(-1.5 * i / T);
    }
    for (int i=0; i<3; i++)
    {
      SineTable[i] = sin((i + 1.0) * 2 * M_PI / 3.0);
      CosTable[i] = cos((i + 1.0) * 2 * M_PI / 3.0);
    }
}

int MCsim::random_state()
{
    double random_pick = dis(rng);
    if (random_pick <= (double)1.0/(double)3.0)
    {
      return (1);
    }
    else if (random_pick <= (double)2.0/(double)3.0)
    {
       return (2);
    }
    else
    {
      return (3);
    }
}

double MCsim::energy_contribution(int s1, int s2)
{
    if (s1 == s2)
    {
      return(-1.0);
    }
    else
    {
      return(0.5);
    }
} 

void MCsim::change_conf()
{
    //random position to change
    changepos[0] = gen0toL(rng);
    changepos[1] = gen0toL(rng);
    
	//save the old state and change position to a new state
    old_s = SpinConf[changepos[0]][changepos[1]]; 
    SpinConf[changepos[0]][changepos[1]] = random_state();
} 

void MCsim::calc_delta_E()
{ 
    double contribution[2] = {0.0, 0.0};                          
    int state[2] = {old_s, SpinConf[changepos[0]][changepos[1]]};    

    //summing nearest neighbour interactions for the old and new state
    for (int i=0; i<2; i++)
    {
      contribution[i] += energy_contribution(state[i], 
	                        SpinConf[F[changepos[0]]][changepos[1]]); 
      contribution[i] += energy_contribution(state[i], 
	                        SpinConf[B[changepos[0]]][changepos[1]]); 
      contribution[i] += energy_contribution(state[i], 
	                        SpinConf[changepos[0]][F[changepos[1]]]); 
      contribution[i] += energy_contribution(state[i], 
	                        SpinConf[changepos[0]][B[changepos[1]]]); 
    }

    //delta E is the difference between the new and the old contribution
    delta_E = contribution[1] - contribution[0];
}

void MCsim::check_pass()
{
    if (delta_E <= 0)
    {
      advance();
    }
    else
    {
      int index = delta_E * 2.0 / 3.0;
      if (ExpTable[index] > dis(rng))
      {
        advance();
      }
      else
      {
        revert();
      }
    } 
}

void MCsim::advance()
{
  E_tot += delta_E;
}

void MCsim::revert()
{
  SpinConf[changepos[0]][changepos[1]] = old_s;
}
