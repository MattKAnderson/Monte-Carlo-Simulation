/* source code file for the Monte Carlo functions used for Exercise 4 of 
 * homework 2 of Physics 4G03 Fall 2018. */

#include "MC_functions.h"

void Exponent_table_calculator(int L, double T, std::vector <double> &table)
{
	int entries = 5;                              //change in Energy can take up to 5 values
	table.resize(entries);                
	for (int i = 0; i < entries; i++)
	{
		table[i] = exp(-(i*4 - 8)/T);
	}
}

double getRandom()
{
	static std::mt19937 generator(getSeed());
	static std::uniform_real_distribution <double> dist(0, 1.0);
	return (dist(generator));
}

unsigned getSeed()
{
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point begin = myclock::now();
	myclock::duration d = myclock::now() - begin;
	return (d.count());
}


int MCSweeps(std::vector < std::vector < int > > &SpinConf, const std::vector <double> &table, int L,double T, int Nsweeps)
{
	int posX, posY;
	int revertSave;
	int newEnergy, EnergyDif;
	int currentEnergy = Energy(SpinConf, L);
	int Nsteps = Nsweeps * L * L;
	double changeProb;
	//std::cout << "Current energy = " << currentEnergy << "\n";	
	for (int i = 0; i < Nsteps; i++)
	{
		//choosing position to flip spin at
		posX = std::round((L - 1) * getRandom());
		posY = std::round((L - 1) * getRandom());
		//std::cout << "posX = " << posX << "\n";
		//std::cout << "posY = " << posY << "\n";
		//flipping spin
		if (SpinConf[posX][posY] == 1)
		{
			SpinConf[posX][posY] = -1;
			revertSave = 1;
		}
		else if (SpinConf[posX][posY] == -1)
		{
			SpinConf[posX][posY] = 1;
			revertSave = -1;
		}
		else
		{
			std::cout << "Error in Spin configurations, spin other than 1 or -1 encountered\n";
		}

		//calculating energy of new spin
		newEnergy = EnergyChange(SpinConf, L, posX, posY, currentEnergy);
		//newEnergy = Energy(SpinConf, L);
		EnergyDif = newEnergy - currentEnergy;
		//std::cout << "Energy change is: " << EnergyDif << std::endl;
		//std::cout << "newEnergy is: " << newEnergy << std::endl;
	
		//getting change probability from pre-calculated table
		changeProb = table[(EnergyDif + 8)/4];
		//changeProb = exp(-EnergyDif/T);
		//checking if pass
		if (changeProb > getRandom())
		{
			currentEnergy = newEnergy;
		}
		else
		{
			SpinConf[posX][posY] = revertSave;
		}
	}	
	return (currentEnergy);
}

int Energy(const std::vector < std::vector <int> > &SpinConf, int L)
{
	int eSum = 0;
	for (int i = 0; i < L; i++)
	{
		//handling boundary condition of all row and columns
		eSum -= SpinConf[L - 1][i] * SpinConf[0][i];  
		eSum -= SpinConf[i][L - 1] * SpinConf[i][0];
		
		//adding all nearest neighbour pairs from rows and columns
		for (int j = 0; j < (L - 1); j++)
		{
			eSum -= SpinConf[j][i] * SpinConf[j + 1][i];   
			eSum -= SpinConf[i][j] * SpinConf[i][j + 1];
		}
	}
	return(eSum);
}

int EnergyChange(const std::vector < std::vector <int> > &SpinConf, int L, int posX, int posY, int cEnergy)
{
	int newContribution;
	newContribution = -SpinConf[posX][posY] * (SpinConf[posX][LB[posY]] + SpinConf[posX][LF[posY]] + SpinConf[LB[posX]][posY] + SpinConf[LF[posX]][posY]);
	//std::cout << "energy change is: " << 2*newContribution << std::endl;
	return(cEnergy + 2 * newContribution);
}

int magnetization(const std::vector < std::vector <int> > &SpinConf, int L)
{
	int magSum = 0;
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{
			magSum += SpinConf[i][j];	
		}
	}
	return(magSum);
}
