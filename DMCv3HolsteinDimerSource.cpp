#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <time.h>
#include <cstdio>
#include <random>
#include <ctime>

using namespace std;

long double Potential(const long double& omega, const long double& position, const long double& F, const int& site)
{
	long double V;
	if (site == 1) //site A V +
	{
		V = 0.5*omega*omega*position*position + F*position; //V_1 = 1/2 * w^2 * x^2 + Fx
	}
	else if (site == 2) //site B V -
	{
		V = 0.5*omega*omega*position*position - F*position; //V_2 = 1/2 * w^2 * x^2 - Fx
	}
	return V;
}

void Initial_distribution(long double walkerpos[], const int& nowalkers,int& dist_type)
{
	normal_distribution<double> distribution(0.0, 1); // mean = 0.5, std = 0.1
	unsigned seed = time(NULL);
	default_random_engine generator(seed);



	if (dist_type == 1) // All set to zero
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 0;
		}
	}
	else if (dist_type == 2) //normally distributed with mean = 0, std = 1
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = distribution(generator);
		}
	}
	else if (dist_type == 3) // +/- 4 uniformly spaced grid
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 8 * (rand() / (double)RAND_MAX) - 4;
		}
	}
	else if (dist_type == 4) // all to the right
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 4;
		}
	}
	else if (dist_type == 5) // all to the left
	{
		for (int i = 0; i < nowalkers; ++i)
		{
			walkerpos[i] = -4;
		}
	}
}

void Initial_walkersite(long double walkersite[], long double& nowalkers)
{
	for (int i = 0; i < nowalkers; i++)
	{
		walkersite[i] = rand() % 2;
	}
}

void move_left(long double& position, long double& step)
{
	position = position + step;
}

void move_right(long double& position, long double& step)
{
	position = position - step;
}

void record_histogram(ofstream& Histogram, long double array1[], long double array2[], long double low_bin, long double high_bin, long double bin_width)
{
	double NUM_BINS = ((high_bin - low_bin) / bin_width) + 1;


	if (Histogram.is_open())
	{
		for (int i = 0; i < NUM_BINS; i++)
		{
			Histogram << low_bin + i*bin_width << "\t" << array1[i] << "\t" << array2[i] << endl;
		}

	}
	else
	{
		cout << "Error Recording Histogram File" << endl;
	}
}

double random_gen()
{
	double randm = rand() / (double)RAND_MAX;
	return randm;
}

void Information_file(int sim_num, ofstream& Info, ofstream& MATLAB, long double& Epi, long double& Epf, long double& Eps, long double& omega, long double& thop, const long double& N0, long double& T,int& dtn, long double dt_array[], long double weighting[])
{
	if (Info.is_open())
	{
		Info << "************* Simulation Number: " << sim_num;
		Info << "*************" << endl;
		Info << "Initial Energy Shift (Epi): ";
		Info << Epi << endl;
		MATLAB << Epi << endl;
		Info << "Final Energy Shift (Epf): ";
		Info << Epf << endl;
		MATLAB << Epf << endl;
		Info << "Energy Shift Step size (Eps): ";
		Info << Eps << endl;
		MATLAB << Eps << endl;
		Info << "Value of t: ";
		Info << thop << endl;
		MATLAB << thop << endl;
		Info << "Initial Number of walkers: ";
		Info << N0 << endl;
		MATLAB << N0 << endl;
		Info << "Total simulation time: ";
		Info << T << endl;
		MATLAB << T << endl;
		Info << "Number of time steps (dtn): ";
		Info << dtn << endl;
		MATLAB << dtn << endl;
		Info << "dt1: ";
		Info << dt_array[0] << endl;
		MATLAB << dt_array[0] << endl;
		Info << "dt2: ";
		Info << dt_array[1] << endl;
		MATLAB << dt_array[1] << endl;
		Info << "dt3: ";
		Info << dt_array[2] << endl;
		MATLAB << dt_array[2] << endl;
		Info << "dt4: ";
		Info << dt_array[3] << endl;
		MATLAB << dt_array[3] << endl;
		Info << "dt1_lim (% of t): ";
		Info << weighting[0] << endl;
		MATLAB << weighting[0] << endl;
		Info << "dt2_lim (% of t): ";
		Info << weighting[1] << endl;
		MATLAB << weighting[1] << endl;
		Info << "dt3_lim (% of t): ";
		Info << weighting[2] << endl;
		MATLAB << weighting[2] << endl;
		Info << endl;
	}
	else { cout << "Error recording Information file" << endl; }
}

void Record_Energy(ofstream& Energy, long double& time, long double& avev)
{
	Energy << time << "\t" << avev << endl;
}

void Record_Walkers(ofstream& Walkers, long double& time, long double& N_site1, long double& N_site2,long double& nowalkers)
{
	Walkers << time << "\t" << N_site1 << "\t" << N_site2 << "\t" << nowalkers << endl;
}

void User_Input(long double& thop, long double& omega, long double& Epi, long double& Epf, long double& Eps, long double& N0, long double& T,int& dtn, long double dt_array[], long double weighting[], long double& convergence_time, int& dist_type)
{
	cout << "Please give the hopping term (t): ";
	cin >> thop;
	cout << endl;

	cout << "Please give the value of omega (w): ";
	cin >> omega;
	cout << endl;

	cout << "Please give the Initial Polaron Energy shift (Epi): ";
	cin >> Epi;
	cout << endl;

	cout << "Please give the Final Polaron Energy Shift (Epf): ";
	cin >> Epf;
	cout << endl;

	cout << "Please give the step size to final energy shift (Eps): ";
	cin >> Eps;
	cout << endl;

	cout << "Please set the initial number of walkers (N0): ";
	cin >> N0;
	cout << endl;

	cout << "Please set the total simulation time (T): ";
	cin >> T;
	cout << endl;

	cout << "Please give the number of time steps (max 4): ";
	cin >> dtn;
	cout << endl;

	if (dtn == 1)
	{
		cout << "dt1: ";
		cin >> dt_array[0];
		cout << endl;
	}
	else if (dtn == 2)
	{
		cout << "dt1: ";
		cin >> dt_array[0];
		cout << endl;
		cout << "dt1 limit (% of t): ";
		cin >> weighting[0];

		cout << "dt2: ";
		cin >> dt_array[1];
		cout << endl;
	}
	else if (dtn == 3)
	{
		cout << "dt1: ";
		cin >> dt_array[0];
		cout << endl;
		cout << "dt1 limit (% of t): ";
		cin >> weighting[0];
		cout << endl;

		cout << "dt2: ";
		cin >> dt_array[1];
		cout << endl;
		cout << "dt2 limit (% of t): ";
		cin >> weighting[1];
		cout << endl;

		cout << "dt3: ";
		cin >> dt_array[2];
		cout << endl;

	}
	else if (dtn == 4)
	{
		cout << "dt1: ";
		cin >> dt_array[0];
		cout << endl;
		cout << "dt1 limit (% of t): ";
		cin >> weighting[0];
		cout << endl;

		cout << "dt2: ";
		cin >> dt_array[1];
		cout << endl;
		cout << "dt2 limit (% of t): ";
		cin >> weighting[1];

		cout << "dt3: ";
		cin >> dt_array[2];
		cout << endl;
		cout << "dt3 limit (% of t): ";
		cin >> weighting[2];

		cout << "dt4: ";
		cin >> dt_array[3];
		cout << endl;

	}

	cout << "Please select the intial distribution type" << endl;
	cout << "\t1 - All walkers at zero" << endl;
	cout << "\t2 - Normal distribution (mean = 0, std = 1)" << endl;
	cout << "\t3 - +/- 4 grid with uniform  distribution" << endl;
	cout << "\t4 - All walkers to the right (+4)" << endl;
	cout << "\t5 - All walkers to the left (-4)" << endl;
	cin >> dist_type;

	cout << "Expected convergence time, for histogram, (% of t) : ";
	cin >> convergence_time;
	cout << endl;
}

void progress_bar(int sim_num, long double total_time, long double current_time)
{
	//int progress = static_cast<int>(current_time / total_time);

	double progress = current_time / total_time;

	int barWidth = 70;
	cout << "Simulation # " << sim_num;
	cout << "[";
	int pos = barWidth * progress;
	for (int i = 0; i < barWidth; ++i)
	{
		if (i < pos) { cout << "="; }
		else if (i == pos) { cout << ">"; }
		else { cout << " "; }
	}
	cout << "] " << int(progress * 100.0) << " %\r";
	cout.flush();
}


int main()
{
	// *******************************************************************************************
	//						VARIABLES INITIALIZATION
	// *******************************************************************************************
	clock_t begin, end;
	srand(time(NULL));//seed for random numbers
	const int pos_size = 20000000; //maximum size of walkerposition array
	
	long double* walkerpos = new long double[pos_size]; //walkerposition pointer to array of double, size (pos_size)
	long double* walkersite = new long double[pos_size]; //walker site pointer to array of double, size (hist_size)

	// ----------- Initialization of DMC parameters
	long double t; // t = total simulation time
	long double avev = 0; // avev = average potential
	long double aveE = 0; // average energy set to 0
	long double refpot = 0; // refpot = reference potential
	long double nowalkersite1; //site A for walkers
	long double nowalkersite2; //site B for walkers
	int dist_type; //Distribution type requested from the user (constrained to 3 for initial testing)
	long double convergence_time;
	long double c; //counter used to count the death or birth of walkers (maybe should be moved within scope) 
	long double N0; // initial number of walkers requested from the user
	  
	// ----------- Initialization of Holstein-Dimer parameters
	// Polaron Energy shifts Ep (initial, final and step)
	long double Epi; //Initial Ep
	long double Epf; //Final Ep
	long double Eps; //Epi to Epf step size
	long double Hop_term; //Polaron hopping term (noted as t in HD model)
	long double omega; //frequency omega 
	long double Vx;  // potential
	const int site1 = 1; //choice parameter for potential calculation site A
	const int site2 = 2; //choice parameter for potential calculation site B
	
	// ---------- Initialization of  histogram parameters
	const int hist_size = 10000; //maximum size of histogram array
	long double histoa[hist_size] = { 0 }; //initialization to 0 & declaration of histogram A 
	long double histob[hist_size] = { 0 }; //initialization to 0 & declaration of histogram B
	const long double low_bin = -10; //low bin of histogram
	const long double high_bin = 10; //high bin of histogram
	const long double bin_width = 0.1; //bin width of histogram
	const long double NUM_BINS = ((high_bin - low_bin) / bin_width) + 1; //given low, high and bin width = NUM_BINS 

	// ---------- Initialization of time step variables
	long double dt_array[4] = { 0 };
	
	long double weighting[4] = { 0 };
	int dtn = 0;

// ---------------------------------------------------------------------------------------------------
//                          USER INTERFACE VIA COMMAND WINDOW  
// ---------------------------------------------------------------------------------------------------

	cout << "Diffusion Monte Carlo Algorithm - Holstein-Dimer" << endl;
	cout << "\t\t Version 3" << endl << endl;
	User_Input(Hop_term, omega, Epi, Epf, Eps, N0, t, dtn, dt_array,weighting,convergence_time,dist_type);  //calls on the user input function

	double Epn = (Epf - Epi) / Eps; // determines the number of simulations to run

	// --------- Assignment of timestep variables and timestep limits
	long double dt1 = dt_array[0];
	long double dt2 = dt_array[1];
	long double dt3 = dt_array[2];
	long double dt4 = dt_array[3];
	long double dt1_lim = weighting[0] / 100;
	long double dt2_lim = weighting[1] / 100;
	long double dt3_lim = weighting[2] / 100;
	long double dt4_lim = weighting[3] / 100;
	
	long double timestep = dt1;
	
	// --------- CREATION OF OUTPUT TEXT FILES --------------------------------
	// -- construct of date and time information 
	time_t now;
	struct tm nowLocal;
	now = time(NULL);
	nowLocal = *localtime(&now);
	int year = nowLocal.tm_year + 1900;
	int month = nowLocal.tm_mon + 1;
	int day = nowLocal.tm_mday;
	int hour = nowLocal.tm_hour;
	int min = nowLocal.tm_min;
	int sec = nowLocal.tm_sec;
	// constructs date-specifing string 
	char datespecifier[64];
	sprintf(datespecifier, "%d-%d-%d-%d-%d-%d", year, month, day, hour, min, sec);
	// constructs filename for each output file
	char filename_E[64];
	sprintf(filename_E, "Archive/%s-Energies.txt", datespecifier);

	char filename_W[64];
	sprintf(filename_W, "Archive/%s-walkers.txt", datespecifier);

	char filename_AVE[64];
	sprintf(filename_AVE, "Archive/%s-Average_Energy.txt", datespecifier);

	char filename_Etest[64];
	sprintf(filename_Etest, "Archive/%s-Energy_test.txt", datespecifier);

	char filename_CPU[64];
	sprintf(filename_CPU, "Archive/%s-performance.txt", datespecifier);

	char filename_H[64];
	sprintf(filename_H, "Archive/%s-Wavefunction.txt", datespecifier);

	char filename_O[64];
	sprintf(filename_O, "Archive/%s-MATLAB.txt", datespecifier);

	char filename_I[64];
	sprintf(filename_I, "Archive/%s-Info.txt", datespecifier);

	char filename_HD[64];
	sprintf(filename_HD, "Archive/%s-HD.txt", datespecifier);
	
	//constructs filestream and creates the file.
	ofstream Energy(filename_E);
	ofstream walkers(filename_W);
	ofstream Average_Energy(filename_AVE);
	ofstream Energy_test(filename_Etest);
	ofstream Histogram(filename_H);
	ofstream Info(filename_I);
	ofstream MATLAB(filename_O);
	ofstream performance(filename_CPU);
	ofstream HD(filename_HD);

	for (int l = 0; l < (Epn + 1); l++) //loop used to iterate through all the coupling strengths 
	{
		long double Ep = Epi + l*Eps;
		const double nowalkers_initial = N0; //fixes the inital number of walkers to const variable 
		long double nowalkers = N0;
		
		const long double groundstate1 = -Ep + (0.5*omega) - (Hop_term*Hop_term / (4 * Ep)); //first E0 approximation
		const long double groundstate2 = -Ep + (0.5*omega) - Hop_term*exp(-2 * Ep / omega) - (Hop_term*Hop_term / (4 * Ep)); //second E0 approximation
		const long double F = omega*sqrt(2 * Ep); // F variable
		const long double sitea = F / (omega*omega); // local minimum of potential
		HD << Ep << "\t" << groundstate1 << "\t" << groundstate2 << "\t" << F << "\t" << sitea << endl;
		

		// ----------- INITIALIZE WALKER DISTRIBUTION --------------------------------
		nowalkersite1 = nowalkers_initial; //walkers set to one side
		nowalkersite2 = 0; //Walkers 

		Initial_distribution(walkerpos, nowalkers_initial, dist_type); //Generates the initial distribution
		Initial_walkersite(walkersite, nowalkers); //generates walkers site 

		long double timestep = dt1;
		long double Phop = Hop_term*timestep; //probability of hopping given a time step
		double steplength = sqrt(3 * dt1); //calculates the steplength of the walkers
		// *******************************************************************************************
		//						Simulation Begins Here
		// *******************************************************************************************
		for (long double time = 0; time < t; time = time + timestep)//no of iterations
		{

			progress_bar(l+1, t, time);

			if (dtn == 1)
			{
				timestep = dt1;
			}
			else if (dtn == 2)
			{
				if (time > dt1_lim*t)
				{
					timestep = dt2;
				}
			}
			else if (dtn == 3)
			{
				if (time < dt1_lim*t)
				{
					timestep = dt1;
				}
				else if ((time > dt1_lim*t) && (time < dt2_lim*t))
				{
					timestep = dt2;
				}
				else if (time > dt2_lim*t)
				{
					timestep = dt3;
				}
			}
			else if (dtn == 4)
			{
				if (time < dt1_lim*t)
				{
					timestep = dt1;
				}
				else if ((time > dt1_lim*t) && (time < dt2_lim*t))
				{
					timestep = dt2;
				}
				else if ((time > dt2_lim*t) && (time < dt3_lim*t))
				{
					timestep = dt3;
				}
				else if (time > dt3_lim*t)
				{
					timestep = dt4;
				}
			}

			Phop = Hop_term * timestep;
			steplength = sqrt(3 * timestep);

			Record_Walkers(walkers, time, nowalkersite1, nowalkersite2, nowalkers); //records the number of walkers
			nowalkersite1 = 0;
			nowalkersite2 = 0;
			long double* walkersite1 = new long double[pos_size]; // pointer to copy array initialization
			long double* walkerpos1 = new long double[pos_size];//Initialization oh second walker postion array used later
			refpot = 0; //reference potential reset to zero
			avev = 0; //Average potential reset to zero

					  //---------------- MOVES WALKERS ------------------------------------------------- 
			for (int i = 0; i < nowalkers; i++)//moves walkers left or right
			{
				long double step = steplength*random_gen(); //random step generation
				int randn = rand() % 2; //random movement left or right

				if (randn == 1) //Move walker to the left by step
				{
					move_left(walkerpos[i], step);
				}
				else if (randn == 0) //Move walker to the right by step
				{
					move_right(walkerpos[i], step);
				}

				// ------ Potential recalculation 
				if (walkersite[i] == 0) //Potential if walker is in site A (site1)
				{
					avev = avev + (Potential(omega, walkerpos[i], F, site1)) / nowalkers;
				}
				else if (walkersite[i] == 1) //Potential if walker is in site B (site2)
				{
					avev = avev + (Potential(omega, walkerpos[i], F, site2)) / nowalkers;
				}
			}

			refpot = avev - Hop_term - ((nowalkers / nowalkers_initial) - 1) / timestep;//calculation of reference potential
			int j = 0;

			c = 0; //walker change

			//--------------- WALKER POPULATION CONTROL --------------------------------
			for (int i = 0; i < nowalkers; i++) // birth and death loop
			{
				if (walkersite[i] == 0) //calculates potential depending on the site 
				{
					Vx = Potential(omega, walkerpos[i], F, site1);
				}
				else if (walkersite[i] == 1)
				{
					Vx = Potential(omega, walkerpos[i], F, site2);
				}

				double a;//term used to move walker positon between arrays, reset to zero each time
				double b;//term used to move walkers site between arrays, reset to zero each time
				long double k = (Vx - refpot)*timestep;
				long double W = exp(-k);
				double m = W + random_gen();

				if (m < 1)//death occurs
				{
					c--; //decrement walkers change by one
				}
				else if (m>2)//birth occurs
				{
					a = walkerpos[i]; //copy position to a
					b = walkersite[i]; //copy site to b
					walkerpos1[j] = a; //assign position a to new array
					walkerpos1[(j + 1)] = a; //duplicate position to next element
					walkersite1[j] = b; //assign site b to new array
					walkersite1[(j + 1)] = b; //assign site to next element
					j = j + 2; //increase iteration by 2
					c++; // increment walker change by one
				}
				else if (m > 1 && m < 2)
				{
					a = walkerpos[i]; //copy position to a
					b = walkersite[i]; //copy site to b
					walkerpos1[j] = a; //asign copy to new array
					walkersite1[j] = b; //asign copy site to new array
					j = j + 1; //increment by one
				}
			}
			nowalkers = nowalkers + c; //calculates the new number of walkers
			avev = 0; //reset average energy to zero

			//----------------- REMOVE DEAD WALKERS ----------------------------------------
			for (int i = 0; i < nowalkers; i++)//removes all dead walkers
			{
				double a; //term used to move walkers between arrays, reset to zero each time
				double b; //term used to move walkers site between arrays, reset to zero each time
				a = walkerpos1[i]; //assign position in copy array to transfer variable a
				walkerpos[i] = a; // assign transfer variable to original position array
				b = walkersite1[i]; // assign walker site to transfer variable
				walkersite[i] = b; //assign transfer variable to original site.

								   //Recalculates average potential energy after walker is removed
				if (walkersite[i] == 0) //Site A
				{
					avev = avev + (Potential(omega, walkerpos[i], F, site1)) / nowalkers; //forgot to divide by /nowalkers
				}
				else if (walkersite[i] == 1) //Site B
				{
					avev = avev + (Potential(omega, walkerpos[i], F, site2)) / nowalkers; //forgot to divide by /nowalkers
				}
				// ---------- WALKER HOPPING ------------------------------------
				double randomt = random_gen(); //generates random number between 0 and 1 randomt number which will compared against Phop to decide if the walker should hop
				if ((randomt <= Phop) && (walkersite[i] == 0))
				{
					walkersite[i] = 1;
				} //causes walker to hop site
				else if ((randomt <= Phop) && (walkersite[i] == 1))
				{
					walkersite[i] = 0;//causes walker to hop site
				}
				if (walkersite[i] == 0)
				{
					nowalkersite1++;
				}
				else if (walkersite[i] == 1)
				{
					nowalkersite2++;
				}


			}

			avev = avev - Hop_term; //corrects average energy

			//Removing copy arrays to prevent memory leask
			delete[] walkersite1;//information in second array deleted to prevent memory leak
			delete[] walkerpos1;//information in second array deleted to prevent memory leak

			Record_Energy(Energy, time, avev); //records average energy as time progresses
			Energy_test << time << "\t" << avev-0.5 << endl;
			//---------------- RECORD HISTOGRAM -----------------------------------------------

			if (time > (convergence_time/100)*t)//loop which creates histograms
			{
				aveE = aveE + (avev*(timestep / (t - convergence_time*t)));
				for (int i = 0; i < NUM_BINS; i++)//creates data for histogram
				{
					double lowerlimit = low_bin - (0.5*bin_width) + (i*bin_width);//lower limit for the bin
					double upperlimit = low_bin + (0.5*bin_width) + (i*bin_width);//upper limit for the bin
					for (int j = 0; j < nowalkers; j++)//Loop which places the walkers in the correct bin
					{
						if ((walkerpos[j] > lowerlimit) && (walkerpos[j] <= upperlimit))
						{
							if (walkersite[j] == 0)
							{
								histoa[i] = histoa[i] + (timestep / (t - convergence_time*t));//records each walker for site 1
							}
							if (walkersite[j] == 1)
							{
								histob[i] = histob[i] + (timestep / (t - convergence_time*t));//records each walker for site 2
							}
						}
					}
				}
			}
		}

		// *******************************************************************************************
		//						Simulation Ends Here
		// *******************************************************************************************

		Average_Energy << aveE << endl; //Append average energy to energy file
		Information_file(l+1,Info, MATLAB, Epi, Epf, Eps, omega, Hop_term, nowalkers_initial, t, dtn, dt_array, weighting); //record simulation parameters
		record_histogram(Histogram, histoa, histob, low_bin, high_bin, bin_width); //record histogram	
	}
	delete[] walkerpos;//deletes walker position histogram to prevent memory leak
	return 0;
}