#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
//#include <time.h>
#include <ctime>
#include <cstdio>
#include "functions.h"

using namespace std;


int main()
{
	clock_t begin, end; //used for measuring execution time

	srand(time(NULL)); //seed the random number generator

					   // ---- Initializing walkerposition and histogram
	const int hist_size = 10000, pos_size = 20000000;
	double* walkerpos = new double[pos_size];
	double c;
	double nowalkers;
	// double nowalkers_initial; //initial number of walkers (should be const)

	double t; //simulation time \tau
	double Er = 0; // Reference energy

						// --- Initializing Potential Modifier Parameters
	double m1 = 0; //modifier 1
	double m2 = 0; //modifier 2
	int potential_choice; //potential choice variable

	int dist_type; //initial distribution choice type
	const double ij = -10; //bin low point
	const double ji = 10; //bin high point
	const double l = 20; //bin difference
	const double dl = 0.1; //bin_width
	const int NUM_BINS = ((ji - ij) / dl) + 1;
	double convergence_time = 0;

	int sim_num;
	//int dt_n;

	//double dummy_variable;

	// --- Time steps and weighting ---
	int dt_n;
	double dt1;
	double dt2;
	double dt3;
	double dt4;
	double dt1_lim;
	double dt2_lim;
	double dt3_lim;
	double dt_array[3] = { 0 };
	double weighting[3] = { 0 };

// ---------------------------------------------------------------------------------------------------
//                          USER INTERFACE VIA COMMAND WINDOW
// ---------------------------------------------------------------------------------------------------

	cout << "Diffusion Monte Carlo Algorithm - Simple Test Potentials" << endl;
	cout << "\t\t Version 3" << endl << endl;

	// ---------- Prompts the user to select a Potential type -----------
	cout << "Please choose the potential V(x) you want to test: " << endl;
	cout << "\t1 - Simple Harmonic Oscillator Potential " << endl;
	cout << "\t2 - Morse Potential " << endl;
	cout << "\t3 - Poschl-Teller Potential " << endl;
	cout << "\t4 - Double-Well Potential " << endl;
	cin >> potential_choice;
	cout << endl;

	// ---------- Prompts the user to select the specific Potential parameter (Modifiers) -----------
	cout << "Please provide the needed modifiers for the chosen potential: " << endl;
	if (potential_choice == 1)
	{
		cout << "k: ";
		cin >> m1;
		cout << endl;
	}
	else if (potential_choice == 2)
	{
		cout << "De: ";
		cin >> m1;
		cout << endl;
		cout << "a: ";
		cin >> m2;
		cout << endl;
	}
	else if (potential_choice == 3)
	{
		cout << "V0: ";
		cin >> m1;
		cout << endl;
		cout << "a: ";
		cin >> m2;
		cout << endl;
	}
	else if (potential_choice == 4)
	{
		cout << "A: ";
		cin >> m1;
		cout << endl;
		cout << "B: ";
		cin >> m2;
		cout << endl;
	}

	// ---------- Prompts the user for the number of sequential simulations to run ------------
	cout << "How Many simulations: ";
	cin >> sim_num;
	cout << endl;

	// ---------- Prompts the user for the initial population size ---------------------
	cout << "Please enter the initial walker population size: "; //Initial number of walkers
	cin >> nowalkers;
	const double nowalkers_initial = nowalkers;

	// ---------- Propts the user for the total simulation time ------------
	cout << endl << "Pleas give the total simulation time (in imaginary units) t: "; //User input simulation time
	cin >> t;
	const double total_time = t;

	// ---------- Prompts the user for the number of time steps to be used -------
	cout << endl << "Please give the number of time steps to be used (max 4): ";
	cin >> dt_n;
	const int time_steps = dt_n;

	if (time_steps == 1)
	{
		cout << "dt1: ";
		cin >> dt1;
		dt_array[0] = dt1;
		cout << endl;
	}
	else if (time_steps == 2)
	{
		cout << "dt1: ";
		cin >> dt1;
		dt_array[0] = dt1;
		cout << "dt1 limit (% of t): ";
		cin >> dt1_lim;
		weighting[0] = dt1_lim;
		dt1_lim = dt1_lim / 100;
		cout << endl;

		cout << "dt2: ";
		cin >> dt2;
		dt_array[1] = dt2;
		cout << endl;
	}
	else if (time_steps == 3)
	{
		cout << "dt1: ";
		cin >> dt1;
		dt_array[0] = dt1;
		cout << "dt1 limit (% of t): ";
		cin >> dt1_lim;
		weighting[0] = dt1_lim;
		dt1_lim = dt1_lim / 100;
		cout << endl;

		cout << "dt2: ";
		cin >> dt2;
		dt_array[1] = dt2;
		cout << "dt2 limit (% of t): ";
		cin >> dt2_lim;
		weighting[1] = dt2_lim;
		dt2_lim = dt2_lim / 100;
		cout << endl;

		cout << "dt3: ";
		cin >> dt3;
		dt_array[2] = dt3;
		cout << endl;
	}
	else if (time_steps == 4)
	{
		cout << "dt1: ";
		cin >> dt1;
		dt_array[0] = dt1;
		cout << "dt1 limit (% of t): ";
		cin >> dt1_lim;
		weighting[0] = dt1_lim;
		dt1_lim = dt1_lim / 100;
		cout << endl;

		cout << "dt2: ";
		cin >> dt2;
		dt_array[1] = dt2;
		cout << "dt2 limit (% of t): ";
		cin >> dt2_lim;
		weighting[1] = dt2_lim;
		dt2_lim = dt2_lim / 100;
		cout << endl;

		cout << "dt3: ";
		cin >> dt3;
		dt_array[2] = dt3;
		cout << "dt3 limit (% of t): ";
		cin >> dt3_lim;
		weighting[2] = dt3_lim;
		dt3_lim = dt3_lim / 100;
		cout << endl;

		cout << "dt4: ";
		cin >> dt4;
		//		dt_array[3] = dt4;
		cout << endl;
	}

	// ---------- Prompts the user to select an initial distribution type -----------
	cout << "Please select the intial distribution type" << endl;
	cout << "\t1 - All walkers at zero" << endl;
	cout << "\t2 - Normal distribution (mean = 0, std = 1)" << endl;
	cout << "\t3 - +/- 4 grid with uniform  distribution" << endl;
	cout << "\t4 - All walkers to the right (+4)" << endl;
	cout << "\t5 - All walkers to the left (-4)" << endl;
	cin >> dist_type;

	// ---------- Prompts the user when the histogram recording for groundstate wavefunction is to be done
	cout << "When do you want to start building Histogram statistics (% of t): ";
	cin >> convergence_time;
	convergence_time = convergence_time / 100;
	cout << endl;

	// building of datespecifier for creation of folders and output files

	char datespecifier[64];
	build_datestring(datespecifier);

	for (int sim_iteration = 1; sim_iteration<(sim_num + 1); sim_iteration++)
	{

		//Creation of specific filenames within the given folder
		// filename = bin/'date-time'+filename

		char filename_E[64];
		sprintf(filename_E, "bin/%s-Energies%d.txt", datespecifier, sim_iteration);

		char filename_W[64];
		sprintf(filename_W, "bin/%s-walkers%d.txt", datespecifier, sim_iteration);

		char filename_T[64];
		sprintf(filename_T, "bin/%s-performance%d.txt", datespecifier, sim_iteration);

		char Filename_I[64];
		sprintf(Filename_I, "bin/%s-Info%d.txt", datespecifier, sim_iteration);

		char Filename_O[64];
		sprintf(Filename_O, "bin/%s-Other%d.txt", datespecifier, sim_iteration);

		char Filename_H[64];
		sprintf(Filename_H, "bin/%s-Wavefunction%d.txt", datespecifier, sim_iteration);

		//creates the output files and their respective pointers
		ofstream Energy(filename_E);
		ofstream walkers(filename_W);
		ofstream cpu_time(filename_T);
		ofstream Info(Filename_I);
		ofstream Other(Filename_O);
		ofstream Histogram(Filename_H);

		//------- initialized histogram to zero
		double wavefunction[hist_size] = { 0 }; //histogram initialization (may be implementation dependent)

		//------- creates the initial distribution of walkers
		//Initial_distribution(walkerpos, nowalkers_initial, dist_type);

		double aveE = 0; //initializes average energy
		begin = clock(); // obtains the clock time to begin
		double timestep = dt1; // initial timestep
		double steplength = sqrt(3 * dt1); //initial steplength
// -----------------------------------------------------------------------------------------------------------------------
// 				Simulation Begins Here
// -----------------------------------------------------------------------------------------------------------------------

		for (double time = 0; time < t; time = time + timestep)
		{
			//-------- Outputs the progress bar to the console
			progress_bar(sim_iteration, t, time);

			if (time_steps == 1)
			{
				timestep = dt1;
			}
			else if (time_steps == 2)
			{
				if (time > dt1_lim*t)
				{
					timestep = dt2;
				}
			}
			else if (time_steps == 3)
			{
				if (time < dt1_lim*t)
				{
					timestep = dt1;
				}
				if ((time > dt1_lim*t) && (time < dt2_lim*t))
				{
					timestep = dt2;
				}
				if (time > dt2_lim*t)
				{
					timestep = dt3;
				}

			}
			else if (time_steps == 4)
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

			steplength = sqrt(3 * timestep); //refreshes the steplength with new timestep

			double* walkerpos1 = new double[pos_size]; //initialize pointer to copy array
			Er = 0; //initialize reference energy

			double avev = 0; //initialize average potential

			// ----- Moves walkers -------------
			for (int i = 0; i < nowalkers; i++)
			{
				double step = steplength*random_gen();
				double randn = rand() % 2;

				if (randn == 1)
				{
					move_left(walkerpos[i], step);
				}
				else
				{
					move_right(walkerpos[i], step);
				}
				// ----- Recalculate average potential
				avev = avev + (Potential(walkerpos[i], m1, m2, potential_choice) / nowalkers);
			}

			// ----- calculate average energy
			Er = avev - ((nowalkers / nowalkers_initial) - 1) / timestep;

			int j = 0; //copy array iteration variable initialization
			c = 0; // initialize death counter

			// ----- Walker population control loop
			for (int i = 0; i < nowalkers; i++)
			{
				double a; //copy variable
				double k = (Potential(walkerpos[i], m1, m2, potential_choice) - Er)*timestep;
				double m = exp(-k) + random_gen();

				if (m < 1) //death condition
				{
					c--;
				}
				else if (m > 2) // duplication (birth) condition
				{
					a = walkerpos[i];
					walkerpos1[j] = a;
					walkerpos1[(j + 1)] = a;
					j = j + 2;
					c++;
				}
				else if (m >1 && m < 2) //staying alive condition
				{
					a = walkerpos[i];
					walkerpos1[j] = a;
					j = j + 1;
				}

			}
			//records walkers
			walkers << time << "\t" << nowalkers << endl;
			nowalkers = nowalkers + c;
			avev = 0;

			//removes dead walkers and recalculates average potential
			for (int i = 0; i < nowalkers; i++)
			{
				double a;
				a = walkerpos1[i];
				walkerpos[i] = a;
				avev = avev + (Potential(walkerpos[i], m1, m2, potential_choice) / nowalkers);
			}

			delete[] walkerpos1; //delete copy array to prevent memory leak
			//records energy
			Energy << time << " " << avev << " " << Er << endl;

			// ------- Records histogram -----------
			if (time > convergence_time*t)
			{

				//average energy calculation
				aveE = aveE + (avev*(timestep / (t - convergence_time*t)));

				build_histogram(wavefunction, walkerpos, nowalkers, NUM_BINS, ij, ji, l, dl, convergence_time, dt1, t);

			}

		}//time loop ends
		end = clock();
		double CPU_Time = ((double)end - (double)begin) / (double)CLOCKS_PER_SEC;
		cpu_time << "Execution Time = " << CPU_Time << "\tseconds" << endl;

		record_histogram(Histogram, wavefunction, ij, ji, dl); //function writes histogram to file

		cout << endl << "Average Energy = " << aveE << endl;

		Information_file(Info, Other, potential_choice, m1,m2,nowalkers_initial, time_steps, t, dt_array,weighting, dist_type, aveE, CPU_Time); //writes data to information file


		Energy.close();	walkers.close();
	}
	delete[] walkerpos;
	//system("pause");
	return EXIT_SUCCESS;
}
