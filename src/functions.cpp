#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

//User Libraries
#include "functions.h"


using namespace std;

void progress_bar(int sim_num, long double total_time, long double current_time)
{
	double progress = current_time / total_time;

	int barWidth = 70;
	cout << "Simulation #" << sim_num;
	cout << "[";
	int pos = barWidth * progress;

	for (int i = 0; i < barWidth; ++i)
	{
		if (i < pos) { cout << "="; }
		else if (i == pos) { cout << ">"; }
		else { cout << " "; }
	}
	cout << "] " << (progress * 100.0) << " %\r";
	cout.flush();

}

long double Potential(long double& x, long double& m1, long double& m2, int& potential_choice)
{
	// m1 = potential modifier 1
	// m2 = potential modifier 2
	long double V; //Potential Return variable (long double)

	if (potential_choice == 1) // Simple Harmonic Potential selection
	{
		// m1 = k (SHO spring constant)
		V = 0.5 * m1 * (x*x);
	}
	else if (potential_choice == 2) // Morse Potential selection
	{
		// m1 = De (Morse Well depth)
		// m2 = a (Morse Well width)
		V = m1*((1 - exp(-m2*x))*(1 - exp(-m2*x)));
	}
	else if (potential_choice == 3)
	{
		// m1 = V0 (Poschl-Teller potential height)
		// m2 = a (Poschl-Teller width)
		V = -m1 * (1 / (cosh(x / m2)*cosh(x / m2)));
	}
	else if (potential_choice == 4)
	{
		// m1 = A (double-well parameter A)
		// m2 = B (double-well parameter B)
		V = -m1*(x*x) + m2*(x*x*x*x);
	}
	else
	{
		cerr << "Potential Function Not Selected - Error terminated simulation " << endl;
	}

	return V;
}

void move_left(long double& position, double& step)
{
	position = position + step;
}

void move_right(long double& position, double& step)
{
	position = position - step;
}

double random_gen()
{
	double randm = rand() / (double)RAND_MAX;
	return randm;
}

void Initial_distribution(long double walkerpos[], const long double& nowalkers, int& dist_type)
{
	normal_distribution<double> distribution(0, 1); // mean = 0 std =1
	unsigned seed = time(NULL); // random time seed (32 to 64 bit warning)
	default_random_engine generator(seed);

	// dist_type 1 : sets all walkers to zero
	if (dist_type == 1)
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 0;
		}
	}
	// dist_type 2 : sets all walkers with a gaussian distribution
	else if (dist_type == 2)
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = distribution(generator);
		}
	}
	// dist_type 3 : creates a unifomly distributed grid between +/-4
	else if (dist_type == 3)
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 8 * (rand() / (double)RAND_MAX) - 4;
		}
	}
	// dist_type 4 : puts all walkers to the right (+4)
	else if (dist_type == 4)
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = 4;
		}
	}
	// dist_type 5 : puts all walkers to the left (-4)
	else if (dist_type == 5)
	{
		for (int i = 0; i < nowalkers; i++)
		{
			walkerpos[i] = -4;
		}
	}

	else
	{
		cerr << "Undefined Initial distribution - program will exit" << endl;
	}

}


void Bias_distribution(long double walkerpos[], long double& dt, long double& nowalkers_initial, int& sim_num)
{

	double xi; // x position variable containment
	double xmin = -10; //min x limit
	double xmax = 10; //max x limit
	double step = (xmax - xmin) / nowalkers_initial; //step size
	double y0, y1, y2, y3, y4; //wavefunctions/ eigenstates
	double G; //gaussian exponential
	const long double PI = acos(-1.0L); //Defining pi with CPU precision limit
	const long double C0 = 1 / (pow(PI, 0.25)); // coefficient  0
	const double C1 = 1 / (sqrt(2)); // coefficient 1
	const double C2 = 1 / (2 * sqrt(2)); //coefficient 2
	const double C3 = 1 / (4 * sqrt(3)); // coefficient 3
	const double C4 = 1 / (8 * sqrt(6)); // coefficient 4

	char filename_Bias[64];
	sprintf(filename_Bias, "Bias_distribution%d.txt", sim_num);
	ofstream Bias(filename_Bias);

	for (int i = 0; i < nowalkers_initial; i++)
	{
		xi = (xmin + i*step) - 5; // x - iteration variable
		G = exp(-0.5*xi*xi); //Gaussian function
		y0 = C0*G; // Ground state
		y1 = C1*C0*G * 2 * xi; // first excited state
		y2 = C2*C0*G*(4 * xi*xi - 2); //second excited state
		y3 = C3*C0*G*(8 * xi*xi*xi - 12 * xi); //third excited state
		y4 = C4*C0*G*(16 * xi*xi*xi*xi - 48 * xi*xi + 12); //fourth excited state
		walkerpos[i] = y0; //combination of eigenstates
		Bias << xi << "\t" << walkerpos[i] << endl; //record bias distribution
	}
}


void record_histogram(std::ofstream& Histogram, long double array[], double low_bin, double high_bin, double bin_width)
{
	//Calculate the number of bins in the histogram
	double NUM_BINS = ((high_bin - low_bin) / bin_width) + 1;


	if (Histogram.is_open()) //check if file is open
	{
		for (int i = 0; i < NUM_BINS; i++)
		{
			Histogram << low_bin + i*bin_width << "\t" << array[i] << endl;
		}
		Histogram.close();
	}
	else
	{
		cerr << "Error Recording Histogram File" << endl;
	}
}

void Information_file(ofstream& Info, ofstream& Other, int& V_choice, long double m1, long double m2, const long double& N0, const int& dtn, long double& simulation_time, long double dt_array[], long double weighting[], int& dist_type, long double& average_energy, double& CPU_Time)
{

	//Function creates the information file associated with each simulation
	// Other records the string identifier of the variable stored
	// Info records the actual value of the information stored

	if ((Info.is_open() && Other.is_open()))
	{
		Other << "Potential Choice" << endl;
		Info << V_choice << endl;

		Other << "Potential modifier m1" << endl;
		Info << m1 << endl;

		Other << "Potential modifier m2" << endl;
		Info << m2 << endl;

		Other << "Initial Number of Walkers" << endl;
		Info << N0 << endl;

		Other << "Total Simulation Time" << endl;
		Info << simulation_time << endl;

		Other << "Number of Timesteps used" << endl;
		Info << dtn << endl;

		Other << "Time step 1 (dt1)" << endl;
		Info << dt_array[0] << endl;

		Other << "Time step 2 (dt2)" << endl;
		Info << dt_array[1] << endl;

		Other << "Time step 3 (dt3)" << endl;
		Info << dt_array[2] << endl;

		Other << "Time step 4 (dt4)" << endl;
		Info << dt_array[3] << endl;

		Other << "Weighting dt1 (%)" << endl;
		Info << weighting[0] << endl;

		Other << "Weighting dt2 (%)" << endl;
		Info << weighting[1] << endl;

		Other << "Weighting dt3 (%)" << endl;
		Info << weighting[2] << endl;

		Other << "Weighting dt4 (%)" << endl;
		Info << weighting[3] << endl;


		Other << "Distribution Type" << endl;
		Info << dist_type << endl;

		Other << "Average Energy" << endl;
		Info << average_energy << endl;

		Other << "CPU Execution Time" << endl;
		Info << CPU_Time << endl;
	}

	else
	{
		cerr << "Unable to record Information and Other output files - program will exit" << endl;
	}

}

void record_snapshot(double array[], double& nowalkers, int& snap_number)
{
	char filename_Snap[64];
	sprintf(filename_Snap, "Snapshot%d.txt", snap_number);
	ofstream snapshot(filename_Snap);

	if (snapshot.is_open())
	{
		for (int i = 0; i < nowalkers; i++)
		{
			snapshot << array[i] << endl;
		}
	}
	else
	{
		cout << "Unable to record snapshot" << endl;
	}
}
