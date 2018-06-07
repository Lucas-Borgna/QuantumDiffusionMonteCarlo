#include <iostream>
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