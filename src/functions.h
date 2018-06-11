#ifndef functions_H
#define functions_H
#include <fstream>

void progress_bar(int, long double, long double);

long double Potential(long double&, long double&, long double&, int&);

void move_right(long double&, double&);

void move_left(long double&, double&);

double random_gen();

void Initial_distribution(long double, const long double&, int&);

void Bias_distribution(long double, long double&, long double&, int&);

void record_histogram(std::ofstream&, long double[], double, double, double);

void Information_file(std::ofstream& Info, std::ofstream&, int&, long double, long double, const long double&, const int&, long double&, long double[], long double[], int&, long double&, double&);

void record_snapshot(double[], double&, int&);

#endif
