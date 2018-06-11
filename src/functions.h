#ifndef functions_H
#define functions_H
#include <fstream>

void progress_bar(int, double, double);

double Potential(double&, double&, double&, int&);

void move_right(double&, double&);

void move_left(double&, double&);

double random_gen();

void Initial_distribution(double, const double&, int&);

void Bias_distribution(double, double&, double&, int&);

void record_histogram(std::ofstream&, double[], double, double, double);

void Information_file(std::ofstream& Info, std::ofstream&, int&, double, double, const double&, const int&, double&, double[], double[], int&, double&, double&);

void record_snapshot(double[], double&, int&);

#endif
