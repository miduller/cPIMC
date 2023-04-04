#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include "pimc.h"

/*compile using: gcc -o step2 step2.c `gsl-config --cflags --libs` -lm */

int main(int argc, char const *argv[])
{
	initialize_rng(31415926);
	struct setup setup;
	setup.steps = 10e7, setup.gap = 100, setup.equilibrationsteps = setup.steps / 10,
	setup.calibrationsteps = 10000, setup.runs = 5, setup.energy = 0., setup.deltax = 1.;
	printf("%s\n", "Input the inverse temperature as a floating point value:");
	scanf("%lf", &setup.beta);
	FILE *output;
	setup.exact = false;
	double energy, epsilon_sq, accept;
	char filename[30];
	snprintf(filename, 30, "energy_beta%lf.dat", setup.beta);
	output = fopen(filename, "w");
	for (int i = 0; i < 5; ++i)
	{
		setup.P = 40 + i * 5, setup.epsilon = setup.beta / setup.P;
		struct var var[setup.P];
		energy = 0.;
		printf("Trotter number of current calculation %d\n", 40 + i * 5);
		for (int j = 0; j < setup.runs; ++j)
		{
			initialize_rng(31415926+100*i+5*j);
			initialize(var, &setup);
			calibration(var, &setup);
			initialize(var, &setup);
			equilibration(var, setup);
			pimc_energy(var, &setup);
			energy += setup.energy;
			accept = acceptanceratio(var, setup.steps, setup.P);
			printf("Acceptance ratio is: %lf\n", accept);
		}
		energy /= setup.runs;
		epsilon_sq = setup.epsilon * setup.epsilon;
		fprintf(output, "%lf %lf\n", epsilon_sq, energy);
	}
	fclose(output);
	return 0;
}