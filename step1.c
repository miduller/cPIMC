#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include "pimc.h"

/*compile using: gcc -o step1 step1.c `gsl-config --cflags --libs` -lm */

int main(int argc, char const *argv[])
{
	initialize_rng(31415926);
	srand(time(NULL));
	struct setup setup;
	setup.steps = 100000, setup.gap = 100, setup.equilibrationsteps = setup.steps / 10, setup.P = 50,
	setup.calibrationsteps = 10000, setup.energy = 0., setup.deltax = 1., setup.runs = 5;
	printf("%s\n", "Input the inverse temperature as a floating point value:");
	scanf("%lf", &setup.beta);
	setup.epsilon = setup.beta / setup.P;
	setup.exact = false;
	double energy = 0., acceptance;
	struct var var[setup.P];
	initialize(var, &setup);
	calibration(var, &setup);
	for (int i = 0; i < setup.runs; ++i)
	{
		initialize_rng(31415926+10*i);
		setup.energy = 0.;
		initialize(var, &setup);
		equilibration(var, setup);
		pimc_energy(var, &setup);
		//printf("%lf\n", setup.energy);
		energy += setup.energy / setup.runs;
	}
	printf("%lf\n", energy);

	return 0;
}