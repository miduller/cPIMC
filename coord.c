#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "pimc.h"

/*compile using: gcc -o coord coord.c `gsl-config --cflags --libs` -lm */

int main(int argc, char const *argv[])
{
	initialize_rng(31415926);
	struct setup setup;
	setup.steps = 10e7, setup.gap = 100, setup.equilibrationsteps = setup.steps / 10, setup.calibrationsteps = 100000,
	setup.energy = 0., setup.deltax = 1., setup.P = 50;
	printf("%s\n", "Input the inverse temperature as a floating point value:");
	scanf("%lf", &setup.beta);
	setup.exact = false;
	setup.epsilon = setup.beta / setup.P;
	struct var var[setup.P];
	initialize(var, &setup);
	calibration(var, &setup);
	initialize(var, &setup);
	equilibration(var, setup);
	pimc_coord(var, &setup);
	printf("%lf\n", setup.energy);
	char filename[50];
	snprintf(filename, 50, "histogram_beta_%lf.dat", setup.beta);
	FILE *output;
	output = fopen(filename, "w");
	double sum = 0., varhist;
	for (int i = 1; i < 100; ++i)
	{
		sum += (setup.histogram[i-1] + setup.histogram[i]) / 2 * 0.08;
	}
	for (int i = 0; i < 100; ++i)
	{
		setup.histogram[i] /= sum;
		varhist = -4.0 + 0.04 + 0.08*i;
		fprintf(output, "%lf %lf\n",varhist , setup.histogram[i]);
	}
	fclose(output);
	return 0;
}