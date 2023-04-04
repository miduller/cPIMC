#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


struct var
{
	double x;
	double x_sq;
	int moves;
};

struct setup
{
	int P;
	long int steps;
	int equilibrationsteps;
	int calibrationsteps;
	int gap;
	int runs;
	bool exact;
	double beta;
	double epsilon;
	double deltax;
	double energy;
	double energy_error;
	double histogram[100];
};

gsl_rng * r;

void initialize_rng(long int seed)
{
	const gsl_rng_type * T;
  	gsl_rng_env_setup();

  	T = gsl_rng_ranlxs2;
  	r = gsl_rng_alloc (T);
  	gsl_rng_set(r, seed);
}

double uniform(double a, double b)
{
  	double u = gsl_rng_uniform(r);
  	return u * (b - a) + a;
}

long int randint(long int a)
{
	long int u = gsl_rng_uniform_int(r, a);
	return u;
}

void update(struct var var[], double deltax, int P, double epsilon, bool exact)
{
	int pick = randint(P);
	double xbeforepick, xafterpick, xtrial, probability;
	if (pick == 0)
	{
		xbeforepick = var[P - 1].x;
	} 
	else
	{
		xbeforepick = var[pick - 1].x;
	}
	if (pick == P - 1)
	{
		xafterpick = var[0].x;
	} 
	else
	{
		xafterpick = var[pick + 1].x;
	}
	xtrial = var[pick].x + uniform(- deltax, deltax);
	if (exact == false)
	{
		double var1 = - epsilon * 0.5 * (xtrial * xtrial - var[pick].x * var[pick].x);
		double var2 = - 0.5 * ((xbeforepick - xtrial) * (xbeforepick - xtrial) + (xtrial - xafterpick) * (xtrial - xafterpick)
			- (xbeforepick - var[pick].x) * (xbeforepick - var[pick].x) - (var[pick].x - xafterpick) * (var[pick].x - xafterpick)) / epsilon;
		probability = exp(var1 + var2);
	}
	if (exact == true)
	{
		double var1 = cosh(epsilon) * (xtrial * xtrial - var[pick].x * var[pick].x);
		double var2 = (xbeforepick + xafterpick) * (var[pick].x - xtrial);
		probability = exp(- 1 / sinh(epsilon) * (var1 + var2));
	}
	if (probability > uniform(0, 1))
	{
		var[pick].x = xtrial;
		var[pick].x_sq = xtrial * xtrial;
		var[pick].moves += 1;
	}

}

double energy(struct var var[], int P)
{
	double temp = 0.;
	for (int i = 0; i < P; ++i)
	{
		temp += var[i].x_sq;
	}
	return 1. / P * temp;
}

double acceptanceratio(struct var var[], int steps, int P)
{
	double temp = 0.;
	for (int i = 0; i < P; ++i)
	{
		temp += var[i].moves;
	}
	return temp / steps;
}

void initialize(struct var var[], struct setup *setup)
{
	for (int i = 0; i < setup->P; ++i)
	{
		var[i].x = gsl_ran_gaussian(r, 1.);
		var[i].x_sq = var[i].x * var[i].x;
		var[i].moves = 0;
	}
	for (int i = 0; i < 100; ++i)
	{
		setup->histogram[i] = 0.;
	}
}

void initialize_moves(struct var var[], struct setup *setup)
{
	for (int i = 0; i < setup->P; ++i)
	{
		var[i].moves = 0;
	}
}

void calibration(struct var var[], struct setup *setup)
{
	bool calibrated = false;
	double acceptance;
	while (calibrated == false)
	{
		for (int i = 0; i < setup->calibrationsteps; ++i)
		{
			update(var, setup->deltax, setup->P, setup->epsilon, setup->exact);
		}
		acceptance = acceptanceratio(var, setup->calibrationsteps, setup->P);
		initialize_moves(var, setup);
		if (acceptance < 0.39)
		{
			setup->deltax *= 0.9;
		}
		if (acceptance > 0.45)
		{
			setup->deltax *= 1.1;
		}
		if (acceptance >= 0.39 && acceptance <= 0.45)
		{
			calibrated = true;
		}

	}
}

void equilibration(struct var var[], struct setup setup)
{
	for (int i = 0; i < setup.equilibrationsteps; ++i)
	{
		update(var, setup.deltax, setup.P, setup.epsilon, setup.exact);
	}
	for (int i = 0; i < setup.P; ++i)
	{
		var[i].moves = 0;
	}
}

void pimc_coord(struct var var[], struct setup *setup)
{
	for (long int i = 0; i < setup->steps; ++i)
	{
		update(var, setup->deltax, setup->P, setup->epsilon, setup->exact);
		//printf("%lf\n", var[10].x);
		if (i%setup->gap == 0)
		{
			for (int j = 0; j < 100; ++j)
			{
				for (int k = 0; k < setup->P; ++k)
				{
					if (-4.0 + 0.08*j < var[k].x && -4.0 + 0.08 * (j+1) > var[k].x)
					{
						setup->histogram[j] += 1.;
					}
				}
			}
		}
	}
}

void pimc_energy(struct var var[], struct setup *setup)
{
	setup->energy = 0.;
	setup->energy_error = 0.;
	for (long int i = 0; i < setup->steps; ++i)
	{
		update(var, setup->deltax, setup->P, setup->epsilon, setup->exact);
		//printf("%lf\n", var[10].x);
		if (i%setup->gap == 0)
		{
			setup->energy += 1. / (setup->steps / setup->gap) * energy(var, setup->P);
			setup->energy_error += 1. / (setup->steps / setup->gap) * energy(var, setup->P) * energy(var, setup->P);
		}
	}
	setup->energy_error = sqrt(1. / (setup->steps / setup->gap) * (setup->energy_error - setup->energy * setup->energy));
}