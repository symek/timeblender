/* This code  uses ALGLIB implementation of:

Barycentric Rational Interpolation with
no Poles and High Rates of
Approximation
Michael S. Floater
and Kai Hormann

Since VEX_Op hasn't arrays interface at the moment, 
I implement 3 and 5 knots float version separately.
skk.

*/



#include <VEX/VEX_VexOp.h>
#include <UT/UT_DSOVersion.h>
#include <interpolation.h>

static void * bri_init()
{
	
	static float *x = NULL;
	return  x;
}


static void bri_cleanup(void * data)

{

	//static float *x = NULL;
	//return  x;

}

using namespace alglib;

/* Barycentric Rational Interpolator with 3 floats knots*/
static void
bri_knots3(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out    = (float *) argv[0];
	float *u      = (float *) argv[1];
	float *prev_p = (float *) argv[2];
	float *pp     = (float *) argv[3];	
	float *post_p = (float *) argv[4];
	double uu     = (double ) *u;

	//Alglib structures:
	barycentricinterpolant p;
	real_1d_array aa;
	real_1d_array bb;

	// Alglib custom arrays: :(
	double a[] = {0.0f,0.5f,1.0f};
	double b[] = {(double)*prev_p, (double)*pp, (double)*post_p};
	aa.setcontent(3, a);
	bb.setcontent(3, b);
	
	// Build a polynomial and evalute:
	polynomialbuild(aa, bb, p);
	out[0] = barycentriccalc(p, uu);
	
}

/* Barycentric Rational Interpolator with 5 float knots*/
static void
bri_knots5(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out    = (float *) argv[0];
	float *u      = (float *) argv[1];
	float *prev_p2= (float *) argv[2];
	float *prev_p = (float *) argv[3];
	float *pp     = (float *) argv[4];	
	float *post_p = (float *) argv[5];
	float *post_p2= (float *) argv[6];
	double uu     = (double ) *u;

	//Alglib structures:
	barycentricinterpolant p;
	real_1d_array aa;
	real_1d_array bb;

	// Alglib custom arrays: :(
	double a[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
	double b[] = {(double)*prev_p2, (double)*prev_p, (double)*pp, (double)*post_p, (double)*post_p2};
	aa.setcontent(5, a);
	bb.setcontent(5, b);
	
	// Build a polynomial and evalute:
	polynomialbuild(aa, bb, p);
	out[0] = barycentriccalc(p, uu);
	
}


void
newVEXOp(void *)
{
	new VEX_VexOp("brispline@&FFFFF", bri_knots3, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);
	new VEX_VexOp("brispline@&FFFFFFF", bri_knots5, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);

}
