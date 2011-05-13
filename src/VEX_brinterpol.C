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
#include <UT/UT_Spline.h>
#include <UT/UT_Color.h>
#include <UT/UT_Vector3.h>
#include <SYS/SYS_Math.h>


class BRInterpolate
{
public:
    BRInterpolate(float *ii, float *x, int n, int d)
    {
        /// initialize and compute weights.
        size = n; order= d;
        idx  = new float[n];
        val  = new float[n];
        wei  = new float[n];
        
        /// copy arrays:
        for (int i=0; i<n; i++)
        {
            idx[i] = ii[i];
            val[i] = x[i];
        }
    
        /// weights:
        int imax, imin;
        float summ;
        for (int k=0; k<n; k++)
        {
            imin = SYSmax(k-d,0);
            imax = (k >= (n-d)) ? n-(d+1): k;
    
            /// Summ:
            for(int i=imin; i<=imax; i++)
            {
                summ = 1.0;
                for (int j=i; j<=i+d; j++)
                {
                    if(j != k) summ *= SYSfabs(idx[k] - idx[j]);
                }
                summ    = 1.0f/summ;
                wei[k] += summ; 
            }
            if (SYSabs(k-d)%2 == 1) wei[k] *= -1;
        }
        
    
    }
   
    /// Evaluation of bri on initialized array of values with unknown weights
    float evaluate(float u)
    {
        float temp, p, q;
        p = 0.0; q = 0.0;
        for (int i=0; i<size; i++)
        {
            if (u == idx[i])
                return val[i];  
           temp = wei[i] / (u-idx[i]);
           p   += val[i]*temp;
           q   += temp;
        }
        return p/q;
    }

    int getSize(void) {return size;}
    
private:
    int   size;
    int   order;
    float *idx;
    float *val;
    float *wei;

};

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

/* UT_Spline test re-implementation 3 knots*/
static void
hspline_knots3(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out    = (float *) argv[0];
	float *u      = (float *) argv[1];
	int   *type   = (int   *) argv[2];
	float *tens   = (float *) argv[3];
	float *prev_p = (float *) argv[4];
	float *pp     = (float *) argv[5];	
	float *post_p = (float *) argv[6];
	

	//UT_Spline structure:
	UT_Spline *spline = new UT_Spline();
	spline->setGlobalBasis((UT_SPLINE_BASIS)*type);
	spline->setTension((fpreal64)*tens);
	spline->setSize(5, 1);
	

	
	int   a[] = {0, 1, 2, 3, 4};
	const fpreal32 b[] = {(const fpreal32)*prev_p, (const fpreal32)*prev_p, 
                          (const fpreal32)*pp, (const fpreal32)*post_p, (const fpreal32)*post_p};
	int i = 0;	

	/// Set values:
	for (i=0; i<5; i++)
	{
		spline->setValue(a[i],  &b[i], 1);
	}
	
	// eval:
	spline->evaluate(*u, out, 1, (UT_ColorType)3);
	
}


/// UT_Spline test re-implementation 3 knots vectors
static void
hspline_knots3_v(int narg, void *argv[], void *data)
{
	// Repack arguments:
	UT_Vector3 *out    = (UT_Vector3 *) argv[0];
	float      *u      = (float      *) argv[1];
	int        *type   = (int        *) argv[2];
	float      *tens   = (float      *) argv[3];
	UT_Vector3 *prev_p = (UT_Vector3 *) argv[4];
	UT_Vector3 *pp     = (UT_Vector3 *) argv[5];	
	UT_Vector3 *post_p = (UT_Vector3 *) argv[6];
	fpreal32 raw_out[] = {0.0f,0.0f,0.0f}; 
	

	//UT_Spline structure:
	UT_Spline *spline = new UT_Spline();
	spline->setGlobalBasis((UT_SPLINE_BASIS)*type);
	spline->setTension((fpreal64)*tens);
	spline->setSize(5, 3);
	

	
	//int   a[] = {0, 1, 2, 3, 4};
	const fpreal32 b1[] = {prev_p->x(), prev_p->y(),prev_p->z()};
	const fpreal32 b2[] = {pp->x(),     pp->y(),   pp->z()};
	const fpreal32 b3[] = {post_p->x(), post_p->y(), post_p->z()};
	

	/// Set values:
	spline->setValue(0,  b1, 3);
	spline->setValue(1,  b1, 3);
	spline->setValue(2,  b2, 3);
	spline->setValue(3,  b3, 3);
	spline->setValue(4,  b3, 3);

	// eval:
	spline->evaluate(*u, raw_out, 3, (UT_ColorType)2);
	out->assign(raw_out[0], raw_out[1], raw_out[2]);
}



/*Own BRI implementation based on "Numerical Recipes" 3 knots*/
static void
nrbri_knots3(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out    = (float *) argv[0];
	float *u      = (float *) argv[1];
	float *prev_p = (float *) argv[2];
	float *pp     = (float *) argv[3];	
	float *post_p = (float *) argv[4];
	

	//UT_Spline structure:
	float a[] = {0.0f, 0.5f, 1.0f};
	float b[] = {*prev_p, *pp, *post_p};
	BRInterpolate *spline = new BRInterpolate(a, b, 3, 2);
	
	// eval:
	float x;
	x = spline->evaluate(*u);
	out[0] = x;
	
}

/*Own BRI implementation based on "Numerical Recipes" 3 knots*/
static void
nrbri_knots5(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out    = (float *) argv[0];
	float *u      = (float *) argv[1];
	float *prev_p = (float *) argv[2];
	float *prev_p2 = (float *) argv[3];
	float *pp     = (float *) argv[4];	
	float *post_p = (float *) argv[5];
	float *post_p2 = (float *) argv[6];
	

	//UT_Spline structure:
	float a[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
	float b[] = {*prev_p2, *prev_p, *pp, *post_p, *post_p2};
	BRInterpolate *spline = new BRInterpolate(a, b, 5, 4);
	
	// eval:
	float x;
	x = spline->evaluate(*u);
	out[0] = x;
	
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
	new VEX_VexOp("hspline@&VFIFVVV", hspline_knots3_v, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);
    new VEX_VexOp("nrbrispline@&FFFFF", nrbri_knots3, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);
	  new VEX_VexOp("nrbrispline@&FFFFFFF", nrbri_knots5, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);
	new VEX_VexOp("hspline@&FFIFFFF", hspline_knots3, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);

}
