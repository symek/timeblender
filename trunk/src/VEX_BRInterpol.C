/*
VEX implementation of BRI:
Barycentric Rational Interpolation with
no Poles and High Rates of Approximation,
by Michael S. Floater and Kai Hormann

TB_Bri is based on "Numerical Recipes" 
(third edition).
skk.
*/

#include <VEX/VEX_VexOp.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Vector3.h>
#include <SYS/SYS_Math.h>
#include <UT/UT_RefArray.h>

#include "TB_GeoInterpolants.h"

using namespace TimeBlender;

/// Persistant interpolant class...
class PersistBRI
{
public:
    PersistBRI()
    {
        counter = 1;
        bri = new TB_Bri();
        cout << " bri inited." << endl;
    }
    ~PersistBRI() { delete bri; cout << "bri deleted" << endl;} 
    int counter;
    TB_Bri * bri;
};

/// ... and the pointer
static PersistBRI * persistBRI = NULL;

/// Init function:
static void * bri_init()
{
	if(!persistBRI)
	    persistBRI = new PersistBRI();
	else
	    persistBRI->counter++;
	return  persistBRI;
}

/// Clean up:
static void bri_cleanup(void * data)
{
    PersistBRI * bri = (PersistBRI *) data; 
    bri->counter--;
    if (!bri->counter)
    {
        delete bri;
        persistBRI = NULL;
        cout << "persist bri deleted" << endl;
    }

}

/*Own BRI implementation based on "Numerical Recipes" */
static void
brinterpol(int narg, void *argv[], void *data)
{
	// Repack arguments:
	float *out                         = (float *) argv[0];
	const UT_RefArray<fpreal32> *knots = (const UT_RefArray<fpreal32> *) argv[1];
	float *u                           = (float *) argv[2];
	int   *order                       = (int   *) argv[3];
	PersistBRI * persistbri            = (PersistBRI *) data;
    TB_Bri     * bri                   = persistbri->bri;
    int alloc                          = 0;
    
    if (persistbri->counter)
    {
	    //TB_Bri structure:
	    int size = knots->entries();
	    float  a[size];
	    const float* b = knots->getRawArray();
	    for (int i = 0; i < size; i++) 
	    {
	        a[i] = 1.0f/(size-1) * i;
	    }
	
	    alloc = bri->initialize(a, (float*)b, size, *order);
	}
	
	// eval:
	
	    float x;
	    x = bri->evaluate(*u);
	    out[0] = x;
}

void
newVEXOp(void *)
{
	new VEX_VexOp("brinterpol@&F[FFI", brinterpol, 
			VEX_ALL_CONTEXT, 
			bri_init,
			bri_cleanup, 
			VEX_OPTIMIZE_2, true);
}
