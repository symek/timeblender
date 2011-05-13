/* 
	Timeblender::VRAY_InterpolatedGeometry v.01, 11.05.2011, 

	This is a VRAY Procedural DSO which computes a series of interpolated geometry
	from time samples at rendertime. Main reason for doing this is when we don't have
	an enough time steps in BGEO sequence to compute smooth non-linear motion blur.
	This can be crappy animation pipeline which doesn't handle subframes,
	or baked simulation (i.e. particles) with too high supersampling cost.

	It's meant to be a replacement for Mantra Deleyed Load Shader, albeit it can't currently 
	deal with materials archives stored in *.ifd files (HDK limitation).
	
	TimeBlender implements barycentric rational interpolation described here:
		"Barycentric Rational Interpolation with no Poles and High Rates of Approximation"
		by Michael S. Floater and Kai Hormann, and explain nicely in "Numerical Recipes".

	This kind of function interpolation performs specially well with small number of samples.
	It also quarantees to have second derivative ~= 0, what gives a desired smoothness. 

	skk.

	TODO: 
	- Fix Bounds!
	- Materials assigment (limited).
	- Own implementation of BRI (done).
	- UT_Splines (native HDK) interpolation (done).
	- Five knots interpolation (done).
	- Shutter retimer.
		-- Extrapolate motion.
		-- Nonlinear shutter retime a'la Pixar (?)
	- Heavy geometry tests (speed and memory).
	- Threading. (HDK native method).
	- Interpolation of non-constant topological geometry 
		-- Via attribute match (particles ID)
		-- Neigbhourhood search. 
	
	- Interpolate attributes: N, uv? 
    - threading (native HDK point loops)
		-- we need thread-safe vector for that.
	
*/

//#ifndef __VRAY_InterpolatedGeometry_C__
//#define __VRAY_InterpolatedGeometry_C__

#include "VRAY_InterpolatedGeometry.h"
#include "TB_PointMatch.h"

#if DEBUG==1
#define DEBUG
#endif

using namespace TimeBlender;

float TB_Bri::evaluate(float u) const
{
    if (!alloc) return 0.0f;
    float temp, p, q;
    p = 0.0; q = 0.0;
    for (int i=0; i<size; i++)
    {
        float t = u-idx[i];
        if ( t == 0.0)
        {
            /// we are exacty on index:
            return val[i];
        }
        else
        {
            temp = wei[i] / t;
            p   += val[i]*temp;
            q   += temp;
        }
    }
    return p/q;
}

int
TB_Bri::initialize(float *ii, float *x, int n, int d)
{
    /// initialize and compute weights.
    size = n; order=d;
    idx  = new float[n];
    val  = new float[n];
    wei  = new float[n];
    
    /// copy arrays:
    for (int i=0; i<n; i++)
    {
        idx[i] = ii[i];
        val[i] = x[i];
    }

    /// Compute weights:
    int imax, imin;
    float summ, temp;
    for (int k=0; k<n; k++)
    {
        imin = SYSmax(k-d, 0);
        imax = (k >= (n-d)) ? n-(d+1): k;
        temp = imin & 1 ? -1.0: 1.0;
        summ = 0.0;
        
        /// Summ:
        for(int i=imin; i<=imax; i++)
        {
            int   jmax = SYSmin(i+d, n-1);
            float term = 1.0;
            for (int j=i; j<=jmax; j++)
            {
                if(j == k) continue; 
                term *= (idx[k] - idx[j]);
            }
            term  = temp/term;
            temp  =-temp;
            summ += term; 
        }
        /// Weights computed:
        wei[k] = summ;
    }
    return 1;
}

void
BRInterpolant::interpolate(float u, GU_Detail * const gdp) const
{

	GEO_Point * ppt;
	int   i = 0;
	float x = 0; float y = 0; float z = 0;
	
	FOR_ALL_GPOINTS(gdp, ppt)
	{
		x = interpolants.at(0)->at(i)->evaluate(u);
		y = interpolants.at(1)->at(i)->evaluate(u);
		z = interpolants.at(2)->at(i)->evaluate(u);
	    ppt->setPos(x, y, z);
		i++;
	}
}

void
BRInterpolant::build(const GU_Detail * prev, const GU_Detail * curr, const GU_Detail * next)
{
    int i = 0;
    const GEO_Point  *currppt, *prevppt, *nextppt;
    float idx[] = {0.0f, 0.5f, 1.0f};
    float val[] = {0.0, 0.0, 0.0};
    
    /// Build correspondence:	

	/// BRI supports only floats type.
	/// TODO: Wy could try multithreading on this loop.
	FOR_ALL_GPOINTS(curr, currppt)
	{   
		prevppt = prev->points()[i];
		nextppt = next->points()[i];

		for (int j = 0; j < 3; j++)
		{
		    val[0] = prevppt->getPos()[j];
		    val[1] = currppt->getPos()[j];
		    val[2] = nextppt->getPos()[j];
			TB_Bri * bri =  new TB_Bri(idx, val, 3, 2);
			interpolants.at(j)->at(i) = bri;
		}
		i++;
	}
    this->valid = true;	
}

void
BRInterpolant::build(const GU_Detail * prev2,
                     const GU_Detail * prev, 
					 const GU_Detail * curr, 
					 const GU_Detail * next,
                     const GU_Detail * next2)
                     
{
    int i = 0;
    const GEO_Point  *currppt, *prevppt, *nextppt,  *prevppt2, *nextppt2;
    float idx[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
    float val[] = {0.0, 0.0, 0.0, 0.0, 0.0};	

	/// BRI supports only floats type.
	/// TODO: Wy could try multithreading on this loop.
	FOR_ALL_GPOINTS(curr, currppt)
	{   
		prevppt = prev->points()[i];
		nextppt = next->points()[i];
		prevppt2 = prev2->points()[i];
		nextppt2 = next2->points()[i];

		for (int j = 0; j < 3; j++)
		{
		    val[0] = prevppt2->getPos()[j];
		    val[1] = prevppt->getPos()[j];
		    val[2] = currppt->getPos()[j];
		    val[3] = nextppt->getPos()[j];
		    val[4] = nextppt2->getPos()[j];
			TB_Bri * bri =  new TB_Bri(idx, val, 5, 4);
			interpolants.at(j)->at(i) = bri;
		}
		i++;
	}
    this->valid = true;	
}

void
SplineInterpolant::interpolate(float u, GU_Detail * const gdp) const
{
	GEO_Point * ppt;
	int   i = 0;
	fpreal32 x[] = {0.0f, 0.0f, 0.0f};
	
	FOR_ALL_GPOINTS(gdp, ppt)
	{
		interpolants.at(i)->evaluate(u, x, 3, (UT_ColorType)2);
	    ppt->setPos(x[0],x[1],x[2]);
		i++;
	}
}

void
SplineInterpolant::build(const GU_Detail * prev, 
					     const GU_Detail * curr, 
					     const GU_Detail * next)

{
	/// Call 5 knots builder:
	build(prev, prev, curr, next, next);

}

void
SplineInterpolant::build(const GU_Detail * prev2,
                         const GU_Detail * prev, 
					     const GU_Detail * curr, 
					     const GU_Detail * next,
                         const GU_Detail * next2)
{
	UT_Spline  *spline;
	fpreal32 v1[3]; 
	int i = 0;
	const GEO_Point * currppt;
	const GU_Detail * gdps[] = {prev2, prev, curr, next, next2};

	/// TODO: Wy could try multithreading on this loop.
	FOR_ALL_GPOINTS(curr, currppt)
	{
		spline  = new UT_Spline(); 
		spline->setGlobalBasis((UT_SPLINE_BASIS)itype);
		spline->setSize(5, 3);

		for (int j=0; j<5;j++)
		{
			v1  = {gdps[j]->points()[i]->getPos().x(),
                   gdps[j]->points()[i]->getPos().y(),
                   gdps[j]->points()[i]->getPos().z()};

			spline->setValue(j, v1, 3);
		}

		interpolants.at(i) = spline;	
		
		i++;
	} 
	/// Interpolant is valid for evaluation.
	/// This is quite optimistic assumption, 
	/// as no checkes were performed.
	this->valid = true;
}

// Arguments:
static VRAY_ProceduralArg theArgs[] =
{
    VRAY_ProceduralArg("threeknots","int",   "1"),
    VRAY_ProceduralArg("prefile2",  "string", ""),
    VRAY_ProceduralArg("prefile",   "string", ""),
    VRAY_ProceduralArg("file",      "string", ""),
    VRAY_ProceduralArg("nextfile",  "string", ""),
    VRAY_ProceduralArg("nextfile2", "string", ""),
    VRAY_ProceduralArg("dointerpolate","int",   "1"),
    VRAY_ProceduralArg("nsamples",     "int",   "6"),
    VRAY_ProceduralArg("itype",        "int",   "4"),
    VRAY_ProceduralArg("shutter",      "real",  "1"),
    VRAY_ProceduralArg("velocityblur", "int",   "0"),
    VRAY_ProceduralArg("shutterretime", "int", "0"),
    /// These two are spare, as proc. get bounds in initialize(*box),
    /// Otherwise they need to be computed by us.
    VRAY_ProceduralArg("minbound", "real", "-1 -1 -1"),
    VRAY_ProceduralArg("maxbound", "real", "1 1 1"),
    VRAY_ProceduralArg()
};

// VRAY allocator:
VRAY_Procedural * 
allocProcedural(const char *)
{
	return new VRAY_IGeometry();
}

// Return argumentes
const VRAY_ProceduralArg * 
getProceduralArgs(const char *)
{
	return theArgs;
}

// Initialiser:
VRAY_IGeometry::VRAY_IGeometry()
{
	myBox.initBounds(0,0,0);
}

//Deallocator:
VRAY_IGeometry::~VRAY_IGeometry() {}

// Classname:
const char * 
VRAY_IGeometry::getClassName()
{
	return "VRAY_IGeometry";
}

// Initilize and set bounds:
int 
VRAY_IGeometry::initialize(const UT_BoundingBox *box)
{
	/// Main file:
	if (!import("file", myfile)) 
	{
		fprintf(stderr, "At least current frame must be specified.");
		return 0;	
	}

    /// Params:
	if (!import("dointerpolate", &mydointerpolate, 1)) 
		mydointerpolate = 1;	
    if (!import("threeknots", &mythreeknots, 1)) 
        mythreeknots = 1;
    if (!import("shutter", &myshutter, 1))
        myshutter = 1;
    if (!import("nsamples", &mynsamples, 1))
        mynsamples = 6;
    if (!import("itype", &myitype, 1))
        myitype = 4;

    /// Time samples:
	if (!import("prefile", myprefile))  cout << "No pre sample to interpolate." << endl;
	if (!import("nextfile", mynextfile)) cout << "No post sample to interpolate." << endl;

    /// Knots:
	if (mythreeknots == 0)
	{
		import("prefile2", myprefile2);
		import("nextfile2", mynextfile2); 
	}

	/// Bounding box (optionally from a file).
	/// TODO: Bounds should be enlarged with all gdps involved
	/// in interpolation:
	if (!box)
	{
        debug("Warning! No bounding box specified. Computing it from a sources.");
        GU_Detail gdp;
        UT_BoundingBox * gdpbox = new UT_BoundingBox();
        gdp.load(myfile, 0);
        gdp.getPointBBox(gdpbox);
        myBox = *gdpbox;
     } 
     else 
     {
        myBox = *box;
     }
	return 1;
}

// Return bounding box of a procedural:
void
VRAY_IGeometry::getBoundingBox(UT_BoundingBox &box)
{ 
	box = myBox;
}

// Actual render:
void 
VRAY_IGeometry::render()
{
	GU_Detail *gdp, *gdpp, *gdpn, *gdpp2, *gdpn2;

	gdpp = allocateGeometry();
	gdp  = allocateGeometry();
	gdpn = allocateGeometry();

	if (gdpp->load(myprefile,0) < 0) 
	{
		cout << "Can't open pre frame geometry: " << myprefile << endl;
		freeGeometry(gdpp);
		gdpp = 0;
	}

	if (gdp->load(myfile, 0) < 0)
	{
		cout << "Can't open current frame geometry: " << myfile << endl;
		freeGeometry(gdp);
		gdp = 0;
		return;
	}

	if (gdpn->load(mynextfile, 0) < 0)
	{
		cout << "Can't open next frame geometry: " << mynextfile << endl;
		freeGeometry(gdpn);
		gdpn = 0;
	}
	
	/// 5 knots mode:
	if (!mythreeknots)
	{
		gdpp2 = allocateGeometry();
		gdpn2 = allocateGeometry();

		if (gdpp2->load(myprefile2, 0) < 0)
		{
			cout << "Can't open pre pre frame geometry: " << myprefile2 << endl;
			freeGeometry(gdpp2);
			gdpp2 = 0;
		}

		if (gdpn2->load(mynextfile2, 0) < 0)
		{
			cout << "Can't open second next frame geometry: " << mynextfile2 << endl;
			freeGeometry(gdpn2);
			gdpn2 = 0;
		}
    }
		
    /// Main part goes here:
    /// TODO: assign shaders
    #ifdef DEBUG 
        TB_PointMatch * matcher = new TB_PointMatch();
        matcher->initialize((const GU_Detail *)gdp);
        cout << "TB_PointMatch entries: "<< matcher->entries() << endl; 
    #endif
    openGeometryObject();
    changeSetting("surface", "plastic diff (1.0 0.8 0.8)", "object");

    #ifdef DEBUG
        debug("openGeometryObject();");
    #endif
	/// Perform geometry interpolation.
    if (mydointerpolate)	
    { 
        GeoInterpolant * gi;
        static const fpreal min = 0.0, max = 1.0, nmin = 0.5; 
		fpreal nmax = 1.0;
		
        if (!mythreeknots) 
            nmax = 0.75; 

        fpreal32 fshutter = 0, shutter= 0;
        GU_Detail *bgdp;
		
		/// Allocate interpolant for npoints, and choose type:
        if (myitype	== INTER_BARYCENTRIC)
        {
            gi = new BRInterpolant(gdp->points().entries());
		}
        else
        {
            gi = new SplineInterpolant(gdp->points().entries(), myitype);
        }
        if (!gi->isAlloc())
        {
            debug("Couldn't allocate storate for the interpolant.");
            return; 
        }

		/// Build the interpolant from provided gdps
		/// with 3 or 5 time samples:
        if (mythreeknots) 
            gi->build(gdpp, gdp, gdpn);
        else 
            gi->build(gdpp2, gdpp, gdp, gdpn, gdpn2);

        if (!gi->isValid())
        {
            debug("Couldn't build the interpolant.");
            return;
        }

		/// Loop over samples generating interpolated geometry and add them to Mantra
		for (int i =0; i <= mynsamples; i++)
		{
			fshutter = shutter = (1.0f*i/mynsamples);

			/// We remap shutter to a proper range, since we have 
			/// effectively motion paths beyond Houdini's shutter:
			/// 0.5->1.0 for 3 knots, 0.5->0.75 for 5 knots.
			/// TODO: This could be controled by artist.
			fshutter = SYSfit(shutter, min, max, nmin, nmax);
            
			/// Allocate blur file and copy initial from current frame gdp.
			bgdp = allocateGeometry();
			bgdp->copy((const GU_Detail ) gdp, 0, false, true);
		
			/// Call interpolator, which replaces points' positions 
			/// (and only positions!(?)), finnaly add new blur file.
			if (bgdp && gi->isValid())
			{
				#ifdef DEBUG
					debug("Interpolating");
				#endif
				/// We evaluate interpolant on remapped time
				/// but add it to a scene on user time (?)
				gi->interpolate(fshutter*myshutter, bgdp);
				addGeometry(bgdp, shutter*myshutter);
				//referenceGeometry(bgdp);
			} 
			else 
			{
				closeObject();
			}
		}
	} 
	else 
	/// Proceed with standard blur file:
	{
		addGeometry(gdp, 0); 
		if (gdpn) 
		{
			addGeometry(gdpn, myshutter);	
		}
	}
	
	closeObject();
}
//#endif
