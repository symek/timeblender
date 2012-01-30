/* 
	Timeblender::VRAY_TimeBlender v.0.1.0, 11.05.2011, 

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
	- Fix Bounds! (done).
	- Materials assigment (limited).
	- Own implementation of BRI (done).
	- UT_Splines (native HDK) interpolation (done).
	- Five knots interpolation (done).
	
	- Interpolation of non-constant topological geometry 
        -- Via attribute match (particles ID) (in progress).
		-- Neigbhourhood search. 
		
	- Shutter retimer.
		-- Extrapolate motion.
		-- Nonlinear shutter retime a'la Pixar (?)
	- Heavy geometry tests (speed and memory).
	
	- Threading. (HDK native method).
	- Interpolate attributes: N, uv? 
    - threading (native HDK point loops)
		-- we need thread-safe vector for that.
	
*/

#include "VRAY_TimeBlender.h"
#include "TB_GeoInterpolants.h"
#include "TB_PointMatch.h"

#if DEBUG==1
#define DEBUG
#endif

using namespace TimeBlender;

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
    VRAY_ProceduralArg("matchbyid",    "int",   "0"),
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
	return new VRAY_TimeBlender();
}

// Return argumentes
const VRAY_ProceduralArg * 
getProceduralArgs(const char *)
{
	return theArgs;
}

// Initialiser:
VRAY_TimeBlender::VRAY_TimeBlender()
{
	myBox.initBounds(0,0,0);
}

//Deallocator:
VRAY_TimeBlender::~VRAY_TimeBlender() {}

// Classname:
const char * 
VRAY_TimeBlender::getClassName()
{
	return "VRAY_IGeometry";
}

// Initilize and set bounds:
int 
VRAY_TimeBlender::initialize(const UT_BoundingBox *box)
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
    if (!import("matchbyid", &mymatchbyid, 1))
        mymatchbyid = 0;

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
VRAY_TimeBlender::getBoundingBox(UT_BoundingBox &box)
{ 
	box = myBox;
}

// Actual render:
void 
VRAY_TimeBlender::render()
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
        matcher->initialize((const GU_Detail *)gdp, 0);
        cout << "TB_PointMatch entries: "<< matcher->entries() << endl;
        GEO_Point * ppt = matcher->find(10);
        if (ppt)
            cout << ppt->getPos()[0] << ppt->getPos()[1] << ppt->getPos()[2] << endl;
        else
            cout << "No point!" << endl;
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
        TB_PointMatch  * prevmatch, * nextmatch;
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
        
        /// If not allocated...
        if (!gi->isAlloc())
        {
            debug("Couldn't allocate storate for the interpolant.");
            return; 
        }

		/// Build the interpolant from provided gdps (3 or 5). 
        if (mythreeknots)
        {
            /// If "id" has been found, build the interpolant based on TB_PointMatch:
            if (mymatchbyid && gdp->getPointAttribute("id").isAttributeValid())
            {
                /// TODO: Properly handle exeption: "no id attrribute found".
                prevmatch = new TB_PointMatch(gdpp, CORR_POINT_ID);
                nextmatch = new TB_PointMatch(gdpn, CORR_POINT_ID);
                gi->build(prevmatch, (const GU_Detail *) gdp, nextmatch);
            }
            else
            {
                /// or use plain GU_Details':
                gi->build(gdpp, gdp, gdpn);
            }
        }
        else
        {
            /// Currently this doesn't work.
            if (mymatchbyid) debug("5 knobs interpolation doesn't work with match by id."); 
            gi->build(gdpp2, gdpp, gdp, gdpn, gdpn2);
        }
        
        /// Is interpolator valid?
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
