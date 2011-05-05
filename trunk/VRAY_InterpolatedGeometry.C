/* 
	Timeblender::VRAY_InterpolatedGeometry v.01, 4.05.2011, 

	This is a VRAY Procedural DSO which computes a series of interpolated geometry
	from time samples at rendertime. Main reason for doing this is when we don't have
	an enough time steps in BGEO sequence to compute smooth non-linear motion blur.
	This can be crappy animation pipeline which doesn't handle subframes,
	or baked simulation (i.e. particles) with too high supersampling cost.

	It's meant to be a replacement for Mantra Deleyed Load Shader, albeit it can't currently 
	deal with materials archives stored in *.ifd files (HDK limitation).
	
	Currently TimeBlender uses ALGLIB library to compute interpolation according to:
		"Barycentric Rational Interpolation with no Poles and High Rates of Approximation"
		by Michael S. Floater and Kai Hormann

	This kind of function interpolation performs specially well with small number of samples.
	It also quarantees to have second derivative ~= 0, what gives a desired smoothness. 
	
	ALGLIB is GPL'ed for non-commercial OS usage, but requires commercial license otherwise.

	skk.

	TODO: 
	- Fix Bounds!
	- Separete files .C/.h for GInterpolant and VRAY_IGeometry (really? nightmare with name-spaces)
	- UT_Splines (native HDK) interpolation.
	- Five knots interpolation.
	- Shutter retimer.
		-- Extrapolate motion.
		-- Nonlinear shutter retime a'la Pixar (?)
	- GInterpolant:
		-- vector proccessing
		-- Native HDK types?
		-- Own implementation of BRI(?)
	- Heavy geometry tests (speed and memory).
	- Threading. (HDK native method).
	- Interpolation of non-constant topological geometry 
		-- Via attribute match (particles ID)
		-- Neigbhourhood search. 
	- Materials assigment.
	- Interpolate attributes: N, uv? 
	
*/

#ifndef __VRAY_InterpolatedGeometry_C__
#define __VRAY_InterpolatedGeometry_C__

#include "VRAY_InterpolatedGeometry.h"


using namespace alglib;
using namespace TimeBlender;



inline void
GeoInterpolant::calc(const int * axis, 
					 const int * i,
					 double    * u, /*TODO Fix mixing types (floats, fpreals, doubles)*/
				     float      &w)
{
	w = alglib::barycentriccalc(*interpolants.at(*axis)->at(*i), *u);
}

/*
inline void
GeoInterpolant::buildPoint(const GEO_Point *prev,
			 	  	       const GEO_Point *curr,
			  	  	       const GEO_Point *next,
					       const int        i, 
					             alglib::barycentricinterpolant &inpt)

{

	real_1d_array idx;
	real_1d_array samples;

	for (int j = 0; j < 3; j++)
	{
		buildArrays(prev, curr, next, j, idx, samples);	
		polynomialbuild(idx, samples, *inpt);
	}

}

*/
void
GeoInterpolant::build(const GU_Detail * prev, 
					  const GU_Detail * curr, 
					  const GU_Detail * next)
{
	int i = 0;
	int j;
	const GEO_Point  *currppt, *prevppt, *nextppt;
	
	/// interpolation type?
	if (this->itype == INTER_BARYCEN_RAT)
	{
		real_1d_array idx;
		real_1d_array samples;

		/// Unfortunatelly for BRI we create a single interpolant
		/// for every component (float).
		/// TODO: Wy could try multithreading on this loop.
		FOR_ALL_GPOINTS(curr, currppt)
		{   
			prevppt = prev->points()[i];
			nextppt = next->points()[i];

			for (j = 0; j < 3; j++)
			{
				build_ALGLIB_Arrays(prevppt, currppt, nextppt, j, idx, samples);
				barycentricinterpolant * bc =  new barycentricinterpolant();
				polynomialbuild(idx, samples, *bc);
				interpolants.at(j)->at(i) = bc;
			}
			i++;
		}

	  	this->valid = true;		

	} else 
	{
		UT_Spline  *spline;
		fpreal32 v[3];

		
		/// Retrieve points' positions and cunstruct
		/// splines from it. Store results in this->hinterpolate <vector>,
		/// so it can be evaluated later. 
		/// TODO: Wy could try multithreading on this loop.
 
		FOR_ALL_GPOINTS(curr, currppt)
		{
			spline  = new UT_Spline(); 
			spline->setGlobalBasis(UT_SPLINE_CATMULL_ROM);
			spline->setSize(3, 3);

			v[0] = prev->points()[i]->getPos()[0];
			v[1] = prev->points()[i]->getPos()[1];
			v[2] = prev->points()[i]->getPos()[2];
			spline->setValue(0, v, 3);

			v[0] = currppt->getPos()[0];
			v[1] = currppt->getPos()[1];
			v[2] = currppt->getPos()[2];
			spline->setValue(1, v, 3);

			v[0] = next->points()[i]->getPos()[0];
			v[1] = next->points()[i]->getPos()[1];
			v[2] = next->points()[i]->getPos()[2];
			spline->setValue(2, v, 3);

			hinterpolants.at(i) = spline;	
			
			i++;
		}
		/// Interpolant is valid for evaluation.
		/// This is quite optimistic assumption, 
		/// no checkes were performed.
		this->valid = true;	
	}
}



inline void
GeoInterpolant::build_ALGLIB_Arrays(const GEO_Point * prev, 
					   				const GEO_Point * curr, 
					   				const GEO_Point * next,
									const int         i,
					   				real_1d_array &idx,
					   				real_1d_array &v)
{

	double rawsamples[3];
	double ii[]  = {0.0f, 0.5f, 1.0f};
	idx.setcontent(3, ii);

	UT_Vector4  currpos,  prevpos,  nextpos;
	prevpos = prev->getPos();
	currpos = curr->getPos();
	nextpos = next->getPos();

	rawsamples[0] = (double)prevpos[i]; 
	rawsamples[1] = (double)currpos[i]; 
	rawsamples[2] = (double)nextpos[i];
	v.setcontent(3, rawsamples);

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
	VRAY_ProceduralArg("nsamples",    "int",    "6"),
	VRAY_ProceduralArg("imethod",     "string", "bri"),
	VRAY_ProceduralArg("shutter",     "real",   "1"),
    VRAY_ProceduralArg("velocityblur","int",    "0"),

	VRAY_ProceduralArg("shutterretime", "int", "0"),

	VRAY_ProceduralArg("minbound", "real", "-10 -10 -10"),
	VRAY_ProceduralArg("maxbound", "real", "10 10 10"),
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
VRAY_IGeometry::initialize(const UT_BoundingBox *)
{
	/* Get all parameters here: */ 
	if (!import("file", myfile)) 
	{
		fprintf(stderr, "At least current frame must be specified.");
		return 0;	
	} else {
		if (!myfile.isstring()) return 0;
	}

	if (!import("dointerpolate", &mydointerpolate, 1)) 
		mydointerpolate = 1;
	
		if (!import("threeknots", &threeknots, 1)) 
			threeknots = 1;
		if (!import("shutter", &myshutter, 1))
			myshutter = 1;
		if (!import("nsamples", &mynsamples, 1))
			mynsamples = 6;
		if (!import("imethod", myimethod))
			myimethod = "bri";

	if (!import("prefile", myprefile))  cout << "No pre sample to interpolate" << endl;
	if (!import("nextfile", mynextfile)) cout << "No post sample to interpolate" << endl;

	if (threeknots == 0)
	{
		import("prefile2", myprefile2);
		import("nextfile2", mynextfile2); 
	}


	// TODO: Bounds don't work!
	fpreal val[3];
	val[0] = val[1] = val[2] = -1;
	import("minbound", val, 3);
	myBox.initBounds(val[0], val[1], val[2]);
	val[0] = val[1] = val[2] = 1;
	import("maxbound", val, 3);
	myBox.enlargeBounds(val[0], val[1], val[2]);
	
	return 1;

}

// Return bounding box of a procedural:
void
 VRAY_IGeometry::getBoundingBox(UT_BoundingBox &box)
{ 
	box = myBox;
}


void
VRAY_IGeometry::interpolate(GeoInterpolant * gi, 
							double         * u, 
							GU_Detail      * gdp)
{

	GEO_Point * ppt;
	int   i = 0;
	float x = 0; float y = 0; float z = 0;
	static const int jx = 0, jy = 1, jz = 2;
	
	FOR_ALL_GPOINTS(gdp, ppt)
	{
		gi->calc(&jx, &i, u, x);
		gi->calc(&jy, &i, u, y);
		gi->calc(&jz, &i, u, z);
	    ppt->setPos(x,y,z);
		i++;
	}

}


int
VRAY_IGeometry::saveGeometry(const GU_Detail * gdp,
							 const UT_String * path)
{
	
		//char tmpname[50];
		//const char * name = "/tmp/IG_test";
		//sprintf(tmpname, "%s.%i.bgeo", name, i+25);   
		//bgdp->save((const char*)tmpname, 0,0);
		return 1;

}




// Actual render:
void 
VRAY_IGeometry::render()
{

	GU_Detail *gdp0, *gdp1, *gdp2;

	gdp0 = allocateGeometry();
	gdp1 = allocateGeometry();
	gdp2 = allocateGeometry();

	if (gdp0->load(myprefile,0) < 0) 
	{
		cout << "Can't open pre frame geometry: " << myprefile << endl;
		freeGeometry(gdp0);
		gdp0 = 0;
	}

	if (gdp1->load(myfile, 0) < 0)
	{
		cout << "Can't open current frame geometry: " << myfile << endl;
		freeGeometry(gdp1);
		gdp1 = 0;
		return;
	}

	if (gdp2->load(mynextfile, 0) < 0)
	{
		cout << "Can't open next frame geometry: " << mynextfile << endl;
		freeGeometry(gdp2);
		gdp2 = 0;
	}


	/// This is main part:
	openGeometryObject();
	changeSetting("surface", "plastic diff (1.0 0.8 0.8)", "object");

	/// Perform geometry interpolation.
	if (mydointerpolate)	
	{
		static const fpreal min = 0, max = 1, nmin = 0.5, nmax = 1.0; 
		double fshutter = 0, shutter= 0;
		GU_Detail  *bgdp = NULL;


		/// Allocate interpolant for npoints,  and  build it with 3 or 5 time sampels. 
		GeoInterpolant * gi = new GeoInterpolant(gdp1->points().entries());
    	gi->build(gdp0, gdp1, gdp2);
		
		//if (0 == 0) return;
		
		/// Loop over samples generating interpolated geometry and add them to Mantra
		for (int i =1; i < mynsamples+1; i++)
		{
			fshutter = shutter = (1.0f*i/mynsamples);

			/// We remap shutter to a proper range, since we have 
			/// effectively motion paths beyond Houdini's shutter:
			/// 0.5->1.0 for 3 knots, 0.5->0.75 for 5 knots.
			/// TODO: This could be controled by artist.
			fit(shutter, min, max, nmin, nmax, fshutter);

			/// Allocate blur file and copy initial from current frame gdp.
  			bgdp = allocateGeometry();
			bgdp->copy((const GU_Detail ) gdp1, 0, false, true);

			/// Call interpolator, which replaces points' positions 
			/// (and only positions!(?)), finnaly add new blur file.
			if (bgdp && gi->isValid())
			{
				/// We evaluate geo on remapped time
				/// but add it to a scene on user time (?)
				interpolate(gi, &fshutter, bgdp);
				addGeometry(bgdp, shutter*myshutter);
			} else 
			{
				closeObject();
			}
		}
	} 
	else 
	/// Proceed with standard blur file:
	{
		addGeometry(gdp1, 0); 
		if (gdp2) 
		{
			addGeometry(gdp2, myshutter);	
		}
	}
	
	closeObject();

}

#endif

