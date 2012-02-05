/* 
    Timeblender::VRAY_TimeBlender v.0.1.1, 05.02.2012, 
    Timeblender::VRAY_TimeBlender v.0.1.0, 11.05.2011, 

    This is a VRAY Procedural DSO which computes a series of interpolated geometry
    from time samples at rendertime. 
	
    TimeBlender implements barycentric rational interpolation described here:
    "Barycentric Rational Interpolation with no Poles and High Rates of Approximation"
    by Michael S. Floater and Kai Hormann, and explain nicely in "Numerical Recipes".

	skk.

	TODO: 
	- Fix Bounds! (enlerment).
	- Materials assigment (limited).
	- Shutter retimer.
		-- Extrapolate motion.
		-- Nonlinear shutter retime a'la Pixar (?)
	- Threading. (HDK native method).
	- Interpolate attributes: N, uv?,
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
  
    VRAY_ProceduralArg("files","int", "0"),
    VRAY_ProceduralArg("filename_string","string", ""),
    VRAY_ProceduralArg("nsamples",     "int",   "6"),
    VRAY_ProceduralArg("itype",        "int",   "0"),
    VRAY_ProceduralArg("shutter",      "real",  "1"),
    VRAY_ProceduralArg("matchbyid",    "int",   "0"),
    VRAY_ProceduralArg("shutterretime", "int", "0"),
    VRAY_ProceduralArg("shutter_start", "real", "0"),
    VRAY_ProceduralArg("shutter_end", "real", "1"),
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
    return "VRAY_TimeBlender";
}

// Initilize and set bounds:
int 
VRAY_TimeBlender::initialize(const UT_BoundingBox *box)
{
	
    /// Import filename_string, ie. whitespace separated list of files: 
	/// "file.bgeo file2.bgeo...". 
    if (!import("filename_string", myfilenamestring))
    {
        debug("Can't read filenames!");
        return 0;
    }
	
	/// Tokenize names and save in myfilenamelist.
    myfilenamestring.tokenize(myfilenamelist, " ");

    /// Params:
    if (!import("shutter_start", &myshutterstart, 1))
        myshutterstart = 0;    
        
    if (!import("shutter_end", &myshutterend, 1))
        myshutterend = 1;  
    
    if (!import("files", &myfiles, 1))
        myfiles = 0;
        
    if (!import("shutter", &myshutter, 1))
        myshutter = 1;
        
    if (!import("nsamples", &mynsamples, 1))
        mynsamples = 6;
        
    if (!import("itype", &myitype, 1))
        myitype = 0;
        
    if (!import("matchbyid", &mymatchbyid, 1))
        mymatchbyid = 0;
    
    
        /// TODO: Do we need this, or not?
        mycurrentframe = 0;
        
        
    /// Bounding box (optionally from a file).
    /// TODO: Bounds should be enlarged with all gdps involved
    /// in interpolation:
    if (!box)
    {
        debug("Warning! No bounding box specified. Computing from an input.");
        GU_Detail gdp;
        UT_BoundingBox * gdpbox = new UT_BoundingBox();
        gdp.load(myfilenamelist(mycurrentframe), 0);
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
    /// Create gdps and store it into array:    
    UT_PtrArray <GU_Detail *> gdps;
    for (int i = 0; i < myfilenamelist.getArgc(); i++)
    {
        GU_Detail *gdp;
        gdp = allocateGeometry();
        if (gdp->load(myfilenamelist(i), 0) < 0)
        {
            cout << "Can't open geometry: " << myfilenamelist(i) << endl;
            freeGeometry(gdp);
            gdp = NULL;
        }
        gdps.append(gdp);
    } 
    
    /// No geometries what-so-ever? return!
    if (gdps.entries() == 0) return;
   
    /// TODO: assign shaders
    openGeometryObject();
    //changeSetting("surface", 1, *shop_materialpath.buffer());
    changeSetting("surface", "plastic diff ( .8 .2 0 )", "object");
    
    /// Perform geometry interpolation.
    if (myitype != TB_INTER_NONE)	
    {   
        GeoInterpolant *gi;
        UT_PtrArray<TB_PointMatch*> match_array;
        
        if (myitype == TB_INTER_BARYCENTRIC) 
            gi = new BRInterpolant(gdps(mycurrentframe)->points().entries());
        else 
            gi = new SplineInterpolant(gdps(mycurrentframe)->points().entries(), myitype);
        
        /// Match points by id instead of numbering:         
        if (mymatchbyid && gdps(mycurrentframe)->getPointAttribute("id").isAttributeValid())
        {
            for (int i = 0; i < gdps.entries(); i++)
            {
                TB_PointMatch * match;
                match = new TB_PointMatch(gdps(i), CORR_POINT_ID);
                match_array.append(match);
            }
            gi->build(match_array, gdps(mycurrentframe));
        }
        else
        {
            /// FIXME: We should  catch cases of non-matching geometries.
            /// Poin-to-point by number:
            gi->build(gdps,  mycurrentframe);  
        }   
        
        /// Loop over samples generating interpolated geometry and add them to Mantra
		for (int i=0; i <= mynsamples-1; i++)
        {
            fpreal fshutter = 0; 
            fpreal shutter  = 0;      
            fshutter = shutter = (1.0f*i/mynsamples);
            fshutter = SYSfit(shutter, 0.0f, 1.0f, myshutterstart, myshutterend);
            cout << "fshutter: " << fshutter << endl;
            
            /// Allocate blur file and copy initialy from current frame gdp.
            GU_Detail  *bgdp;
            bgdp = allocateGeometry();
            bgdp->copy((const GU_Detail ) gdps(mycurrentframe), 0, false, true);
		
            /// Call interpolator, which replaces points' positions 
            if (bgdp && gi->isValid())
            {
                gi->interpolate(fshutter, bgdp);
                addGeometry(bgdp, shutter*myshutter);
            } 
            else 
            {
                closeObject();
            }
        }
        delete gi;
    } 
    else
    {
        debug("Proceeding with standard blur files.");
        for (int i = 0; i < gdps.entries(); i++)
        {
            cout << "shutter: " << 1.0f*i/gdps.entries() * myshutter << endl;
            addGeometry(gdps(i), 1.0f*i/gdps.entries() * myshutter); 
        }
    }

    closeObject();	
}
