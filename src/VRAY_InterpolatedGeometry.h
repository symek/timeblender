#ifndef __VRAY_InterpolatedGeometry_h__
#define __VRAY_InterpolatedGeometry_h__


#include <VRAY/VRAY_Procedural.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_BoundingBox.h>
#include <GEO/GEO_Point.h>
#include <GU/GU_Detail.h>
#include <GB/GB_Macros.h>
#include <UT/UT_Vector4.h>
#include <UT/UT_Vector3.h>
#include <UT/UT_Spline.h>
#include <UT/UT_Color.h>

#include <interpolation.h>
#include <vector>



namespace TimeBlender
{

/// Interpolation type enumarator. 
typedef enum {

	INTER_CONSTANT, /*Native HDK UT_Splines*/
	INTER_LINEAR,					 
	INTER_CATMULLROM,
	INTER_CUBIC,
	INTER_BARYCENTRIC,  /*Barycentric Rational Interpolation*/          

} IG_INTER_TYPES;



////-------- Utility functions. ------------

inline void fit(const double x, 
                const fpreal a, 
                const fpreal b, 
                const fpreal c, 
                const fpreal d, 
                double &o) { o = c+((x-a)/(b-a))*(d-c);} 


inline void debug(const char * info) 
		{ cout << info << endl;}



/* Abstract class for storing interpolants. */
class GeoInterpolant
{
public:
	/// Call it in a default constractor.
	virtual int   initialize(int size, int type) = 0;

	/// Build structures necessery to do work.
	virtual void  build(const GU_Detail * prev, 
                        const GU_Detail * curr, 
                        const GU_Detail * next)  = 0;

	virtual void  build(const GU_Detail *,
                        const GU_Detail *,
                        const GU_Detail *, 
                        const GU_Detail *, 
                        const GU_Detail *)  = 0;


	/// This interpolates previously built stractures, and modifies
	/// GU_Detail accoring to it.
	virtual void interpolate(const float, GU_Detail  * const) const = 0;

	/// Object details.
	/// TODO: implement getMemoryUsage()
	/// TODO: implement allocation error checker.
	/// TODO: implement load(), save().

	virtual int getitype() const = 0;
	virtual bool isValid() const = 0;
	virtual bool isAlloc() const = 0;

private:
     /// TODO: This probably shouldn't be in an abstract class. 
     int  mySize;
     bool valid;
     bool alloc;
     int  itype;
};


///------ALGLIB specific GeoInterpolant child.---- ////
/// TODO: ALGLIB should be replaced by the own implementation.

class ALGLIB_Interpolant : public GeoInterpolant
{
public:
	// Allocate vectors for storing interpolants structures
	ALGLIB_Interpolant(int size, int type = INTER_BARYCENTRIC)
	{
		int success = initialize(size, type);
		if (!success) alloc = false;
	}

	/// This requires initialization.
	ALGLIB_Interpolant()
	{
		valid              = false;
		alloc              = false;
	}	

	/// Allocate memory, set flags.
	int initialize(int size, int type)
	{
		X.resize(size); Y.resize(size);
		Z.resize(size);
		interpolants.resize(3);  
		interpolants.at(0) = &X;
		interpolants.at(1) = &Y;
		interpolants.at(2) = &Z;
	
		itype              = INTER_BARYCENTRIC; // NOTE: type is ignored.
		mySize             = size;
		valid              = false;
		alloc              = true;
		return               1;
	}

	/// Build interpolants vector un mass for a whole geometry (3 and 5 knots).
	/// Can eat some memory.
	void  build(const GU_Detail * prev, 
                const GU_Detail * curr, 
                const GU_Detail * next);

	void  build(const GU_Detail * prev2,
                const GU_Detail * prev,  
                const GU_Detail * curr,  
                const GU_Detail * next,
                const GU_Detail * next2);


	/// Perform interpolation on GU_Detail *.
	void interpolate(const float, GU_Detail * const) const;


	/// Check out interpolation type. This can only be set on class allocation.
	int getitype() const { return itype; };
	bool isValid() const { return valid; };
	bool isAlloc() const { return alloc; }; 

private:

	/// Calculate interpolat at x,y: x {0, 1, 2}, y = {0,...,Npoints}, 
	/// at u parametric length.
	inline void evaluate(const int *x, int *y, 
                         float *u, float &w) const;

	/// Get to the specific interpolant.
	const alglib::barycentricinterpolant * getInterpolant(int axis, int i)
	{ 
		return interpolants.at(axis)->at(i);
	}

	/// Builds arrays of scalars used later by a single interpolant.
	inline void build_Arrays(const GEO_Point * prev, const GEO_Point * curr, 
                             const GEO_Point * next, const int i, 
                             alglib::real_1d_array &idx,
                             alglib::real_1d_array &v);

	inline void build_Arrays(const GEO_Point * prev2, const GEO_Point * prev, 
                             const GEO_Point * curr, const GEO_Point * next,
                             const GEO_Point * next2, const int i, 
                             alglib::real_1d_array &idx,
                             alglib::real_1d_array &v);

	/// Internal data structures.
	int  mySize;  
	bool valid; 
	int  itype;
	bool alloc;

	/// 3 vectors of 1-dim interpolants.
	vector <alglib::barycentricinterpolant *> X;
	vector <alglib::barycentricinterpolant *> Y;
	vector <alglib::barycentricinterpolant *> Z;
	vector <vector <alglib::barycentricinterpolant *> *> interpolants;
};




class SplineInterpolant : public GeoInterpolant
{
public:
	// Allocate vectors for storing interpolants structures
	SplineInterpolant(int size, int type = INTER_CUBIC)
	{
		int success = initialize(size, type);
		if (!success) alloc = false;
	}	
	/// This requires initialization.
	SplineInterpolant()
	{
		valid              = false;
		alloc              = false;
	}

	/// Allocate memory, set flags.
	int initialize(int size, int type)
	{
		interpolants.resize(size);  
		itype              = type;
		mySize             = size;
		valid              = false;
		alloc              = true;
		return               1;
	}

	void  build(const GU_Detail * prev, 
			    const GU_Detail * curr, 
			    const GU_Detail * next);

	void  build(const GU_Detail * prev2,
			    const GU_Detail * prev,  
			    const GU_Detail * curr, 
			    const GU_Detail * next,
				const GU_Detail * next2);


	void interpolate(const float,
					  GU_Detail * const) const;

	int getitype() const { return itype; };
	bool isValid() const { return valid; };
	bool isAlloc() const { return alloc; }; 

private:

	inline void evaluate(const int *x, int *y, 
                          float *u, float *w) const;

	/// Get to the specific interpolant.
	const UT_Spline * getInterpolant(int i)
	{ 
		return interpolants.at(i);
	}

	int  mySize;  
	bool valid;
	int  itype;
	bool alloc;

	/// Vector of 3-dim interpolants.
	vector <UT_Spline  *>  interpolants;
	
};





class VRAY_IGeometry: public VRAY_Procedural 
{
public:
		 VRAY_IGeometry();
virtual ~VRAY_IGeometry();

	virtual const char *getClassName();
	virtual int         initialize(const UT_BoundingBox *);
	virtual void        getBoundingBox(UT_BoundingBox &box);
	virtual void        render();

private:
	int        saveGeometry(const GU_Detail *,
                                const UT_String *);

	UT_BoundingBox  myBox;

	int             mydointerpolate;
	int             mythreeknots;
	fpreal          myshutter;
	int             mynsamples;
	int             myitype;

	UT_String       myfile;
	UT_String       myprefile;
	UT_String       myprefile2;
	UT_String       mynextfile;
	UT_String       mynextfile2;
};

}//End of timeblender namescape

#endif
