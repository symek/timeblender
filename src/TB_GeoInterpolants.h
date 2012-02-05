#ifndef __TB_GeoInterpolants_h__
#define __TB_GeoInterpolants_h__

#include <vector>

#include <GEO/GEO_Point.h>
#include <GU/GU_Detail.h>
#include <GB/GB_Macros.h>
#include <UT/UT_Spline.h>
#include <UT/UT_Color.h>

#include "TB_PointMatch.h"

namespace TimeBlender
{
/// Interpolation type enumarator. 
typedef enum {

	TB_INTER_NONE,   /* None */
	TB_INTER_LINEAR, /*Native HDK UT_Splines*/				 
	TB_INTER_CATMULLROM,
	TB_INTER_CUBIC,
	TB_INTER_BARYCENTRIC,  /*Barycentric Rational Interpolation*/          

} TB_INTER_TYPES;

/// Utility functions.
inline void fit(const double x, 
                const fpreal a, 
                const fpreal b, 
                const fpreal c, 
                const fpreal d, 
                double &o) { o = c+((x-a)/(b-a))*(d-c);} 


inline void debug(const char * info) 
		{ cout << info << endl;}


/*********************************************************
 Abstract class for storing interpolants. 
 ********************************************************/
 
class GeoInterpolant
{
public:
	/// Call it in a default constractor.
	virtual int initialize(int size, int typem) = 0;
	
	/// Build interpolant with plain GU_Details:
    virtual void build(UT_PtrArray<GU_Detail*> gdps, int current_frame) = 0;
    
    /// Build interpolant with Point Match: red/black tree 
    /// (std::map) for point-point by id matching
    virtual void build(UT_PtrArray<TB_PointMatch *>, const GU_Detail *) = 0;    
       
	/// This computes interpolattion and modifies GU_Detail's 'P' accoring to it.
	virtual void interpolate(const float, GU_Detail  * const) const = 0;

	/// Avarage mem usage:
	virtual int getMemoryUsage() const = 0;
	
	/// Utilities:
	virtual int init_arrays(float *a, float *b,
	                float *c,  float *d, int n) = 0;
	                
	virtual int getitype() const = 0;
	virtual bool isValid() const = 0;
	virtual bool isAlloc() const = 0;

private:
     ///  This probably shouldn't be in an abstract class?
     int  mySize;
     bool valid;
     bool alloc;
     int  itype;
};


/* ***************************************************************************************
/// An implementation of Barycentric Rational Interpolation based on: "Numerical Recipes",
/// and "Barycentric Rational Interpolation with no Poles and High Rates of Approximation",
/// by Michael S. Floater  and Kai Hormann.
******************************************************************************************/

class TB_Bri
{
public:
    /// Needs initialisation:
    TB_Bri() {};
    /// ii: index array; x: value array; n: arrays' size; 
    /// d: interpolation order <= n-1.
    TB_Bri(float *ii, float *x, int n, int d)
    { 
        alloc = initialize(ii, x, n, d); 
    }
    
    ~TB_Bri() 
    {
        delete [] idx;
        delete [] val;
        delete [] wei;
        cout << "idx, val, wei deleted" << endl; 
    }
    
    int initialize(float *ii, float *x, int n, int d);
    
    float evaluate(float u) const;
    
    int isAlloc() {return alloc;}
    int getMemoryUsage() { return size * sizeof(float) * 3;}
    
    /// TODO: make possible to reinitialize interpolant 
    /// with different order value (keeping arrays)
    /// void setOrder(int d) { order = d;} 
    
private:
    int alloc;
    int size;
    int order;
    float *idx;
    float *val;
    float *wei;
};

/****************************************************************
/ The interpolator based on TB_Bri (see above).
*****************************************************************/

class BRInterpolant : public GeoInterpolant
{
public:
	// Allocate vectors for storing interpolants structures
	BRInterpolant(int size, int type = TB_INTER_BARYCENTRIC)
	{
		if(!initialize(size, type)) alloc = false;
	};
	
	/// This requires initialization.
	BRInterpolant()
	{
		valid = false;
		alloc = false;
	};

	/// Allocate memory, set flags.
	int initialize(int size, int type)
	{
		interpolants.resize(3);
		X.resize(size, NULL); Y.resize(size, NULL);
		Z.resize(size, NULL);
		interpolants.at(0) = &X;
		interpolants.at(1) = &Y;
		interpolants.at(2) = &Z;
		
		itype      = TB_INTER_BARYCENTRIC;
		mySize     = size;
		valid      = false;
		alloc      = true;
		return 1;
	};
    
    /// Three main methods, builds with gdps, with red/black trees, 
    /// and interpolate positions in gdp:
    void build(UT_PtrArray<GU_Detail*> gdps, int current_frame);        
    void build(UT_PtrArray<TB_PointMatch *>, const GU_Detail *);
	void interpolate(const float, GU_Detail * const) const;
	
	/// The summ of ocupied memory:
	int getMemoryUsage() const 
	{ 
	    int mem = 0;
	    if (alloc) 
	        mem = X.at(0)->getMemoryUsage() * mySize * 3;
	    return mem;
	};
	
	/// Utilities:
	int init_arrays(float *a, float *b, float *c,  float *d, int n);
	int getitype() const { return itype; };
	bool isValid() const { return valid; };
	bool isAlloc() const { return alloc; }; 

private:
	int  mySize;  
	bool valid;
	int  itype;
	bool alloc;

	/// 3 vectors of 1-dim interpolants.
	vector <TB_Bri *> X;
	vector <TB_Bri *> Y;
	vector <TB_Bri *> Z;
	vector <vector<TB_Bri  *> *>  interpolants;
};


/*************************************************************
/ Interpolator based on HDK UT_Spline.
*************************************************************/

class SplineInterpolant : public GeoInterpolant
{
public:
	// Allocate vectors for storing interpolants structures
	SplineInterpolant(int size, int type = TB_INTER_CUBIC)
	{
		if (!initialize(size, type)) alloc = false;
	}
		
	/// This requires initialization.
	SplineInterpolant()
	{
		valid              = false;
		alloc              = false;
	};

	/// Allocate memory, set flags.
	int initialize(int size, int type)
	{
		interpolants.resize(size);
		itype              = type;
		mySize             = size;
		valid              = false;
		alloc              = true;
		return               1;
	};
	            
    void build(UT_PtrArray<GU_Detail*> gdps, int current_frame);
    void build(UT_PtrArray<TB_PointMatch *>, const GU_Detail *);
	void interpolate(const float, GU_Detail * const) const;

	int getitype() const { return itype; };
	bool isValid() const { return valid; };
	bool isAlloc() const { return alloc; };
	
	
	int getMemoryUsage() const 
	{ 
	    // No idea...
	    int mem = 0;
	    return mem;
	};
	
	int init_arrays(float *a, float *b,
	                float *c,  float *d, int n);

private:
	inline void evaluate(const int *x, int *y, float *u, float *w) const;

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
}//End of timeblender namescape
#endif
