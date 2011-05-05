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

#include <interpolation.h>
#include <vector>

#define IG_TYPE_BRI 0
#define IG_TYPE_CAT 1

namespace TimeBlender
{


typedef enum {

	INTER_LINEAR,           /* */
	INTER_CATMULL_ROM,      /*Native HDK UT_Spline*/
	INTER_BARYCEN_RAT,      /*Barycentric Rational Interpolation*/

} IG_INTER_TYPES;






inline void fit( const double x, 
				 const fpreal a, 
				 const fpreal b, 
				 const fpreal c, 
				 const fpreal d, 
				 double &o) { o = c+((x-a)/(b-a))*(d-c);} 


class GeoInterpolant
{
public:
	// Allocate vectors for storing interpolants structures
	GeoInterpolant(int size, int type = INTER_BARYCEN_RAT) {
	if (type == INTER_BARYCEN_RAT)
	{
		interpolantsX.resize(size);
		interpolantsY.resize(size); 
		interpolantsZ.resize(size); 
		interpolants.resize(3);  
		interpolants.at(0) = &interpolantsX;
		interpolants.at(1) = &interpolantsY;
		interpolants.at(2) = &interpolantsZ;
	} else 
	{
		hinterpolants.resize(size);
	}

		itype              = type;
		mySize             = size;
		valid              = false;}

	~GeoInterpolant();

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


	/// Build interpolant per point (3 and 5 knots).
	inline void buildPoint(const GEO_Point *prev,
			 	  	       const GEO_Point *curr,
			  	  	       const GEO_Point *next,
					       const int        i,
					             alglib::barycentricinterpolant &);

	inline void buildPoint(const GEO_Point *prev,
			 	  	       const GEO_Point *curr,
			  	  	       const GEO_Point *next,
					       const int       *i,
					             UT_Spline  &);


	inline void buildPoint(const GEO_Point *prev2,
				  	       const GEO_Point *prev,
			 	  	       const GEO_Point *curr,
			  	  	       const GEO_Point *next,
				  	       const GEO_Point *next2,
						   const int       *i,
								 alglib::barycentricinterpolant &);

	inline void buildPoint(const GEO_Point *prev2,
				  	       const GEO_Point *prev,
			 	  	       const GEO_Point *curr,
			  	  	       const GEO_Point *next,
				  	       const GEO_Point *next2,
						   const int       *i,
								 UT_Spline   &);






	/// Calculate interpolat at idices: axis {0..2}, i = {0,....Npoints}, 
	/// at u parametric length.
	inline void calc( const int    *axis, 
					  const int    *i, 
					  		double *u, 
					  		float   &w);


	/// Check out interpolation type. This can only be set on class allocation.
	int getitype() {return itype; };

	bool isValid() { return valid;};

private:

	/// Get to the specific interpolant.
	const alglib::barycentricinterpolant * getInterpolant(int axis, int i)
	{ 
		return interpolants.at(axis)->at(i);
	}

	///Builds arrays of scalars per single interpolant.
	inline void build_ALGLIB_Arrays(const GEO_Point * prev, 
					   				const GEO_Point * curr, 
					   				const GEO_Point * next,
									const int         i,
							   		alglib::real_1d_array &idx,
							   		alglib::real_1d_array &v);

	/// Internal data structures.
	int  mySize;  
	bool valid; // Is set true, when interpolants were created succesfully.
	int  itype;

	/// Not perfectly, but we are going to store separete shots
	/// for BRI and HDK UT_Vector3
	vector <UT_Spline  *>                     hinterpolants;
	vector <alglib::barycentricinterpolant *> interpolantsX;
	vector <alglib::barycentricinterpolant *> interpolantsY;
	vector <alglib::barycentricinterpolant *> interpolantsZ;
	vector <vector <alglib::barycentricinterpolant *> *> interpolants;
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

	void       interpolate(GeoInterpolant *, 
						   double         *, 
						   GU_Detail      *);

	int        saveGeometry(const GU_Detail *,
							const UT_String *);

	UT_BoundingBox  myBox;

	int             mydointerpolate;
	int             threeknots;
	fpreal          myshutter;
	int             mynsamples;
	UT_String       myimethod;

	UT_String       myfile;
	UT_String       myprefile;
	UT_String       myprefile2;
	UT_String       mynextfile;
	UT_String       mynextfile2;
};

}//End of timeblender namescape

#endif
