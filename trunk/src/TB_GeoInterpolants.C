#include "TB_GeoInterpolants.h"

using namespace TimeBlender;

/// This code comes from a "Numerical Recipies" (third edition),
/// and implements: "Barycentric Rational Interpolation with
/// no Poles and High Rates of Approximation",
/// by Michael S. Floater and Kai Hormann.
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



void
BRInterpolant::build(const GU_Detail * prev, const GU_Detail * curr, const GU_Detail * next)
{
    int i = 0;
    const GEO_Point  *currppt, *prevppt, *nextppt;
    float idx[] = {0.0f, 0.5f, 1.0f};
    float val[] = {0.0, 0.0, 0.0};

	/// BRI currently supports only floats.
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

	/// BRI currently supports only floats.
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

/// This version of build deals with TB_PointMatch.
void
BRInterpolant::build(TB_PointMatch * prev, const GU_Detail * curr, TB_PointMatch * next)
{
    int i = 0; int id;
    const GEO_Point  *currppt;
    GEO_Point *prevppt, *nextppt;
    float idx[] = {0.0f, 0.5f, 1.0f};
    float val[] = {0.0, 0.0, 0.0};
    
    /// Get id of a current point.
    GEO_AttributeHandle handle = curr->getAttribute(GEO_POINT_DICT, "id");
   
	FOR_ALL_GPOINTS(curr, currppt)
	{
	    handle.setElement(currppt);
        id = handle.getI();
        /// TODO: Fix missing ids!
        /// Points with no forward or backward
        /// correspondences should be deleted from the gdp.    
		prevppt = prev->find(id);
		if (!prevppt) prevppt = (GEO_Point*)currppt;
		nextppt = next->find(id);
		if (!nextppt) nextppt = (GEO_Point*)currppt;
		
        /// BRI currently supports only floats.     
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

/// Interpolant based on UT_Spline class (native HDK splines).
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

void
SplineInterpolant::build(TB_PointMatch   * prev, 
					     const GU_Detail * curr, 
					     TB_PointMatch   * next)
{
	/// Call 5 knots builder:
	/// TODO: placeholder
	cout << "THIS IS ONLY PLACEHOLDER" << endl;
	build(curr, curr, curr, curr, curr);
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

