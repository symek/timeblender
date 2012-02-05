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
    alloc =1;
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

int
BRInterpolant::init_arrays(float *a, float *b, float *c, float *d, int n)
{
    debug("Initing!");
    a = new float[n]; 
    b = new float[n];
    c = new float[n];
    d = new float[n];
    
    if (!a || !b || !c || !d) 
        return 0;
        
    debug("after a or b or c");
    
    for (int i = 0; i < n; i++) 
        a[i] = 1.0f*i/n;
        
    memset(b, 0.0, sizeof(float) * n);
    memset(c, 0.0, sizeof(float) * n);
    memset(d, 0.0, sizeof(float) * n);

    debug("After memset");
    return 1;
}


void
BRInterpolant::build(UT_PtrArray<GU_Detail*> gdps, int current_frame)
{
    /// We work on current_gdp, taken from gdps.
    GU_Detail *current_gdp, *gdp;
    const GEO_Point  *ppt;
    current_gdp = gdps(current_frame);
    int entries = gdps.entries();
    
    /// Initialize index and values arrays:
    float *idx, *valx, *valy, *valz;
    idx  = new float[entries]; 
    valx = new float[entries];
    valy = new float[entries];
    valz = new float[entries];
    
    for (int i = 0; i < entries; i++) 
        idx[i] = 1.0f*i/entries;
        
    memset(valx, 0.0, sizeof(float) * entries);
    memset(valy, 0.0, sizeof(float) * entries);
    memset(valz, 0.0, sizeof(float) * entries);
    
    /// Loop over points -> then gdps. 
    /// Is it goint to be much slower, than opposite?
    for (int i=0; i < current_gdp->points().entries(); i++)
    {   
        for (int g = 0; g < entries; g++)
        {
            gdp = gdps(g);
            // FIXME: This stays valid for constant point number, add exception.
            ppt = gdp->points()(i);
	        valx[g] = ppt->getPos().x();
		    valy[g] = ppt->getPos().y();
		    valz[g] = ppt->getPos().z(); 
	     }
	    /// Three splines are created, values copied, and assigned to 'this':
        TB_Bri * brix =  new TB_Bri(idx, valx, entries, entries-1);
        TB_Bri * briy =  new TB_Bri(idx, valy, entries, entries-1);
        TB_Bri * briz =  new TB_Bri(idx, valz, entries, entries-1);
        interpolants.at(0)->at(i) = brix;
        interpolants.at(1)->at(i) = briy;
        interpolants.at(2)->at(i) = briz;
    }
    delete idx; 
    delete valx; 
    delete valy;
    delete valz;
    valid = true;	
}


void 
BRInterpolant::build(UT_PtrArray<TB_PointMatch *> matches, const GU_Detail * gdp)
{
    //GU_Detail *current_gdp,;
    TB_PointMatch *current;
    const GEO_Point  *ppt;
    int entries = matches.entries();
    int id;
    
     /// Initialize index and values arrays:
    float *idx, *valx, *valy, *valz;
    idx  = new float[entries]; 
    valx = new float[entries];
    valy = new float[entries];
    valz = new float[entries];
    
    for (int i = 0; i < entries; i++) 
        idx[i] = 1.0f*i/entries;
        
    memset(valx, 0.0, sizeof(float) * entries);
    memset(valy, 0.0, sizeof(float) * entries);
    memset(valz, 0.0, sizeof(float) * entries);
    
    GEO_AttributeHandle handle = gdp->getAttribute(GEO_POINT_DICT, "id");
    
    for (int i=0; i < gdp->points().entries(); i++)
    {   
        handle.setElement(gdp->points()(i));
        id = handle.getI();
        
        for (int g = 0; g < entries; g++)
        {
            current = matches(g);
            ppt     = current->find(id);
            if (!ppt) ppt = gdp->points()(i);
	        valx[g] = ppt->getPos().x();
		    valy[g] = ppt->getPos().y();
		    valz[g] = ppt->getPos().z(); 
	     }
	       
	    /// Three splines are created, values copied, and assigned to 'this':
        TB_Bri * brix =  new TB_Bri(idx, valx, entries, entries-1);
        TB_Bri * briy =  new TB_Bri(idx, valy, entries, entries-1);
        TB_Bri * briz =  new TB_Bri(idx, valz, entries, entries-1);
        interpolants.at(0)->at(i) = brix;
        interpolants.at(1)->at(i) = briy;
        interpolants.at(2)->at(i) = briz;
    }
    
    delete idx; 
    delete valx; 
    delete valy;
    delete valz;
    valid = true;	
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


int
SplineInterpolant::init_arrays(float *a, float *b, float *c,  float *d, int n) { return 0; }


void 
SplineInterpolant::build(UT_PtrArray<GU_Detail*> gdps, int current_frame)
{
    GU_Detail *current_gdp, *gdp;
    const GEO_Point  *ppt;
    current_gdp = gdps(current_frame);
    int entries = gdps.entries();
    UT_Spline  *spline;
    fpreal32 v1[3]; 
    
    for(int i = 0; i < current_gdp->points().entries(); i++)
    {
        spline  = new UT_Spline(); 
        spline->setGlobalBasis((UT_SPLINE_BASIS)itype);
        spline->setSize(entries, 3);
        
        for (int g = 0; g < entries; g++)
        {
            gdp    = gdps(g);
            ppt    = gdp->points()(i);
            v1[0]  = ppt->getPos().x();
            v1[1]  = ppt->getPos().y();
            v1[2]  = ppt->getPos().z();
            spline->setValue(g, v1, 3);
        }
        interpolants.at(i) = spline;
    }
    valid = true;
}


void 
SplineInterpolant::build(UT_PtrArray<TB_PointMatch *> matches, const GU_Detail * gdp)
{
    TB_PointMatch *current;
    const GEO_Point  *ppt;
    int entries = matches.entries();
    UT_Spline  *spline;
    fpreal32 v1[3]; 
    int id;
    
    GEO_AttributeHandle handle = gdp->getAttribute(GEO_POINT_DICT, "id");
    
    for(int i = 0; i < gdp->points().entries(); i++)
    {
        handle.setElement(gdp->points()(i));
        id = handle.getI();
        
        spline  = new UT_Spline(); 
        spline->setGlobalBasis((UT_SPLINE_BASIS)itype);
        spline->setSize(entries, 3);
        
        for (int g = 0; g < entries; g++)
        {
            current = matches(g);
            ppt     = current->find(id);
            if (!ppt) ppt = gdp->points()(i);
            v1[0]  = ppt->getPos().x();
            v1[1]  = ppt->getPos().y();
            v1[2]  = ppt->getPos().z();
            spline->setValue(g, v1, 3);
        }
        interpolants.at(i) = spline;
    }
    valid = true;
    
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

