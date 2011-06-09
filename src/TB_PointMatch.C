#include "TB_PointMatch.h"

using namespace TimeBlender;

int 
TB_PointMatch::initialize(GU_Detail * gdp, int correspond = 0)
{
    /// Init accelerator:
    tree = new Tree();
    GEO_AttributeHandle handle = gdp->getPointAttribute("id");
    std::map<int, GEO_Point*>::iterator it;
    GEO_Point * ppt;
    int id;
    
	/// Create red/black tree accelerator:
	FOR_ALL_GPOINTS(gdp, ppt)
	{
	    it = tree->end();
	    handle.setElement(ppt);
	    id = handle.getI();
	    tree->insert(it, std::pair<int, GEO_Point*>(id, ppt));
	}
	
	detail = gdp;
	alloc  = true; 
	return 1;
}

GEO_Point * TB_PointMatch::find(int id)
{
    if (alloc) 
        return tree->find(id)->second;
    else
        return NULL;
}

