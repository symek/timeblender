#include "TB_PointMatch.h"

using namespace TimeBlender;

int
compareSplayNodes(const void * a, const void * b)
{
    if ( ((TB_SplayNode*)a)->id < ((TB_SplayNode*)b)->id )
    {
        return -1;
    }
    else if ( ((TB_SplayNode *)a)->id > ((TB_SplayNode *)b)->id ) 
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int 
TB_PointMatch::initialize(const GU_Detail * gdp, int correspond = 0)
{
    /// Init accelerator:
    tree = new UT_SplayTree(&compareSplayNodes);
    GEO_AttributeHandle handle = gdp->getAttribute(GEO_POINT_DICT, "id");
    const GEO_Point * ppt;
	TB_SplayNode * snode;
	
	/// Create accelerator:
	FOR_ALL_GPOINTS(gdp, ppt)
	{
	    handle.setElement(ppt);
	    snode = new TB_SplayNode((int)handle.getI());
	    snode->ppt = (const GEO_Point*)ppt;
	    tree->add(snode);
	}
	detail = gdp;
	alloc  = true; 
	return 1;
}

TB_SplayNode * TB_PointMatch::find(int id)
{
    if (alloc)
    {
        TB_SplayNode node = TB_SplayNode(id);
        return (TB_SplayNode *) tree->find(&node);
    }
}

GEO_Point * TB_PointMatch::findPoint(int id)
{
    if (alloc)
    {
        TB_SplayNode node = TB_SplayNode(id);
        return ((GEO_Point *)tree->find(&node)->ppt);
    }
}

GEO_Point * TB_PointMatch::get(int id)
{
    if (detail)
    {
        return detail->points()[id];
    }

}

