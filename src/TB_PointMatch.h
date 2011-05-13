#ifndef __TB_PointMatch_h__
#define __TB_PointMatch_h__

#include <GU/GU_Detail.h>
#include <UT/UT_SplayTree.h>
#include <GEO/GEO_Point.h>
#include <UT/UT_Vector3.h>

/// This is a basic infrastructure for finding correspondence between geometries. 
/// At first we'd like to match points with variadic indices via points' "id".
/// I'm thinking also about adding another options...
/// The main idea is to replace point access in Interpolants classes:
///     FOR_ALL_GPOINTS() {...
///             nextppt = gdp->points()p[i] ...}
///
/// with method of TB_PointMatch:
///     TB_PointMatch * pm = TB_PointMatch(gdp);
///     pm->matchAttribute(id, *ppt);
/// 
/// where idx could be particle id for example.
///
/// Even better TB_PointMatch interface should allow to return
/// positions indead of point numbers.   

namespace TimeBlender
{

/// I'm not inherating after UT_SplayTree because
/// I really don't know where this class will go ahead. 
/// Perhpas later on this will be fixed. 
class TB_SplayNode 
{
public:
    TB_SplayNode();
    TB_SplayNode(int i) {id=i;};
    TB_SplayNode(int i, const GEO_Point *p) {id=i; ppt=p;};  
    int id;
    const GEO_Point * ppt;
};

class TB_PointMatch
{
public:
    TB_PointMatch() {};
    ~TB_PointMatch() {};
    int initialize(const GU_Detail *, int correspondence);
    UT_SplayTree * getTree() {return  tree;}
    TB_SplayNode * find(int id);
    GEO_Point    * findPoint(int id);
    GEO_Point    * get(int id);
    UT_Vector3   * getPos(int id);
    int entries() {return tree->entries();}
    bool isAlloc() {return alloc;}
 
private:
    bool alloc;
    UT_SplayTree    * tree;
    const GU_Detail * detail;
};

} // End of Timeblender namespace
#endif
