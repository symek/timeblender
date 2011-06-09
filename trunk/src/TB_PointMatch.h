#ifndef __TB_PointMatch_h__
#define __TB_PointMatch_h__

#include <map>
#include <GU/GU_Detail.h>
#include <UT/UT_SplayTree.h>
#include <GEO/GEO_Point.h>
#include <UT/UT_Vector3.h>

/// This is a basic infrastructure for finding correspondence between geometries. 
/// At first we'd like to match points via points' ids, but I'm thinking also 
/// about adding another options (see bellow TB_CORR_TYPES for a hint).
/// The main idea is to replace point access in Interpolants classes:
///     FOR_ALL_GPOINTS() {...
///             nextppt = gdp->points()p[i] ...}
///
/// with a method of TB_PointMatch:
///     TB_PointMatch * pm = TB_PointMatch(gdp);
///     pm->matchAttribute(id, *ppt);
/// 
/// where idx could be particle id for example.
///
/// Even better TB_PointMatch interface should allow to return
/// positions of a virtual point instead of an actual point, since
/// some of the methods may base on fitting displacements of non-equal
/// geometries.   

namespace TimeBlender
{
/// Types for building geometry correspondence.
/// Currently only CORR_POINT_ID is implemented.
typedef enum {
    CORR_POINT_ID,
    CORR_UV_SPACE,
    CORR_SPHERICAL_SPACE,
    CORR_TEATRA_MATCH,
} TB_CORR_TYPES;

/// Node typedef for STL map (red/black tree):
typedef std::map<int, GEO_Point*> Tree;

/// TODO: This should be plain structre?
class TB_SplayNode 
{
public:
    TB_SplayNode();
    TB_SplayNode(int i) {id=i;};
    TB_SplayNode(int i, const GEO_Point *p) {id=i; ppt=p;};  
    int id;
    const GEO_Point * ppt;
};

/// I'm removing Splay tree from here, since they're not optimial for
/// this kind of queries (with no coherence between sequantial searches)
/// and additionally not thread safe at access (important as we'd like to expose 
/// point match in VEX). Generally speaking using Splay Tree here had
/// no sense at all.

class TB_PointMatch
{
public:
    /// Requires initialization:
    TB_PointMatch() {};
    
    /// Deallocate splay tree on exit:
    ~TB_PointMatch() { delete tree;}
    
    /// Alloc wiht parameters:
    TB_PointMatch(GU_Detail * gdp, int correspond)
    {
        if (!initialize(gdp, correspond)) alloc = false;
    }
    /// Initilialize splay tree:
    int initialize(GU_Detail *, int correspond);
    
    /// Find corresponding point based on provided ID,
    /// and assing finding to GEO_Point *.
    GEO_Point * find(int d);
    
    /// Get to the raw tree.
    Tree * getTree() {return  tree;}
    
    /// Tree entries == npoints.
    int entries() const {return (int)tree->size();}
    
    /// Are we ready for searach.
    bool isAlloc() {return alloc;}
 
private:
    bool  alloc;
          Tree      * tree;
    const GU_Detail * detail;
};
} // End of Timeblender namespace
#endif
