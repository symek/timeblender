#ifndef __VRAY_TimeBlender_h__
#define __VRAY_TimeBlender_h__

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
#include <SYS/SYS_Math.h>
#include <UT/UT_WorkArgs.h>
#include <UT/UT_PtrArray.h>
#include <vector>

namespace TimeBlender
{
/// VRAY_IGeometry, the main worker.
class VRAY_TimeBlender: public VRAY_Procedural 
{
public:
        VRAY_TimeBlender();
virtual ~VRAY_TimeBlender();

    virtual const char *getClassName();
    virtual int         initialize(const UT_BoundingBox *);
    virtual void        getBoundingBox(UT_BoundingBox &box);
    virtual void        render();

private:
	int saveGeometry(const GU_Detail *, const UT_String *);

    UT_BoundingBox  myBox;
    fpreal          myshutter;
    int             mynsamples;
    int             myitype;
    fpreal          myshutterstart;
    fpreal          myshutterend;
    int             mycurrentframe;
    int             mymatchbyid;
    int             myfiles;
    UT_String       shop_materialpath;
    UT_String       myfilenamestring;
    UT_WorkArgs     myfilenamelist;
};
}//End of timeblender namescape

#endif
