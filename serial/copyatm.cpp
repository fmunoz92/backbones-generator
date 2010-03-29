#include "copyatm.h"



void copyatm(const ATOM * source, ATOM * dest) 
{
   dest->x = source->x;
   dest->y = source->y;
   dest->z = source->z;
   dest->vdw = source->vdw;
}
