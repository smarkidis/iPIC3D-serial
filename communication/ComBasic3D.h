/***************************************************************************
    ComBasic.h  -  Library to handle Basic Communication
                       ------------------
 ***************************************************************************/

#ifndef ComBasic_H
#define ComBasic_H

#include "ComParser3D.h"


/** communicate particles along a direction **/
inline void communicateParticlesDIR(int buffer_size, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *b_right, double *b_left){
  
}
/** communicate ghost along a direction **/
inline void communicateGhostFace(int b_len, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightFace, double *ghostLeftFace){
  // just swap the buffer if you have just a1 processor in 1 direction
  if (right_neighbor == myrank && left_neighbor == myrank)
    swapBuffer(b_len,ghostLeftFace,ghostRightFace);
}
/** communicate ghost edge along a direction; there are 6 Diagonal directions through which we exchange Ghost Edges :
          0 = from   XrightYrightZsame to YleftZleftZsame; we exchange Z edge
          1 = from   XrightYleftZsame to XleftYrightZsame; we exchange Z edge
          2 = from   XrightYsameZright to XleftYsameZsame; we exchange Y edge
          3 = from   XrightYsameZleft to XleftYsameZright; we exchange Y edge
          4 = from   XsameYrightZright to XsameYleftZleft; we exchange X edge
          5 = from   XsameYrightZleft to XsameYleftZright; we exchange X edge
*/
inline void communicateGhostEdge(int b_len, int myrank, int right_neighborD, int left_neighborD, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightEdge, double *ghostLeftEdge){
 // swap
 if (right_neighborD == myrank && left_neighborD == myrank){
     swapBuffer(b_len,ghostLeftEdge,ghostRightEdge);
	 
  }
}
/** Communicate ghost corners along a direction; there are 4 Diagonal directions through which we exchange Ghost Corners:
          0 =  from XrightYrightZright to XleftYleftZleft
          1 =  from XrightYleftZright  to XleftYrightZleft
          2 =  from XleftYrightZright  to XrightYleftZleft
          3 =  from XleftYleftZright   to XrightYrightZleft
*/
inline void communicateGhostCorner(int myrank, int right_neighborC, int left_neighborC, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightCorner, double *ghostLeftCorner){
   // swap
  if (right_neighborC == myrank && left_neighborC == myrank){
     swapBuffer(1,ghostLeftCorner,ghostRightCorner);
  }
}

















#endif
