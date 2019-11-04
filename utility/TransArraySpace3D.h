#ifndef TransArraySpace3D_H
#define TransArraySpace3D_H

/** method to convert a 1D field in a 3D field not considering guard cells*/
inline void solver2phys(double*** vectPhys, double* vectSolver, int nx, int ny, int nz){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++)
           vectPhys[i][j][k] = *vectSolver++;
          
}
/** method to convert a 1D field in a 3D field not considering guard cells*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double* vectSolver, int nx, int ny, int nz){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++){
          vectPhys1[i][j][k] = *vectSolver++;
          vectPhys2[i][j][k] = *vectSolver++;
          vectPhys3[i][j][k] = *vectSolver++;
        }
}
/** method to convert a 1D field in a 3D field not considering guard cells for FLUID simulations*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double**** vectJx, double**** vectJy, double**** vectJz,double**** rhoDATA, double* vectSolver, int nx, int ny, int nz, int ns){
    for (register int i=1; i < nx-1; i++)
        for (register int j=1; j < ny-1; j++)
            for (register int k=1; k < nz-1; k++){
                vectPhys1[i][j][k] = *vectSolver++;
                vectPhys2[i][j][k] = *vectSolver++;
                vectPhys3[i][j][k] = *vectSolver++;
                for (register int kk=0; kk < ns; kk++){
                    vectJx[kk][i][j][k] = *vectSolver++;
                    vectJy[kk][i][j][k] = *vectSolver++;
                    vectJz[kk][i][j][k] = *vectSolver++;
                    rhoDATA[kk][i][j][k] = *vectSolver++;
                }
            }
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys, int nx, int ny, int nz){
     for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++)
           *vectSolver++ = vectPhys[i][j][k];
           

}
/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,double*** vectPhys3, int nx, int ny, int nz){
     for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++){
           *vectSolver++ =  vectPhys1[i][j][k];
           *vectSolver++ =  vectPhys2[i][j][k];
           *vectSolver++ =  vectPhys3[i][j][k];
        }

}
/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double**** vectJx, double**** vectJy, double**** vectJz,double**** rhoDATA, int nx, int ny, int nz, int ns){
    for (register int i=1; i < nx-1; i++)
        for (register int j=1; j < ny-1; j++)
            for (register int k=1; k < nz-1; k++){
                *vectSolver++ =  vectPhys1[i][j][k];
                *vectSolver++ =  vectPhys2[i][j][k];
                *vectSolver++ =  vectPhys3[i][j][k];
                for (register int kk=0; kk < ns; kk++){
                   *vectSolver++ =  vectJx[kk][i][j][k];
                   *vectSolver++ =  vectJy[kk][i][j][k];
                   *vectSolver++ =  vectJz[kk][i][j][k];
                   *vectSolver++ =  rhoDATA[kk][i][j][k];
                }
            }
    
}
#endif
