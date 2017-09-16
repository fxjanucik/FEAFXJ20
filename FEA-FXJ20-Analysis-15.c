
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//
//double MINV3X3(int ielm,int ig,double mj[3][3],double mji[3][3])
//void   MMULTAB(int ni,int nj,int nk,double c[ni][nj],double a[ni][nk],double b[nk][nj])
//void   MMULTATB(int ni,int nj,int nk,double c[ni][nj],double a[nk][ni],double b[nk][nj])
//void   E20PARTS(int ip,double gp[27][3],double spqr[3][20])
//void   E20SHAPES(double pqr[],double sn[])
//void   GAUSSQUAD(int norder,double wt[],double gp[27][3])
//int    ij2ipt(int i,int j,int nsize)
//void   BTDB(int etype,double ke[60][60],double mat[3],double sxyz[3][20],double dv)
//void   S08PARTS(int ip,double gp[],double s8[8],double s8pq[2][8])
//void   E20PRESS(int ielm,int iface,double press,int element[],double mxyz[],int nf[],double ff[])
//void   ASSEMBLEK(int ielm,int element[],int etype,double ks[],double ke[60][60],int ssize)
//void   CHOLESKYLDLT(int nnode,double ks[])
//
// **********************************************************
double
MINV3X3(int ielm, int ig, double mj[3][3], double mji[3][3])
// **********************************************************
//  3 x 3 matrix inverse and determinant
//  mj = mj[0][0] mj[0][1] mj[0][2]
//       mj[1][0] mj[1][1] mj[1][2]
//       mj[2][0] mj[2][1] mj[2][2]
// **********************************************************
{
int i,j;
double detj,rdetj;
//
mji[0][0] = mj[1][1]*mj[2][2]-mj[2][1]*mj[1][2];
mji[0][1] = mj[2][0]*mj[1][2]-mj[1][0]*mj[2][2];
mji[0][2] = mj[1][0]*mj[2][1]-mj[2][0]*mj[1][1];

detj =  mj[0][0]*mji[0][0]+mj[0][1]*mji[0][1]+mj[0][2]*mji[0][2];

if(detj > 0.){
  mji[1][0] = mj[2][1]*mj[0][2]-mj[0][1]*mj[2][2];
  mji[1][1] = mj[0][0]*mj[2][2]-mj[2][0]*mj[0][2];
  mji[1][2] = mj[2][0]*mj[0][1]-mj[0][0]*mj[2][1];
  
  mji[2][0] = mj[0][1]*mj[1][2]-mj[1][1]*mj[0][2];
  mji[2][1] = mj[1][0]*mj[0][2]-mj[0][0]*mj[1][2];
  mji[2][2] = mj[0][0]*mj[1][1]-mj[1][0]*mj[0][1];
  
  rdetj = mji[1][0];
  mji[1][0] = mji[0][1];
  mji[0][1] = rdetj;
  
  rdetj = mji[2][0];
  mji[2][0] = mji[0][2];
  mji[0][2] = rdetj;  
  
  rdetj = mji[2][1];
  mji[2][1] = mji[1][2];
  mji[1][2] = rdetj;  
  
  rdetj = 1./detj;
  for(i=0;i<3; i++)
   for(j=0; j<3; j++)
     mji[i][j]=rdetj*mji[i][j];
  }
else{
  printf(" Warning: Negative or zero determinant for element %d gauss pt %d  = %8.2f \n",ielm,ig,detj);
  }
return detj;
}

// ***********************************************************************************
void
MMULTAB(int ni, int nj, int nk, double c[ni][nj], double a[ni][nk], double b[nk][nj])
// ***********************************************************************************
// matrix multiply c[ni][nj] = a[ni][nk] b[nk][nj]
// MMULTAB(3,3,20,jacob,spqr,xyz);
// ***********************************************************************************
{
int i,j,k;
for(i=0; i<ni; i++){
   for(j=0; j<nj; j++){
      c[i][j] = 0.;
      for(k=0;k<nk; k++){
         c[i][j] = c[i][j] + a[i][k]*b[k][j];
         }
      }
   }
return;
}

// ************************************************************************************
void
MMULTATB(int ni, int nj, int nk, double c[ni][nj], double a[nk][ni], double b[nk][nj])
// ************************************************************************************
// matrix multiply c[ni][nj] = a[nk][ni](Transposed) b[nk][nj]
//                           = a[ni][nk] b[nk][nj]            
// ************************************************************************************
{
int i,j,k;
for(i=0; i<ni; i++){
   for(j=0; j<nj; j++){
      c[i][j] = 0.;
      for(k=0;k<nk; k++){
         c[i][j] = c[i][j] + a[k][i]*b[k][j];
         }
      }
   }
return;
}

// ****************************************************************
void
E20PARTS(int ip, double gp[27][3], double spqr[3][20])
// ****************************************************************
// evaluate partial derivatives of shape functions in curvilinear
// coordinate space at gauss point coordinates pqr[3]
// ****************************************************************
//           7------ 14 ------6
//          /|               /|
//         / |              / |
//       15  |            13  |    +Q
//       /  19            /  18      |
//      /    |           /    |      |
//     4--------12------5     |      |
//     |     |          |     |      *-----> +P
//     |     3------10--|-----2     /
//     |    /           |    /     /
//    16   /           17   /    +R
//     | 11             |  9
//     | /              | /
//     |/               |/
//     0------- 8 ------1
// x = 0-1
// y = 0-3
// z = 0-12
// ****************************************************************
{
int i;
double p,q,r;
double opp,omp,ompp,opq,omq,omqq,opr,omr,omrr;

p    = gp[ip][0];
q    = gp[ip][1];
r    = gp[ip][2];

opp  = 1.+p;
omp  = 1.-p;
ompp = 1.-p*p;

opq  = 1.+q;
omq  = 1.-q;
omqq = 1.-q*q;

opr  = 1.+r;
omr  = 1.-r;
omrr = 1.-r*r;

// Shape Function Partials wrt p
spqr[0][ 0] =  .125*omq*opr*( 2.*p+q-r+1.);
spqr[0][ 1] =  .125*omq*opr*( 2.*p-q+r-1.);
spqr[0][ 2] =  .125*omq*omr*( 2.*p-q-r-1.);
spqr[0][ 3] =  .125*omq*omr*( 2.*p+q+r+1.);
spqr[0][ 4] =  .125*opq*opr*( 2.*p-q-r+1.);
spqr[0][ 5] =  .125*opq*opr*( 2.*p+q+r-1.);
spqr[0][ 6] =  .125*opq*omr*( 2.*p+q-r-1.);
spqr[0][ 7] =  .125*opq*omr*( 2.*p-q+r+1.);
spqr[0][ 8] = -.50*p*omq*opr;
spqr[0][ 9] =  .25*omrr*omq;
spqr[0][10] = -.50*p*omq*omr;
spqr[0][11] = -.25*omrr*omq;
spqr[0][12] = -.50*p*opq*opr;
spqr[0][13] =  .25*omrr*opq;
spqr[0][14] = -.50*p*opq*omr;
spqr[0][15] = -.25*omrr*opq;
spqr[0][16] = -.25*omqq*opr;
spqr[0][17] =  .25*omqq*opr;
spqr[0][18] =  .25*omqq*omr;
spqr[0][19] = -.25*omqq*omr;

// Shape Function Partials wrt q
spqr[1][ 0] =  .125*omp*opr*( p+2.*q-r+1.);
spqr[1][ 1] =  .125*opp*opr*(-p+2.*q-r+1.);
spqr[1][ 2] =  .125*opp*omr*(-p+2.*q+r+1.);
spqr[1][ 3] =  .125*omp*omr*( p+2.*q+r+1.);
spqr[1][ 4] =  .125*omp*opr*(-p+2.*q+r-1.);
spqr[1][ 5] =  .125*opp*opr*( p+2.*q+r-1.);
spqr[1][ 6] =  .125*opp*omr*( p+2.*q-r-1.);
spqr[1][ 7] =  .125*omp*omr*(-p+2.*q-r-1.);
spqr[1][ 8] = -.25*ompp*opr;
spqr[1][ 9] = -.25*omrr*opp;
spqr[1][10] = -.25*ompp*omr;
spqr[1][11] = -.25*omrr*omp;
spqr[1][12] =  .25*ompp*opr;
spqr[1][13] =  .25*omrr*opp;
spqr[1][14] =  .25*ompp*omr;
spqr[1][15] =  .25*omrr*omp;
spqr[1][16] = -.50*q*omp*opr;
spqr[1][17] = -.50*q*opp*opr;
spqr[1][18] = -.50*q*opp*omr;
spqr[1][19] = -.50*q*omp*omr;

// Shape Function Partials wrt r
spqr[2][ 0] =  .125*omp*omq*(-p-q+2.*r-1.);
spqr[2][ 1] =  .125*opp*omq*( p-q+2.*r-1.);
spqr[2][ 2] =  .125*opp*omq*(-p+q+2.*r+1.);
spqr[2][ 3] =  .125*omp*omq*( p+q+2.*r+1.);
spqr[2][ 4] =  .125*omp*opq*(-p+q+2.*r-1.);
spqr[2][ 5] =  .125*opp*opq*( p+q+2.*r-1.);
spqr[2][ 6] =  .125*opp*opq*(-p-q+2.*r+1.);
spqr[2][ 7] =  .125*omp*opq*( p-q+2.*r+1.);
spqr[2][ 8] =  .25*ompp*omq;
spqr[2][ 9] = -.50*r*opp*omq;
spqr[2][10] = -.25*ompp*omq;
spqr[2][11] = -.50*r*omp*omq;
spqr[2][12] =  .25*ompp*opq;
spqr[2][13] = -.50*r*opp*opq;
spqr[2][14] = -.25*ompp*opq;
spqr[2][15] = -.50*r*omp*opq;
spqr[2][16] =  .25*omqq*omp;
spqr[2][17] =  .25*omqq*opp;
spqr[2][18] = -.25*omqq*opp;
spqr[2][19] = -.25*omqq*omp;

return;
}

// ***********************************************
void
E20SHAPES(double pqr[], double sn[])
// ***********************************************
// Evaluate 20 node element shape functions
// at gauss point coordinates pqr[3]
// ***********************************************
//           7------ 14 ------6
//          /|               /|
//         / |              / |
//       15  |            13  |    +Q
//       /  19            /  18      |
//      /    |           /    |      |
//     4--------12------5     |      |
//     |     |          |     |      *-----> +P
//     |     3------10--|-----2     /
//     |    /           |    /     /
//    16   /           17   /    +R
//     | 11             |  9
//     | /              | /
//     |/               |/
//     0------- 8 ------1
// ***********************************************
{
int i;
double p,q,r;
double opp,omp,ompp,opq,omq,omqq,opr,omr,omrr;

p    = pqr[0];
q    = pqr[1];
r    = pqr[2];

opp  = 1.+p;
omp  = 1.-p;
ompp = 1.-p*p;

opq  = 1.+q;
omq  = 1.-q;
omqq = 1.-q*q;

opr  = 1.+r;
omr  = 1.-r;
omrr = 1.-r*r;

// Shape Functions
sn[ 0] = .125*omp*omq*opr*(-p-q+r-2.);  // (-1,-1,+1)
sn[ 1] = .125*opp*omq*opr*( p-q+r-2.);  // (+1,-1,+1)
sn[ 2] = .125*opp*omq*omr*( p-q-r-2.);  // (+1,-1,-1)
sn[ 3] = .125*omp*omq*omr*(-p-q-r-2.);  // (-1,-1,-1)
sn[ 4] = .125*omp*opq*opr*(-p+q+r-2.);  // (-1,+1,+1)
sn[ 5] = .125*opp*opq*opr*( p+q+r-2.);  // (+1,+1,+1)
sn[ 6] = .125*opp*opq*omr*( p+q-r-2.);  // (+1,+1,-1)
sn[ 7] = .125*omp*opq*omr*(-p+q-r-2.);  // (-1,+1,-1)
sn[ 8] = .250*ompp*omq*opr;             // ( 0,-1,+1)
sn[ 9] = .250*opp*omq*omrr;             // (+1,-1, 0)
sn[10] = .250*ompp*omq*omr;             // ( 0,-1,-1)
sn[11] = .250*omp*omq*omrr;             // (-1,-1, 0)
sn[12] = .250*ompp*opq*opr;             // ( 0,+1,+1)
sn[13] = .250*opp*opq*omrr;             // (+1,+1, 0)
sn[14] = .250*ompp*opq*omr;             // ( 0,+1,-1)
sn[15] = .250*omp*opq*omrr;             // (-1,+1, 0)
sn[16] = .250*omp*omqq*opr;             // (-1, 0,+1)
sn[17] = .250*opp*omqq*opr;             // (+1, 0,+1)
sn[18] = .250*opp*omqq*omr;             // (+1, 0,-1)
sn[19] = .250*omp*omqq*omr;             // (-1, 0,-1)

return;
}

// ***********************************************
void
E08SHAPES(double pqr[], double sn[])
// ***********************************************
// Evaluate 8 node element shape functions
// at gauss point coordinates pqr[3]
// ***********************************************
//           7----------------6
//          /|               /|
//         / |              / |
//        /  |             /  |    +Q
//       /   |            /   |      |
//      /    |           /    |      |
//     4----------------5     |      |
//     |     |          |     |      *-----> +P
//     |     3----------|-----2     /
//     |    /           |    /     /
//     |   /            |   /    +R
//     |  /             |  /
//     | /              | /
//     |/               |/
//     0----------------1
// ***********************************************
{
int i;
double p,q,r,opp,omp,opq,omq,opr,omr;

p    = pqr[0];
q    = pqr[1];
r    = pqr[2];

opp  = 1.+p;
omp  = 1.-p;

opq  = 1.+q;
omq  = 1.-q;

opr  = 1.+r;
omr  = 1.-r;

// 8 Node 3D Isoparametric Shape Functions
sn[ 0] = .125*omp*omq*opr;  // (-1,-1,+1)
sn[ 1] = .125*opp*omq*opr;  // (+1,-1,+1)
sn[ 2] = .125*opp*omq*omr;  // (+1,-1,-1)
sn[ 3] = .125*omp*omq*omr;  // (-1,-1,-1)
sn[ 4] = .125*omp*opq*opr;  // (-1,+1,+1)
sn[ 5] = .125*opp*opq*opr;  // (+1,+1,+1)
sn[ 6] = .125*opp*opq*omr;  // (+1,+1,-1)
sn[ 7] = .125*omp*opq*omr;  // (-1,+1,-1)

return;
}

// **********************************************************************************
void
GAUSSQUAD(int norder, double wt[], double gp[27][3])
// **********************************************************************************
//                          *                           *                 
//                          * Q=.7746  7-----14------6  * Q= +.57735 7-----------6
//      7-------14-------6  *          |             |  *            |           |
//     /|    +Q         /|  *         15     21-->  13  *            |     *---> |
//   15 |     |       13 |  *          |      |  +P  |  *            |     |  +P |   
//   /  |     |       /  |  *          4-----12------5  *            4-----|-----5      
//  /  19     |      /   |  *               +R|         *                +R|        
// 4------12--------5  18|  *                           *
// |    |     |     |    |  * Q=0.    19-----24-----18  * Q= -.57735 3-----------2      
// |    |     *-----|> +P|  *          |             |  *            |           |      
// |    3----/-10---|----2  *         25     26-->  23  *            |     *---> |     
//16   /    /      17   /   *          |      |  +P  |  *            |     |  +P |     
// | 11   +R        |  9    *         16-----22-----17  *            0-----|-----1       
// | /              | /     *               +R|         *                +R|
// |/               |/      *                           *   
// 0-------8--------1       * Q=-.7746 3-----10------2  * 
//                          *          |             |  *
//                          *         11     20-->   9  * 
//                          *          |      |  +P  |  * 
//                          *          0------8------1  *  
//                          *               +R|         *
// **********************************************************************************
{
int i,j,k,ic;
double pm, p0, pp;
double wm, w0, wp;
static double w3[3] = { .5555555555555556,.8888888888888889,.5555555555555556};
static double g3[3] = {-.7745966692414834,0.               ,.7745966692414834};

if(norder != 3){
// 2x2x2 gauss quadrature weights and abscissae
    for(i=0; i<8; i++)wt[i] = 1.;
    
    pm = -.5773502691896257;
    pp =  .5773502691896257;
    
    gp[0][0] = pm;
    gp[0][1] = pm;
    gp[0][2] = pp;
    
    gp[1][0] = pp;
    gp[1][1] = pm;
    gp[1][2] = pp;
    
    gp[2][0] = pp;
    gp[2][1] = pm;
    gp[2][2] = pm;
    
    gp[3][0] = pm;
    gp[3][1] = pm;
    gp[3][2] = pm;
               
    gp[4][0] = pm;
    gp[4][1] = pp;
    gp[4][2] = pp;
    
    gp[5][0] = pp;
    gp[5][1] = pp;
    gp[5][2] = pp;
    
    gp[6][0] = pp;
    gp[6][1] = pp;
    gp[6][2] = pm;
    
    gp[7][0] = pm;
    gp[7][1] = pp;
    gp[7][2] = pm;             
    }
else{
    // 3x3x3 gauss quadrature weights and abscissae
    pm = -.7745966692414834;
    p0 = 0.;
    pp =  .7745966692414834;
    
    wm = .5555555555555556;
    w0 = .8888888888888889;
    wp = .5555555555555556;

    wt[ 0]    = wm*wm*wp;
    gp[ 0][0] = pm;
    gp[ 0][1] = pm;
    gp[ 0][2] = pp;
    
    wt[ 1]    = wp*wm*wp;
    gp[ 1][0] = pp;
    gp[ 1][1] = pm;
    gp[ 1][2] = pp;
        
    wt[ 2]    = wp*wm*wm;    
    gp[ 2][0] = pp;
    gp[ 2][1] = pm;
    gp[ 2][2] = pm;
        
    wt[ 3]    = wm*wm*wm;    
    gp[ 3][0] = pm;
    gp[ 3][1] = pm;
    gp[ 3][2] = pm;
        
    wt[ 4]    = wm*wp*wp;    
    gp[ 4][0] = pm;
    gp[ 4][1] = pp;
    gp[ 4][2] = pp;
        
    wt[ 5]    = wp*wp*wp;    
    gp[ 5][0] = pp;
    gp[ 5][1] = pp;
    gp[ 5][2] = pp;
        
    wt[ 6]    = wp*wp*wm;    
    gp[ 6][0] = pp;
    gp[ 6][1] = pp;
    gp[ 6][2] = pm;
            
    wt[ 7]    = wm*wp*wm;    
    gp[ 7][0] = pm;
    gp[ 7][1] = pp;
    gp[ 7][2] = pm;
            
    wt[ 8]    = w0*wm*wp;    
    gp[ 8][0] = p0;
    gp[ 8][1] = pm;
    gp[ 8][2] = pp;
        
    wt[ 9]    = wp*wm*w0;    
    gp[ 9][0] = pp;
    gp[ 9][1] = pm;
    gp[ 9][2] = p0;
        
    wt[10]    = w0*wm*wm;    
    gp[10][0] = p0;
    gp[10][1] = pm;
    gp[10][2] = pm;
        
    wt[11]    = wm*wm*w0;    
    gp[11][0] = pm;
    gp[11][1] = pm;
    gp[11][2] = p0;
        
    wt[12]    = w0*wp*wp;    
    gp[12][0] = p0;
    gp[12][1] = pp;
    gp[12][2] = pp;
        
    wt[13]    = wp*wp*w0;    
    gp[13][0] = pp;
    gp[13][1] = pp;
    gp[13][2] = p0;
        
    wt[14]    = w0*wp*wm;    
    gp[14][0] = p0;
    gp[14][1] = pp;
    gp[14][2] = pm;
        
    wt[15]    = wm*wp*w0;    
    gp[15][0] = pm;
    gp[15][1] = pp;
    gp[15][2] = p0;
            
    wt[16]    = wm*w0*wp;    
    gp[16][0] = pm;
    gp[16][1] = p0;
    gp[16][2] = pp;
        
    wt[17]    = wp*w0*wp;    
    gp[17][0] = pp;
    gp[17][1] = p0;
    gp[17][2] = pp;
        
    wt[18]    = wp*w0*wm;    
    gp[18][0] = pp;
    gp[18][1] = p0;
    gp[18][2] = pm;
        
    wt[19]    = wm*w0*wm;    
    gp[19][0] = pm;
    gp[19][1] = p0;
    gp[19][2] = pm;
        
    wt[20]    = w0*wm*w0;    
    gp[20][0] = p0;
    gp[20][1] = pm;
    gp[20][2] = p0;
        
    wt[21]    = w0*wp*w0;    
    gp[21][0] = p0;
    gp[21][1] = pp;
    gp[21][2] = p0;
        
    wt[22]    = w0*w0*wp;    
    gp[22][0] = p0;
    gp[22][1] = p0;
    gp[22][2] = pp;
        
    wt[23]    = wp*w0*w0;    
    gp[23][0] = pp;
    gp[23][1] = p0;
    gp[23][2] = p0;
        
    wt[24]    = w0*w0*wm;    
    gp[24][0] = p0;
    gp[24][1] = p0;
    gp[24][2] = pm;
            
    wt[25]    = wm*w0*w0;    
    gp[25][0] = pm;
    gp[25][1] = p0;
    gp[25][2] = p0;
        
    wt[26]    = w0*w0*w0;    
    gp[26][0] = p0;
    gp[26][1] = p0;
    gp[26][2] = p0;
    }
return;
}

// ************************************************************************************
int
ij2ipt(int i, int j, int nsize)
// ************************************************************************************
// map (nsize x nsize) square matrix (i,j) to ipt in a lower triangle matrix           
// ************************************************************************************
{
int ipt;
ipt = .5*(j*(2*nsize-j-1)+2*i);
return ipt;
}

// ************************************************************************************
void
printkblock(int ielm, int ir, int ic, double kij[3][3])
// ************************************************************************************
{
printf(" Element stiffness matrix block for element %d ir %d ic %d \n",ielm,ir,ic);
printf(" %8.2e  %8.2e %8.2e \n",kij[0][0],kij[0][1],kij[0][2]);
printf(" %8.2e  %8.2e %8.2e \n",kij[1][0],kij[1][1],kij[1][2]);
printf(" %8.2e  %8.2e %8.2e \n",kij[2][0],kij[2][1],kij[2][2]);
return;
}

// ******************************************************************************
void
BTDB(int etype, double ke[60][60], double mat[3], double sxyz[3][20], double dv)
// ******************************************************************************
// calculate lower triangle of element stiffness matrix
// by computing and summing kij[3][3] stiffness blocks for [row][column] [i][j]
// ******************************************************************************
{
int i,j,ir,jc,kr,kc,ipt,nsize;
double bbi[3][3],dbbbj[3][3],kij[3][3];

// calculate lower triangle 3x3 stiffness matrix blocks kij[3][3]

// loop on element rows (ir = 0 - 19)
for(ir=0; ir<etype; ir++){
   kr = 3*ir;
   // loop on element columns (jc = 0 - i)
  for(jc=0; jc<(ir+1); jc++){
    kc = 3*jc;

   // initialize stiffness block (i,j) to 0.
   for(i=0; i<3; i++)
     for(j=0; j<3; j++)
       kij[i][j] = 0.;
       
    // define bbi[3][3]
    bbi[0][0] = .5*sxyz[1][ir];
    bbi[0][1] = .5*sxyz[0][ir];
    bbi[0][2] =  0.;
    
    bbi[1][0] =  0.;
    bbi[1][1] = .5*sxyz[2][ir];
    bbi[1][2] = .5*sxyz[1][ir];
    
    bbi[2][0] = .5*sxyz[2][ir];
    bbi[2][1] =  0.;
    bbi[2][2] = .5*sxyz[0][ir];
    
    //define dbbbj[3][3]
    dbbbj[0][0] = .5*sxyz[1][jc]*mat[2];
    dbbbj[0][1] = .5*sxyz[0][jc]*mat[2];
    dbbbj[0][2] =  0.;
    
    dbbbj[1][0] =  0.;
    dbbbj[1][1] = .5*sxyz[2][jc]*mat[2];
    dbbbj[1][2] = .5*sxyz[1][jc]*mat[2];
    
    dbbbj[2][0] = .5*sxyz[2][jc]*mat[2];
    dbbbj[2][1] =  0.;
    dbbbj[2][2] = .5*sxyz[0][jc]*mat[2];
 
    // bbi[3][3](T)*dbbbj[3][3]
    MMULTATB(3,3,3,kij,bbi,dbbbj);
    
    // add in bai[3][3](T)DA[3][3]baj[3][3] equivalent
    kij[0][0] = kij[0][0] + mat[0]*sxyz[0][ir]*sxyz[0][jc];
    kij[0][1] = kij[0][1] + mat[1]*sxyz[0][ir]*sxyz[1][jc];
    kij[0][2] = kij[0][2] + mat[1]*sxyz[0][ir]*sxyz[2][jc];
    
    kij[1][0] = kij[1][0] + mat[1]*sxyz[1][ir]*sxyz[0][jc];
    kij[1][1] = kij[1][1] + mat[0]*sxyz[1][ir]*sxyz[1][jc];
    kij[1][2] = kij[1][2] + mat[1]*sxyz[1][ir]*sxyz[2][jc];
    
    kij[2][0] = kij[2][0] + mat[1]*sxyz[2][ir]*sxyz[0][jc];
    kij[2][1] = kij[2][1] + mat[1]*sxyz[2][ir]*sxyz[1][jc];
    kij[2][2] = kij[2][2] + mat[0]*sxyz[2][ir]*sxyz[2][jc];

//    if(ir == jc)printkblock(-1,ir,jc,kij);
            
    // sum [kt] into lower triangle of ke[60][60]
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        ke[kr+i][kc+j] = ke[kr+i][kc+j] + dv*kij[i][j];
    }
  }
return;
}

// ******************************************************************************
void
BTDBEXP(int etype, double ke[60][60], double d[6][6], double sxyz[3][20], double dv)
// ******************************************************************************
// calculate full element stiffness matrix [ke] = [B](T)[D][B]dv
// ******************************************************************************
{
int i,j,jc,kr,kc,ipt,nsize;
double b[6][60],db[6][60],kt[60][60];

// Define [B]
for(i=0; i<6; i++)
   for(j=0; j<60; j++)
      b[i][j] = 0.;

for(j=0; j<etype; j++){
   jc = 3*j;
   b[0][jc  ] =    sxyz[0][j];
   b[1][jc+1] =    sxyz[1][j];
   b[2][jc+2] =    sxyz[2][j];
   b[3][jc  ] = .5*sxyz[1][j];
   b[3][jc+1] = .5*sxyz[0][j];
   b[4][jc+1] = .5*sxyz[2][j];
   b[4][jc+2] = .5*sxyz[1][j];
   b[5][jc  ] = .5*sxyz[2][j];
   b[5][jc+2] = .5*sxyz[0][j];
   }
// Calculate [DB] = [D][B]
// *************************************************************************************
// MMULTAB(int ni, int nj, int nk, double c[ni][nj], double a[ni][nk], double b[nk][nj])
// matrix multiply c[ni][nj] = a[ni][nk] b[nk][nj]
// *************************************************************************************
MMULTAB(6,60,6,db,d,b);

// Calculate [kt] = [B](T)[DB]
// ************************************************************************************
// MMULTATB(int ni, int nj, int nk, double c[ni][nj], double a[nk][ni], double b[nk][nj])
// matrix multiply c[ni][nj] = a[nk][ni](Transposed) b[nk][nj]
//                           = a[ni][nk] b[nk][nj]            
// ************************************************************************************
MMULTATB(60,60,6,kt,b,db);
  
// sum [kt] into lower triangle of [ke]
for(i=0; i<60; i++)
    for(j=0; j<(i+1); j++)
      ke[i][j] = ke[i][j] + kt[i][j]*dv;
return;
}

// ****************************************************************
void
S08PARTS(int ip, double gp[], double s8[8], double s8pq[2][8])
// ****************************************************************
// evaluate shape functions and partial derivatives of 2D 8 node surface
// in curvilinear coordinate space at gauss point coordinates (p,q)
// ****************************************************************
{
double p,q;
double opp,omp,ompp,opq,omq,omqq;

p    = gp[2*ip  ];
q    = gp[2*ip+1];

opp  = 1.+p;
omp  = 1.-p;
ompp = 1.-p*p;

opq  = 1.+q;
omq  = 1.-q;
omqq = 1.-q*q;

// Shape Functions
s8[0] = .25*omp*omq*(-p-q-1);
s8[1] = .25*opp*omq*( p-q-1);
s8[2] = .25*opp*opq*( p+q-1);
s8[3] = .25*omp*opq*(-p+q-1);
s8[4] = .5*ompp*omq;
s8[5] = .5*opp*omqq;
s8[6] = .5*ompp*opq;
s8[7] = .5*omp*omqq;

// Shape Function Partials wrt p
s8pq[0][0] = .25*omq*( 2.*p+q);
s8pq[0][1] = .25*omq*( 2.*p-q);
s8pq[0][2] = .25*opq*( 2.*p+q);
s8pq[0][3] = .25*opq*( 2.*p-q);
s8pq[0][4] =  -p*omq;
s8pq[0][5] =  .5*omqq;
s8pq[0][6] =  -p*opq;
s8pq[0][7] = -.5*omqq;

// Shape Function Partials wrt q
s8pq[1][0] = .25*omp*( p+2.*q);
s8pq[1][1] = .25*opp*(-p+2.*q);
s8pq[1][2] = .25*opp*( p+2.*q);
s8pq[1][3] = .25*omp*(-p+2.*q);
s8pq[1][4] = -.5*ompp;
s8pq[1][5] =  -q*opp;
s8pq[1][6] =  .5*ompp;
s8pq[1][7] =  -q*omp;

return;
}

// **************************************************************************************
void
E20PRESS(int ielm,int iface,double press,int element[],double mxyz[],int nf[],double ff[])
// **************************************************************************************
// calculate forces due to face normal pressure press
// in curvilinear coordinate space at gauss point coordinates pq[2]
// **************************************************************************************
//         7------ 14 ------6               5----13-----6  7----15-----4   7----14-----6
//        /|               /|               |           |  |           |   |           |
//       / |              / |              17    +P    18 19    -P    16  15    +Q    13 
//     15  |            13  |  +Q           |           |  |           |   |           |
//     /  19            /  18    |          |           |  |           |   |           |
//    /    |           /    |    |          1-----9-----2  3----11-----0   4----12-----5
//   4--------12------5     |    |             Face 0          Face 1          Face 2 
//   |     |          |     |    *----> +P  0-----8-----1  4----12-----5   6----14-----7
//   |     3------10--|-----2   /           |           |  |           |   |           |
//   |    /           |    /   /            |           |  |           |   |           |
//  16   /           17   /  +R            11    -Q     9 16    +R    17  18    -R    19
//   | 11             |  9                  |           |  |           |   |           |
//   | /              | /                   |           |  |           |   |           |
//   |/               |/                    3----10-----2  0-----8-----1   2----10-----3
//   0------- 8 ------1                        Face 3          Face 4          Face 5
// **************************************************************************************
{
int    i,j,en,mn,ipt;
double fxyz[3][8],s8[8],s8pq[2][8],v1[3],v2[3],v3[3],dl1,dl2,ds;

static int nface[6][8] = {  1, 2, 6, 5,  9, 18, 13, 17,
                            3, 0, 4, 7, 11, 16, 15, 19,
                            4, 5, 6, 7, 12, 13, 14, 15, 
                            3, 2, 1, 0, 10,  9,  8, 11,
                            0, 1, 5, 4,  8, 17, 12, 16,  
                            2, 3, 7, 6, 10, 19, 14, 18};                    
static double gp[8]    = {-.5773502691896257,-.5773502691896257,
                           .5773502691896257,-.5773502691896257,
                           .5773502691896257, .5773502691896257,
                          -.5773502691896257, .5773502691896257};
// get element face coordinates
for(i=0; i<8; i++){
   en    = nface[iface][i];
   mn    = element[20*ielm + en];
   nf[i] = mn;                       // face nodes having pressure based forces
   for(j=0; j<3; j++)
     fxyz[j][i] = mxyz[3*mn+j];
   }
// initialize face pressure forces
for(i=0; i<24; i++)ff[i] = 0.;

// loop over 2x2 gauss points
for(ipt=0; ipt<4; ipt++){
   // evaluate face shape function partial derivatives at curvilinear coordinates ip
   S08PARTS(ipt,gp,s8,s8pq);
   
   // get tangent vectors v1[3] and v2[3]
   for(i=0; i<3; i++){
     v1[i] = 0.;
     for(j=0; j<8; j++)
       v1[i] = v1[i] + s8pq[0][j]*fxyz[i][j];
     }
   for(i=0; i<3; i++){
     v2[i] = 0.;
     for(j=0; j<8; j++)
       v2[i] = v2[i] + s8pq[1][j]*fxyz[i][j];
     } 
   // v3[3] = v1[3] x v2[3]
   v3[0] = v1[1]*v2[2]-v2[1]*v1[2];  
   v3[1] = v1[2]*v2[0]-v2[2]*v1[0];       
   v3[2] = v1[0]*v2[1]-v2[0]*v1[1];
        
   ds = sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]);
   
   // sum qauss point force into face force vector
   for(i=0; i<8; i++){
      ff[3*i  ] = ff[3*i  ] + press*v3[0]*s8[i];
      ff[3*i+1] = ff[3*i+1] + press*v3[1]*s8[i];
      ff[3*i+2] = ff[3*i+2] + press*v3[2]*s8[i];       
      }
   }
return;
}

// ***************************************************************************************
void
ASSEMBLEK(int ielm, int element[], int etype, double ks[], double ke[60][60], int ssize)
// ***************************************************************************************
// assemble element ke into model lower triangle model stiffness ks
// using lower triangle (i,j) kb[3][3] blocks from ke[60][60] 
// ***************************************************************************************
//          n0         n1         n2         n3         n4     ...    n19
//  n0  kb[ n0,n0]
//  n1  kb[ n1,n0] kb[ n1,n1]
//  n2  kb[ n2,n0] kb[ n2,n1] kb[ n2,n2]
//  n3  kb[ n3,n0] kb[ n3,n1] kb[ n3,n2] kb[ n3,n3]
//  n4  kb[ n4,n0] kb[ n4,n1] kb[ n4,n2] kb[ n4,n3] kb[ n4,n4]
// ...     ...        ...        ...        ...        ...
// n19  kb[n19,n0] kb[n19,n1] kb[n19,n2] kb[n19,n3] kb[n19,n4] ... kb[n19,n19]
//
{
int i,j,ir,jc,kr,kc,nr,nc,ipt;
int row,col;
double kb[3][3],temp;

// loop over 3x3 ke[ir][jc] blocks
                                      
// loop on element ielm node row blocks (ir = 0 - etype)
for(ir=0; ir<etype; ir++){ 
    kr = 3*ir;                        // row start location in ke
    nr = element[etype*ielm+ir];      // row start location in ks
   
    // loop on element ielm node columns (jc = 0 - ir)
    for(jc=0; jc<(ir+1); jc++){  
        kc = 3*jc;                    // column start location in ke
        nc = element[etype*ielm+jc];  // column start location in ks
        
        // put ke 3x3 matrix block at (kr,kc) in kb[3][3]
        for(i=0; i<3; i++)
           for(j=0; j<3; j++)
              kb[i][j] = ke[kr+i][kc+j];       
        // add kb into ks
        if(nr == nc){
          // diagonal case - add lower triangle of block only
          ipt = ij2ipt(3*nr  ,3*nc  ,ssize);
          ks[ipt  ] = ks[ipt  ] + kb[0][0];
          ks[ipt+1] = ks[ipt+1] + kb[1][0];
          ks[ipt+2] = ks[ipt+2] + kb[2][0];
       
          ipt = ij2ipt(3*nr+1,3*nc+1,ssize);    
          ks[ipt  ] = ks[ipt  ] + kb[1][1];
          ks[ipt+1] = ks[ipt+1] + kb[2][1];
    
          ipt = ij2ipt(3*nr+2,3*nc+2,ssize);     
          ks[ipt  ] = ks[ipt  ] + kb[2][2];
          }
        else{
          // off-diagonal case - add full block
          row = nr;
          col = nc;
          if(nr < nc){
            row = nc;
            col = nr;
            }
          ipt = ij2ipt(3*row,3*col  ,ssize);
          for(j=0; j<3; j++)
            ks[ipt+j] = ks[ipt+j] + kb[j][0];
      
          ipt = ij2ipt(3*row,3*col+1,ssize);
          for(j=0; j<3; j++)    
            ks[ipt+j] = ks[ipt+j] + kb[j][1];
          
          ipt = ij2ipt(3*row,3*col+2,ssize); 
          for(j=0; j<3; j++)    
            ks[ipt+j] = ks[ipt+j] + kb[j][2];       
          }          
        }   
    }
return;
}

// **************************************************************************
void
CHOLESKYLDLT(int ncol, int diags[], double ks[])
// **************************************************************************
// Cholesky LDL(T) decomposition of model stiffness ks where
// ks is stored as a linear array representing lower triangle 
// **************************************************************************
// e.g.
// l0                       d[0][0]
// l1  l6                   l[1][0] d[1][1]
// l2  l7 l11               l[2][0] l[2][1] d[2][2]
// l3  l8 l12 l15           l[3][0] l[3][1] l[3][2] d[3][3]  
// l4  l9 l13 l16 l18       l[4][0] l[4][1] l[4][2] l[4][3] d[4][4]
// l5 l10 l14 l17 l19 l20   l[5][0] l[5][1] l[5][2] l[5][3] l[5][4] d[5][5]
{
int i,j,k,irow,icol,ipt,kpt,nlength;
double temp;

// replace ks[] terms with d[i][i] and l[i][j] coefficients
// column 0
for(i=1; i<ncol; i++)
   ks[i] = ks[i]/ks[0];
   
// loop over columns 1 thru ncol
for(j=1; j<ncol; j++){
  ipt = diags[j];
  // d[j][j] coefficient
  for(k=0; k<(j); k++){
    kpt = diags[k];
    ks[ipt] = ks[ipt] - ks[kpt]*ks[kpt+j-k]*ks[kpt+j-k];
    }
  // l[i][j] coefficients for rows (j+1) to nrow/ncol
  for(i=(j+1); i<(ncol); i++){
    for(k=0; k<(j); k++){
      kpt = diags[k];
      ks[ipt+i-j] = ks[ipt+i-j] - ks[kpt]*ks[kpt+j-k]*ks[kpt+i-k];
      }
    ks[ipt+i-j] = ks[ipt+i-j]/ks[ipt];    
    }
  }
return;
}

// ********************************************************************************
void
CHOLESKYSOLVE(int ncol,int diags[],double ks[],double ds[],double fs[],double y[])
// ********************************************************************************
// Cholesky [L][D][L](T){ds} = {fs} solution
// ks is a linear array representing lower triangle [L] and [D]
// diags[] are pointers to dii terms in ks[]
// First solve [L][D]{y} = {fs} then solve [L](T){ds} = {y}
// ********************************************************************************
{
int i,j,k,l,ii;

// forward solution: [L][D]{y} = {fs}
y[0] = fs[0];
for(i=1; i<ncol; i++){
   y[i] = fs[i];
   for(j=0; j<i; j++){
      k = diags[j];
      y[i] = y[i] - ks[k+i-j]*y[j];
      }
   }
// backward solution: [L](T){ds} = {y}
k = diags[ncol-1];
ds[ncol-1] = y[ncol-1]/ks[k];
i = ncol-1;
for(ii=0; ii<(ncol); ii++){
   k = diags[i];
   ds[i] = y[i]/ks[k];
   l = 1;
   for(j=(i+1); j<ncol; j++){
      ds[i] = ds[i] - ks[k+l]*ds[j];
      l++;
      }
   i = i-1;
   }
return;
}

//**************************************************************
void
lookatd(int nnode, double ks[])
//**************************************************************
// print choleski ldl(t) d matrix
//**************************************************************
{
int i,k,ipt;

ipt = 0;
k   = 3*nnode;
printf(" System Stiffness Matrix Diagonal Terms \n");
for(i=0; i<3*nnode; i++){
   printf(" ks %d d %d %d = %e \n",ipt,i,i,ks[ipt]);
   ipt = ipt + k;
   k = k - 1;
   }
   printf(" \n");
return;
}

//***************************************************************************************
void
modelpic()
//***************************************************************************************
// 20 node element shape functions and element pressure faces
// **************************************************************************************
//         7------ 14 -------6              5----13-----6  7----15-----4   7----14-----6
//        /|                /|              |           |  |           |   |           |
//       / |               / |              |           |  |           |   |           |
//      /  |              /  |             17    +P    18 19    -P    16  15    +Q    13 
//    15   |            13   |  +Q          |           |  |           |   |           |
//    /   19            /   18    |         |           |  |           |   |           |
//   /     |           /     |    |         1-----9-----2  3----11-----0   4----12-----5
//  4--------12-------5      |    |            Face 0          Face 1          Face 2 
//  |      |          |      |    *---> +P  0-----8-----1  4----12-----5   6----14-----7
//  |      3------10--|------2   /          |           |  |           |   |           |
//  |     /           |     /   /           |           |  |           |   |           |
// 16    /           17    /  +R           11    -Q     9 16    +R    17  18    -R    19
//  |  11             |   9                 |           |  |           |   |           |
//  |  /              |  /                  |           |  |           |   |           |
//  | /               | /                   3----10-----2  0-----8-----1   2----10-----3
//  0------- 8 -------1/                       Face 3          Face 4          Face 5
// **************************************************************************************
{
printf("  FEAFXJ20 Development Test Model:                                           \n");
printf("                                                                             \n");
printf("              Y|(Q)                                                          \n");
printf("               |                                                             \n");
printf("               |                                                             \n");
printf("               7----------14-----------6       5----13-----6   7----15-----4 \n");
printf("              /|                      /|       |           |   |           | \n");
printf("             / |                     / |       |           |   |           | \n");
printf("            /  |                    /  |      17    +P    18  19    -P    16 \n");
printf("           /   |                   /   |       |           |   |           | \n");
printf("          /    |                  /    |       |           |   |           | \n");
printf("        15    19                13    18       1-----9-----2   3----11-----0 \n");
printf("        /      |                /      |          Face 0          Face 1     \n");
printf("  Z    /       |               /       |                                     \n");
printf("   |  /        |              /        |       7----14-----6   0-----8-----1 \n");
printf("   | /         3----------10-/---------2----   |           |   |           | \n");
printf("   |/         /             /         /    X   |           |   |           | \n");
printf("   4---------/-12----------5         /    (P) 15    +Q    13  11    -Q     9 \n");
printf("   |        /              |        /          |           |   |           | \n");
printf("   |       /               |       /           |           |   |           | \n");
printf("   |     11                |      9            4----12-----5   3----10-----2 \n");
printf("   |     /                 |     /                Face 2          Face 3     \n");
printf("  16    /                 17    /                                            \n");
printf("   |   /                   |   /               4----12-----5   6----14-----7 \n");
printf("   |  /                    |  /                |           |   |           | \n");
printf("   | /                     | /                 |           |   |           | \n");
printf("   |/                      |/                 16    +R    17  18    -R    19 \n");
printf("   0-----------8-----------1                   |           |   |           | \n");
printf("  /                                            |           |   |           | \n");
printf(" / Z (R)                                       0-----8-----1   2----10-----3 \n");
printf("                                                   Face 4          Face 5    \n");
printf("          (1x1x1 cube)                                                       \n");
printf("                                                                             \n");
return;
}

//****************************************************************************************
void
modelin(int nn[8],double mxyz[],int element[],int emats[],  double matprops[],
                                              int fixnode[],double fixdisp[],
                                              int fpressn[],double pnorm[],
                                              int cforcen[],double cforce[])
//****************************************************************************************
// int     nn[8]       // nnode,etype,nelem,ngauss,nmats,ndisp,npress,nconc
// double *mxyz;       // nnode*(x,y,z)
// int    *element;    // nelem*(n1,n2...n20)
// int    *emats;      // nelem*(ematid)
// double *matprops;   // nmats*(E,nu,alpha)
// int    *fixnode;    // ndisp*4  = ndisp*(node, nfx, nfy, nfz)
// double *fixdisp;    // ndisp*3  = ndisp*(dx, dy, dz)
// int    *fpressn;    // npress*2 = npress*(node, face)
// double *pnorm;      // npress   = npress*(normal pressure)
// double *cforce;     // nconc*3  = nconc*(fx, fy, fz)
// int    *cforcen;    // nconc    = nconc nf
//****************************************************************************************
{
int i,j;
static int    nmodel[ 8] = { 20, 20, 1,  2,  1,  8,  1,  1};
static double coords[60] = {0.,0.,1.,  1.,0.,1.,  1.,0.,0.,  0.,0.,0.,
                            0.,1.,1.,  1.,1.,1.,  1.,1.,0.,  0.,1.,0.,
                            .5,0.,1.,  1.,0.,.5,  .5,0.,0.,  0.,0.,.5,
                            .5,1.,1.,  1.,1.,.5,  .5,1.,0.,  0.,1.,.5,
                            0.,.5,1.,  1.,.5,1.,  1.,.5,0.,  0.,.5,0.};
static int fixn[8] = {2, 3, 6, 7,10,14,18,19};
                                                                                  
// print picture of one element model geometry and nodes used to develop FEAFXJ20
modelpic();

// control information
for(i=0; i<8; i++)
   nn[i] = nmodel[i];
   
// node coordinates
for(i=0; i<60; i++)
   mxyz[i] = coords[i];

// element nodes
for(i=0;i<20; i++)
   element[i] = i;

// material properties
emats[0] = 0;

matprops[0] = 30000000.;
matprops[1] = 0.0;
matprops[2] = 0.0001;

// fixed nodes
for(i=0; i<8; i++){
   fixnode[4*i  ] = fixn[i];
   fixnode[4*i+1] = 1;
   fixnode[4*i+2] = 1;   
   fixnode[4*i+3] = 1;   
   }
for(i=0; i<24; i++)
   fixdisp[i] = 0.0000;

// element face normal pressures
fpressn[0] = 0;
fpressn[1] = 4;

pnorm[0]   = -1500.;

// node point forces
cforcen[1] = 0;

cforce[0]  = 0.;
cforce[1]  = 0.;
cforce[2]  = 0.;

return;
}

//**********************************************************************************
void
readmodelin(int nn[8],    double mxyz[],    int element[],
            int emats[],  double matprops[],
            int fixnode[],double fixdisp[],
            int fpressn[],double pnorm[],
            int cforcen[],double cforce[])
//***********************************************************************************
{
int    nnode, etype, nelem, ngauss, nmats, ndisp, npress, nconc;
int i,j;

FILE * pFile;

// open FEA-FXJ20 model file for reading
pFile = fopen ("feafxj20.txt","r");
if(pFile == NULL ){
  printf(" Error opening FEA model file - terminate \n");
  return;
  }
// read and header data
fscanf(pFile,"%d %d %d %d %d %d %d %d",
       &nnode,&etype,&nelem,&ngauss,&nmats,&ndisp,&npress,&nconc);
       
// material properties
for(i=0; i<nmats; i++)
   fscanf(pFile,"%lE %lE %lE",&matprops[3*i],&matprops[3*i+1],&matprops[3*i+2]);
   
// model node coordinates
for(i=0; i<nnode; i++)
   fscanf(pFile,"%lE %lE %lE",&mxyz[3*i],&mxyz[3*i+1],&mxyz[3*i+2]);
   
// model element indices
for(i=0; i<nelem; i++)
   for(j=0; j<20; j++)
       fscanf(pFile,"%d",&element[20*i+j]);
       
// element material ids
for(i=0; i<nelem; i++)
   fscanf(pFile,"%d",&emats[i]);
   
// fixed nodes
for(i=0; i<ndisp; i++){
   fscanf(pFile,"%d %d %d %d",&fixnode[4*i],&fixnode[4*i+1],&fixnode[4*i+2],&fixnode[4*i+3]);
   }
// element face normal pressures
for(i=0; i<npress; i++){
   fscanf(pFile,"%d %d %lE",&fpressn[2*i],&fpressn[2*i+1],&pnorm[i]);      
   }
// node point forces
for(i=0; i<nconc; i++){
   fscanf(pFile,"%d",&cforcen[i]);  
   for(j=0; j<3; j++)
      fscanf(pFile,"%lE",&cforce[3*i+j]);
   }
// close files
fclose (pFile);
   
return;
}

// ******************************************************************************
void
ESTRAIN(int ngauss, int etype, double eeps[8][6], double eneps[20][6])
// ******************************************************************************
// calculate element node strains from gauss point strains
//           7------ 14 ------6
//          /|               /|
//         / |              / |
//       15  |            13  |    +Q
//       /  19            /  18      |
//      /    |           /    |      |
//     4--------12------5     |      |
//     |     |          |     |      *-----> +P
//     |     3------10--|-----2     /
//     |    /           |    /     /
//    16   /           17   /    +R
//     | 11             |  9
//     | /              | /
//     |/               |/
//     0------- 8 ------1
// ******************************************************************************
{
int j;
double loc,pqr[3];
static double extrap[20][3] = {
             -1.,-1., 1.,   1.,-1., 1.,   1.,-1.,-1.,   -1.,-1.,-1.,
             -1., 1., 1.,   1., 1., 1.,   1., 1.,-1.,   -1., 1.,-1.,
              0.,-1., 1.,   1.,-1., 0.,   0.,-1.,-1.,   -1.,-1., 0.,                             
              0., 1., 1.,   1., 1., 0.,   0., 1.,-1.,   -1., 1., 0.,                                
             -1., 0., 1.,   1., 0., 1.,   1., 0.,-1.,   -1., 0.,-1. };
loc = sqrt(3.);
// loop over element nodes
for(i=0; i<20; i++){
  // get node extrapolation coordinates
  for(j=0; j<3; j++)
     pqr[j] = loc*extrap[i][j];
  // evaluate 8 gauss point shape functions
  E08SHAPES(pqr,sn);
  // calculate node point strains
  for(j=0; j<6; j++){
     eneps[i][j] = 0.;
     for(k=0; k<8; k++)
        eneps[i][j] = eneps[i][j] + sn[k]*eeps[k][j];
     }
   }
return;
}

// ******************************************************************************
void
PSTRESS(double sig[6], double psig[3])
// ******************************************************************************
// calculate principal stress values
// ******************************************************************************
{
int j;
double i1,i2,i3,phi,num,sr,denom;
static double pi = 3.141592654;

i1 =  sig[0]+sig[1]+sig[2];
i2 =  sig[0]*sig[1] + sig[1]*sig[2] + sig[2]*sig[0]
    -(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]);
i3 =  sig[0]*sig[1]*sig[2] + 2.*sig[3]*sig[4]*sig[5]
    -(sig[0]*sig[4]*sig[4] + sig[1]*sig[5]*sig[5] + sig[2]*sig[3]*sig[3]);
    
num   = 2.*i1*i1*i1 - 9.*i1*i2 + 27.*i3;
sr    = sqrt((i1*i1 - 3.*i2));
denom = 2.*sr*sr*sr;
phi   = (acos(num/denom))/3.;

num     = 2.*sqrt((i1*i1-3.*i2))/3.;
psig[0] = (i1/3.) + num*cos(phi);
psig[1] = (i1/3.) + num*cos(phi-(2.*pi/3.));
psig[2] = (i1/3.) + num*cos(phi-(4.*pi/3.));

return;
}

//**************************************************************
int
main()
//**************************************************************
{
double *cforce;       // nconc*3  = nconc*(fx, fy, fz)
int    *cforcen;      // nconc    = nconc nf
int    *diags;        // nnodes*3 = pointers to ks(i,i) diagonal terms
double *ds;           // nnode*3 =  nnode*(dx, dy, dz)
int    *element;      // nelem*(n1,n2...n20)
int    *emats;        // nelem*(ematid)
double *fixdisp;      // ndisp*3  = ndisp*(dx, dy, dz)
int    *fixnode;      // ndisp*4  = ndisp*(node, nfx, nfy, nfz)
int    *fpressn;      // npress*2 = npress*(node, face)
double *fs;           // nnode*3 =  nnode*(fx, fy, fz)
double *ks;           // (3*nnode*(3*nnode+1)/2
double *matprops;     // nmats*(E,nu,alpha)
double *mxyz;         // nnode*(x,y,z)
double *pnorm;        // npress   = npress*(normal pressure)
double *y;            // nnode*3
double *avgeps        // nnode*7 = nnode*(count, ex, ey, ez, exy, eyz, ezx)
int    *matbynode     // nnode*nmats = nnode*(mat0 flag, mat1 flag, ...)

int    nn[8],nnode, etype, nelem, ngauss, nmats, ndisp, npress, nconc;
int    pout[7]; // nodexyzs,enodes,bcforces,bcdisps,ndisps,nstrains,nstresses
int    i,j,k,l,ielm,imat,iface,ipt,npt,ndof,nsize;
int    ir,jc;
int    nf[8];

double emod,pois,alpha,mat[3];
double wt[27],gp[27][3],detj,dv;
double exyz[20][3];
double spqr[3][20],sxyz[3][20];
double jacob[3][3],jinv[3][3];
double d[6][6],ke[60][60],press,ff[24];
double edelta[20][3],eeps[8][6],eneps[20][6],factor;


clock_t start, end;
double cpu_time_used[10];

// ************* time to read in model and initialize matrices *************
start = clock();

// explicitly define test case
nnode  = 20;
etype  = 20;
nelem  =  1;
ngauss =  2;
nmats  =  1;
ndisp  =  8;
npress =  1;
nconc  =  1;

// allocate arrays
matprops = (double *)malloc(sizeof(double)*3*nmats);
mxyz     = (double *)malloc(sizeof(double)*3*nnode);
element  =    (int *)malloc(sizeof(int)*20*nelem);
emats    =    (int *)malloc(sizeof(int)*nelem);
fixnode  =    (int *)malloc(sizeof(int)*4*ndisp);
fixdisp  = (double *)malloc(sizeof(double)*3*ndisp);
fpressn  =    (int *)malloc(sizeof(int)*2*npress);
pnorm    = (double *)malloc(sizeof(double)*npress);
cforcen  =    (int *)malloc(sizeof(int)*nconc);
cforce   = (double *)malloc(sizeof(double)*3*nconc);

// define development model
modelin(nn,mxyz,element,emats,matprops,fixnode,fixdisp,fpressn,pnorm,cforcen,cforce);

// normal model readin from file
// readmodelin();

// model lower triangle stiffness matrix ks initialized to zero
i  = 3*nnode*(3*nnode+1)/2;
ks = (double *)calloc(i,sizeof(double));

// model force and displacement vectors initialized to zero
// {ds} = (nnode*(dx, dy, dz)) and {fs} = (nnode*(fx, fy, fz))
i = 3*nnode;
ds = (double *)calloc(i,sizeof(double));
fs = (double *)calloc(i,sizeof(double));

end =clock();
cpu_time_used[0] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ************* time to echo model input *************
start = clock();
// echo input data
printf(" \n");
printf(" Number of nodes...............%6.0d \n",nnode);
printf(" Element type..................%6.0d \n",etype);
printf(" Number of elements............%6.0d \n",nelem);
printf(" Gauss integration order.......%6.0d \n",ngauss);
printf(" Number of materials...........%6.0d \n",nmats);
printf(" Number of fixed nodes.........%6.0d \n",ndisp);
printf(" Number of pressure types......%6.0d \n",npress);
printf(" Number of node forces.........%6.0d \n",nconc);
printf(" Model stiffness matrix size...%6.0d \n",(3*nnode*(3*nnode+1)/2));
printf(" \n");
printf(" Material Properties: \n");
for(i=0; i<nmats; i++)
  printf(" Material %d Youngs Modulus %E Poissons Ratio %E Alpha %E \n",
           i,matprops[3*i],matprops[3*i+1],matprops[3*i+2]);
printf(" \n");
printf(" Node Coordinates: \n");
for(i=0;i<nnode; i++)
   printf(" node %4.0d %8.3f %8.3f %8.3f \n",i,mxyz[3*i],mxyz[3*i+1],mxyz[3*i+2]);
printf(" \n");
printf(" Element Material ID and Indices: \n");
for(i=0;i<nelem; i++){
   printf(" element %d materialid %d Nodes: ",i,emats[i]);
   for(j=0; j<etype; j++)
       printf("%4.0d ",element[20*i+j]);
   printf("\n");
   }
printf(" \n");
printf(" Fixed Nodes: \n");
for(i=0;i<ndisp; i++){
   printf(" node %4.0d fx %4.0d fy %4.0d fz %4.0d %e %e %e \n",fixnode[4*i],fixnode[4*i+1],
          fixnode[4*i+2],fixnode[4*i+3],fixdisp[3*i],fixdisp[3*i+1],fixdisp[3*i+2]);
   }
printf(" \n");
printf(" Element Face Normal Pressures: \n");
for(i=0;i<npress; i++)
  printf(" element: %5.0d face: %3.0d npress: %e \n",fpressn[2*i],fpressn[2*i+1],pnorm[i]);

printf(" \n");
printf(" Concentrated Node Forces: \n");
for(i=0;i<nconc; i++){
  printf(" node: %5.0d ",cforcen[i]);
  printf(" Fx: %8.3f Fy: %8.3f Fz: %8.3f \n",cforce[3*i],cforce[3*i+1],cforce[3*i+2]);
  }

// size of system stiffness matrix
nsize = 3*nnode;

// get element gauss points
GAUSSQUAD(ngauss,wt,gp);
// set number of integration points
npt = 8;
if(ngauss == 3)npt = 27;

end =clock();
cpu_time_used[1] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ************* time to compute and assemble element matrices *************
start = clock();
printf(" Computing and Assembling Element Stiffness Matrices \n");

// loop over elements
for(ielm=0; ielm<nelem; ielm++){
   // get element node coordinates
      for(i=0; i<20; i++){
        j = element[20*ielm+i];
        for(k=0; k<3; k++)
           exyz[i][k] = mxyz[3*j+k];
        }
   // element material properties & elasticity matrix constants and element [D] matrix
      i = emats[ielm];
      emod = matprops[3*i];
      pois = matprops[3*i+1];
      mat[0] = emod*(1.-pois)/((1.+pois)*(1.-2.*pois));
      mat[1] = emod*pois/((1.+pois)*(1.-2.*pois));
      mat[2] = emod/(2.*(1.+pois)); 
   // initialize element stiffness [ke] to zero
      ndof = 3*etype;
      for(i=0; i<ndof; i++)
        for(j=0; j<ndof; j++)
          ke[i][j] = 0.;
   // loop over gauss points (ipt)
   for(ipt=0; ipt < npt; ipt++){
      // evaluate partial derivatives in curvilinear coordinates spqr[3][20]
         E20PARTS(ipt,gp,spqr);
      // calculate jacobian matrix
         MMULTAB(3,3,20,jacob,spqr,exyz);
      // jacobian matrix inverse and determinant jinv[3][3]
         detj = MINV3X3(ielm,ipt,jacob,jinv);     
      // partial derivatives in global coordinates sxyz[3][20] = jinv[3][3] spqr[3][20]
         MMULTAB(3,20,3,sxyz,jinv,spqr);
      // calculate [ke] = [ke] + [B]T[D][B]dv
         dv = wt[ipt]*detj;               
      // compute lower triangle element stiffness matrix  
         BTDB(etype,ke,mat,sxyz,dv);
      }
  // assemble element ke into model lower triangle model stiffness ks
     ASSEMBLEK(ielm,element,etype,ks,ke,nsize);
  }

end =clock();
cpu_time_used[2] = ((double) (end - start)) / CLOCKS_PER_SEC;    

// ************* time to apply fixed displacements and forces *************
start = clock();
// Normal Pressure Loads
for(i=0; i<npress; i++){
   ielm  = fpressn[2*i];
   iface = fpressn[2*i+1];
   press = pnorm[i];
   // calculate equivalent node point force vector
   E20PRESS(ielm,iface,press,element,mxyz,nf,ff);
   // add element node forces into global force vector
   for(j=0; j<8; j++){
      k = nf[j];
      for(l=0; l<3; l++)
         fs[3*k+l] = fs[3*k+l] + ff[3*j+l];
      }
   }
// Concentrated Forces
for(i=0; i<nconc; i++){
   k = cforcen[i];
   for(j=0; j<3; j++)
      fs[3*k+j] = fs[3*k+j] + cforce[3*i+j];
   }
   
// Print Total Node Forces
printf("\n Node Forces \n");
for(i=0; i<nnode; i++)
   printf(" Node %5.0d dx %14.6e dy %14.6e dz %14.6e \n",i,fs[3*i],fs[3*i+1],fs[3*i+2]);

// Prescribed Displacements
for(i=0; i<ndisp; i++){
   k = fixnode[4*i];
   for(j=0; j<3; j++){
     if(fixnode[4*i+j+1] == 1){
        // apply fixed displacement to dof
        ipt = ij2ipt((3*k+j),(3*k+j),nsize);
        ks[ipt] = 1.E+20*ks[ipt];
        fs[3*k+j] = ks[ipt]*fixdisp[3*i+j];
        }
     }  
   }
// define an array of pointers to diagonal elements in array ks
diags = (int *)malloc(sizeof(int)*nsize);
diags[0] = 0;
for(i=1;i<nsize;i++)
   diags[i] = diags[i-1] +(nsize+1-i);
  
end =clock();
cpu_time_used[3] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ******** time to decompose assembled system matrix ks ********
start = clock();

// Cholesky [L](T)[D][L]
printf("\n Cholesky L D L(T) Decomposition \n");
CHOLESKYLDLT(nsize,diags,ks);

end =clock();
cpu_time_used[4] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ******** time to solve assembled system equations [L][D][[L](T){ds} = {fs} ********
start = clock();

// Cholesky Forward [L]{y} = {fs} and Backward [D][L](T){ds} = {fs}
printf("\n Cholesky [L][D][L](T){ds} = {fs} Solver \n");
y = (double *)malloc(sizeof(double)*nsize);
CHOLESKYSOLVE(nsize,diags,ks,ds,fs,y);

// free stiffness matrices and vectors
free(diags);
free(ks);
free(y);

end =clock();
cpu_time_used[5] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ******** time to compute node point strains and stresses ********
start = clock();

// allocate and initialize node point average strains array
i = 7*nnode;
avgeps = (double *)calloc(i,sizeof(double));

// loop over elements
for(ielm=0; ielm<nelem; ielm++){
   // get element node coordinates and displacements
      for(i=0; i<20; i++){
        j = element[20*ielm+i];
        for(k=0; k<3; k++){
           exyz[i][k]   =  mxyz[3*j+k];
           edelta[i][k] =    ds[3*j+k];
           }
        }
   // element material properties & elasticity matrix constants and element [D] matrix
      i = emats[ielm];
      emod = matprops[3*i];
      pois = matprops[3*i+1];
      mat[0] = emod*(1.-pois)/((1.+pois)*(1.-2.*pois));
      mat[1] = emod*pois/((1.+pois)*(1.-2.*pois));
      mat[2] = emod/(2.*(1.+pois)); 
   // initialize element strains at gauss points eeps[8][6] to zero
      for(i=0; i<npt; i++)
        for(j=0; j<6; j++)
          eeps[i][j] = 0.;
   // loop over gauss points (ipt)
   for(ipt=0; ipt < npt; ipt++){
      // evaluate partial derivatives in curvilinear coordinates spqr[3][20]
      E20PARTS(ipt,gp,spqr);
      // calculate jacobian matrix
      MMULTAB(3,3,20,jacob,spqr,exyz);
      // jacobian matrix inverse and determinant jinv[3][3]
      detj = MINV3X3(ielm,ipt,jacob,jinv);     
      // partial derivatives in global coordinates sxyz[3][20] = jinv[3][3] spqr[3][20]
      MMULTAB(3,20,3,sxyz,jinv,spqr);              
      // compute element gauss point strains 
      for(j=0; j<etype; j++){
        eeps[ipt][0] = eeps[ipt][0]+    sxyz[0][j]*edelta[j][0];
        eeps[ipt][1] = eeps[ipt][1]+    sxyz[1][j]*edelta[j][1];
        eeps[ipt][2] = eeps[ipt][2]+    sxyz[2][j]*edelta[j][2];
        eeps[ipt][3] = eeps[ipt][3]+.5*(sxyz[1][j]*edelta[j][0]+sxyz[0][j]*edelta[j][1]);
        eeps[ipt][4] = eeps[ipt][4]+.5*(sxyz[2][j]*edelta[j][1]+sxyz[1][j]*edelta[j][2]);
        eeps[ipt][5] = eeps[ipt][5]+.5*(sxyz[2][j]*edelta[j][0]+sxyz[0][j]*edelta[j][2]);
        }       
      }
  // extrapolate element gauss point strains to element node points
     ESTRAINS(ngauss,etype,eeps,eneps);  
     
  // add element node point strains into avgeps[7*nnode] array
     for(i=0; i<20; i++){
       j = element[20*ielm+i];
       avgeps[7*j] = avgeps[7*j]+1.;
       for(k=0; k<6; k++)
          avgeps[7*j+k+1] = avgeps[7*j+k+1] + eneps[i][k];
       }
  }
// calculate average strains
for(i=0; i<nnode; i++){
   if(avgeps[7*i] > 0.){
      factor = 1./avgeps[7*i];
      for(k=1; k<7; k++)avgeps[7*i+k] = factor*avgeps[7*i+k];
      }
   }

end =clock();
cpu_time_used[6] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ******** time to compute node point stresses ********
start = clock();

// allocate and initialize node by material array matbynode[nnode*nmat]
i = nnode*nmats;
matbynode = (int *)calloc(i,sizeof(int));

if(nmats == 1){
   for(i=0; i<nnode; i++)
      matbynode[i] = 1;
   }
else{
   for(ielm=0; ielm<nelem; ielm++){
       imat = emats[ielm];  
       for(i=0; i<20; i++){
          j = element[20*ielm+i];
          matbynode[j] = 1;
          }
       }
   } 

end =clock();
cpu_time_used[7] = ((double) (end - start)) / CLOCKS_PER_SEC;

// ******** time to output results ********
start = clock();

// Print Node Displacements
printf("\n Node Displacements \n");
for(i=0; i<nnode; i++)
   printf(" Node %5.0d dx %14.6e dy %14.6e dz %14.6e \n",i,ds[3*i],ds[3*i+1],ds[3*i+2]);

// Print Node Point Strains


// Print Node Point Stresses


end =clock();
cpu_time_used[8] = ((double) (end - start)) / CLOCKS_PER_SEC;

// print timing results
printf(" \n");
printf(" Timing Results: \n");
printf(" Time to read in model............................%12.7f \n",cpu_time_used[0]);
printf(" Time to echo model input.........................%12.7f \n",cpu_time_used[1]);
printf(" Time to compute and assemble element stiffness...%12.7f \n",cpu_time_used[2]);
printf(" Time to apply fixed displacements and forces.....%12.7f \n",cpu_time_used[3]);
printf(" Time to decompose assembled stiffness matrix.....%12.7f \n",cpu_time_used[4]);
printf(" Time to solve assembled system equations.........%12.7f \n",cpu_time_used[5]);
printf(" Time to compute node point strains...............%12.7f \n",cpu_time_used[6]);
printf(" Time to compute node point stresses..............%12.7f \n",cpu_time_used[7]);
printf(" Time to printout results.........................%12.7f \n",cpu_time_used[8]);

return 0;
}
