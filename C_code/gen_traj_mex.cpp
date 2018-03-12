#include <iostream>
#include <fstream>
#include <cmath>
#include "mex.h"
using namespace std;

#define M_PI 3.14159265358979323846

double haltonnumber(int index,int base)
{
  double result = 0;
  double f = 1.0;
  int i = index;
  while (i>0)
  {
    f = f/base;
    result += f*fmod(i,base);
    i = i/base;
  }
  return result;
}

void haltonSeq(double *m_adAzimuthalAngle,double *m_adPolarAngle,long num_frames,long num_projPerFrame)
{
  int p1 = 2;
  int p2 = 3;
  double z;
  double phi;
  int linter;
  int lk;                  // index for projections in a frame
  int lFrame;              // index for number of frames

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      z = haltonnumber(lk+1,p1) * 2 - 1;
      phi = 2 * M_PI * haltonnumber(lk+1,p2);
      //calculate polar and azimuthal angle on a unit sphere
      m_adPolarAngle[linter] = acos(z);
      m_adAzimuthalAngle[linter] = phi;
    }
  }
}

void spiralSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   llin;                // linearly increasing index
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double dPreviousAngle;      // previous angle value
  double dH;                  // parameter for the calculation

  long num_totalProjections = num_frames * num_projPerFrame;

  dPreviousAngle = 0;

  for (lk = 0; lk < num_projPerFrame; lk ++ )
  {
    for (lFrame = 0; lFrame < num_frames; lFrame ++ )
    {
      llin    = lFrame + lk * num_frames;
      linter  = lk + lFrame * num_projPerFrame;

      dH = -1 + 2 * llin / double ( num_totalProjections );
      m_polar[linter] = acos (dH);
      if (llin == 0)
      m_azi[linter] = 0;
      else
      m_azi[linter] = fmod ( dPreviousAngle + 3.6/ ( sqrt( num_totalProjections * (1. - dH*dH) ) ) , 2.0 * M_PI );

      dPreviousAngle = m_azi[linter];

    }//endfor lk
  }//endfor lInt
}

void archimedianSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double dZ;
  double dAngle;

  dAngle = (3.0-sqrt(5.0))*M_PI;
  dZ = 2.0/(num_projPerFrame - 1.0);

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      m_polar[linter] = acos(1.0 - dZ*lk);
      m_azi[linter] = lk*dAngle;
    }
  }
}

void dgoldenMSeq(double *m_azi,double *m_polar,long num_frames,long num_projPerFrame)
{
  long   lk;
  long   linter;              // interleaved increasing index
  long   lFrame;              // index for number of frames
  double goldmean1;
  double goldmean2;
  goldmean1 = 0.465571231876768;
  goldmean2 = 0.682327803828019;

  for (lFrame = 0; lFrame < num_frames; lFrame ++ )
  {
    for (lk = 0; lk < num_projPerFrame; lk ++ )
    {
      linter  = lk + lFrame * num_projPerFrame;
      m_polar[linter] = acos(2.0 * fmod(lk * goldmean1, 1) - 1);
      m_azi[linter] = 2 * M_PI * fmod(lk * goldmean2, 1);
    }
  }
}

// Functions for Random Spiral
void swap (double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

long partition(double *ht_polar, double *sp_polar, double *sp_azi, long low, long high)
{
    double pivot = ht_polar[high];
    long i = low - 1;

    for(long j = low; j<=high-1;j++)
    {
        if(ht_polar[j] <= pivot)
        {
            i++;
            swap(&ht_polar[i],&ht_polar[j]);
            swap(&sp_polar[i],&sp_polar[j]);
            swap(&sp_azi[i],&sp_azi[j]);
        }
    }
    swap(&ht_polar[i+1],&ht_polar[high]);
    swap(&sp_polar[i+1],&sp_polar[high]);
    swap(&sp_azi[i+1],&sp_azi[high]);

    return (i+1);
}

void quickSort(double *ht_polar, double *sp_polar, double *sp_azi, long low, long high)
{
    if(low < high)
    {
        long pi = partition(ht_polar,sp_polar,sp_azi,low,high);

        quickSort(ht_polar,sp_polar,sp_azi,low,pi-1);
        quickSort(ht_polar,sp_polar,sp_azi,pi+1,high);
    }
}

// Ziyi function randomizes spiral traj with Halton sequence
void randomSpiral(double *m_adAzimuthalAngle,double *m_adPolarAngle,long num_projPerFrame)
{
  double *ht_adAzimu = new double [num_projPerFrame];
  double *ht_adPolar = new double [num_projPerFrame];
  haltonSeq(ht_adAzimu,ht_adPolar,1,num_projPerFrame);

  quickSort(ht_adPolar,m_adPolarAngle,m_adAzimuthalAngle,0,num_projPerFrame-1);

  delete [] ht_adPolar;
  delete [] ht_adAzimu;
}

// calculate Trajectory
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /********* input check ***********/
  /* check for proper number of arguments */
  if(nrhs!=2) {
      mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
  }
  if(nlhs!=3) {
      mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Three output required.");
  }
  /* make sure the first input argument is scalar */
  if( !mxIsDouble(prhs[0]) ||
       mxIsComplex(prhs[0]) ||
       mxGetNumberOfElements(prhs[0])!=1 ) {
      mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","First Input must be a scalar.");
  }

  /* make sure the second input argument is scalar */
  if( !mxIsDouble(prhs[1]) ||
       mxIsComplex(prhs[1]) ||
       mxGetNumberOfElements(prhs[1])!=1 ) {
      mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Second Input must be a scalar.");
  }
  /********* finished check ***********/

  /********  acquire input and output ************/
  /* get the value of the scalar input  */
  long m_lTrajectoryType = mxGetScalar(prhs[0]);

  /* create a pointer to the real data in the input matrix  */
  long m_lProjectionsPerFrame = mxGetScalar(prhs[1]);

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(1,(mwSize)m_lProjectionsPerFrame,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,(mwSize)m_lProjectionsPerFrame,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,(mwSize)m_lProjectionsPerFrame,mxREAL);

  /* get a pointer to the real data in the output matrix */
  double * cordinateX = mxGetPr(plhs[0]);
  double * cordinateY = mxGetPr(plhs[1]);
  double * cordinateZ = mxGetPr(plhs[2]);
  /********  finished input and output ************/

  double * m_adAzimuthalAngle = new double [m_lProjectionsPerFrame];
  double * m_adPolarAngle     = new double [m_lProjectionsPerFrame];

  long m_lNumberOfFrames = 1;

  switch(m_lTrajectoryType)
  {
    //Ziyi: siemens spiral
    case 1: {spiralSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi, case 2: archimedian Seq
    case 2: {haltonSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi: case 3 haltoned spiral
    case 3:
    {
      spiralSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);
      randomSpiral(m_adAzimuthalAngle,m_adPolarAngle,m_lProjectionsPerFrame);
      break;
    }

    //Ziyi, case 4: double golden means
    case 4: {archimedianSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}

    //Ziyi, case 5: double golden means
    default: {dgoldenMSeq(m_adAzimuthalAngle,m_adPolarAngle,m_lNumberOfFrames,m_lProjectionsPerFrame);break;}
  }

  //convert angles to x, y ,z coordinate
  for(int k = 0; k<m_lProjectionsPerFrame; k++)
  {
    cordinateX[k] = sin(m_adPolarAngle[k])*cos(m_adAzimuthalAngle[k]);
    cordinateY[k] = sin(m_adPolarAngle[k])*sin(m_adAzimuthalAngle[k]);
    cordinateZ[k] = cos(m_adPolarAngle[k]);
  }
}
