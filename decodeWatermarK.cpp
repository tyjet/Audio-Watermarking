#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define L (32)
#define LENGTH (1024)
#define NUM_SEGMENTS (32)
#define NUM_SEGMENTS_OVER_TWO (16)
#define ONE_OVER_NUM_SEGMENTS (0.03125)
#define NUM_COEFFICIENTS (16)
#define INT_HI_MASK (4294901760)
#define INT_LO_MASK (65535)
#define L_OVER_TWO (16)
#define ONE_OVER_L (0.0625)
#define ONE_OVER_SQRT_L (0.25)
#define SQRT_TWO (1.41421356237)
#define R (4)
#define ONE_OVER_R (0.25)
#define M (2) /* M must be a power of 2 */
#define TWO_M (4)
#define ONE_OVER_M (0.5)
#define ALPHA_ONE (1.41421356237)
#define THRESHOLD (6) /* TODO: figure out a good value empirically */
#define ALPHA_TWO (0.70710678118)
#define PI (3.14159265358979323846264338327950288419716939937510)

/* Function for embedding a watermark within a given bitstream using the
 * pathwork algorithm decribed in [Natgunanathan 2012]. The paper can be found
 * at http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6198872 *
 *
 * It is only possible to encode NUM_SEGMENTS bits per stream.
 */
void decodeWatermarK(int** inputStream, int* secretKey, char* &outputStream)
{
  /*
   ****************************************************************************
   * Step 1: DCT coefficient and fragment generation 
   ****************************************************************************
  */

  double coefficientsF[NUM_SEGMENTS][L];
  double coefficientsR[NUM_SEGMENTS][L];
  double dctMultiplier;
  int inputStreamF, inputStreamR;
  double coefficientF, coefficientR;
  /* Generate the DCT coefficients for each segment */
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    coefficientF = 0.0;
    coefficientR = 0.0;
    for (int coefficientIdx = 0; coefficientIdx < L_OVER_TWO; ++coefficientIdx)
    {
      for (int sampleIdx = 0; sampleIdx < L_OVER_TWO; ++sampleIdx)
      {            
        inputStreamF = inputStream[segmentIdx][sampleIdx];
        inputStreamR = inputStream[segmentIdx][sampleIdx + L_OVER_TWO];
        
        dctMultiplier = (coefficientIdx) ? 2 * ONE_OVER_SQRT_L : SQRT_TWO * ONE_OVER_SQRT_L;
        coefficientF += dctMultiplier * inputStreamF * cos(PI * (2 * sampleIdx + 1) * coefficientIdx) * ONE_OVER_L;
        coefficientR += dctMultiplier * inputStreamR * cos(PI * (2 * sampleIdx + 1) * coefficientIdx) * ONE_OVER_L;
      } 
      
      coefficientsF[segmentIdx][coefficientIdx] = coefficientF;
      coefficientsR[segmentIdx][coefficientIdx] = coefficientR;
    }
  }
  
  /* Split the DCT coefficients into R frames of length 2M */
  double framesF[NUM_SEGMENTS][R][TWO_M];
  double framesR[NUM_SEGMENTS][R][TWO_M];
  int frameIdx = 0;
  int kIdx = 0;
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    for (int coefficientIdx = 0; coefficientIdx < NUM_COEFFICIENTS; ++coefficientIdx)
    {
      framesF[segmentIdx][frameIdx][kIdx] = coefficientsF[segmentIdx][coefficientIdx];
      framesR[segmentIdx][frameIdx][kIdx] = coefficientsR[segmentIdx][coefficientIdx];
      kIdx++;
      if (kIdx == TWO_M)
      {
        kIdx = 0;
        frameIdx++;        
      }      
    }
  }
  
  /* Split each frame into M fragments using the PN sequence in the secret key*/
  double fragmentsF1[NUM_SEGMENTS][R][M];
  double fragmentsF2[NUM_SEGMENTS][R][M];
  double fragmentsR1[NUM_SEGMENTS][R][M];
  double fragmentsR2[NUM_SEGMENTS][R][M];

  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    for (int frameIdx = 0; frameIdx < R; ++frameIdx)
    {
      for (int fragmentIdx = 0; fragmentIdx < M; ++fragmentIdx)
      {
        /* TODO: verify this is the correct fragment implementation */
        fragmentsF1[segmentIdx][frameIdx][fragmentIdx] = framesF[segmentIdx][frameIdx][secretKey[fragmentIdx]];
        fragmentsF2[segmentIdx][frameIdx][fragmentIdx] = framesF[segmentIdx][frameIdx][secretKey[M + fragmentIdx]];
        fragmentsR1[segmentIdx][frameIdx][fragmentIdx] = framesR[segmentIdx][frameIdx][secretKey[fragmentIdx]];
        fragmentsR2[segmentIdx][frameIdx][fragmentIdx] = framesR[segmentIdx][frameIdx][secretKey[M + fragmentIdx]];        
      }
    }
  }
  
  /*
   ****************************************************************************
   * Step 2 : Determine which DCT frame pairs contain a watermarK 
   ****************************************************************************
  */
   
  /* Compute the expected value for frames across all segments */
  double p[NUM_SEGMENTS][R];
  double p1[NUM_SEGMENTS][R];
  double p2[NUM_SEGMENTS][R];
  double pTilde[NUM_SEGMENTS][R];
  double q[NUM_SEGMENTS][R];
  double q1[NUM_SEGMENTS][R];
  double q2[NUM_SEGMENTS][R];
  double qTilde[NUM_SEGMENTS][R];
  char   containsWatermarK[NUM_SEGMENTS][R];
   
  double pVal, p1Val, p2Val, pTildeVal, qVal, q1Val, q2Val, qTildeVal;
  char containsWatermarKVal, isSilent;
  
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  { 
    pVal  = 0.0;
    p1Val = 0.0;
    p2Val = 0.0;
    qVal  = 0.0;
    q1Val = 0.0;
    q2Val = 0.0;
    
    for (int frameIdx = 0; frameIdx < R; ++frameIdx)
    {      
      /* Compute the expected value for fragments across all segments */
      for (int fragmentIdx = 0; fragmentIdx < M; ++fragmentIdx)
      {
        p1Val += abs(fragmentsF1[segmentIdx][frameIdx][fragmentIdx]);
        p2Val += abs(fragmentsF2[segmentIdx][frameIdx][fragmentIdx]);
        q1Val += abs(fragmentsR1[segmentIdx][frameIdx][fragmentIdx]);
        q2Val += abs(fragmentsR2[segmentIdx][frameIdx][fragmentIdx]);
      }
      
      p1Val *= ONE_OVER_M;
      p2Val *= ONE_OVER_M;
      q1Val *= ONE_OVER_NUM_SEGMENTS;
      q2Val *= ONE_OVER_NUM_SEGMENTS;
      
      /* Compute the expected value of this frame */
      pVal = 0.5 * (p1Val + p2Val);  
      qVal = 0.5 * (q1Val + q2Val);  
    
      /* Test statistic used to determine if this frame should be used for
         watermarKing. */
      pTildeVal   = abs(p1Val - p2Val) - ALPHA_ONE * pVal; 
      qTildeVal   = abs(q1Val - q2Val) - ALPHA_ONE * qVal;
    
      /* Boolean indicating whether of not this frame should be used for
         watermarKing. */
      isSilent    = (pVal - qVal >= 0) ? qVal < THRESHOLD : pVal < THRESHOLD; /* THRESHOLD will be determined emperically. 
                                                                                isSilent is used to maKe sure we do not
                                                                                insert sound during a silent portion of
                                                                                the host audio segment
                                                                              */
      containsWatermarKVal = pTildeVal <= 0.0 && qTildeVal <= 0.0 && !isSilent;
      
      /* Store the values so they can be accessed later */
      p[segmentIdx][frameIdx]        = pVal;
      p1[segmentIdx][frameIdx]       = p1Val;
      p2[segmentIdx][frameIdx]       = p2Val;
      pTilde[segmentIdx][frameIdx]   = pTildeVal;
      q[segmentIdx][frameIdx]        = qVal;
      q1[segmentIdx][frameIdx]       = q1Val;
      q2[segmentIdx][frameIdx]       = q2Val;
      qTilde[segmentIdx][frameIdx]   = qTildeVal;
      containsWatermarK[segmentIdx][frameIdx] = containsWatermarKVal;
    }                    
  }
  
  /*
   ****************************************************************************
   * Step 3 : For every frame pair that contains a watermarK, extract the 
   *          watermarK. 
   ****************************************************************************
  */
  
  double r1, r2;
  int numOnes = 0;
  int numZeros = 0;
  /* Calculate the modified expected value for each segment */
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    numOnes = 0;
    numZeros = 0;
    /* Each watermarK bit can be encoded in multiple frames, so we encode it in
     * every valid frame pair
     */
    for (int frameIdx = 0; frameIdx < R; ++frameIdx)
    {
      if (containsWatermarK[segmentIdx][frameIdx])
      {
        /* TODO: Verify that equation 28 in Natgunanathan, 2012 is entirely prime'd or not prime'd */
        r1 = p1[segmentIdx][frameIdx] - p[segmentIdx][frameIdx] + q[segmentIdx][frameIdx] - q1[segmentIdx][frameIdx];
        r2 = p2[segmentIdx][frameIdx] - p[segmentIdx][frameIdx] + q[segmentIdx][frameIdx] - q2[segmentIdx][frameIdx];
        
        if (r1 > 0.0 && r2 < 0.0) numZeros++;
        else                      numOnes++;
      }
    }
    if (numZeros > numOnes) 
    {
      outputStream[segmentIdx] = 0;
    }
    else if (numOnes > numZeros)
    {
	  outputStream[segmentIdx] = 1;
    }
    else
    {
      outputStream = NULL;
      printf("You got haxxed son!\n");
      return;
    }
  }
}