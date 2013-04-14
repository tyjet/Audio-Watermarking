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
void embedWatermarK(int** inputStream, char* watermarK, int watermarKLength, int** &outputStream, int* &secretKey)
{
  /*
   ****************************************************************************
   * Step 1: DCT coefficient and fragment generation 
   ****************************************************************************
  */

  double coefficientsF[NUM_SEGMENTS][L_OVER_TWO];
  double coefficientsR[NUM_SEGMENTS][L_OVER_TWO];
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
  
  /* Split each frame into M fragments using a PN sequence */
  int pnSequence[TWO_M];
  for (int idx = 0; idx < TWO_M; ++idx)
  {
    pnSequence[idx] = idx;
  }
  
  /* Generate random indices for the pnSequence */
  int tmp;
  int rnd;
  for (int idx = 0; idx < TWO_M; ++idx)
  {
    rnd = rand() % TWO_M;
    tmp = pnSequence[rnd];
    pnSequence[rnd] = pnSequence[idx];
    pnSequence[idx] = tmp;
  }
  
  /* Assign the PN sequence to the secret key */
  secretKey = pnSequence;
  
  /* Generate the fragments for each frame */
  double fragmentsF1[NUM_SEGMENTS][R][M];
  double fragmentsF2[NUM_SEGMENTS][R][M];
  double fragmentsR1[NUM_SEGMENTS][R][M];
  double fragmentsR2[NUM_SEGMENTS][R][M];

  /* This array is not needed in the updated code, but I feel liKe it should be Kept around */
  /*
  double k[R] // TODO: figure out how to reference k
  for (int idx = 0; idx < R; ++idx)
  {
    k[idx] = 0;
  }
  */
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    for (int frameIdx = 0; frameIdx < R; ++frameIdx)
    {
      for (int fragmentIdx = 0; fragmentIdx < M; ++fragmentIdx)
      {
        /* This method doesn't worK in the updated code, but I feel liKe it should be Kept around */
        /*
        fragmentsF1[segmentIdx][frameIdx][fragmentIdx] = framesF[k[frameIdx] + pnSequence[fragmentIdx]][frameIdx];
        fragmentsF2[segmentIdx][frameIdx][fragmentIdx] = framesF[k[frameIdx] + pnSequence[M + fragmentIdx]][frameIdx];
        fragmentsR1[segmentIdx][frameIdx][fragmentIdx] = framesR[k[frameIdx] + pnSequence[fragmentIdx]][frameIdx];
        fragmentsR2[segmentIdx][frameIdx][fragmentIdx] = framesR[k[frameIdx] + pnSequence[M + fragmentIdx]][frameIdx];
        */
        
        /* TODO: verify this is the correct fragment implementation */
        fragmentsF1[segmentIdx][frameIdx][fragmentIdx] = framesF[segmentIdx][frameIdx][pnSequence[fragmentIdx]];
        fragmentsF2[segmentIdx][frameIdx][fragmentIdx] = framesF[segmentIdx][frameIdx][pnSequence[M + fragmentIdx]];
        fragmentsR1[segmentIdx][frameIdx][fragmentIdx] = framesR[segmentIdx][frameIdx][pnSequence[fragmentIdx]];
        fragmentsR2[segmentIdx][frameIdx][fragmentIdx] = framesR[segmentIdx][frameIdx][pnSequence[M + fragmentIdx]];        
      }
    }
  }
  
  /*
   ****************************************************************************
   * Step 2 : Select DCT frame pairs for watermarking 
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
  char   useFrame[NUM_SEGMENTS][R];
   
  double pVal, p1Val, p2Val, pTildeVal, qVal, q1Val, q2Val, qTildeVal;
  char useFrameVal, isSilent;
  
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
      
      p1Val      *= ONE_OVER_M;
      p2Val      *= ONE_OVER_M;
      q1Val      *= ONE_OVER_NUM_SEGMENTS;
      q2Val      *= ONE_OVER_NUM_SEGMENTS;
      
      /* Compute the expected value of this frame */
      pVal        = 0.5 * (p1Val + p2Val);  
      qVal        = 0.5 * (q1Val + q2Val);  
    
      /* Test statistic used to determine if this frame should be used for
         watermarKing. */
      pTildeVal   = abs(p1Val - p2Val) - ALPHA_ONE * pVal; 
      qTildeVal   = abs(q1Val - q2Val) - ALPHA_ONE * qVal;
    
      /* Boolean indicating whether of not this frame should be used for
         watermarKing. */
      isSilent    = (pVal - qVal >=0) ? qVal < THRESHOLD : pVal < THRESHOLD; /* THRESHOLD will be determined emperically. 
                                                                                isSilent is used to maKe sure we do not
                                                                                insert sound during a silent portion of
                                                                                the host audio segment
                                                                              */
      useFrameVal = pTildeVal <= 0.0 && qTildeVal <= 0.0 && !isSilent;
      
      /* Store the values so they can be accessed later */
      p[segmentIdx][frameIdx]        = pVal;
      p1[segmentIdx][frameIdx]       = p1Val;
      p2[segmentIdx][frameIdx]       = p2Val;
      pTilde[segmentIdx][frameIdx]   = pTildeVal;
      q[segmentIdx][frameIdx]        = qVal;
      q1[segmentIdx][frameIdx]       = q1Val;
      q2[segmentIdx][frameIdx]       = q2Val;
      qTilde[segmentIdx][frameIdx]   = qTildeVal;
      useFrame[segmentIdx][frameIdx] = useFrameVal;
    }                    
  }
  
  /*
   ****************************************************************************
   * Step 3 : For every selected frame pair, embed a bit of the watermarK 
   ****************************************************************************
  */
  
  int curWatermarKBit = 0;
  double pPrime, p1Prime, p2Prime, qPrime, q1Prime, q2Prime;
  /* Calculate the modified expected value for each segment */
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS && curWatermarKBit < watermarKLength; ++segmentIdx)
  {
  
    int sampleIdx = 0;
    /* Each watermarK bit can be encoded in multiple frames, so we encode it in
     * every valid frame pair
     */
    for (int frameIdx = 0; frameIdx < R; ++frameIdx)
    {
      if (useFrame[segmentIdx][frameIdx])
      {
        /* Encode a '0' bit */
        if (!watermarK[curWatermarKBit])
        {
          /* Check what p1' and p2' should be */
          if (p1[segmentIdx][frameIdx] - p2[segmentIdx][frameIdx] >= ALPHA_TWO * p[segmentIdx][frameIdx])
          {
            p1Prime = p1[segmentIdx][frameIdx];
            p2Prime = p2[segmentIdx][frameIdx];
          }
          else
          {
            p1Prime = (1.0 + 0.5 * ALPHA_TWO) * p[segmentIdx][frameIdx];
            p2Prime = (1.0 - 0.5 * ALPHA_TWO) * p[segmentIdx][frameIdx];
          }
          
          /* Check what q1' and q2' should be */
          if (q2[segmentIdx][frameIdx] - q1[segmentIdx][frameIdx] >= ALPHA_TWO * q[segmentIdx][frameIdx])
          {
            q1Prime = q1[segmentIdx][frameIdx];
            q2Prime = q2[segmentIdx][frameIdx];
          }
          else
          {
            q1Prime = (1.0 - 0.5 * ALPHA_TWO) * q[segmentIdx][frameIdx];
            q2Prime = (1.0 + 0.5 * ALPHA_TWO) * q[segmentIdx][frameIdx];
          }
        }
        
        /* Encode a '1' bit */
        else if (watermarK[curWatermarKBit])
        {
          /* Check what p1' and p2' should be */
          if (p2[segmentIdx][frameIdx] - p1[segmentIdx][frameIdx] >= ALPHA_TWO * p[segmentIdx][frameIdx])
          {
            p1Prime = p1[segmentIdx][frameIdx];
            p2Prime = p2[segmentIdx][frameIdx];
          }
          else
          {
            p1Prime = (1.0 - 0.5 * ALPHA_TWO) * p[segmentIdx][frameIdx];
            p2Prime = (1.0 + 0.5 * ALPHA_TWO) * p[segmentIdx][frameIdx];
          }
          
          /* Check what q1' and q2' should be */
          if (q1[segmentIdx][frameIdx] - q2[segmentIdx][frameIdx] >= ALPHA_TWO * q[segmentIdx][frameIdx])
          {
            q1Prime = q1[segmentIdx][frameIdx];
            q2Prime = q2[segmentIdx][frameIdx];
          }
          else
          {
            q1Prime = (1.0 + 0.5 * ALPHA_TWO) * q[segmentIdx][frameIdx];
            q2Prime = (1.0 - 0.5 * ALPHA_TWO) * q[segmentIdx][frameIdx];
          }
        }
        
        pPrime = 0.5 * (p1Prime + p2Prime);
        qPrime = 0.5 * (q1Prime + q2Prime);
        
        for (int coefficientIdx = 0; coefficientIdx < TWO_M; ++coefficientIdx)
        {
          coefficientsF[segmentIdx][sampleIdx] = pPrime / p[segmentIdx][frameIdx];
          coefficientsR[segmentIdx][sampleIdx] = qPrime / q[segmentIdx][frameIdx];
          sampleIdx++;
        }          
      }
    }         
    curWatermarKBit++;
  }
  
  double idctMultiplier;
  double outputStreamF = 0.0;
  double outputStreamR = 0.0;
  for (int segmentIdx = 0; segmentIdx < NUM_SEGMENTS; ++segmentIdx)
  {
    for (int sampleIdx = 0; sampleIdx < L_OVER_TWO; ++sampleIdx)
    {
      for (int coefficientIdx = 0; coefficientIdx < NUM_COEFFICIENTS; ++coefficientIdx)
      {
        idctMultiplier = (sampleIdx) ? 4 * ONE_OVER_SQRT_L : SQRT_TWO * ONE_OVER_SQRT_L;
        outputStreamF += idctMultiplier * coefficientsF[segmentIdx][coefficientIdx] * cos(PI * (2 * segmentIdx + 1) * coefficientIdx) * ONE_OVER_L;
        outputStreamR += idctMultiplier * coefficientsR[segmentIdx][coefficientIdx] * cos(PI * (2 * segmentIdx + 1) * coefficientIdx) * ONE_OVER_L;            
      }
    
      outputStream[segmentIdx][sampleIdx] = (int)outputStreamF;
      outputStream[segmentIdx][sampleIdx + L_OVER_TWO] = (int)outputStreamR;
    }
  }
}

int main(char* argv, int argn)
{
	return -1;
}