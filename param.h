
#include <math.h>

#ifndef PARAM_H
#define PARAM_H

/************** Defines some global paramaters -- *************/
	
//  units:
//  (nS) = 
//  (nA) = (nanoamp)
//  (mV) = (millivolt)
//  (umho) = (micromho)
//  (mM) = (milli/liter)
//  t = msec

#define WINVC 1
#define UNIX  0

//==============================================================================
//     [simulation parameters and settings]
//==============================================================================

//------------------------------------------------------------------------------
// general settings
//------------------------------------------------------------------------------

const int NTRIAL = 2;
const int ST = 2;         // time without stimuli

const double T = 333;  // [5000 2500 1667 1250 1000 833 714 625 555 500] (ms) period of stimulus 
const double SHORTTIME = T;
const double TMAX =  (NTRIAL+ST) * SHORTTIME;       //simulation time (m sec)

const int SEED = 33;  // seed for the random number sequence
const int RANDIN = 1;

const double TONICI = 900;  // 

const double THR = 0;   //  threshold to cutoff the amplitude of velocity

const int STIMON = 1;  // 1: apply stimuli; 0: no external stim
const double MFRATE  = 26 * 1.0 * 1e-3;     // spontanous firing rate
const double A = 5./3.;      // ( o/s) amplitude of stimulus velocity

const double GCRATE = 5;     // average GC firing rate
const double ALPHA = 0.01;   // time constant of controling

//------------------------------------------------------------------------------
//  network simulation parameter
//------------------------------------------------------------------------------
const int NMF = 5*10;
const int NGC = 5*90;
const int NUBC = 50;

//------------------------------------------------------------------------------
//  stimulus parameter
//------------------------------------------------------------------------------

// Number of input synapses for post cell 
const int N_iG = 4;     // iMF -> GC 4
const int N_eG = N_iG;  // eMF -> GC 4
const int N_MU = 1;     //  MF -> UBC

const int NMFG  = NGC * N_eG ;
const int NMFU  = NUBC * N_MU;

//------------------------------------------------------------------------------
//  synaptic simulation parameter
//------------------------------------------------------------------------------
// synapse weight in units of nS
const double AMPA_MAXiG = 30;     // 1.6  for single MF-GC
const double AMPA_MAXeG = 30;     // 1.6  for single MF-GC
const double weg = 0.4; //1.2  1.6;     
const double wig = 1.6;                 
   
const double NAR_MG   = 2.4;       // MF-GC   = NMDA  / AMPA  
const double NAR_MG2  = 2.0;       // MF-GC   = AMPAslow  / AMPA  

// for gaussian
const double var = 0.3;     //variance for weights with var = avgW*SDW
const double varUBC = 0.3;  //variance for weights with var = avgW*SDW

//------------------------------------------------------------------------------
//  cell properties
//------------------------------------------------------------------------------
const double THR_GC = -49;
const double THR_UBC = -50;    
const int SPKTIME = 2000;

const double DT = 0.1;  
const double TERR = DT; 
//------------------------------------------------------------------------------

#endif
