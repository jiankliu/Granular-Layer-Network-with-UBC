
/* 
 This is the basic code for the neural network model used in the following paper:

 Zampini,V.*, Liu, J.K.*, Diana, M.A., Maldonado, P.P., Brunel, N., Dieudonné, S. 
 (*These authors contributed equally to this work), 
 Mechanisms and functional roles of glutamatergic synapse diversity in a cerebellar circuit. eLfie, 5:e15872 (2016)

 Please fell free to let me know if you have any questions.
 https://sites.google.com/site/jiankliu/home

 Jian K. Liu
 15/01/2017
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <cstdlib> 

using namespace std;

#include "param.h"    // file for all parameters 
#include "synapse.h"  // synapse models

//------------Number of ODE for each cell -------------------------------
const int N_EQ_G = 4;
const int N_EQ_M = 2; 
const int N_EQ_U = 2; 

#define N_EQ   ( N_EQ_G*NGC + N_EQ_U*NUBC + N_EQ_M*NMF )    //  Complete number of ODE

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++

//----------external variables ---------------------------------------------
int no_GC[NGC][N_EQ_G], no_UBC[NUBC][N_EQ_U], no_eMF[NMF][N_EQ_M];  
MatDoub g_eG(NGC,N_eG),    g_iG(NGC,N_iG) ;
MatInt pre_eG(NGC,N_eG, 0),  pre_iG(NGC,N_iG, 0); 

//----------external classes (beginning of initialization)------------------
GCmodel  GC[NGC];  

//mossy fibers
MFIAF  eMF[NMF];
MFIAF  UBC[NUBC];

// postsynapse currernt 
IMFGC *syn_iG[NGC];     // iMF -> GC
IMFGC *syn_eG[NGC];     // eMF -> GC

//   -----external functions----------------------------------------------
void fun(double x, double *y_ini, double *f_ini);

void euler(int n, void fun(double, double*, double*),
		double h, double x, double* y, double* f);


////////////////////////////////////////////////////////////////////////////////////
//++++++Main program+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////////////////////////////////////////////////////////////////////////////////////

int main(int argc,char **argv) 
{
	//MPI
	int num_procs;

	// running time counter
	clock_t start,finish;
	double totaltime;
	start = omp_get_wtime();

	FILE *f1, *f7;

	//---------allocate place for ALL variables and functions-----------
	double y_ini[N_EQ], f_ini[N_EQ];

	//---------allocate place for TWO temporal arrays using by RK.C solver 
	double y1[N_EQ], y2[N_EQ];

	//---------general parameters----------------------------------------------
	double t = 0, h = DT;
	int i, j, k, ii = 0;

	
   //
   //  How many processors are available?
   //
   num_procs = omp_get_num_procs ( );

   cout << "\n";
   cout << "  The number of processors available:\n";
   cout << "  OMP_GET_NUM_PROCS () = " << num_procs << "\n";
   
  int nthreads, thread_id;
  nthreads = 4;
  omp_set_num_threads ( nthreads );

// Fork the threads and give them their own copy of the variables 
#pragma omp parallel private(nthreads, thread_id)
  {
    
    //* Get the thread number 
    thread_id = omp_get_thread_num();
    printf("Thread %d says: Hello World\n", thread_id);

    if (thread_id == 0)     // This is only done by the master 
	{
		nthreads = omp_get_num_threads();
		printf("Thread %d reports: the number of threads are %d\n", thread_id, nthreads);
	}
	
  }    // Now all threads join master thread and they are disbanded 

   /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.
    */
   srand( SEED );

   //----------arrays initialization----------------------------------------------
#pragma omp parallel private(i,j,k)
{

#pragma omp for nowait
   for(i=0; i<N_EQ; i++){
	   y_ini[i] = 0,  f_ini[i] = 0; 
	   y1[i] = 0,     y2[i] = 0;
   }

   //----------classes initialization (continue)----------------------------
#pragma omp for nowait
   for(i=0; i < NGC; i++) {
	   syn_iG[i] = new IMFGC[N_iG];
	   syn_eG[i] = new IMFGC[N_eG];
   }

   //----------creating the integer arrays containing the addresses--------
   //----------of  ALL internal variables for ALL objects        ----------
   //----------and classes (e.g., no_GC[i][j][k] is the address --------	
   //----------of the variable y[k] for the object GC[i][j]) ---------
   //----------NOTE: this is the relative addresses and you should ---------
   //----------add the REAL address of the first element of the -----------
   //----------original 1D array to use these arrays----------------------- 
#pragma omp for nowait
   for(i=0; i < NGC; i++)
	   for(k=0; k < N_EQ_G; k++)
		   no_GC[i][k] = k + i * N_EQ_G;

#pragma omp for nowait
   for(i=0; i < NUBC; i++)
	   for(k=0; k < N_EQ_U; k++)
		   no_UBC[i][k] = NGC*N_EQ_G + k + i * N_EQ_U;  

#pragma omp for nowait
   for(i=0; i < NMF; i++)
	   for(k=0; k < N_EQ_M; k++)
		   no_eMF[i][k] = NGC*N_EQ_G + NUBC*N_EQ_U + k + i * N_EQ_M;  

} //end omp

   //--------------open ALL files-------------------------------------
   if(    ( f1 = fopen("Raster.dat", "w") ) == NULL ){
		   printf("can't open files\n");
		   exit(0);
   }

   //----------Connection matrix-------------------------------
   printf("\n Begin Connect Matrix");

   //-------------- Network topology ------------------------------------
   int pre, exists;
   
   int MFinput = floor(NUBC*0.5);

	  // iMF -> GC
	   for(i=0; i<NGC; i++) {
		   for(j=0; j<N_iG; j++) {
			   do{
				   exists = 0;		// avoid multiple synapses
				   pre = GetRandom(NUBC);
				   for (k=0;k<j;k++) if (pre_iG[i][k]==pre) exists = 1;	// synapse already exists  
			   }while (exists == 1);
			   pre_iG[i][j]= pre;
		   }
	   }
	   		   
	  // eMF -> GC
	  for(i=0; i<NGC; i++) {
		  for(j=0; j<N_eG; j++) {
			   do{
				   exists = 0;		// avoid multiple synapses
				   pre = GetRandom(NMF);
				   for (k=0;k<j;k++) if (pre_eG[i][k]==pre) exists = 1;	// synapse already exists  
			   }while (exists == 1);
			   pre_eG[i][j]= pre;
		  }
	  }

   //--------the initial conductances----------------------------------- 
   // normal distribution of synaptic weights
   double WiG = wig;    
   Normaldev grandWiG( WiG, WiG * var, SEED);
  
   double WeG = weg;   
   Normaldev grandWeG( WeG, WeG * var, SEED);
     
   for(i=0; i<NGC; i++) {
	   for(j=0; j<N_iG; j++) {
		   g_iG[i][j] = grandWiG.dev(); 
		   if (g_iG[i][j] > AMPA_MAXiG )  g_iG[i][j] = AMPA_MAXiG;
		   else if ( g_iG[i][j] < 0. )    g_iG[i][j] = GetRandomDoub(2.*WiG);
	   }

	   for(j=0; j<N_eG; j++) {
		   g_eG[i][j]   = grandWeG.dev(); 
		   if (g_eG[i][j] > AMPA_MAXeG )  g_eG[i][j] = AMPA_MAXeG;
		   else if ( g_eG[i][j] < 0. )    g_eG[i][j] = GetRandomDoub(2.*WeG);
	   }
   }

   double c1 = 0;
   int c2 = 0;
   // select eMF or iMF for each GC
   for(i=0; i<NGC; i++) {
	   c1 = GetRandomDoub(1.0);
	   if (c1 < 0.5) {
		   // 1 eMF 0 iMF 
		   c2 = 1;
	   }else {
		   // 0 eMF 1 iMF
		   c2 = 0;
	   }		   			

	   for(j=0; j<c2; j++) {
		   g_iG[i][j]   = 0; 
	   }
	   for(j=c2; j<N_eG; j++) {
		   g_eG[i][j]   = 0; 
	   }
   }

   //----------here we are changing some variables----------------------
   // get exact stimuli, run rand.repick only at the begining
   for(i=0; i<NMF; i++) {
	   eMF[i].TYPE =  GetRandom(2) + 1;  //randomly select TYPE I or II
	   //if (i< NMF/2 ) {
		//   eMF[i].TYPE =  1;  // only TYPE I
	   //} else {
		 //  eMF[i].TYPE =  2;  // only TYPE II
	   //}
	   eMF[i].k = GetRandomDoub(1);
	   eMF[i].meanRate =  MFRATE;

	   //input
	   eMF[i].slow_invl  = T;
	   eMF[i].start       = TMAX;
	   eMF[i].end        = TMAX;
	   eMF[i].noise = 1;
	   eMF[i].init(y_ini+no_eMF[i][0]);
	   eMF[i].reset();
	   eMF[i].reset_input();

	   eMF[i].ampa_fast_trec = 600; 
	   eMF[i].ampa_slow_trec = 600; 
	   eMF[i].para(); 
   }

   Normaldev r3( 1, varUBC*1, SEED);
   Normaldev r4( 1, varUBC*1, SEED);
   Normaldev grand2( THR_UBC, fabs(THR_UBC * 0.05), SEED);
   Expondev  frate(1./7., SEED);
   Normaldev sig1( -0.5, fabs(0.3), SEED);
   Normaldev sig2( 3.0,  fabs(0.3), SEED);
   Normaldev phaseON( 2./3.*PI,   fabs(2./3.*PI  * 0.1), SEED);
   Normaldev phaseOFF( 4./3.*PI,  fabs(4./3.*PI  * 0.1), SEED);
  
   double UBCratio = 0.0;
   for(i=0; i<NUBC; i++) {
	   UBC[i].TYPE =  GetRandom(2) + 1;  //randomly select TYPE I or II
		   UBC[i].k = GetRandomDoub(1);
		   if (GetRandomDoub(1.0) < 0.5) { //ON type
			   UBC[i].v1 = frate.dev()*1e-3;
			   UBC[i].v2 = frate.dev()*1e-3;
			   if (GetRandomDoub(1.0) < 0.5) {
				   UBC[i].sigma = pow(10., sig1.dev());
			   }
			   else {
				   UBC[i].sigma = pow(10., sig2.dev());
			   }
			   UBC[i].phase = phaseON.dev() - PI / 2;
		   }
		   else {
			   UBC[i].v1 = frate.dev()*1e-3;
			   UBC[i].v2 = frate.dev()*1e-3;
			   if (GetRandomDoub(1.0) > 0.2) {
				   UBC[i].sigma = pow(10., sig1.dev());
			   }
			   else {
				   UBC[i].sigma = pow(10., sig2.dev());
			   }
			   UBC[i].phase = phaseOFF.dev() - PI / 2;
		   }

	   //input
	   UBC[i].slow_invl  = T;
	   UBC[i].start       = TMAX;
	   UBC[i].end        = TMAX;
	   UBC[i].noise = 1;

	   UBC[i].init(y_ini+no_UBC[i][0]);
	   UBC[i].reset();
	   UBC[i].reset_input();

	   //short-term depression
	   UBC[i].ampa_fast_trec = 12; 
	   UBC[i].ampa_slow_trec = 12; 
	   UBC[i].para(); 
   }

   Normaldev grand1( THR_GC, fabs(THR_GC * 0.05), SEED);
   for(i=0; i<NGC; i++) {
	   GC[i].Vthr = grand1.dev();
	   GC[i].init(y_ini+no_GC[i][0]);
	   GC[i].reset();
	   for(j = 0; j < N_eG; ++j){
		   syn_eG[i][j].AMPA_NMDA_RATIO = NAR_MG;
		   syn_eG[i][j].AMPA_AMPAsl_RATIO = NAR_MG2;
	   }
	   for(j = 0; j < N_iG; ++j){
		   syn_iG[i][j].AMPA_NMDA_RATIO = NAR_MG;
		   syn_iG[i][j].AMPA_AMPAsl_RATIO = NAR_MG2;
	   }
   }

   printf("\n Network setting up is done ! ");
   //--------------end variability------------------------------------


   int it = 0, trial = 0;
   double PSC=0.0, GCFR = 0.0, GCFRmean = 0.0, GCfired = 0.0, Gcontrol =0.0, dummy=0. ;
   double currentTIME = 0.0;

   //----------------CALCULATION----------------------------------------
   printf("\n CALCULATION IN PROGICSS!!!: TMAX= %lf (min)", TMAX/ (1000*60));

   while( t <= TMAX){ 
	   if ( t >= ST*SHORTTIME ) {	   
		   //---------- modulate stimulating patterns ----------------------
		   for(i=0; i<NMF; i++) {
			   eMF[i].vstim =  eMF[i].meanRate *A* 1./T*1e3 * sin( 2*PI * t/T  ) ;
			   eMF[i].mu =  eMF[i].k * eMF[i].vstim ;
		   }		   
		   for(i=0; i<NUBC; i++) {
			   dummy = UBC[i].sigma*UBC[i].sigma;
			   UBC[i].vstim = (UBC[i].v2-UBC[i].v1) * ( exp( (cos(2*PI*t/T-UBC[i].phase))/dummy ) - exp(-1/dummy) ) / (exp(1/dummy)-exp(-1/dummy) );
			   UBC[i].mu =  UBC[i].k * UBC[i].vstim ;
		   }
	   } else { 
		   //---------- modulate stimulating patterns ----------------------
		   for(i=0; i<NMF; i++) {
			   eMF[i].mu =  0.;
		   }
		   for(i=0; i<NUBC; i++) {
			   UBC[i].mu =  0.;
		   }
	   }
	   // bidirectional varying stimuli
	   for(i=0; i<NMF; i++) {
			   if ( ( eMF[i].TYPE==1) && (GetRandomDoub(1.0) < DT * (eMF[i].meanRate+eMF[i].mu) ) ) {
				   eMF[i].Inoise =  -200;
			   } else if ( ( eMF[i].TYPE==2) && (GetRandomDoub(1.0) < DT * (eMF[i].meanRate-eMF[i].mu) ) ){
				   eMF[i].Inoise =  -200;
			   } else {
				   eMF[i].Inoise = GetRandomDoub(1.0)-0.5;    
			   } 
	   }
		
	   for(i=0; i<NUBC; i++) {
		   if ( ( UBC[i].TYPE==1) && (GetRandomDoub(1.0) < DT * (UBC[i].v1+UBC[i].mu) ) ) {
			   UBC[i].Inoise =  -200;
		   } else if ( ( UBC[i].TYPE==2) && (GetRandomDoub(1.0) < DT * (UBC[i].v1-UBC[i].mu) ) ){
			   UBC[i].Inoise =  -200;
		   } else {
			   UBC[i].Inoise = GetRandomDoub(1.0)-0.5;    
		   } 
	   }

	   euler(N_EQ, fun, h, t, y_ini, f_ini);

	   t += h;
	   ii++ ;

	   //--------  at the end of one trial -----------------------------------
	   if( ii % int(SHORTTIME/h) ==0 ) {
		   it ++;
		   trial = it;

		   printf("Trial = %d \n", it); 	  
		  
		   //------save spike data -----------------------------------------------------------------
		   fprintf(f1, " %4d %4d %4d\n",9999, it, 1 );
		   
		   GCFRmean = 0.0;
		   GCfired = 0.0;
		   for(i = 0; i < NGC; i++) {
			   for (j =0; j<GC[i].Ca; j++)  {
				   fprintf(f1,"%4d %4d %lf \n", trial, i, GC[i].spiketimes[j]- (it-1)*SHORTTIME );
			   }
			   GCfired += 1.0;
			   GCFR = GC[i].Ca / (T/1000);
			   GCFRmean += GCFR;			   
			   GC[i].Ginh += ALPHA * (GCRATE - GCFR );
		   }
		   GCFRmean /= NGC;

		   printf("Trial=%4d, GCfired=%lf, GCrate=%lf\n", trial, GCfired, GCFRmean);

		   for(i = 0; i < NUBC; i++) {
			   for (j =0; j < UBC[i].Ca; j++)  {	
				   fprintf(f1,"%4d %4d %lf \n", trial, NGC+i, UBC[i].spiketimes[j] - (it-1)*SHORTTIME );
			   }
		   }
		   
		   // reset variables
	       for(i=0; i<NMF; i++)    eMF[i].reset();
		   for(i=0; i<NGC; i++)    GC[i].reset();
		   for(i=0; i<NUBC; i++)   UBC[i].reset();
	  }// end of one trial 
  }//--------------------END CALCULATION-------------------------------
      
  //-----------------close ALL files-----------------------------------
  fprintf(f1, " %4d %4d %4d\n",9999, NTRIAL, 1 );
  fclose(f1);


  // free memory
  for(i=0; i < NGC; i++) {
	  delete [] syn_iG[i];
	  delete [] syn_eG[i];
  }

   finish = omp_get_wtime();
   totaltime = (double)(finish-start);
   cout<<"\n running time is "<<totaltime<<" sec"<<endl;

   return 0;
}



//----------external functions----------------------------------------------
void fun(double x, double *y_ini, double *f_ini){
	int i, j; 

        //========here the MAIN loop to calculate intrinsic conductances===========
    	//--------(f_ini IS changed, y_ini IS NOT changed)-------------------------
#pragma omp parallel private(i,j)
{
#pragma omp for nowait
	for(i=0; i<NGC; ++i)    GC[i].calc( x, y_ini+no_GC[i][0],  f_ini+no_GC[i][0]); 
#pragma omp for nowait
	for(i=0; i<NUBC; ++i) { 
		UBC[i].calc(x, y_ini+no_UBC[i][0], f_ini+no_UBC[i][0]); 
	}
#pragma omp for nowait
	for(i=0; i<NMF; ++i) {
		eMF[i].calc(x, y_ini+no_eMF[i][0], f_ini+no_eMF[i][0]); 
	}
#pragma omp barrier

		//========here the MAIN loop to calculate synaptic conductances=============
		//--------(f_ini IS changed, y_ini IS NOT changed) -------------------------
	    // i-post; j-pre
		// update GC post synapse
#pragma omp for nowait
		for(i = 0; i < NGC; ++i)	{						

			//-----------AMPA and NMDA from iMF to GC cells--------------------------------------
			for(j = 0; j < N_iG; ++j){
				syn_iG[i][j].calc( g_iG[i][j], UBC[ pre_iG[i][j] ].ampaMG, UBC[ pre_iG[i][j] ].nmdaMG, UBC[ pre_iG[i][j] ].ampaMGsl, 
					y_ini[no_GC[i][0]], GC[i].postB );
				f_ini[no_GC[i][0]] -= GC[i].Ginh * syn_iG[i][j].I / GC[i].Cm ;
			}

			//-----------AMPA and NMDA from eMF to GC cells--------------------------------------
			for(j = 0; j < N_eG; ++j) {
				syn_eG[i][j].calc( g_eG[i][j], eMF[ pre_eG[i][j] ].ampaMG, eMF[ pre_eG[i][j] ].nmdaMG, eMF[ pre_eG[i][j] ].ampaMGsl, 
					y_ini[no_GC[i][0]], GC[i].postB );
				f_ini[no_GC[i][0]] -= GC[i].Ginh * syn_eG[i][j].I / GC[i].Cm;
			}
		}
}
		//=============END of MAIN loop==============================================
}


// 1st order Rounge-Kutta ((forward) Euler method) solver for ODE -------------------------------------------------
void euler(int n, void fun(double, double*, double*), 
        double h, double x, double* y, double* f)
{
	int i;

	//k1
	fun(x, y, f);    

#pragma omp parallel for
	for(i = 0; i < n; ++i) 	{
		y[i] += h * f[i]; 
	}
}

