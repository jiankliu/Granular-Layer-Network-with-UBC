

//use nr codes
#include "../nr3/nr3.h"
#include "../nr3/ran.h"
#include "../nr3/gamma.h"
#include "../nr3/deviates.h"

#include "param.h"          // file for all parameters 

//---------------CONSTANTs initialization----------------------------
const double  PI = 3.1415926;

#define Max     ((NE > Mcx) ? NE : Mcx)
#define GetRandom(max) ((rand()%(int)((max)))) // random integer between 0 and max-1
#define GetRandomDoub(max) ((double)(rand()/(double)(RAND_MAX)*((max)))) // random double between 0 and max

//------------first order kiner model for slow AMPA synapse---------------------
class AMPAsl {
	static double Cdur, Deadtime; 
	double alpha, C, S_nmda;
	double lastrelease, releaseat, TimeCount;

	//SHORT-TERM PLASTICITY
    double  q, R, u; 

public:
	double nmda, tauR, tauD;
    double  U, trec, tfac;
    double  Delay;
	AMPAsl() {
		Delay = 1.0;
		q = 0;

		C = 0;
		nmda = 0;
		S_nmda = 0;
		tauR = 163;
		tauD = 504;
		alpha = 1./tauR * pow(tauD/tauR, tauR/(tauD-tauR));

		lastrelease = -9999;
		releaseat = 0;
		TimeCount = -1.0;	//(ms)		: time counter

		tfac = 0;     // (ms)	this value should be close to 0 for no facilitation 
		trec = 500;   // (ms) 	recovery from depression time constant 
		U    = 0.5;   //	    percent of transmitter released on first pulse
		R    = 1;     //        Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
	}

	void calc(double lastspike, double x);
};

double AMPAsl::Cdur = 1, AMPAsl::Deadtime = 1;

void AMPAsl::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeCount -= DT;                 //time since last release ended

	if (TimeCount < -Deadtime+TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeCount = Cdur;
			
			// Synaptic connections displaying depression are characterized by
			// negligible values of tfacil and hence 
			R=1+(R*(1-u)-1)*exp(-q/trec);
			u= U + (1-U) * u * exp(-q/tfac);			
			
			C = R*u;	// start new release, turn on 
			S_nmda += C;
		}
    } else if (TimeCount > 0) {		//: still releasing?
									// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	
	nmda += DT * (- nmda /tauD + 0.015 * S_nmda * (1-nmda) );
	S_nmda += -DT * S_nmda/tauR;
}

//------------first order kiner model for AMPA and NMDA synapse---------------------
class AMPA_NMDA2 {
	static double Cdur, Deadtime; 
	double Alpha1, Alpha2, C, S_ampa, S_nmda;
	double lastrelease, releaseat, TimeCount;

	//SHORT-TERM PLASTICITY
    double  q, R, u, R2, u2; 

public:
	double ampa, nmda, tauR1, tauD1, tauR2, tauD2;
    double  U, trec, tfac, U2, trec2, tfac2;
    double  Delay;
	AMPA_NMDA2() {
		Delay = 1.0;
		q = 0;

		ampa = 0, C = 0;
		nmda = 0;
		S_ampa = 0;
		S_nmda = 0;
		tauR1 = 0.4;
		tauD1 = 2;
		Alpha1 = 1./tauR1 * pow(tauD1/tauR1, tauR1/(tauD1-tauR1));
		tauR2 = 5;
		tauD2 = 100;
		Alpha2 = 1./tauR2 * pow(tauD2/tauR2, tauR2/(tauD2-tauR2));

		lastrelease = -9999;
		releaseat = 0;
		TimeCount = -1.0;	//(ms)		: time counter

		tfac = 0;     // (ms)	this value should be close to 0 for no facilitation 
		trec = 500;   // (ms) 	recovery from depression time constant 
		U    = 0.5;   //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 

		tfac2 = 400;     // (ms)	this value should be close to 0 for no facilitation 
		trec2 = 20;      // (ms) 	recovery from depression time constant 
		U2    = 0.05;    //	        percent of transmitter released on first pulse
		R2    = 1;       //  Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
		u2   = U2;    //  for running value of U	
	}

	void calc(double lastspike, double x);
};

double AMPA_NMDA2::Cdur = 1, AMPA_NMDA2::Deadtime = 1;

void AMPA_NMDA2::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeCount -= DT;                 //time since last release ended

	if ( TimeCount < -Deadtime+TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeCount = Cdur;
			
			// Synaptic connections displaying depression are characterized by
			// negligible values of tfacil and hence 
			R=1+(R*(1-u)-1)*exp(-q/trec);
			u= U + (1-U) * u * exp(-q/tfac);			
			C = R*u;	//			: start new release, turn on 
			S_ampa += R*u;

			R2 = 1 + (R2*(1-u2)-1) * exp(-q/trec2);
			u2 = U2 + (1-U2) * u2 * exp(-q/tfac2);			
			S_nmda += R2*u2;
		}
    } else if (TimeCount > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	
	ampa += DT * (- ampa /tauD1 + 3. * S_ampa * (1-ampa) );
	S_ampa += -DT * S_ampa/tauR1;
	nmda += DT * (- nmda /tauD2 + 0.35 * S_nmda * (1-nmda) );
	S_nmda += -DT * S_nmda/tauR2;
}


//------------ e(i)MF-GC inputs plastic soma ---------------------
class MGsoma: public AMPA_NMDA2, public AMPAsl {

public:
	double ampaMG,nmdaMG,ampaMGsl;
	double ampa_fast_trec, ampa_slow_trec;
	MGsoma() {
		ampaMG=0;
		nmdaMG=0;
		AMPA_NMDA2::U = 0.5;
		AMPA_NMDA2::tfac = 12;
		AMPA_NMDA2::trec = 12; 
		AMPA_NMDA2::tauR1 = 0.3;
		AMPA_NMDA2::tauD1 = 0.8;

		AMPA_NMDA2::U2    = 0.05;
		AMPA_NMDA2::tfac2 = 0.001;
		AMPA_NMDA2::trec2 = 0.001; 
		AMPA_NMDA2::tauR2 = 8.;
		AMPA_NMDA2::tauD2 = 30.;

		ampaMGsl=0;
		AMPAsl::U = 0.5;
		AMPAsl::tfac = 12;
		AMPAsl::trec = 12; //12
		AMPAsl::tauR = 0.5;
		AMPAsl::tauD = 5.0;
	}

 	void para(){
		AMPA_NMDA2::trec = ampa_fast_trec; //12
		AMPAsl::trec = ampa_slow_trec; //12
	} 

	void calc(double lastspk, double x){
		AMPA_NMDA2::calc(lastspk, x); 
		ampaMG = AMPA_NMDA2::ampa;
		nmdaMG = AMPA_NMDA2::nmda;
		AMPAsl::calc(lastspk, x); 
		ampaMGsl = AMPAsl::nmda;
	} 
};



//------------ Ex MF current ---------------------
class IMFGC {
	static double E1, E2;

public:
	double I, AMPA_NMDA_RATIO, AMPA_AMPAsl_RATIO;
	double I1;
	IMFGC() {
		I = 0;
		I1 =0;
		AMPA_NMDA_RATIO = NAR_MG;
		AMPA_AMPAsl_RATIO = NAR_MG2;
	}
	void calc(double g_AMPA, double ampa, double nmda, double ampasl, double y_post,  double postB);
};

double IMFGC::E1 = 0., IMFGC::E2 = 2.404;

void IMFGC::calc(double g_AMPA, double ampa, double nmda, double ampasl, double y_post,  double postB) {
	I = 1 * g_AMPA * ampa * (y_post - E1) + 1* g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2)
		+ 1* g_AMPA * AMPA_AMPAsl_RATIO * ampasl * (y_post - E1);
}


/*
Merged  INPUT.mod and SpkGen.mod
INPUT.mod is used as the input, and drives the synapses to the first layer Ex cells
SpkGen.mod is used to drive the input with patterns	
Dean Buonomano 4/1/01

modified from pregen.mod - Mainen and Destexhe 
modified to incorporate arbitrary sequences of pulses with 
assigned inter-burst intervals
spkgenmode = 0 is the same as the previous version
spkgenmode = 1 permits the user to assign inter-burst intervals
in an aperiodic fashion. In this mode the on_times variable
must be initialized to +9e9 in the .oc or .hoc file.
Potential bug was fixed by adding dt/2
Dean Buonomano 5/20/99
*/


class INPUT: public Poissondev {
  static double spikedur, refact;
  double on_times[10];
  int spkgenmode;
  double burst;
  int    pulsecntr;

  double lcntr, scntr, bcntr;

public:
  double Istim, Spike;
  double noise, fast_invl, slow_invl, burst_len, start, end;
  INPUT(): Poissondev(1,SEED) {
	  spkgenmode = 0;
	  fast_invl  = 1;      //: time between spikes in a burst (msec)
	  slow_invl  = 50;     //: burst period (msec)
	  burst_len  = 1.;     //: burst length (# spikes)
	  start      = 50.;    //: location for first burst
	  end        = 200.;   //: time to stop bursting
	  noise      = 0;

	  Spike      = 0.;  
	  pulsecntr  = -1;
	  burst      = 0;
  }

  void reset_input() {
	  scntr = fast_invl;
	  bcntr = burst_len;

	  if (noise != 0) {
		  scntr	= noise * Poissondev::dev(fast_invl ) + (1 - noise) * fast_invl;
		  bcntr	= noise * Poissondev::dev(burst_len ) + (1 - noise) * burst_len;
	  }
	  lcntr = start;
  }

  void calc( double y_post, double x);
};

void INPUT::calc(double y_post, double x) {
	if ( x < end ) {
		Spike = 0.;
		scntr -= DT;
		lcntr -= DT;

		if ( burst ) {
			// in a burst
			if ( scntr < DT + TERR ) {        
				// a spike
				//if ( fabs(bcntr-spikedur) < TERR ) {
				if ( bcntr <= 1.0 ) {
					// last spike in burst?
					burst = 0;
					if ( spkgenmode ==0 ) {
						if ( noise ==0 ) {
							lcntr = slow_invl - (burst_len-1) * fast_invl;
						} else {
							lcntr = noise * Poissondev::dev(slow_invl) + (1 - noise) * slow_invl;
						}					
					} else if ( spkgenmode ==1 ) {
						lcntr = on_times[pulsecntr] + DT;
					}
				}

				Spike = 10;
				bcntr -= 1;

				if (noise==0) {
					scntr = fast_invl + DT;
				} else {
					scntr = noise * Poissondev::dev(fast_invl) + (1 - noise) * fast_invl; 
				}
			}
		} else {                  
			//  between burst
			if (lcntr <= DT + TERR) {	//	: there is a round off problem here
				if (noise==0) { 
					bcntr = burst_len;
				} else {
					bcntr = noise * Poissondev::dev(burst_len) + (1 - noise) * burst_len;
				}
				burst = 1;
				pulsecntr = pulsecntr + 1;	
			}
		} 
	} // x end
	Istim = Spike * ( y_post  + 20 );
}


//-----------------------------------------------------------------------
//            Define Cells
//-----------------------------------------------------------------------

//-------------------extral MF -----------------------------------------------
class MFIAF: public MGsoma, public INPUT {
	static double V0; 

	static double VON, VOFF; 
	static double spkdur, refact; 
	double Imax, TimeCounter;

public:
	double Cm, G_l, E_l, lastspk;
	int    TYPE, Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
    
    double Vthr, q;
	double tauAHP, gAHPbar, eAHP;
	double k, mu, Inoise, meanRate, phase, A, vstim, v1, v2, sigma;

	MFIAF(){
		TYPE = 0;
		k = 0;
		mu = 0;
		meanRate = 0;
		phase = 0;
		A = 0;
		vstim = 0;
		v1 = 0;
		v2 = 0.0;
		sigma = 0.0;
		
		Cm = 1;
		E_l  = -60.;      // mV    10^-3 
        G_l  = 0.025;     // mS/cm^2  The maximum specific leakage conductance
		lastspk = -9999;
        TimeCounter = -1.0;
		Ca   = 0;
		postB = 0;

		Imax = 1;     //  uA/cm^2
		Inoise = 0;   //  uA/cm^2
		Vthr = -40;

  	    tauAHP   = 10;
		gAHPbar = 0.00007 * 1e3;   // mS/cm2; peak of AHP current
		eAHP = -90;

		q= 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset(){
		Ca = 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
        y[0] = V0;
        y[1] = 0; //gAHP;
	}
	void calc(double x, double *y, double *f); 
};  

double MFIAF::V0 = -60;

double MFIAF::VON=40, MFIAF::VOFF=-60;
double MFIAF::spkdur=1, MFIAF::refact=2;


void MFIAF::calc(double x, double *y, double *f)
{
    q = (x- lastspk) - spkdur; //time since last spike ended
    TimeCounter -= DT;              

	if (TimeCounter < -refact+TERR) {
    // refactory period over?
		if ( y[0]>Vthr) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike counted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca ++;
		}
    } else if (TimeCounter > 0-TERR) {	// spike still on
			y[0] = VON;
    } 	else {                  // spike off, refactory period on
		if ( y[0] > 0) {	 	// turn spike off
			y[0] = VOFF;
			y[1] += gAHPbar;
		}
	}

	if ( Inoise < -1. ) {
		y[0] = VON;
		lastspk = x;
	} else {
		y[0] = VOFF;
	}
	f[0] = 0.0;
	//*/
	//f[0] = (1./Cm) * ( -G_l*(y[0]-E_l) -y[1]*(y[0]-eAHP) - Inoise ); 

    f[1] = - y[1]/tauAHP;

	MGsoma::calc(lastspk, x);
}

//-------------------Grangular cell model from ExpIAF by N. Brunel -----------------------------------------
class GCmodel: public Normaldev {
	double V0, VOFF; 
	static double VON; 
	static double spkdur, refact; 
	double lastspk, TimeCounter;

public:
	double Cm, G_l, E_l, taum, Ginh;
	int   Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
    
    double Vthr, q;
	double tauAHP, gAHPbar, eAHP;

	double DeltaT, tauxAHP, xAHP, xAHPbar, G_ka, tauN, sigma;

	GCmodel(): Normaldev ( 0.0, 1.0, rand()*rand() ){
		Ginh = 0.0;
		Cm = 4.9;     // pF
		taum = 7.0;   // ms			

		// I_Na
		DeltaT = 1.0;  // mV
		Vthr = -50;
		E_l  = -90.0;  // mV
		V0 = -85.0;
		VOFF = -65.0;
        G_l  = 1.5;    // S  The maximum specific leakage conductance
		lastspk = -9999;
        TimeCounter = -1.0 * spkdur;
		Ca   = 0;
		postB = 0;

		tauAHP  = 3.0;
		gAHPbar = 1.0; // nS; peak of AHP current
		eAHP = -90;
		tauxAHP = 1.0;
		xAHP = 0.0;
		xAHPbar = 1.0;
		G_ka = 0; //1.2; // nS

		//noise
		tauN = 1000.0;
		sigma = 0.12;


		q= 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset(){
		Ca = 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
        y[0] = V0;
        y[1] = 0.0; // xAHP;
		y[2] = 0.0; // zAHP;
		y[3] = 0.0; // noise
	}

	void calc(double x, double *y, double *f); 
};  

double GCmodel::VON=40;
double GCmodel::spkdur=0.6, GCmodel::refact=2;

void GCmodel::calc(double x, double *y, double *f)
{
    q = (x- lastspk) - spkdur; //time since last spike ended
    TimeCounter -= DT;              
	
	if (TimeCounter < -refact+TERR) {
		// refactory period over?
		if ( y[0]>=Vthr) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike counted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca ++;
			y[1] += xAHPbar;
		}
    } else if (TimeCounter > 0-TERR) {	// spike still on
			y[0] = VON;
    } 	else {                  // spike off, refactory period on
		if ( y[0] > 0.0) {	 	// turn spike off
			y[0] = VOFF;
		}
	}

	//cout<<" vm is "<< y[0] <<" mv"<<endl;
	
	f[1] = - y[1]/tauxAHP;

	f[2] =  (1 - y[2]) * y[1] - y[2]/tauAHP;

	f[3] =  (- y[3] + sigma * sqrt(tauN) * Normaldev::dev()/sqrt(DT) )/tauN;

#if V_CLAMP
	y[0] = -65; // -75//mGluR2  -72//ex -70//inh
#else
	f[0] = ( -G_l*(y[0]-E_l)*exp(-(y[0]-E_l)/5.) - (G_ka+y[3])*y[0] -gAHPbar*y[2]*(y[0]-eAHP)  -(TONICI/1000.0*(y[0]+75)) ) / Cm;
#endif

	// from Stephane's data
    postB = (exp((y[0]+119.51)/38.427) + exp(-(y[0]+45.895)/28.357))
		/ (exp((y[0]+119.51)/38.427) + exp(-(y[0]+45.895)/28.357) + exp(-(y[0]-84.784)/38.427));

}

