#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mersenne.h"

/* Mersenne Twister, a pseudorandom number generater */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
static unsigned long state[N]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long *next;

/******************************************/
/***     Mersenne Twister functions     ***/
/******************************************/
/* 
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed) 
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/
/* initializes state[N] with a seed */
void init_genrand(unsigned long s)
{
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    left = 1; initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(init_key, key_length)
unsigned long init_key[], key_length;
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    left = 1; initf = 1;
}

static void next_state(void)
{
    unsigned long *p=state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (initf==0) init_genrand(5489UL);

    left = N;
    next = state;
    
    for (j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* int main(void) */
/* { */
/*     int i; */
/*     unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4; */
/*     init_by_array(init, length); */
    /* This is an example of initializing by an array.       */
    /* You may use init_genrand(seed) with any 32bit integer */
    /* as a seed for a simpler initialization                */
/*     printf("1000 outputs of genrand_int32()\n"); */
/*     for (i=0; i<1000; i++) { */
/*       printf("%10lu ", genrand_int32()); */
/*       if (i%5==4) printf("\n"); */
/*     } */
/*     printf("\n1000 outputs of genrand_real2()\n"); */
/*     for (i=0; i<1000; i++) { */
/*       printf("%10.8f ", genrand_real2()); */
/*       if (i%5==4) printf("\n"); */
/*     } */

/*     return 0; */
/* } */

/************** End of Mersenne Twister ******************/


/************ Some usefull prob. distributions **************/
#define PI 3.141592654

double bnldev(double pp, int n)
{
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (genrand_real2() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= genrand_real2();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*genrand_real2();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (genrand_real2() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software %6<5@)1.#109. */



/*********** Some functions to calculate distributions ******/
/*
  double gammln(double xx)
  Parameter: put xx, you will get ln(Gamma(xx)).
  Description: Numerical Recipe C
  Return value: ln(Gamma(xx));
*/
double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software %6<5@)1.#109. */

/*
  double factln(int n)
  Parameter: put n, you will get ln(n!)
  Description: Numerical Recipe C
  Return value: ln(n!).
*/
double factln(int n)
{
  double gammln(double xx);
  void nrerror(char error_text[]);
  static double a[101];
  
  if (n < 0) fprintf(stderr,"Negative factorial in routine factln\n");
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}
/* (C) Copr. 1986-92 Numerical Recipes Software %6<5@)1.#109. */

/*
  double combinat(int n, int r)
  Parameter: The number of ways to chose r out of n. n!/r!/(n-r)!
  Description:
  Return value: (double) n!/r!/(n-r)!
*/
double combinat(int n, int r)
{
  static double fact[13]={
    1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3628800.,39916800.,479001600.
  };
  
  if(r>n){
    fprintf(stderr,"combinat: invalid parameters. r>n in n!/r!/(n-r)!\n");
    return 0.;
  }
  else if(r==n)
    return 1.;
  else if(n<13)
    return fact[n]/fact[r]/fact[n-r];
  else 
    return floor(0.5+exp(factln(n)-factln(r)-factln(n-r)));
}

/* random intager between [min,max] */
int genrand_int(int min, int max){
  int result, relMax;
  relMax = max-min+1;
  result = genrand_real2()*relMax;
  return result+min;
}


/*
  void  rand_choice(int length,int nchoice,int result[])

  Description: the function will chose nchoice numbers from the
  length. If length=100 and nchoice=13, then the function will
  chose 13 numbers randomly from the range [0,99]. The function
  will put the chosen numbers in result[]. Thus, result[] must be
  long enough to hold the values.

  Caution: the result is not in order.
*/
void rand_choice(int length, int nchoice, int *result)
{
  int i,temp,counter;
  static int *array=NULL;
  static int array_size=0;

  if(length<1){
    fprintf(stderr,"rand_choice(): length<1 is illegal.\n");
    exit(-1);
  }else if(nchoice<1){
    fprintf(stderr,"rand_choice(): nchoice<1 is illegal.\n");
    exit(-1);
  }else if(nchoice>length){
    fprintf(stderr,"rand_choice(): nchoice>length is illegal.\n");
    exit(-1);
  }else if(result==NULL){
    fprintf(stderr,"rand_choice(): result=NULL is illegal.\n");
    exit(-1);
  }
  /* Initilizaion */
  if(array_size<length){
    if((array=realloc(array,sizeof(int)*length))==NULL){
      fprintf(stderr,"rand_choice(): No memory\n");
      exit(-1);
    }
    array_size = length;
  }

  counter = (length+7)/8;
  i=0;
  switch( length%8 )
    {
    case 0: do{ array[i] = i; i++;
      case 7: array[i] = i; i++;
      case 6: array[i] = i; i++;
      case 5: array[i] = i; i++;
      case 4: array[i] = i; i++;
      case 3: array[i] = i; i++;
      case 2: array[i] = i; i++;
      case 1: array[i] = i; i++;
      }while(--counter>0);
    }

  /* select random numbers */
  counter = (nchoice+7)/8;
  i=0;
  switch( nchoice%8 )
    {
    case 0: do {     
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 7:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 6:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 5:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 4:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 3:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 2:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      case 1:
	temp = (length-i)*genrand_real2();
	result[i] = array[temp];
	if(temp!=length-i-1)
	  array[temp] = array[length-i-1];
	i++;
      }while(--counter>0);
    }
}

/* 
   Given we have probability x0,x1,x2,... xn. We chose one of
   them, or may not choose anything If Sum(xi) is more than 1, we
   always choose one. Note that the function first calculates the
   sum of probabilities Sum(xi), and the function will not order
   the array of probabilities. Thus, if there is a huge
   order-wise variation in the values in the array, the user
   should order the array in ascending order.

   Parameters: "nitem" should be the number of probabilities
   (i.e., size of "array"). "array" should contains the
   probabilities. "array" will be transformed into the
   follwoign. x[i] = Sum(x[k]) where k = 0,..,i.

   Return: -1 if non was chosen. Otherwise, the corresponding
   index of probability starting from 0.
*/
int choose_claim(int nitem,double *array)
{
  double rand,sum=0.;
  int counter;
  int index=0;
  int result=0;


  if(nitem==0){
    return -1;
  }

  rand = genrand_real2();

  /* calculate the sum of probability */
  counter = (nitem+7)/8;
  switch( nitem%8 )
    {
    case 0: do{ 
	sum+=array[index];
	array[index]=sum;
	/* We first assume that real_sum<1 and compute the result
	   according to the assumption. If rand>=partial_sum
	   (i.e., the sum up to the current index), then result
	   must be at least one larger. Thus we do result++. If
	   rand<partial_sum, then the result is correct. Thus do
	   not incliment result. */
	if(rand>=sum)
	  result++;
	index++;
      case 7: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 6: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 5: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 4: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 3: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 2: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      case 1: sum+=array[index];array[index]=sum;
	if(rand>=sum) result++;
	index++;
      }while(--counter>0);
    }

  /* if real_sum<=1. then our assumption is correct. Return the
     result. */
  if(sum<=1.){
    if(result < nitem)
      return result;
    else if (result == nitem)
      return -1;
    else{
      fprintf(stderr,"choose_claim(): Bug, result > nitem when sum<=1.\n");
      exit(-1);
    }
  }

  /* if real_sum >1. then we have to re-compute the result */
  /* If "rand < 0.5", we search the array from [0] to [nitem-1]
     for the result. If "rand >= 0.5", then we search the array
     from [nitem-1] to [0]. */
  if(rand<0.5){
    /* we expand result */
    rand *= sum;

    counter = (nitem+7)/8;
    index = 0;
    result = 0;
    switch( nitem%8 )
      {
      case 0: do{ 
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 7:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 6:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 5:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 4:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 3:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 2:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	case 1:
	  if(rand>=array[index]){result++;}
	  else {break;}
	  index++;
	}while(--counter>0);
      }
    if(result<nitem)
      return result;
    else{
      fprintf(stderr,"choose_claim(): Bug, result >= nitem when sum>1.\n");
      exit(-1);
    }
  }else{
    /* we expand result */
    rand *= sum;
    
    counter = (nitem+7)/8;
    index = nitem-1;
    result = nitem;
    switch( nitem%8 )
      {
      case 0: do{ 
	  /* we search the array backwards */
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 7:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 6:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 5:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 4:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 3:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 2:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	case 1:
	  if(rand<array[index]){result--;}
	  else {break;}
	  index--;
	}while(--counter>0);
      }
    if(result>=0 && result < nitem)
      return result;
    else{
      fprintf(stderr,"choose_claim(): Bug, nothing was chosen when sum>1.\n");
      exit(-1);
    }
  }
}
