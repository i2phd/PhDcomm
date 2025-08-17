#include "reedsolomon.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// galois arithmetic tables
int gexp[512];
int glog[256];
int polydeg;
int pBytes[MAXDEG];
int synBytes[MAXDEG];
int genPoly[MAXDEG*2];
int Lambda[MAXDEG];
int Omega[MAXDEG];
int ErrorLocs[512];
int NErrors;
int ErasureLocs[512];
int NErasures;


int gmult(int a, int b)
{
  int i,j;

  if (a==0 || b == 0) return (0);
  i = glog[a];
  j = glog[b];
  return (gexp[i+j]);
}

int ginv (int elt)
{
  return (gexp[255-glog[elt]]);
}

void add_polys (int dst[], int src[])
{
  int i;
  for (i = 0; i < polydeg; i++) dst[i] ^= src[i];
}

void copy_poly (int dst[], int src[])
{
  int i;
  for (i = 0; i < polydeg; i++) dst[i] = src[i];
}

void scale_poly (int k, int poly[])
{
  int i;
  for (i = 0; i < polydeg; i++) poly[i] = gmult(k, poly[i]);
}

void zero_poly (int poly[])
{
  int i;
  for (i = 0; i < polydeg; i++) poly[i] = 0;
}

void mul_z_poly (int src[])
{
  int i;
  for (i = polydeg-1; i > 0; i--) src[i] = src[i-1];
  src[0] = 0;
}

void compute_genpoly (int nbytes, int genpoly[]);

void initialize_ecc(int npar)
{
	polydeg = npar*2;
	compute_genpoly(npar, genPoly);
}

void build_codeword (unsigned char msg[], int nbytes, int npar, unsigned char dst[])
{
  int i;

  for (i = 0; i < nbytes; i++) dst[i] = msg[i];

  for (i = 0; i < npar; i++)   dst[nbytes+i] = pBytes[npar-1-i];
}

void decode_data(unsigned char data[], int nbytes)
{
  int i, j, sum;
  for (j = 0; j < NPAR;  j++)
  {
	sum	= 0;
	for (i = 0; i < nbytes; i++)
	{
	  sum = data[i] ^ gmult(gexp[j+1], sum);
	}
	synBytes[j]  = sum;
  }
}

int check_syndrome ()
{
 int i;
 for (i =0 ; i < NPAR; i++)
 {
  if (synBytes[i] != 0) return 1;
 }

 return 0;
}

void mult_polys (int dst[], int p1[], int p2[]);

void compute_genpoly (int nbytes, int genpoly[])
{
  int i, tp[512], tp1[512];

  zero_poly(tp1);
  tp1[0] = 1;

  for (i = 1; i <= nbytes; i++)
  {
	zero_poly(tp);
	tp[0] = gexp[i];
	tp[1] = 1;

	mult_polys(genpoly, tp, tp1);
	copy_poly(tp1, genpoly);
  }
}

void encode_data (unsigned char msg[], int nbytes, unsigned char dst[])
{
  int i, LFSR[MAXNPAR+1],dbyte, j;

  for(i=0; i < NPAR+1; i++) LFSR[i]=0;

  for (i = 0; i < nbytes; i++)
  {
	dbyte = msg[i] ^ LFSR[NPAR-1];
	for (j = NPAR-1; j > 0; j--) {
	  LFSR[j] = LFSR[j-1] ^ gmult(genPoly[j], dbyte);
	}
	LFSR[0] = gmult(genPoly[0], dbyte);
  }

  for (i = 0; i < NPAR; i++)
	pBytes[i] = LFSR[i];

  build_codeword(msg, nbytes, NPAR, dst);
}

void Modified_Berlekamp_Massey (void);
void Find_Roots (void);
void compute_modified_omega (void);

int correct_errors_erasures (unsigned char codeword[],
			 int csize,
			 int nerasures,
			 int erasures[])
{

  int r, i, j, err;

  NErasures = nerasures;
  for (i = 0; i < NErasures; i++) ErasureLocs[i] = erasures[i];

  Modified_Berlekamp_Massey();
  Find_Roots();


  if ((NErrors <= NPAR) && NErrors > 0)
  {

	for (r = 0; r < NErrors; r++)
	{
	  if (ErrorLocs[r] >= csize)
	  {
		return(0);
	  }
	}

	for (r = 0; r < NErrors; r++)
	{
	  int num, denom;
	  i = ErrorLocs[r];

	  num = 0;
	  for (j = 0; j < MAXDEG; j++)
		  num ^= gmult(Omega[j], gexp[((255-i)*j)%255]);

	  denom = 0;
	  for (j = 1; j < MAXDEG; j += 2)
	  {
		denom ^= gmult(Lambda[j], gexp[((255-i)*(j-1)) % 255]);
	  }

	  err = gmult(num, ginv(denom));

	  codeword[csize-i-1] ^= err;
	}
	return(1);
  }
  else
  {
	return(0);
  }
}

void init_gamma (int gamma[]);
int compute_discrepancy (int lambda[], int S[], int L, int n);

void Modified_Berlekamp_Massey (void)
{
  int n, L, L2, k, d, i;
  int psi[MAXDEG], psi2[MAXDEG], D[MAXDEG];
  int gamma[MAXDEG];

  init_gamma(gamma);

  copy_poly(D, gamma);
  mul_z_poly(D);

  copy_poly(psi, gamma);
  k = -1; L = NErasures;

  for (n = NErasures; n < NPAR; n++)
  {

	d = compute_discrepancy(psi, synBytes, L, n);

	if (d != 0)
	{

	  for (i = 0; i < polydeg; i++) psi2[i] = psi[i] ^ gmult(d, D[i]);


	  if (L < (n-k))
	  {
		L2 = n-k;
		k = n-L;
		for (i = 0; i < polydeg; i++) D[i] = gmult(psi[i], ginv(d));
		L = L2;
	  }

	  for (i = 0; i < polydeg; i++) psi[i] = psi2[i];
	}

	mul_z_poly(D);
  }

  for(i = 0; i < polydeg; i++) Lambda[i] = psi[i];
  compute_modified_omega();
}


void Find_Roots (void)
{
  int sum, r, k;
  NErrors = 0;

  for (r = 1; r < 256; r++)
  {
	sum = 0;
	for (k = 0; k < NPAR+1; k++) {
	  sum ^= gmult(gexp[(k*r)%255], Lambda[k]);
	}
	if (sum == 0)
	{
	  ErrorLocs[NErrors] = (255-r); NErrors++;
	}
  }
}

void compute_modified_omega ()
{
  int i;
  int product[MAXDEG*2];

  mult_polys(product, Lambda, synBytes);
  zero_poly(Omega);
  for(i = 0; i < NPAR; i++) Omega[i] = product[i];

}

void mult_polys (int dst[], int p1[], int p2[])
{
  int i, j;
  int tmp1[MAXDEG*2];

  for (i=0; i < (polydeg*2); i++) dst[i] = 0;

  for (i = 0; i < polydeg; i++)
  {
	for(j=polydeg; j<(polydeg*2); j++)
	   tmp1[j]=0;

	for(j=0; j<polydeg; j++) tmp1[j]=gmult(p2[j], p1[i]);
	for (j = (polydeg*2)-1; j >= i; j--) tmp1[j] = tmp1[j-i];
	for (j = 0; j < i; j++) tmp1[j] = 0;

	for(j=0; j < (polydeg*2); j++) dst[j] ^= tmp1[j];
  }
}

void init_gamma (int gamma[])
{
  int e, tmp[MAXDEG];

  zero_poly(gamma);
  zero_poly(tmp);
  gamma[0] = 1;

  for (e = 0; e < NErasures; e++) {
	copy_poly(tmp, gamma);
	scale_poly(gexp[ErasureLocs[e]], tmp);
	mul_z_poly(tmp);
	add_polys(gamma, tmp);
  }
}

void compute_next_omega (int d, int A[], int dst[], int src[])
{
  int i;
  for ( i = 0; i < polydeg;  i++) {
	dst[i] = src[i] ^ gmult(d, A[i]);
  }
}

int compute_discrepancy (int lambda[], int S[], int L, int n)
{
  int i, sum=0;

  for (i = 0; i <= L; i++)
	sum ^= gmult(lambda[i], S[n-i]);
  return (sum);
}

void init_galois_tables (void)
{
  int i, z;
  int pinit,p1,p2,p3,p4,p5,p6,p7,p8;

  pinit = p2 = p3 = p4 = p5 = p6 = p7 = p8 = 0;
  p1 = 1;

  gexp[0] = 1;
  gexp[255] = gexp[0];
  glog[0] = 0;

  for (i = 1; i < 256; i++)
  {
    pinit = p8;
    p8 = p7;
	p7 = p6;
    p6 = p5;
    p5 = p4 ^ pinit;
    p4 = p3 ^ pinit;
    p3 = p2 ^ pinit;
    p2 = p1;
    p1 = pinit;
    gexp[i] = p1 + p2*2 + p3*4 + p4*8 + p5*16 + p6*32 + p7*64 + p8*128;
    gexp[i+255] = gexp[i];
  }

  for (i = 1; i < 256; i++)
  {
    for (z = 0; z < 256; z++)
    {
      if (gexp[z] == i)
      {
	      glog[i] = z;
	      break;
      }
    }
  }
}
