//------------------------------------------------------------------------------
//
// TVit211 component - Viterbi k=11, rate=1/2 decoder (convertito da k=9)
//
//------------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

//------------------------------------------------------------------------------

#include "Vit211.h"

//------------------------------------------------------------------------------
#pragma resource "Vit211.res"
#pragma package(smart_init)

//------------------------------------------------------------------------------
//
// VALIDATOR
//
//------------------------------------------------------------------------------

static inline void ValidCtrCheck(TVit211*)
{
   new TVit211(NULL);
}

//------------------------------------------------------------------------------
// TVit211: constructor

__fastcall TVit211::TVit211(TComponent *aOwner)
   : TComponent(aOwner), Init(false), P_init(false)
{
  partab_init();
}

//------------------------------------------------------------------------------
// TVit211: destructor

__fastcall TVit211::~TVit211()
{
  // La deallocazione dovrebbe essere gestita esternamente tramite delete_viterbi211
}

//------------------------------------------------------------------------------
/* Create 256-entry odd-parity lookup table */
void TVit211::partab_init(void)
{
  int i,cnt,ti;

  if(P_init) return;

  for(i=0;i<256;i++)
  {
	cnt = 0;
	ti = i;
	while(ti)
	{
	  if(ti & 1) cnt++;
	  ti >>= 1;
	}
	Partab[i] = cnt & 1;
  }
  P_init=true;
}
//----------------------------------------------------------------------------------------

inline int TVit211::parity(int x)
{
  // Questa funzione assume un int di almeno 16 bit
  unsigned char y, z;
  y = Partab[(x >> 8) & 0xff];
  z = Partab[x & 0xff];
  return y ^ z;
}

//----------------------------------------------------------------------------------------
// Pre-calcola le uscite per ogni transizione di stato
void TVit211::set_viterbi211_polynomial(int polys[2])
{
  int state;

  // Per K=11, abbiamo 1024 stati, quindi 512 coppie di stati da calcolare
  for(state=0; state < 512; state++)
  {
	Branchtab211[0].c[state] = (polys[0] < 0) ^ parity((2*state) & abs(polys[0])) ? 255 : 0;
	Branchtab211[1].c[state] = (polys[1] < 0) ^ parity((2*state) & abs(polys[1])) ? 255 : 0;
  }
  Init = true;
}

//----------------------------------------------------------------------------------------
/* Inizializza il decodificatore Viterbi per un nuovo frame */
int TVit211::init_viterbi211(void *p,int starting_state)
{
  struct v211 *vp = (v211 *)p;
  int i;

  if(p == NULL)
	return -1;
  // Per K=11, ci sono 1024 stati
  for(i=0;i<1024;i++)
	vp->metrics1.w[i] = 63;

  vp->old_metrics = &vp->metrics1;
  vp->new_metrics = &vp->metrics2;
  vp->dp = vp->decisions;
  // Azzera la metrica per lo stato di partenza, mascherando per 1024 stati
  vp->old_metrics->w[starting_state & 1023] = 0;
  return 0;
}

//----------------------------------------------------------------------------------------
/* Crea una nuova istanza del decodificatore Viterbi */
void *TVit211::create_viterbi211(int len)
{
  struct v211 *vp;

  if(!Init)
  {
	int polys[2] = { V211POLYA, V211POLYB };
	set_viterbi211_polynomial(polys);
  }
  if((vp = (v211 *)malloc(sizeof(struct v211))) == NULL)
	 return NULL;

  // Alloca memoria per le decisioni, K-1=10
  if((vp->decisions = (decision211_t *)malloc((len+10)*sizeof(decision211_t))) == NULL)
  {
	free(vp);
	return NULL;
  }
  init_viterbi211(vp,0);

  return vp;
}

//----------------------------------------------------------------------------------------
/* Codifica un frame di dati */
unsigned char * TVit211::encode211(int indata[], int framebits, bool init, bool flush)
{
	unsigned int i, flushlen=0;
	static unsigned int sr=0; // Usiamo unsigned int per contenere 10 bit

	if(init)
	  sr = 0;

	if(flush)
	  flushlen = 10; // K-1 per K=11

	for(i=0;i<framebits+flushlen;i++)
	{
	  int bit = (i < framebits) ? indata[i] : 0;
	  sr = (sr << 1) | bit;
	  // Non è necessario mascherare sr perché i polinomi lo fanno implicitamente
	  symbols[2*i+0] = parity(sr & V211POLYA);
	  symbols[2*i+1] = parity(sr & V211POLYB);
	}
	return symbols;
}
//----------------------------------------------------------------------------------------
/* Viterbi chainback */
int TVit211::chainback_viterbi211(void *p, unsigned char *data, unsigned int nbits, unsigned int endstate)
{
  struct v211 *vp = (v211 *)p;
  decision211_t *d;

  if(p == NULL) return -1;

  d = vp->decisions;
  // Maschera lo stato finale per K=11 (1024 stati)
  endstate %= 1024;

  // Aggiunge un offset per guardare oltre la "coda" di flushing
  d += 10; // K-1 = 10

  while(nbits-- != 0)
  {
	int k;
    // Estrae il bit di decisione per lo stato corrente
    k = (d[nbits].w[endstate/32] >> (endstate%32)) & 1;
    // Aggiorna lo stato precedente e salva il bit decodificato
    data[nbits>>3] = endstate = (endstate >> 1) | (k << 9); // k << (K-2)
  }
  return 0;
}
//----------------------------------------------------------------------------------------
/* Cancella un'istanza del decodificatore Viterbi */
void TVit211::delete_viterbi211(void *p)
{
  struct v211 *vp = (struct v211 *) p;

  if(vp != NULL)
  {
    free(vp->decisions);
    free(vp);
  }
  Init = false;
  p = NULL;
}

//----------------------------------------------------------------------------------------
/* Macro "butterfly" per il calcolo di Viterbi, aggiornata per K=11 */
#define BFLY(i) \
{\
	unsigned int metric,m0,m1,decision;\
\
	metric      = (Branchtab211[0].c[i] ^ sym0) + (Branchtab211[1].c[i] ^ sym1);\
	m0          = vp->old_metrics->w[i] + metric;\
	m1          = vp->old_metrics->w[i+512] + (510 - metric); /* i + num_states/2 */\
	decision    = (signed int)(m0-m1) > 0;\
	vp->new_metrics->w[2*i] = decision ? m1 : m0;\
	d->w[i/16] |= decision << ((2*i)&31);\
	m0         -= (metric+metric-510);\
	m1         += (metric+metric-510);\
	decision    = (signed int)(m0-m1) > 0;\
	vp->new_metrics->w[2*i+1] = decision ? m1 : m0;\
    d->w[i/16] |= decision << ((2*i+1)&31);\
}

//----------------------------------------------------------------------------------------
/* Aggiorna il decodificatore con un blocco di simboli demodulati */
int TVit211::update_viterbi211_blk(void *p,unsigned char *syms,int nbits)
{
  struct v211 *vp = (v211 *)p;
  void *tmp;
  decision211_t *d;

  if(p == NULL)	return -1;

  d = (decision211_t *)vp->dp;
  while(nbits--)
  {
    unsigned char sym0,sym1;
    int i;

    // Azzera i bit di decisione per il passo corrente
    for(i=0;i<32;i++) d->w[i] = 0; // 32 * 32bit = 1024 bit

    sym0 = *syms++;
    sym1 = *syms++;

    // Esegui il calcolo per tutti gli stati
    for(i=0;i<512;i++) // num_states / 2
      BFLY(i);

    d++;
    // Scambia i buffer delle metriche
    tmp = vp->old_metrics;
    vp->old_metrics = vp->new_metrics;
    vp->new_metrics = (metric211_t *)tmp;
  }
  vp->dp = d;
  return 0;
}

//----------------------------------------------------------------------------------------
// Registrazione del componente VCL

namespace Vit211
{
   void __fastcall PACKAGE Register()
   {
	  TComponentClass classes[1] = {__classid(TVit211)};
	  RegisterComponents("I2PHD", classes, 0);
   }
}

// End of file
//----------------------------------------------------------------------------------------