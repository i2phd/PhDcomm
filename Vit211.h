#ifndef Vit211H
#define Vit211H

//------------------------------------------------------------------------------

#include <sysutils.hpp>
#include <controls.hpp>
#include <classes.hpp>

/* r=1/2 k=11 convolutional encoder polynomials (ottimali) */
#define V211POLYA	03345
#define V211POLYB	03613

#define MAXBYTES    4096 // Mantenuto dal codice originale

// Strutture dati aggiornate per K=11 (1024 stati)
typedef union { unsigned int w[1024]; } metric211_t;
typedef union { unsigned long w[32];} decision211_t; // 1024 stati / 32 bits per long = 32
typedef union { unsigned char c[512]; } branchtab211; // 1024 stati / 2 = 512

//------------------------------------------------------------------------------
/* State info for instance of Viterbi decoder */
struct v211
{
  metric211_t metrics1; /* path metric buffer 1 */
  metric211_t metrics2; /* path metric buffer 2 */
  decision211_t *dp;          /* Pointer to current decision */
  metric211_t *old_metrics,*new_metrics; /* Pointers to path metrics, swapped on every bit */
  decision211_t *decisions;   /* Beginning of decisions for block */
};

//------------------------------------------------------------------------------
class TVit211;

class PACKAGE TVit211: public TComponent
{

   public:
	  virtual __fastcall TVit211(TComponent *aOwner);
	  virtual __fastcall ~TVit211();

	  unsigned char *encode211(int indata[], int framebits, bool init, bool flush);
	  void *create_viterbi211(int len);
	  void set_viterbi211_polynomial(int polys[2]);
	  int  init_viterbi211(void *vp,int starting_state);
	  int  update_viterbi211_blk(void *vp,unsigned char sym[],int npairs);
	  int  chainback_viterbi211(void *vp, unsigned char *data,unsigned int nbits,
							   unsigned int endstate);
	  void delete_viterbi211(void *vp);


   private:
	  bool Init, P_init;
	  unsigned char Partab[256];
	  unsigned char symbols[8*2*(MAXBYTES+10)]; // Aumentato per flush K-1=10

	  void partab_init(void);
	  inline int parity(int x);

	  #pragma pack(16)
	  branchtab211 Branchtab211[2];

   __published:
};

//------------------------------------------------------------------------------

#endif