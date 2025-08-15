#ifndef RS_LLVM_H
#define RS_LLVM_H

// From rs_llvm.c
extern int alpha_to[];
extern int index_of[];
extern int gg[];
extern int recd[];
extern int data[];
extern int bb[];

void generate_gf();
void gen_poly();
void encode_rs();
void decode_rs();

#endif
