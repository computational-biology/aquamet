//
// Created by parthajit on 11/7/20.
//

#ifndef CPPMET_RNABP_H
#define CPPMET_RNABP_H

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "biodefs.h"


struct basepair {
      int corid;
      int cifid;
      char resname[4];
      char chain[4];
      char ins[5];
      char name[5];
      char type[3];
      double eval;
      struct basepair* bp[4];
      int numbp;

};

void bp_free(struct basepair* self);
void bp_fprint_short(struct basepair* self, FILE* fp, char is_target, struct atom* met, struct atom* water, char sec_seq, char* loc);



void bp_fprint(struct basepair* self, FILE* fp);
void bp_fprint_pymol(FILE* fp, struct basepair* self);




struct rnabp{ 
      struct basepair* bp;
      int nres;
      int atm;
      int hetatm;
};

void rnabp_fprint(struct rnabp* self, FILE* fp);


void rnabp_free(struct rnabp* self);

void rnabp_scanf(struct rnabp* self, char* outfile);


//void rnabp_free(struct rnabp* self) {
//    self->bp = NULL;
//}


#endif //CPPMET_RNABP_H
