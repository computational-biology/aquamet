//
// Created by parthajit on 10/7/20.
//



#include "residue.h"


void resi_init(struct residue* self){
      self->size = 0;
}

void resi_add(struct residue* self, struct atom* atm) {
    if(self->size>=40){
        fprintf(stderr, "Error in function %s().... too many atoms in residue1 %s for atom id  %d\n",__func__, self->atom[0].resname, self->atom[0].id);
        exit(1);
    }

    self->atom[self->size] = *atm;
    /*self->atom[self->size].type = atm->type;
    self->atom[self->size].center.x = atm->center.x;
    self->atom[self->size].center.y = atm->center.y;
    self->atom[self->size].center.z = atm->center.z;
    self->atom[self->size].id = atm->id;
    self->atom[self->size].resid = atm->resid;
    self->atom[self->size].occu = atm->occu;
    strcpy(self->atom[self->size].chain, atm->chain);
    strcpy(self->atom[self->size].resname, atm->resname);
    strcpy(self->atom[self->size].loc, atm->loc);*/

    self->size ++;
}

//Residue::Residue() {
//    size = 0;
//}


void resi_copy(struct residue* dst, struct residue *src) {
    dst->size = src->size;
    for(long i=0; i<src->size; ++i){
        dst->atom[i] = src->atom[i];
    }
}



