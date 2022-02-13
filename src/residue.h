//
// Created by parthajit on 10/7/20.
//

#ifndef CPPMET_RESIDUE_H
#define CPPMET_RESIDUE_H

#include "biodefs.h"

struct residue{
    struct atom atom[40];
    int size;
};

void resi_init(struct residue* self);

void resi_add(struct residue* self, struct atom* atm); 

//Residue::Residue() {
//    size = 0;
//}


void resi_copy(struct residue* dst, struct residue *src);


#endif //CPPMET_RESIDUE_H
