//
// Created by parthajit on 10/7/20.
//

#ifndef CPPMET_MOLECULE_H
#define CPPMET_MOLECULE_H


#include "assert.h"
#include "biodefs.h"
#include "residue.h"
#include "rnabp.h"



struct molecule {

    int max;
    struct residue *residue;
    int size;
    int atm;
   int hetatm;
};
//    void reset(bool atm, bool hetatm);
//    Molecule(char* corfile, rnabp* rnabp);
//    Molecule(int maxsize, bool atm, bool hetatm);
//    void scan_cif(char* ciffile, int (*pf)(char*));
//    Molecule(Molecule* mol);

void mol_init(struct molecule* self);
    int mol_atom_count(struct molecule* self);
    void mol_get_all_atoms(struct molecule* self, struct atom* atmarray[], int* atmcount);

void mol_scan_cif(struct molecule* self, char *ciffile,  int (*pf)(char *)) ;

void mol_polulate(struct molecule* self, struct atom* atoms, int size);

void mol_scan_rna(struct molecule* self, char *corfile, struct rnabp *rnabp) ;


void mol_free(struct molecule* self); 

void mol_copy(struct molecule* dst, struct molecule* src); 

void mol_create(struct molecule* self, int maxsize, int atm, int hetatm) ;
     

void mol_reset(struct molecule* self, int atm, int hetatm);


#endif //CPPMET_MOLECULE_H
