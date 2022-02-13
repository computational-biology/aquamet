//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_SECSTRUCT_H
#define CPPMET_SECSTRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>



struct structure {

    char* secseq;
    char* dbn;
    char* priseq;
    char chain[200][6];
    int num_chain;
    int num_res;
    };

    void structure_init(struct structure* self, char* dat_file, int size);
    void structure_free(struct structure* self);





#endif //CPPMET_SECSTRUCT_H
