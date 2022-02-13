//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_METALDEFS_H
#define CPPMET_METALDEFS_H

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "geom3d.h"


#define NUM_METAL 128




struct parameters{
    int n;
    char name[NUM_METAL][5];
    int atomic_no[NUM_METAL];
    int coord_no[NUM_METAL];
    double hb_dst[NUM_METAL];
    double vdw_rad[NUM_METAL];
    double ion_rad[NUM_METAL];
    double o_dst[NUM_METAL];
    double n_dst[NUM_METAL];
    double hoh_dst[NUM_METAL];
    double c_dst[NUM_METAL];
    double p_dst[NUM_METAL];
    double s_dst[NUM_METAL];  // sulpher dist
    double met_dst[NUM_METAL];
};
    void parameters_init(struct parameters* self);

void parameters_create(struct parameters* self, const char *param_file) ;
int is_metal(char* metal);

int is_metal1(char* metal);

int get_atomic_no(char* metal);

double calc_angle_energy(char* metal, double angle);


struct current_params{
    char name[5];
    int atomic_no;
    int coord_no;
    double hb_dst;
    double vdw_rad;
    double ion_rad;
    double o_dst2;
    double n_dst2;
    double hoh_dst2;
    double c_dst2;
    double s_dst2;
    double p_dst2;
    double met_dst2;
};
    void currprm_print(struct current_params* self);
    void currprm_adjst(struct current_params* self, struct parameters* prm, int atomic_no, char rule, double adjust);

#endif //CPPMET_METALDEFS_H
