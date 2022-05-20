//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_METALBIND_H
#define CPPMET_METALBIND_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "biodefs.h"
#include "geom3d.h"
#include "residue.h"
#include "molecule.h"
#include "struct.h"
#include "sysparams.h"
#include "metaldefs.h"
#include "rnabp.h"
#include "util.h"
#include "libnuparm.h"

#define MAX_SIZE 90
#define MAX_PROX 400
//class Map;




struct proximity{
    struct molecule* mol[MAX_PROX];
    int      resindx[MAX_PROX];
    char restype[MAX_PROX];
    int size;
};
    void prox_init(struct proximity* self);
    void prox_add(struct proximity* self, struct molecule* mol, int resindx, char type);
    struct residue* prox_get_residue(struct proximity* self, int index);
    struct molecule* prox_get_molecule(struct proximity* self, int index);



struct ligand{
    struct molecule* mol[MAX_SIZE];
    int resindx[MAX_SIZE];
    int offset[MAX_SIZE];
    char loc[MAX_SIZE];      // S for sugar, N for base, P for php etc.
    char restype[MAX_SIZE]; //p for protein n for nuc w for water M for met
    char rule[MAX_SIZE];
    double dist[MAX_SIZE];
    int size;
};
    void ligand_init(struct ligand* self);
    void ligand_add(struct ligand* self, struct proximity* prox, int index, int offset, char loc, char rule, double dist);

struct mediated{
    int waterindx[MAX_SIZE];
    int link_count[MAX_SIZE];
    double angle[MAX_SIZE];
    struct ligand ligand;
    int size;
};
    void mediated_init(struct mediated* self);
    void mediated_add(struct mediated* self, struct proximity* prox,  int index, int offset, char loc, char rule, int water_indx, double dist, double
    angle);

struct site{
    struct atom* metal;
    int atomic_no;

    struct proximity prox;



    struct ligand ligand;
    struct mediated wmed;

};


    void site_init(struct site* self);

    void site_fprint_summary(struct site* self, FILE* fp);
    void site_fprint_wmed_basepair(struct site* self, FILE *fp, struct rnabp *bp, int *flag, int detail_flag, int* bpflag, int allbaseflag);

    void site_fprint_basepair(struct site* self, FILE* fp, struct rnabp* bp, int* flag, int detail_flag, int* bpflag, int allbaseflag);

    void site_fprint_basepair_motifs(struct site* self, FILE* fp, struct rnabp* bp, struct structure* struc, int* flag, int detail_flag);
    void site_fprint_wmed_basepair_motifs(struct site* self, FILE* fp, struct rnabp* bp, struct structure* struc, int* flag, int detail_flag);

    void site_fprint_wmed(struct site* self, FILE* fp);



    void site_fprint_angle(struct site* self, FILE* fp);

    void site_fprint_pml(struct site* self, struct rnabp* rnabp, FILE* fp);
    
    void site_fprint_inligand(struct site* self, struct runparams* runpar);


    void site_add_metal(struct site* self, struct atom* met);

void site_fill_proximity(struct site* self, struct molecule *mol, char moltype) ;



//Map::Map(struct molecule* mol) {
//self->mol = mol;
//int nres = mol->size;
//self->metal1 = (Site**) malloc(nres * sizeof(Site*));
//self->metal2 = (Site**) malloc(nres * sizeof(struct Site*));
//self->metal3 = (Site**) malloc(nres * sizeof(struct Site*));
//for(int i=0; i<nres; ++i){
//    self->metal1[i] = NULL;
//    self->metal2[i] = NULL;
//    self->metal3[i] = NULL;
//}
//}
//
//void Map::add_site(Site* site, Ligand* lig, int ligindx) {
//      assert(self->mol == lig->mol[ligindx]);
//      int resindx = lig->resindx[ligindx];
//      if(self->metal1[resindx] == NULL){
//	    self->metal1[resindx] = site;
//      }else if(self->metal2[resindx] == NULL){
//	    self->metal2[resindx] = site;
//      }else if(self->metal3[resindx] == NULL){
//	    self->metal3[resindx] = site;
//      }else {
//	    ; //fprintf(stderr, "Error... too many sites bind to same molecule\n");
//	    //exit(1);
//      }
//}
//

void site_comp_nucleic(struct site* self, char rule,
                     struct current_params* mprm);
void site_comp_protein(struct site* self, char rule,
                     struct current_params* mprm);

void site_comp_water(struct site* self, char rule, struct current_params *mprm);

void site_comp_water_mediated(struct site* self, char rule, struct current_params *mprm); 
//    void fill_proximity(struct molecule *mol, char moltype);
//    void comp_nucleic(char rule,
//                 Current_params* mprm);
//    void comp_protein(char rule,
//                 Current_params* mprm);
//    void comp_water(char rule, Current_params* mprm);

//    void comp_water_mediated(char rule, Current_params* mprm);

//    void site_populate(struct molecule* rna,
//                       struct molecule* pro,
//                       struct molecule* hoh,
//                       struct molecule* met,
//                       char rule,
//                       Paremeters* prm);



//class Map{
//public:
//    struct molecule* mol;
//    Site** metal1;
//    Site** metal2;
//    Site** metal3;
//public:
//    Map(struct molecule* mol);
//    void add_site(Site* site, Ligand* lig, int ligindx);
//    void gen_verna(char *file){
//        //FILE* fp = fopen(file, "w");
//        //assert(fp != NULL);
//        //fprintf(fp, "<applet  code=\"VARNA.class\" codebase=\"/usr/local/\" archive=\"VARNAv3-93-src.jar\" "
//          //          "height=6000 width=4000>\n");
//        int count = 0;
//        for(int i=0; i< mol->size; ++i){
//            if(self->metal1[i] != NULL){
//                   //fprintf(fp, " WORKING\n");
//                count ++;
//            }
//
//        }
//        //fclose(fp);
//        if(count >0){
//	      ;
//        }
//    }
//};
//


void site_populate(struct site* self,
                   struct molecule* rna,
                   struct molecule* pro,
                   struct molecule* hoh,
                   struct molecule* met,
                   char rule,
                   struct parameters* prm);


void comp_metal_sites(struct molecule* met,
                      struct molecule* hoh,
                      struct molecule* rna,
                      struct rnabp* rnabp,
                      struct nuparm* nuparm,
                      struct molecule* pro,
                      struct parameters* prm,
                      struct structure* sec,
                      char rule,
                      char* file_path,
                      char* file_name,
		      struct runparams* runpar,
		      struct sysparams* syspar); 




#endif //CPPMET_METALBIND_H
