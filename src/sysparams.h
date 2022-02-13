//
// Created by parthajit on 23/2/20.
//

#ifndef CPPMET_SYSPARAMS_H
#define CPPMET_SYSPARAMS_H
#include <stdio.h>
#include <string.h>

struct runparams{
      FILE* summaryfp;
      FILE* metalfp;
      FILE* metdetailfp;
      FILE* hohfp;
      FILE* hohdetailfp;
      int detailflag;
      int allbaseflag;
};



struct sysparams{

	    char cifparam[512];
	    char accnparam[512];
	    char htparam[512];
	    char hdparam[512];
	    char hdvalparam[512];
	    char chainparam[512];
	    char chainvalparam[512];
	    char angparam[512];
	    char angvalparam[512];
	    char chparam[512];
	    char sgparam[512];
	    char evaltypeparam[512];
	    char corparam[50];
	    char nmrparam[50];
        char nmrvalparam[50];
	    char type[10];
	    int cleaned_res;
	    char file_dir[512];
	    char accn[100];
	    char ext[20];
};
	    void sysparams_init(struct sysparams* self);
	    void syspar_print_params(struct sysparams* self, FILE* fp);




#endif //CPPMET_SYSPARAMS_H
