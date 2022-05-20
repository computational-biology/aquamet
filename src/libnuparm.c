/*
 * =====================================================================================
 *
 *       Filename:  libnuparm.c
 *
 *    Description:  Nuparm library.
 *
 *        Version:  1.0
 *        Created:  Thursday 10 February 2022 11:10:17  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */




#include "libnuparm.h"



void nuparm_init(struct nuparm* self, int size, const struct nuparm* opt){
      self->flag.is_orient = opt->flag.is_orient;
      self->flag.is_step   = opt->flag.is_step;
      self->size = size;
      if(self->flag.is_orient == TRUE){
	    self->orient	= (struct bporient*) malloc ( size * sizeof(struct bporient) );
	    if ( self->orient==NULL ) {
		  fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
		  exit (EXIT_FAILURE);
	    }
	    for(int i=0; i<size; ++i){
		  self->orient[i].res1 = -999;
	    }
      }

      if(self->flag.is_step == TRUE){
	    self->step = (struct bpstep*) malloc ( size * sizeof(struct bpstep) );
	    if ( self->step==NULL ) {
		  fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
		  exit (EXIT_FAILURE);
	    }
	    for(int i=0; i<size; ++i){
		  self->step[i].res1 = -999;
	    }
      }

}

static void nuparm_read_orient(struct nuparm* self, FILE* fp)
{
      if(self->flag.is_orient == FALSE){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... \n", __func__);
	    exit(EXIT_FAILURE);
      }

      
      char line[PRM_LINE_MAX];
      char sep[] = "\t \n";
      char* token;

      int flag = 0;
      int i =0;
      while( fgets(line, PRM_LINE_MAX, fp) != NULL ){
	    if(strncmp(line, "BL ", 3) != 0 ){
		  if(flag == 0){
			continue;
		  }else{
			break;
		  }
	    }

	    //printf("The Line: %s\n", line);

	    token = strtok(line, sep);

	    token = strtok(NULL, sep);
	    int res1 = atoi(token);
	    int indx = res1 -1;
	    self->orient[indx].res1 = res1;

	    token = strtok(NULL, sep);
	    
	    token = strtok(NULL, sep);
	    self->orient[indx].res2 = atoi(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].buckle = atof(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].open = atof(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].propel = atof(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].stagger = atof(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].shear = atof(token);

	    token = strtok(NULL, sep);
	    self->orient[indx].stretch = atof(token);
	    
	    i++;
      }

}

void nuparm_read(struct nuparm* self, char* prmfile){

      FILE	*fp;										/* input-file pointer */
      char	*fp_file_name = prmfile;		/* input-file name    */

      fp	= fopen( fp_file_name, "r" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			fp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      

      nuparm_read_orient(self, fp);

      if( fclose(fp) == EOF ) {			/* close input file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			fp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }

}




void nuparm_free(struct nuparm* self){
      if(self->flag.is_orient == TRUE){
	    free ( self->orient );
	    self->orient = NULL;
      }
      if(self->flag.is_step == TRUE){
	    free ( self->step );
	    self->step = NULL;
      }

}





