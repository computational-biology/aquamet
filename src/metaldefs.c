//
//
// Created by parthajit on 13/7/20.
//

#include "metaldefs.h"


char _global_metal[NUM_METAL][5];
int _global_metal_atomic[NUM_METAL];
int _global_metal_size=0;


    void parameters_init(struct parameters* self){
        self->n = 0;
    }

void parameters_fprint(FILE* fp, struct parameters* self){
    fprintf(fp, "       MET   NO     OXY     NITRO    SULPH     MET\n");
    fprintf(fp, "       ---------------------------------------------\n");
    for(int i=0; i<128; ++i){
        if(self->name[i][0] != '-'){
            
			fprintf(fp, "PARAMS %3s  %3d  %6.2lf    %6.2lf   %6.2lf   %6.2lf\n", 
				  	self->name[i],
					self->atomic_no[i],
					self->o_dst[i],
					self->n_dst[i],
					self->s_dst[i],
					self->met_dst[i]);
		}
	}
	fprintf(fp, "\nNote: The values that are -ve are not used for computations.\n\n");
}
void parameters_create_default(struct parameters* self) {

      char strprm[128][1024];
/*                      TAG   SYMB ATM_NO COORD   VDW_RAD   ION_RAD   O-DST  N-DST  H2O-DST  C-DST  S-DST   P-DST    MET-DST */
      
      strcpy(strprm[0], "DIST    LI    3      6      2.25      ?    2.40   2.40   2.40       ?      ?       ?       2.40    \n"); 
      strcpy(strprm[1], "DIST    NA   11      6      2.40      ?    2.90   3.10   3.00       ?      ?       ?       2.90    \n"); 
      strcpy(strprm[2], "DIST    MG   12      6      2.22      ?    2.58   2.70   2.60       ?      ?       ?       2.70    \n"); 
      strcpy(strprm[3], "DIST    AL   13      6      2.19      ?    3.80   3.80   3.80       ?      ?       ?       3.80    \n"); 
      strcpy(strprm[4], "DIST    K    19      6      2.67      ?    3.38   3.50   3.38       ?      ?       ?       3.38    \n"); 
      strcpy(strprm[5], "DIST    CA   20      6      2.43      ?    2.84   2.84   2.84       ?      ?       ?       2.84    \n");
      strcpy(strprm[6], "DIST    SC   21      6      2.37      ?    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[7], "DIST    MN   25      6      2.15      ?    2.19   2.29   2.19       ?    2.64      ?       2.58    \n"); 
      strcpy(strprm[8], "DIST    FE   26      6      2.19      ?    2.04   2.08   2.10       ?    2.28      ?       2.58    \n"); 
      strcpy(strprm[9], "DIST    CO   27      6      2.16      ?    2.10   2.14   2.10       ?    2.26      ?       2.58    \n"); 
      strcpy(strprm[10],"DIST    NI   28      6      2.14      ?    2.07   2.09   2.08       ?    2.46      ?       2.58    \n"); 
      strcpy(strprm[11],"DIST    CU   29      6      2.16      ?    2.12   2.03   2.37       ?    2.33      ?       2.58    \n"); 
      strcpy(strprm[12],"DIST    ZN   30      6      2.09      ?    2.15   2.10   2.09       ?    2.38      ?       2.58    \n"); 
      strcpy(strprm[13],"DIST    SR   38      6      2.54      ?    3.80   3.80   3.80       ?      ?       ?       3.80    \n"); 
      strcpy(strprm[14],"DIST    CD   48      6      2.18      ?    3.25   3.25   3.25       ?      ?       ?       3.25    \n"); 
      strcpy(strprm[15],"DIST    AU   79      6      2.18      ?    2.58   2.70   2.60       ?      ?       ?       2.58    \n");
       strcpy(strprm[16],"#\n");
      
      
      
      
      
     /* strcpy(strprm[0], "DIST    LI    3      6      2.25      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[1], "DIST    BE    4      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[2], "DIST    NA   11      6      3.00      3.00    2.90   3.10   3.00       ?      ?       ?       2.90    \n"); 
      strcpy(strprm[3], "DIST    MG   12      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[4], "DIST    AL   13      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[5], "DIST    K    19      6      3.00      3.00    3.38   3.50   3.38       ?      ?       ?       3.38    \n"); 
      strcpy(strprm[6], "DIST    CA   20      6      3.00      3.00    2.37   3.10   2.42       ?      ?       ?       2.90    \n"); 
      strcpy(strprm[7], "DIST    SC   21      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[8], "DIST    MN   25      6      3.00      3.00    2.19   2.29   2.19       ?    2.64      ?       2.58    \n"); 
      strcpy(strprm[9], "DIST    FE   26      6      3.00      3.00    2.04   2.08   2.10       ?    2.28      ?       2.58    \n"); 
      strcpy(strprm[10],"DIST    CO   27      6      3.00      3.00    2.10   2.14   2.10       ?    2.26      ?       2.58    \n"); 
      strcpy(strprm[11],"DIST    NI   28      6      3.00      3.00    2.07   2.09   2.08       ?    2.46      ?       2.58    \n"); 
      strcpy(strprm[12],"DIST    CU   29      6      3.00      3.00    2.12   2.03   2.37       ?    2.33      ?       2.58    \n"); 
      strcpy(strprm[13],"DIST    ZN   30      6      3.00      3.00    2.15   2.10   2.09       ?    2.38      ?       2.58    \n"); 
      strcpy(strprm[14],"DIST    GA   31      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[15],"DIST    RB   37      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[16],"DIST    SR   38      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[17],"DIST    CD   48      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[18],"DIST    AU   79      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      strcpy(strprm[19],"DIST    PB   82      6      3.00      3.00    2.58   2.70   2.60       ?      ?       ?       2.58    \n"); 
      */
  
      const char* sep="\t \n";
      char line[1024];
      char* token;
      char* token1;
      for(int i=0; i< NUM_METAL; ++i){
	    strcpy(self->name[i],"-");
      }
      self->n = 0;
      _global_metal_size = 0;
      strcpy(line, strprm[self->n]);
      do{
	    // printf("%s", line);
	    if(line[0] == '#') break;

	    token1 = strtok(line, sep);  // DIST TAG

	    token1 = strtok(NULL, sep);


	    token = strtok(NULL, sep); // atomic_no; // atomic_no is the index;
	    int atomic_no = atoi(token);
	    strcpy(self->name[atomic_no], token1);
	    strcpy(_global_metal[_global_metal_size], token1);
	    _global_metal_atomic[_global_metal_size] = atomic_no;
	    _global_metal_size++;

	    self->atomic_no[atomic_no] = atomic_no;

	    self->hb_dst[atomic_no] = 3.8;

	    token = strtok(NULL, sep);
	    self->coord_no[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atoi(token);

	    token = strtok(NULL, sep);
	    self->vdw_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



	    token = strtok(NULL, sep);
	    self->ion_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




	    token = strtok(NULL, sep);
	    self->o_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



	    token = strtok(NULL, sep);
	    self->n_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);


	    token = strtok(NULL, sep);
	    self->hoh_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



	    token = strtok(NULL, sep);
	    self->c_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);


	    token = strtok(NULL, sep);
	    self->s_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



	    token = strtok(NULL, sep);
	    self->p_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




	    token = strtok(NULL, sep);
	    self->met_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);

	    self->n ++;
	    strcpy(line, strprm[self->n]);
      }while(1);
      
      //parameters_fprint(stdout, self);

}
void parameters_create(struct parameters* self, const char *param_file) {
   FILE* fp = fopen(param_file,"r");
   if(fp == NULL){    /* Exception Handling */ 
	 fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Cannot open metal param from %s. make sure that it exists and readable.\n", __func__, __FILE__, __LINE__, param_file);
	 exit(EXIT_FAILURE);
   }
   const char* sep="\t \n";
   char line[1024];
   char* token;
   char* token1;
   while(fgets(line, 1024, fp) != NULL){
      if(strncmp(line, "DIST", 4) != 0) continue;
         break;
   }
   for(int i=0; i< NUM_METAL; ++i){
      strcpy(self->name[i],"-");
   }
   self->n = 0;
   _global_metal_size = 0;
   do{
      // printf("%s", line);
      if(line[0] == '#') break;
      
      token1 = strtok(line, sep);  // DIST TAG
      
      token1 = strtok(NULL, sep);


      token = strtok(NULL, sep); // atomic_no; // atomic_no is the index;
      int atomic_no = atoi(token);
      strcpy(self->name[atomic_no], token1);
      strcpy(_global_metal[_global_metal_size], token1);
      _global_metal_atomic[_global_metal_size] = atomic_no;
      _global_metal_size++;

      self->atomic_no[atomic_no] = atomic_no;

      self->hb_dst[atomic_no] = 3.8;

      token = strtok(NULL, sep);
      self->coord_no[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atoi(token);

      token = strtok(NULL, sep);
      self->vdw_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      self->ion_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




      token = strtok(NULL, sep);
      self->o_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      self->n_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);


      token = strtok(NULL, sep);
      self->hoh_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      self->c_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);
      
      
      token = strtok(NULL, sep);
      self->s_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      self->p_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




      token = strtok(NULL, sep);
      self->met_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);

      self->n ++;
   }while(fgets(line, 1024, fp) != NULL);
   fclose(fp);
}

int is_metal(char* metal){
    for(int i=0; i<_global_metal_size; ++i){
        if(strcmp(metal, _global_metal[i]) == 0) return 1;
    }
    return 0;
}



int get_atomic_no(char* metal){
    for(int i=0; i<_global_metal_size; ++i){
        if(strcmp(metal, _global_metal[i]) == 0){
            return _global_metal_atomic[i];
        }
    }
    return -1;
}

double calc_angle_energy(char* metal, double angle){
   double ktheta;

   if(strcmp(metal, "MG") == 0){
      ktheta = 50.0;
      if(angle > (3.0*PI)/4.0){
         return ktheta * (PI - angle) * (PI - angle);
      }else{
         return ktheta * ((PI/2.0) - angle) * ((PI/2.0) - angle);
      }
   }
   return -1.0;
}


    void currprm_print(struct current_params* self){
        printf("o2= %lf, n2=%lf, c2=%lf, p2=%lf\n", self->o_dst2,
        self->n_dst2, self->c_dst2, self->p_dst2);
    }
    void currprm_adjst(struct current_params* self, struct parameters* prm, int atomic_no, char rule, double adjust){
        //printf("%d %d\n", prm->atomic_no[atomic_no], atomic_no);
        assert(prm->atomic_no[atomic_no] == atomic_no);
        strcpy(self->name, prm->name[atomic_no]);
        self->atomic_no = prm->atomic_no[atomic_no];
        self->coord_no = prm->coord_no[atomic_no];
        self->hb_dst = prm->hb_dst[atomic_no];
        self->vdw_rad = prm->vdw_rad[atomic_no];
        self->ion_rad = prm->vdw_rad[atomic_no];
        self->o_dst2 = prm->o_dst[atomic_no];
        self->n_dst2 = prm->n_dst[atomic_no];
        self->hoh_dst2 = prm->hoh_dst[atomic_no];
        self->c_dst2 = prm->c_dst[atomic_no];
        self->p_dst2 = prm->p_dst[atomic_no];
        self->s_dst2= prm->s_dst[atomic_no];
        self->met_dst2 = prm->met_dst[atomic_no];
        /*if(rule != '0'){
            self->c_dst2 = -9.9;
            self->p_dst2 = -9.9;
        }*/

        if(self->o_dst2 >= 0.0)
            self->o_dst2 = (self->o_dst2 + adjust) * (self->o_dst2 + adjust);
        if(self->n_dst2 >= 0.0)
            self->n_dst2 = (self->n_dst2 + adjust) * (self->n_dst2 + adjust);
        if(self->hoh_dst2 >= 0.0)
            self->hoh_dst2 = (self->hoh_dst2 + adjust) * (self->hoh_dst2 + adjust);
        if(self->c_dst2 >= 0.0)
            self->c_dst2 = (self->c_dst2 + adjust) * (self->c_dst2 + adjust);
        if(self->p_dst2 >= 0.0)
            self->p_dst2 = (self->p_dst2 + adjust) * (self->p_dst2 + adjust);
        if(self->s_dst2 >= 0.0)
            self->s_dst2 = (self->s_dst2 + adjust) * (self->s_dst2 + adjust);
        if(self->met_dst2 >= 0.0)
            self->met_dst2 = (self->met_dst2 + adjust) * (self->met_dst2 + adjust);
    }




