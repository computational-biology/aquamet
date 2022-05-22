//
// Created by parthajit on 13/7/20.
//

#include "metalbind.h"

void prox_init(struct proximity* self){
      self->size = 0;
}
void prox_add(struct proximity* self, struct molecule* mol, int resindx, char type){
      assert(self->size < MAX_PROX);
      self->mol[self->size] = mol;
      self->resindx[self->size] = resindx;
      self->restype[self->size] = type;
      self->size = self->size + 1;
}
struct residue* prox_get_residue(struct proximity* self, int index){
      return self->mol[index]->residue+self->resindx[index];
}
struct molecule* prox_get_molecule(struct proximity* self, int index){
      return self->mol[index];
}


void ligand_init(struct ligand* self){
      self->size = 0;
}
void ligand_add(struct ligand* self, struct proximity* prox, int index, int offset, char loc, char rule, double dist){
      assert(self->size < MAX_SIZE);
      self->mol[self->size] = prox->mol[index];
      self->resindx[self->size] = prox->resindx[index];
      self->loc[self->size] = loc;
      self->restype[self->size] = prox->restype[index];
      self->rule[self->size] = rule;
      self->offset[self->size] = offset;
      self->dist[self->size] = dist;
      self->size = self->size + 1;
}


void mediated_init(struct mediated* self){
      self->size = 0;
      ligand_init(&(self->ligand));
      for(int i=0; i<MAX_SIZE; ++i){
	    self->link_count[i] = 0;
      }
}
void mediated_add(struct mediated* self, struct proximity* prox,  int index, int offset, char loc, char rule, int water_indx, double dist, double
	    angle){
      assert(self->size < MAX_SIZE);
      self->waterindx[self->size] = water_indx;
      self->angle[self->size] = angle;


      ligand_add(&(self->ligand), prox, index, offset, loc, rule, dist);
      int max_link = self->link_count[self->size];
      self->size =self->size + 1;
      for(int i=0; i<self->size; ++i){
	    if(self->waterindx[i] == water_indx){
		  if(self->link_count[i] > max_link){
			max_link = self->link_count[i];
		  }
	    }
      }
      for(int i=0; i< self->size; ++i){
	    if(self->waterindx[i] == water_indx){
		  self->link_count[i] = max_link +1;
	    }
      }
}



void site_init(struct site* self) {
      self->metal = NULL;
      self->atomic_no = -999;
      //        self->prox = Proximity();
      prox_init(&(self->prox));
      //        self->ligand = Ligand();
      ligand_init(&(self->ligand));
      mediated_init(&(self->wmed));
}

void site_fprint_summary(struct site* self, FILE* fp){
      int pcount = 0;
      int ncount = 0;
      int wcount = 0;
      int mcount = 0;
      int wmpcount = 0;
      int wmncount = 0;
      int wmwcount = 0;
      int wmmcount = 0;
      for(int i=0; i<self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N') ncount++;
	    else if(self->ligand.restype[i] == 'P') pcount++;
	    else if(self->ligand.restype[i] == 'W') wcount ++;
	    else if(self->ligand.restype[i] == 'M') mcount ++;
      }
      for(int i=0; i<self->wmed.size; ++i){
	    if(self->wmed.ligand.restype[i] == 'N') wmncount++;
	    else if(self->wmed.ligand.restype[i] == 'P') wmpcount++;
	    else if(self->wmed.ligand.restype[i] == 'W') wmwcount ++;
	    else if(self->wmed.ligand.restype[i] == 'M') wmmcount ++;
	    else{
		  fprintf(stderr,"Error...  Wrong type %c\n", self->wmed.ligand.restype[i]);
		  exit(EXIT_FAILURE);
	    }
      }
      fprintf(fp, "SUMM  %6d  %-3s  %-3s  %3d =  %3d  %3d  %3d  %3d    %3d  =  %3d  %3d  %3d  %3d\n",
		  self->metal->resid,
		  self->metal->chain,
		  self->metal->resname,
		  self->ligand.size, ncount, pcount, wcount, mcount,
		  self->wmed.size,
		  wmncount,
		  wmpcount,
		  wmwcount,
		  wmmcount
		  );
		  
		  fprintf(fp, "SUMM  %6d  %-3s  %-3s  %3d =  %3d  %3d  %3d  %3d\n",
		  self->metal->resid,
		  self->metal->chain,
		  self->metal->resname,
		  self->ligand.size, ncount, pcount, wcount, mcount
		  );
}

void site_fprint_wmed_bporient(struct site* self, FILE *fp, struct rnabp *bp, struct nuparm* nuparm, int *flag, int detail_flag, int* bpflag,
	    int allbaseflag){
      if(bp == NULL) return;

      //int resindx1;
      //int link = 80;



      /*for(int i=0; i<ligand.size; ++i){
	if(ligand.restype[i] == 'W'){
	resindx1 = ligand.resindx[i];
	if(bp->bp[resindx1].numbp > 0){
	link++;
	}
	}
	}*/
      int tag  = 0;
      int previndex = -999;
      for(int i=0; i<self->wmed.size; ++i){
	    int windx = self->wmed.waterindx[i];

	    struct atom* water = self->ligand.mol[windx]->residue[self->ligand.resindx[windx]].atom+ self->ligand.offset[windx];
	    if(self->wmed.ligand.restype[i] != 'N') continue;

	    struct molecule* mol = self->wmed.ligand.mol[i];
	    int resindx = self->wmed.ligand.resindx[i];
	    int offset  = self->wmed.ligand.offset[i];
	    struct atom* atm = mol->residue[resindx].atom + offset;
	    if(atm->resid != bp->bp[resindx].cifid){
		  fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
		  exit(EXIT_FAILURE);
	    }
	    if(allbaseflag == 0){
		  if(bp->bp[resindx].numbp <= 0) continue;
	    }
	    char location[5] = "---";
	    if(self->wmed.ligand.loc[i] == 'N'){
		  strcpy(location,"NUC");
	    }else if(self->wmed.ligand.loc[i] == 'P'){
		  strcpy(location, "PHP");
	    }else if(self->wmed.ligand.loc[i] == 'S'){
		  strcpy(location, "SUG");
	    }
	    else{
		  fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    if(detail_flag == 0){
		  if(strcmp(location, "PHP") ==0) continue;
		  if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;

	    }

	    if(*flag == 0){
		  *flag = 1;
		    fprintf(fp, "  |     Metal    | Water      |      Base Pair Details        |                 outcome             |\n");
		    fprintf(fp, "  |    Detail    |  Details   |\n");
  fprintf(fp, "  ---------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "      id ch  mtl   id  ch  res lc res atm  id ch    id chn   base-pair E-Val   M-W   W-N   ang\n");
  fprintf(fp, "  ---------------------------------------------------------------------------------------------------\n");
//		  fprintf(fp, "\n\n      |  Metal Detail   | Water Details   |      Base Pair Details                  |"
//			      "  outcome      "
//			      "    |\n");
//		  fprintf(fp,      "      "
//			      "---------------------------------------------------------------------------------------------------\n");
//		  fprintf(fp, "       resid  chn  mtl    resid  chn  res  loc  res  atm   resid  chn   resid  chn     pair------>> "
//			      "\n");
//		  fprintf(fp, "      "
//			      "---------------------------------------------------------------------------------------------------\n");


	    }
	    tag = 1;
	    *bpflag = 1;
	    if(resindx != previndex){
		  previndex = resindx;
		  if(bp->bp[resindx].numbp > 0){
		  //fprintf(fp, "ST %5d %-3s %-3s %5d %-3s %-3s %c %3s  %-3s",
		  fprintf(fp, "ST %-3s -- %-3s %3s (%-3s) ",
			      //self->metal->resid,
			      //self->metal->chain,
			      self->metal->resname,
			      //water->resid,
			      //water->chain,
			      water->resname,
			      //location[0],
			      atm->resname,
			      atm->loc);

		  if(allbaseflag == 0){
			bp_fprint(bp->bp+resindx,fp);
			
		  }else{
			if(bp->bp[resindx].numbp <= 0){
			      ;//fprintf(fp, "    ~~  ~        ~~  ~       ~~~~~~~~    ~~~~  ");
			}else{
			      bp_fprint(bp->bp+resindx,fp);
			}
		  }
		  fprintf(fp, "%6.3lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf", 
			      nuparm->orient[resindx].buckle,
			      nuparm->orient[resindx].open,
			      nuparm->orient[resindx].propel,
			      nuparm->orient[resindx].stagger,
			      nuparm->orient[resindx].shear,
			      nuparm->orient[resindx].stretch);
		  fprintf(fp, "\n");
		  }
	    }
	    //fprintf(fp,"\n");

      }

      //fprintf(fp, "WMBP  %6d  %-3s  %-3s\n",
      //        water->resid,
      //        water->chain,
      //        water->resname);
      if(tag == 1)
	    fprintf(fp,"\n");

}
void site_fprint_wmed_basepair(struct site* self, FILE *fp, struct rnabp *bp, int *flag, int detail_flag, int* bpflag,
	    int allbaseflag){
      if(bp == NULL) return;

      //int resindx1;
      //int link = 80;



      /*for(int i=0; i<ligand.size; ++i){
	if(ligand.restype[i] == 'W'){
	resindx1 = ligand.resindx[i];
	if(bp->bp[resindx1].numbp > 0){
	link++;
	}
	}
	}*/
      int tag  = 0;
      for(int i=0; i< self->wmed.size; ++i){
	    int windx = self->wmed.waterindx[i];

	    struct atom* water = self->ligand.mol[windx]->residue[self->ligand.resindx[windx]].atom+ self->ligand.offset[windx];
	    if(self->wmed.ligand.restype[i] != 'N') continue;

	    struct molecule* mol = self->wmed.ligand.mol[i];
	    int resindx = self->wmed.ligand.resindx[i];
	    int offset  = self->wmed.ligand.offset[i];
	    struct atom* atm = mol->residue[resindx].atom + offset;
	    if(atm->resid != bp->bp[resindx].cifid){
		  fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
		  exit(EXIT_FAILURE);
	    }
	    if(allbaseflag == 0){
		  if(bp->bp[resindx].numbp <= 0) continue;
	    }
	    char location[5] = "---";
	    if(self->wmed.ligand.loc[i] == 'N'){
		  strcpy(location,"NUC");
	    }else if(self->wmed.ligand.loc[i] == 'P'){
		  strcpy(location, "PHP");
	    }else if(self->wmed.ligand.loc[i] == 'S'){
		  strcpy(location, "SUG");
	    }
	    else{
		  fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    if(detail_flag == 0){
		  if(strcmp(location, "PHP") ==0) continue;
		  if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;

	    }

	    if(*flag == 0){
		  *flag = 1;
		    fprintf(fp, "  |     Metal    | Water |   Base1      | Atom, Met-H2O & H2O-Atom  |       Base2 pair                  |\n");
		    fprintf(fp, "  |    Detail    | Detail|   Details    | dist & Met-H2O-Atom Algle |       Details                     |\n");
  fprintf(fp, "  -------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "      id ch  mtl  id  ch   res  ch   res atm lc  M-W   W-A   ang     bp      id ch     base-pair Q-Indx  |\n");
  fprintf(fp, "  -------------------------------------------------------------------------------------------------------\n");
//		  fprintf(fp, "\n\n      |  Metal Detail   | Water Details   |      Base Pair Details                  |"
//			      "  outcome      "
//			      "    |\n");
//		  fprintf(fp,      "      "
//			      "---------------------------------------------------------------------------------------------------\n");
//		  fprintf(fp, "       resid  chn  mtl    resid  chn  res  loc  res  atm   resid  chn   resid  chn     pair------>> "
//			      "\n");
//		  fprintf(fp, "      "
//			      "---------------------------------------------------------------------------------------------------\n");


	    }
	    tag = 1;
	    *bpflag = 1;
	    fprintf(fp, "WM %5d %-3s %-3s %5d %-3s  %5d %-3s %3s %-3s %c",
			self->metal->resid,
			self->metal->chain,
			self->metal->resname,
			water->resid,
			water->chain,
//			water->resname,
			
			atm->resid,
			atm->chain,
			atm->resname,
			atm->loc,
			location[0]);
			
			fprintf(fp, "  %4.2lf", dist(self->metal->center, water->center));
		    fprintf(fp, "  %4.2lf", dist(water->center, atm->center));

	        fprintf(fp, "  %6.2lf", todeg(angle3d(self->metal->center, water->center, atm->center)));
	          
		     if(bp->bp[resindx].numbp <= 0){
		       fprintf(fp, "  NP\n");
		       continue;
	        }else{
		       fprintf(fp, "  BP  ");
		    }

	    if(allbaseflag == 0){
		  bp_fprint(bp->bp+resindx,fp);
	    }else{
		  if(bp->bp[resindx].numbp <= 0){
			fprintf(fp, "    ~~  ~        ~~  ~       ~~~~~~~~    ~~~~  ");
		  }else{
			bp_fprint(bp->bp+resindx,fp);
		  }
	    }
		
	    fprintf(fp, "\n");
	    //fprintf(fp,"\n");

      }

      //fprintf(fp, "WMBP  %6d  %-3s  %-3s\n",
      //        water->resid,
      //        water->chain,
      //        water->resname);
      if(tag == 1)
	    fprintf(fp,"\n");

}

void site_bp_flag(struct site* self, struct rnabp* bp,int detail_flag, int* bpflag,
	    int allbaseflag){
      if(bp == NULL) return;
      int resindx;
      for(int i=0; i< self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N'){
		  resindx = self->ligand.resindx[i];
		  struct molecule* mol = self->ligand.mol[i];
		  struct atom* atm = mol->residue[resindx].atom +self->ligand.offset[i];
		  if(atm->resid != bp->bp[resindx].cifid){
			fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
			exit(EXIT_FAILURE);
		  }
		  if(allbaseflag ==0){
			if(bp->bp[resindx].numbp <= 0) continue;
		  }
		  char location[5] = "---";
		  if(self->ligand.loc[i] == 'N'){
			strcpy(location,"NUC");
		  }else if(self->ligand.loc[i] == 'P'){
			strcpy(location, "PHP");
		  }else if(self->ligand.loc[i] == 'S'){
			strcpy(location, "SUG");
		  }
		  else{
			fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
			exit(EXIT_FAILURE);
		  }
		  if(detail_flag == 0){
			if(strcmp(location, "PHP") ==0) continue;
			if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;

		  }


		  *bpflag = 1;
	    }
      }
}
void site_fprint_basepair(struct site* self, FILE* fp, struct rnabp* bp, int* flag, int detail_flag, int* bpflag,
	    int allbaseflag){
      if(bp == NULL) return;
      int resindx;
      int link = 0;
      for(int i=0; i<self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N'){
		  resindx = self->ligand.resindx[i];
		  if(bp->bp[resindx].numbp > 0){
			link++;
		  }
	    }
      }
      for(int i=0; i< self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N'){
		  resindx = self->ligand.resindx[i];
		  struct molecule* mol = self->ligand.mol[i];
		  struct atom* atm = mol->residue[resindx].atom +self->ligand.offset[i];
		  if(atm->resid != bp->bp[resindx].cifid){
			fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
			exit(EXIT_FAILURE);
		  }
		  if(allbaseflag ==0){
			if(bp->bp[resindx].numbp <= 0) continue;
		  }
		  char location[5] = "---";
		  if(self->ligand.loc[i] == 'N'){
			strcpy(location,"NUC");
		  }else if(self->ligand.loc[i] == 'P'){
			strcpy(location, "PHP");
		  }else if(self->ligand.loc[i] == 'S'){
			strcpy(location, "SUG");
		  }
		  else{
			fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
			exit(EXIT_FAILURE);
		  }
		  if(detail_flag == 0){
			if(strcmp(location, "PHP") ==0) continue;
			if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;

		  }

		  if(*flag == 0){
			*flag = 1;
			fprintf(fp, "\n\n      |  Metal Detail       |   Base Pair Details                    |   outcome      "
				    "         "
				    "          |\n");
			fprintf(fp,      "      "
				    "---------------------------------------------------------------------------------------------------\n");
			fprintf(fp, "       resid  chn  mtl  lnk  loc  res  atm   resid  chn   resid  chn     pair------>> \n");
			fprintf(fp, "      "
				    "---------------------------------------------------------------------------------------------------\n");


		  }


		  *bpflag = 1;
		  fprintf(fp, "BP    %6d  %-3s  %-3s  %3d  %-3s  %3s  %-3s  ",
			      self->metal->resid,
			      self->metal->chain,
			      self->metal->resname,
			      link,
			      location,
			      atm->resname,
			      atm->loc);

		  if(allbaseflag == 0){
			bp_fprint(bp->bp+resindx,fp);
		  }else{
			if(bp->bp[resindx].numbp <= 0){
			      fprintf(fp, "    ~~  ~        ~~  ~       ~~~~~~~~    ~~~~  ");

			}else{
			      bp_fprint(bp->bp+resindx,fp);
			}
		  }
		  fprintf(fp, "  %6.3lf", dist(self->metal->center, atm->center));
		  fprintf(fp, "\n");
		  // bp->bp[resindx-1].fprint_bp_short(fp,'F', self->metal->resname);
		  //bp->bp[resindx].fprint_bp_short(fp, 'T', self->metal->resname);
		  //bp->bp[resindx+1].fprint_bp_short(fp, 'F', self->metal->resname);
		  //fprintf(fp,"\n\n");

	    }
      }
      //if(*flag == 1){
      //    fprintf(fp, "#");
      //}
      //fprintf(fp, "\n");
}

void site_fprint_basepair_motifs(struct site* self, FILE* fp, struct rnabp* bp, struct structure* struc, int* flag, int detail_flag){
      if(bp == NULL) return;
      int resindx;
      int link = 0;
      for(int i=0; i<self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N'){
		  resindx = self->ligand.resindx[i];
		  if(bp->bp[resindx].numbp > 0){
			link++;
		  }
	    }
      }
      for(int i=0; i< self->ligand.size; ++i){
	    if(self->ligand.restype[i] == 'N'){
		  resindx = self->ligand.resindx[i];
		  struct molecule* mol = self->ligand.mol[i];
		  struct atom* atm = mol->residue[resindx].atom +self->ligand.offset[i];
		  if(atm->resid != bp->bp[resindx].cifid){
			fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
			exit(EXIT_FAILURE);
		  }
		  //if(bp->bp[resindx].numbp <= 0) continue;
		  char location[5] = "---";
		  if(self->ligand.loc[i] == 'N'){
			strcpy(location,"NUC");
		  }else if(self->ligand.loc[i] == 'P'){
			strcpy(location, "PHP");
		  }else if(self->ligand.loc[i] == 'S'){
			strcpy(location, "SUG");
		  }
		  else{
			fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
			exit(EXIT_FAILURE);
		  }
		  //		if(detail_flag == 0){
		  //		      if(bp->bp[resindx].numbp <= 0) continue;
		  //		      if(strcmp(location, "PHP") ==0) continue;
		  //		      if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;
		  //
		  //		}

		  if(*flag == 0){
			*flag = 1;
			fprintf(fp, "\n\n      |  Metal Detail  |   Base Pair Details                                   "
				    "         "
				    "           |\n");
			fprintf(fp,      "      "
				    "----------------------------------------------------------------------------------------------\n");
			fprintf(fp, "       resid  chn  mtl sec  bp    type   atm         resid \n");
			fprintf(fp, "      "
				    "----------------------------------------------------------------------------------------------\n");


		  }


		  /*fprintf(fp, "BP    %6d  %-3s  %-3s  %3d  %-3s  %-3s  ",
		    metal->resid,
		    metal->chain,
		    metal->resname,
		    link,
		    location,
		    atm->loc);*/

		  //bp->bp[resindx].fprint_bp(fp);
		  if(resindx !=0){
			bp_fprint_short(bp->bp+resindx-1,fp,'F', self->metal, NULL, struc->secseq[resindx-1], atm->loc);
		  }
		  bp_fprint_short(bp->bp+resindx,fp, 'T', self->metal, NULL, struc->secseq[resindx], atm->loc);
		  if(resindx != bp->nres-1){
			bp_fprint_short(bp->bp+resindx+1,fp, 'F', self->metal, NULL, struc->secseq[resindx+1], atm->loc);
		  }
		  fprintf(fp,"\n\n");

	    }
      }
      //if(*flag == 1){
      //    fprintf(fp, "#");
      //}
      //fprintf(fp, "\n");
}
void site_fprint_wmed_basepair_motifs(struct site* self, FILE* fp, struct rnabp* bp, struct structure* struc, int* flag, int detail_flag){
      if(bp == NULL) return;
      int resindx;
      int link = 0;
      for(int i=0; i<self->ligand.size; ++i){
	    if(self->wmed.ligand.restype[i] == 'N'){
		  resindx = self->wmed.ligand.resindx[i];
		  if(bp->bp[resindx].numbp > 0){
			link++;
		  }
	    }
      }
      for(int i=0; i< self->wmed.ligand.size; ++i){
	    if(self->wmed.ligand.restype[i] == 'N'){
		  resindx = self->wmed.ligand.resindx[i];
		  struct molecule* mol = self->wmed.ligand.mol[i];
		  struct atom* atm = mol->residue[resindx].atom +self->wmed.ligand.offset[i];
		  if(atm->resid != bp->bp[resindx].cifid){
			fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
			exit(EXIT_FAILURE);
		  }
		  //if(bp->bp[resindx].numbp <= 0) continue;
		  char location[5] = "---";
		  if(self->wmed.ligand.loc[i] == 'N'){
			strcpy(location,"NUC");
		  }else if(self->wmed.ligand.loc[i] == 'P'){
			strcpy(location, "PHP");
		  }else if(self->wmed.ligand.loc[i] == 'S'){
			strcpy(location, "SUG");
		  }
		  else{
			fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
			exit(EXIT_FAILURE);
		  }

		  if(detail_flag == 0){
			if(bp->bp[resindx].numbp <= 0) continue;
			if(strcmp(location, "PHP") ==0) continue;
			if(strcmp(location, "SUG") == 0 && strcmp(atm->loc,"O2*") != 0) continue;

		  }
		  if(*flag == 0){
			*flag = 1;

			fprintf(fp, "\n\n      |  Metal Detail  |   Base Pair Detsils                                   "
				    "         "
				    "           |\n");
			fprintf(fp,      "      "
				    "----------------------------------------------------------------------------------------------\n");
			fprintf(fp, "       resid  chn  mtl sec  bp    type   atm         resid \n");
			fprintf(fp, "      "
				    "----------------------------------------------------------------------------------------------\n");


		  }

		  /*fprintf(fp, "BP    %6d  %-3s  %-3s  %3d  %-3s  %-3s  ",
		    metal->resid,
		    metal->chain,
		    metal->resname,
		    link,
		    location,
		    atm->loc);*/

		  //bp->bp[resindx].fprint_bp(fp);
		  int waterindx = self->wmed.waterindx[i];
		  struct atom* water = self->ligand.mol[waterindx]->residue[self->ligand.resindx[waterindx]].atom+self->ligand.offset[waterindx];
		  if(resindx != 0){
			bp_fprint_short(bp->bp+resindx-1,fp,'F', self->metal, water,  struc->secseq[resindx-1], atm->loc);
		  }

		  bp_fprint_short(bp->bp+resindx,fp, 'T', self->metal, water, struc->secseq[resindx], atm->loc);
		  if(resindx != bp->nres-1){
			bp_fprint_short(bp->bp+resindx+1,fp, 'F', self->metal,water, struc->secseq[resindx+1], atm->loc);
		  }
		  fprintf(fp,"\n\n");

	    }
      }
      //if(*flag == 1){
      //    fprintf(fp, "#");
      //}
      //fprintf(fp, "\n");
}

void site_fprint_wmed(struct site* self, FILE* fp){
      if(self->wmed.size <= 0) return;
      fprintf(fp, "      |  Metal Detail  |   Water atom detail    |  Mediated atom detail  |      outcome   "
		  "         |\n");
      fprintf(fp,      "      "
		  "----------------------------------------------------------------------------------------------\n");
      fprintf(fp, "       resid  chn  mtl   resid  chn  res  lnk      resid  chn  res  atm     loc    dist    "
		  "angle\n");
      fprintf(fp, "      "
		  "----------------------------------------------------------------------------------------------\n");
      for(int i=0; i<self->wmed.size; ++i){
	    char location[5] = "---";
	    char restype1 = self->wmed.ligand.restype[i];
	    char loc1 = self->wmed.ligand.loc[i];
	    if(restype1 == 'P'){
		  strcpy(location, "PRO");
	    }else if(restype1 == 'W'){
		  strcpy(location, "H2O");
	    }else if(restype1 == 'M'){
		  strcpy(location, "MET");
	    }else if(restype1 == 'N'){
		  if(loc1 == 'N'){
			strcpy(location, "NUC");
		  }else if(loc1 == 'P'){
			strcpy(location, "PHP");
		  }else if(loc1 == 'S'){
			strcpy(location, "SUG");
		  }
		  else{
			fprintf(stderr, "Error in function %s()... invalid ligand type found.\n", __func__);
			exit(EXIT_FAILURE);
		  }
	    }
	    int k = self->wmed.waterindx[i];
	    struct atom* water = self->ligand.mol[k]->residue[self->ligand.resindx[k]].atom + self->ligand.offset[k];
	    assert(strcmp(water->resname, "HOH")==0);
	    struct residue* res1 = self->wmed.ligand.mol[i]->residue+self->wmed.ligand.resindx[i];
	    struct atom* medatm = res1->atom + self->wmed.ligand.offset[i];
	    double dist = sqrt(self->wmed.ligand.dist[i]);
	    double angle = self->wmed.angle[i];
	    int link = self->wmed.link_count[i];

	    fprintf(fp, "WMED  %6d  %-3s  %-3s  %6d  %-3s  %-3s  %3d     %6d  %-3s  %-3s  %-3s     %-3s  %6.3lf   "
			"%6.3lf\n",
			self->metal->resid,
			self->metal->chain,
			self->metal->resname,
			water->resid,
			water->chain,
			water->resname,
			link,
			medatm->resid,
			medatm->chain,
			medatm->resname,
			medatm->loc,
			location,
			dist,
			angle

		   );
      }
      fprintf(fp, "#\n");
}



void site_fprint_angle(struct site* self, FILE* fp){
      if(self->ligand.size < 2) return;
      fprintf(fp, "      |  Metal Detail  |  C-1 atom detail      |   C-2 atom detail      |      outcome   "
		  "         |\n");
      fprintf(fp,      "      "
		  "---------------------------------------------------------------------------------------------\n");
      fprintf(fp, "       resid  chn  mtl   resid  chn  res  atm      resid  chn  res  atm     angle    \n");
      fprintf(fp, "      "
		  "---------------------------------------------------------------------------------------------\n");
      for(int i=0; i<self->ligand.size; ++i){
	    for(int j=i+1; j<self->ligand.size; ++j){
		  struct atom* a = self->ligand.mol[i]->residue[self->ligand.resindx[i]].atom+self->ligand.offset[i];
		  struct atom* c = self->ligand.mol[j]->residue[self->ligand.resindx[j]].atom+self->ligand.offset[j];
		  struct atom* b = self->metal;
		  //                double angdeg = b->angl_deg(a, c);
		  double angrad = angle3d(a->center, b->center, c->center);
		  double angdeg = todeg(angrad);
		  fprintf(fp, "ANGL  %6d  %-3s  %-3s  %6d  %-3s  %-3s  %-3s  :  %6d  %-3s  %-3s  %-3s    %6.2lf\n",
			      self->metal->resid,
			      self->metal->chain,
			      self->metal->resname,
			      a->resid,
			      a->chain,
			      a->resname,
			      a->loc,
			      c->resid,
			      c->chain,
			      c->resname,
			      c->loc,
			      angdeg);

	    }
      }
      fprintf(fp, "#\n");
}

void wmed_fprint_pml(struct site* self, struct rnabp* bp, int ligindex, FILE* fp){
      int i = ligindex;
      char loca[10];
      struct atom* a = self->wmed.ligand.mol[i]->residue[self->wmed.ligand.resindx[i]].atom+ self->wmed.ligand.offset[i];
      int resindx = self->wmed.ligand.resindx[i];
      
      if(self->wmed.ligand.restype[i] == 'N'){
	    if(self->wmed.ligand.loc[i] == 'N'){
		  if(bp->bp[resindx].numbp > 0){ 
			fprintf(fp, "select tmp, ");
			bp_fprint_pymol(fp, bp->bp + resindx);
			fprintf(fp, "color white, tmp\n");
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
//			fprintf(fp, "color red, tmp\n");
			fprintf(fp, "util.cbag tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }else{
			fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			fprintf(fp, "color blue, tmp\n");
//			fprintf(fp, "util.cbay tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }

	    }else{
		  return;
	    }//else if(self->wmed.ligand.loc[i] == 'P'){
//		  fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
//		  fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
//		  fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
//		  fprintf(fp, "color green, tmp\n");
//		  fprintf(fp, "util.cbag tmp\n");
//		  fprintf(fp, "show_as cartoon, tmp\n");
//		  fprintf(fp, "show line, tmp\n");
//	    }else if(self->wmed.ligand.loc[i] == 'S'){
//		  fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
//		  fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
//		  fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
//		  fprintf(fp, "color red, tmp\n");
//		  fprintf(fp, "util.cbao tmp\n");
//		  fprintf(fp, "show_as cartoon, tmp\n");
//		  fprintf(fp, "show line, tmp\n");
//	    }else{
//		  fprintf(stderr,"Error... Wrong type given\n");
//		  exit(EXIT_FAILURE);
//	    }
      }else if(self->wmed.ligand.restype[i] == 'P'){
	    fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
	    fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
	    fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
	    fprintf(fp, "color blue, tmp\n");
	    fprintf(fp, "set cartoon_color, blue, tmp\n");
	    fprintf(fp, "show_as cartoon, tmp\n");
	    fprintf(fp, "show line, tmp\n");
      }else if(self->wmed.ligand.restype[i] == 'W'){
	    strcpy(loca, "H2O");
	    return;
      }else{
	    fprintf(stderr,"Error... Wrong restype given\n");
	    exit(EXIT_FAILURE);
      }

      int k = self->wmed.waterindx[i];
      struct atom* water = self->ligand.mol[k]->residue[self->ligand.resindx[k]].atom + self->ligand.offset[k];
      assert(strcmp(water->resname, "HOH")==0);

      if(self->wmed.ligand.mol[i]->residue[self->wmed.ligand.resindx[i]].size >= 1){
          char loc[5];
	    strcpy(loc, a->loc);
          if(strlen(loc) == 3 && loc[2] == '*'){
	        loc[2] = '\'';
	    }
	    fprintf(fp, "select wmed, (resi %d and chain %s and name %s)\n", a->resid, a->chain, loc);
	    fprintf(fp, "set sphere_scale, 0.40, wmed\n");
	    fprintf(fp, "show spheres,  wmed\n");
	
	
	    fprintf(fp, "distance dwmed,  (resi %d and chain %s),  (resi %d and chain %s and name %s), 15, 4\n",
			water->resid, water->chain, a->resid, a->chain, loc);
	    fprintf(fp, "set dash_width, 1.0, dwmed\n");
	    fprintf(fp, "set dash_gap, 0.2, dwmed\n");
	    fprintf(fp, "set dash_color, red, dwmed\n");
	    fprintf(fp, "hide label, dwmed\n");
      }else{

	    fprintf(fp, "select wmed, (resi %d and chain %s)\n", a->resid, a->chain);
	    fprintf(fp, "color white, wmed\n");
	    fprintf(fp, "set cartoon_ring_mode, 2, wmed\n");
	    fprintf(fp, "set cartoon_ladder_mode, 0, wmed\n");
	    fprintf(fp, "set cartoon_color, white, wmed\n");
	    fprintf(fp, "show_as cartoon, wmed\n");
	    fprintf(fp, "show line, wmed\n");
	    fprintf(fp, "select wmed, (resi %d and chain %s and name %s)\n", a->resid, a->chain, a->loc);
	    fprintf(fp, "set sphere_scale, 0.40, wmed\n");
	    fprintf(fp, "show spheres,  wmed\n");
	    fprintf(fp, "distance dwmed,  (resi %d and chain %s),  (resi %d and chain %s and name %s), 15, 4\n",
			water->resid, water->chain, a->resid, a->chain, a->loc);
	    fprintf(fp, "set dash_width, 1.0, dwmed\n");
	    fprintf(fp, "set dash_gap, 0.2, dwmed\n");
	    fprintf(fp, "set dash_color, red, dwmed\n");
	    fprintf(fp, "hide label, dwmed\n");
      }


}
void ligand_fprint_pml(struct site* self, struct rnabp* bp, int ligindex, FILE* fp){
      int i = ligindex;
      char loca[10];
      struct atom* a = self->ligand.mol[i]->residue[self->ligand.resindx[i]].atom+ self->ligand.offset[i];
      int resindx = self->ligand.resindx[i];
      
      if(self->ligand.restype[i] == 'N'){
	    if(self->ligand.loc[i] == 'N'){
		  if(bp->bp[resindx].numbp > 0){ 
			fprintf(fp, "select tmp, ");
			bp_fprint_pymol(fp, bp->bp + resindx);
			fprintf(fp, "color white, tmp\n");
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			fprintf(fp, "color red, tmp\n");
//			fprintf(fp, "util.cbaw tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }else{
			fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			fprintf(fp, "color yellow, tmp\n");
//			fprintf(fp, "util.cbay tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }

	    }
	    //else if(self->ligand.loc[i] == 'P'){
//		  fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
//		  fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
//		  fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
//		  fprintf(fp, "color green, tmp\n");
//		  fprintf(fp, "util.cbag tmp\n");
//		  fprintf(fp, "show_as cartoon, tmp\n");
//		  fprintf(fp, "show line, tmp\n");
//	    }else if(self->ligand.loc[i] == 'S'){
//		  fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
//		  fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
//		  fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
//		  fprintf(fp, "color red, tmp\n");
//		  fprintf(fp, "util.cbao tmp\n");
//		  fprintf(fp, "show_as cartoon, tmp\n");
//		  fprintf(fp, "show line, tmp\n");
//	    }else{
//		  fprintf(stderr,"Error... Wrong type given\n");
//		  exit(EXIT_FAILURE);
//	    }
      }else if(self->ligand.restype[i] == 'P'){
	    fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
	    fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
	    fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
	    fprintf(fp, "color blue, tmp\n");
	    fprintf(fp, "set cartoon_color, blue, tmp\n");
	    fprintf(fp, "show_as cartoon, tmp\n");
	    fprintf(fp, "show line, tmp\n");
      }else if(self->ligand.restype[i] == 'W'){
	    strcpy(loca, "H2O");
      }else{
	    fprintf(stderr,"Error... Wrong restype given\n");
	    exit(EXIT_FAILURE);
      }
      
      if(self->ligand.mol[i]->residue[self->ligand.resindx[i]].size >= 1){
          char loc[5];
	    strcpy(loc, a->loc);
          if(strlen(loc) == 3 && loc[2] == '*'){
	        loc[2] = '\'';
	    }
	    fprintf(fp, "select ligand, (resi %d and chain %s and name %s)\n", a->resid, a->chain, loc);
	    fprintf(fp, "set sphere_scale, 0.40, ligand\n");
	    fprintf(fp, "show spheres,  ligand\n");
	
	
	    fprintf(fp, "distance d1,  (resi %d and chain %s),  (resi %d and chain %s and name %s), 15, 4\n",
			self->metal->resid, self->metal->chain, a->resid, a->chain, loc);
	    fprintf(fp, "set dash_width, 1.0, d1\n");
	    fprintf(fp, "set dash_gap, 0.2, d1\n");
	    //fprintf(fp, "set dash_color, white, d1\n");
	    fprintf(fp, "hide label, d1\n");
      }else{

	    fprintf(fp, "select ligand, (resi %d and chain %s)\n", a->resid, a->chain);
	    fprintf(fp, "color white, ligand\n");
	    fprintf(fp, "set cartoon_ring_mode, 2, ligand\n");
	    fprintf(fp, "set cartoon_ladder_mode, 0, ligand\n");
	    fprintf(fp, "set cartoon_color, white, ligand\n");
	    fprintf(fp, "show_as cartoon, ligand\n");
	    fprintf(fp, "show line, ligand\n");
	    fprintf(fp, "select ligand, (resi %d and chain %s and name %s)\n",
			a->resid, a->chain, a->loc);
	    fprintf(fp, "set sphere_scale, 0.40, ligand\n");
	    fprintf(fp, "show spheres,  ligand\n");
	    fprintf(fp, "distance d1,  (resi %d and chain %s),  (resi %d and chain %s and name %s), 15, 4\n",  self->metal->resid, self->metal->chain, a->resid, a->chain, a->loc);
	    fprintf(fp, "set dash_width, 1.0, d1\n");
	    fprintf(fp, "set dash_gap, 0.2, d1\n");
	    fprintf(fp, "set dash_color, white, d1\n");
	    fprintf(fp, "hide label, d1\n");
      }


}


void site_fprint_pml_wmed(struct site* self, struct rnabp* bp, FILE* fp)
{
      if(self->ligand.size == 0) return;
      fprintf(fp, "select %s%d%s, ",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "(resi %d and chain %s)\n", self->metal->resid, self->metal->chain);
      fprintf(fp, "set sphere_scale, 0.20, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "show spheres, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      //fprintf(fp, "select %s%d%sligand, ",self->metal->resname, self->metal->resid, self->metal->chain);
      for(int i=0; i<self->ligand.size; ++i){
	    ligand_fprint_pml(self, bp, i, fp);

      }
      for(int i=0; i<self->wmed.size; ++i){
	    wmed_fprint_pml(self, bp, i, fp);
      }

      //fprintf(fp, "TOTAL LIGANDS FOR %3s: %4ld",self->metal->resname, ligand.size);
}

void site_fprint_pml(struct site* self, struct rnabp* bp, FILE* fp)
{
      if(self->ligand.size == 0) return;
      fprintf(fp, "select %s%d%s, ",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "(resi %d and chain %s)\n", self->metal->resid, self->metal->chain);
      fprintf(fp, "set sphere_scale, 0.40, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "show spheres, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      //fprintf(fp, "select %s%d%sligand, ",self->metal->resname, self->metal->resid, self->metal->chain);
      for(int i=0; i<self->ligand.size; ++i){
	    ligand_fprint_pml(self, bp, i, fp);

      }

      //fprintf(fp, "TOTAL LIGANDS FOR %3s: %4ld",self->metal->resname, ligand.size);
}

void site_fprint_pml1(struct site* self, struct rnabp* bp, FILE* fp)
{
      if(self->ligand.size == 0) return;
      fprintf(fp, "select %s%d%s, ",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "(resi %d and chain %s)\n", self->metal->resid, self->metal->chain);
      fprintf(fp, "set sphere_scale, 0.40, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      fprintf(fp, "show spheres, %s%d%s\n",self->metal->resname, self->metal->resid, self->metal->chain);
      //fprintf(fp, "select %s%d%sligand, ",self->metal->resname, self->metal->resid, self->metal->chain);
      for(int i=0; i<self->ligand.size; ++i){
	    struct atom* a = self->ligand.mol[i]->residue[self->ligand.resindx[i]].atom+ self->ligand.offset[i];
	    if(self->ligand.mol[i]->residue[self->ligand.resindx[i]].size == 1){
		  fprintf(fp, "select ligand, (resi %d and chain %s)\n", a->resid, a->chain);
		  fprintf(fp, "set sphere_scale, 0.40, ligand\n");
		  fprintf(fp, "show spheres,  ligand\n");
		  fprintf(fp, "distance d1,  (resi %d and chain %s),  (resi %d and chain %s), 15, 4\n",  self->metal->resid, self->metal->chain, a->resid, a->chain);
		  fprintf(fp, "set dash_width, 1.0, d1\n");
		  fprintf(fp, "set dash_gap, 0.2, d1\n");
		  fprintf(fp, "set dash_color, white, d1\n");
		  fprintf(fp, "hide label, d1\n");
	    }else{
		  fprintf(fp, "select ligand, (resi %d and chain %s)\n", a->resid, a->chain);
		  fprintf(fp, "color white, ligand\n");
		  fprintf(fp, "set cartoon_ring_mode, 2, ligand\n");
		  fprintf(fp, "set cartoon_ladder_mode, 0, ligand\n");
		  fprintf(fp, "set cartoon_color, white, ligand\n");
		  fprintf(fp, "show_as cartoon, ligand\n");
		  fprintf(fp, "show line, ligand\n");
		  fprintf(fp, "distance d1,  (resi %d and chain %s),  (resi %d and chain %s and name %s), 15, 4\n",  self->metal->resid, self->metal->chain, a->resid, a->chain, a->loc);
		  fprintf(fp, "set dash_width, 1.0, d1\n");
		  fprintf(fp, "set dash_gap, 0.2, d1\n");
		  fprintf(fp, "set dash_color, white, d1\n");
		  fprintf(fp, "hide label, d1\n");
	    }
      }
      fprintf(fp, "\n");
      char loca[10];
      for(int i=0; i<self->ligand.size; ++i){
	    int resindx = self->ligand.resindx[i];
	    struct atom* a = self->ligand.mol[i]->residue[self->ligand.resindx[i]].atom+ self->ligand.offset[i];
	    if(self->ligand.restype[i] == 'N'){

		  if(self->ligand.loc[i] == 'N'){
			if(bp->bp[resindx].numbp > 0){ 
			      fprintf(fp, "select tmp, ");
			      bp_fprint_pymol(fp, bp->bp + resindx);
			      fprintf(fp, "color white, tmp\n");
			      fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			      fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			      fprintf(fp, "set cartoon_color, white, tmp\n");
			      fprintf(fp, "show_as cartoon, tmp\n");
			      fprintf(fp, "show line, tmp\n");
			}else{
			      fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
			      fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			      fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			      fprintf(fp, "color yellow, tmp\n");
			      fprintf(fp, "set cartoon_color, yellow, tmp\n");
			      fprintf(fp, "show_as cartoon, tmp\n");
			      fprintf(fp, "show line, tmp\n");
			}

		  }else if(self->ligand.loc[i] == 'P'){
			fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			fprintf(fp, "color green, tmp\n");
			fprintf(fp, "set cartoon_color, green, tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }else if(self->ligand.loc[i] == 'S'){
			fprintf(fp, "select tmp,  (resi %d and chain %s)\n", a->resid, a->chain);
			fprintf(fp, "set cartoon_ring_mode, 2, tmp\n");
			fprintf(fp, "set cartoon_ladder_mode, 0, tmp\n");
			fprintf(fp, "color red, tmp\n");
			fprintf(fp, "set cartoon_color, red, tmp\n");
			fprintf(fp, "show_as cartoon, tmp\n");
			fprintf(fp, "show line, tmp\n");
		  }else{
			fprintf(stderr,"Error... Wrong type given\n");
			exit(EXIT_FAILURE);
		  }
	    }else if(self->ligand.restype[i] == 'P'){
		  strcpy(loca, "PRO");
	    }else if(self->ligand.restype[i] == 'W'){
		  strcpy(loca, "H2O");
	    }else{
		  fprintf(stderr,"Error... Wrong restype given\n");
		  exit(EXIT_FAILURE);
	    }
      }

      //fprintf(fp, "TOTAL LIGANDS FOR %3s: %4ld",self->metal->resname, ligand.size);
}
void site_fprint_inligand(struct site* self, struct runparams* runpar){
      if(self->ligand.size == 0) return;
      FILE* fp = runpar->metdetailfp;
      fprintf(fp,"\n\n");
      fprintf(fp, "        +---------------------S I T E    S U M M A R Y---------------------+\n\n");
      fprintf(fp, "      |  Metal Detail  |      Coordination atom detail    |\n");
      fprintf(fp, "      -----------------------------------------------------\n");
      fprintf(fp, "       resid  chn  mtl   loc  res  atm   resid  chn   dist\n");
      fprintf(fp, "      -----------------------------------------------------\n");
      for(int i=0; i<self->ligand.size; ++i){
	    char loca[6];
	    if(self->ligand.restype[i] == 'N'){
		  if(self->ligand.loc[i] == 'N'){
			strcpy(loca, "NUC");
		  }else if(self->ligand.loc[i] == 'P'){
			strcpy(loca, "PHP");
		  }else if(self->ligand.loc[i] == 'S'){
			strcpy(loca, "SUG");
		  }else{
			fprintf(stderr,"Error... Wrong type given\n");
			exit(EXIT_FAILURE);
		  }
	    }else if(self->ligand.restype[i] == 'P'){
		  strcpy(loca, "PRO");
	    }else if(self->ligand.restype[i] == 'W'){
		  strcpy(loca, "H2O");
	    }else{
		  fprintf(stderr,"Error... Wrong restype given\n");
		  exit(EXIT_FAILURE);
	    }
	    struct atom* a = self->ligand.mol[i]->residue[self->ligand.resindx[i]].atom+ self->ligand.offset[i];
	    fprintf(fp, "BIND  %6d  %-3s  %-3s   %-3s  %-3s  %-3s  %6d  %-3s %6.3lf\n",
			self->metal->resid,
			self->metal->chain,
			self->metal->resname,
			loca,
			a->resname,
			a->loc,
			a->resid,
			a->chain,
			self->ligand.dist[i]);

      }
      fprintf(runpar->metdetailfp,"#\n");

      //fprintf(fp, "TOTAL LIGANDS FOR %3s: %4ld",self->metal->resname, ligand.size);
}


void site_add_metal(struct site* self, struct atom* met){
      self->metal = met;
      self->atomic_no = get_atomic_no(met->resname);
      if(self->atomic_no <0){
	    fprintf(stderr, "Error... not  valid metal: %s\n", met->resname);
      }
}

void site_fill_proximity(struct site* self, struct molecule *mol, char moltype) {
      for(int i=0; i< mol->size; ++i){
	    struct atom* a = mol->residue[i].atom+0;
	    double val = distsqr(self->metal->center, a->center);
	    double radius = 20.0*20.0;
	    if(moltype == 'W'){
		  radius = 8.5 * 8.5;
	    }
	    if(val< radius){  // skip resi which are far than 20A
		  prox_add(&(self->prox),mol, i, moltype);
	    }
      }
}



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
	    struct current_params* mprm){
      struct atom* currmetal = self->metal;

      //mprm->print();
      double dst2;
      for(int resindx=0; resindx<self->prox.size; ++resindx){
	    // if(prox.restype[resindx] != 'N') continue;


	    if(self->prox.restype[resindx] != 'N') continue; // Not necleic acid

	    //int start = rna->resbeg[self->proxres[resindx]];
	    //int end   = rna->resend[self->proxres[resindx]];
	    struct residue* resi = prox_get_residue(&(self->prox),resindx); //  ->residue+resindx;
	    for(int offset=0; offset<resi->size; ++offset){
		  struct atom* curatom = resi->atom+offset;
		  //currmetal->print();
		  //curatom->print();
		  dst2 = distsqr(currmetal->center, curatom->center)+0.00001;
		  //printf("dst=%lf\n", dst2);

		  int len = strlen(curatom->loc);
		  char symb = curatom->loc[0];
		  if(len == 2){ //Then it is Nucleobase
			if((symb == 'O' && dst2 <= mprm->o_dst2)
				    ||( symb == 'N' && dst2 <= mprm->n_dst2)
				    ||( symb == 'C' && dst2 <= mprm->c_dst2)
				    ||( symb == 'S' && dst2 <= mprm->s_dst2)){

			      ligand_add(&(self->ligand), &(self->prox), resindx, offset,
					  'N',
					  rule,
					  sqrt(dst2));
			      //                        map->add_site(this, &self->ligand, ligand.size-1);
			      //printf("Added NUC\n");
			}
		  }else if(len == 3 && (curatom->loc[2] == '*' || curatom->loc[2] == '\'')){ // sugar backbone for cor
			// file
			if((symb == 'O' && dst2 <= mprm->o_dst2)
				    ||( symb == 'C' && dst2 <= mprm->c_dst2)
				    ||( symb == 'S' && dst2 <= mprm->s_dst2)){

			      ligand_add(&(self->ligand),&(self->prox), resindx, offset, 'S', rule, sqrt(dst2));
			      //                        map->add_site(this, &self->ligand, ligand.size-1);
			      //printf("Added SUG\n");
			      //molecule_metal_add_site(map, rna, site, self->proxres[resindx]);
			}
		  }else{ // phosphate case
			if((len == 3 && symb == 'O' && dst2 <= mprm->o_dst2)
				    ||(len == 3 && symb == 'P' && dst2 <= mprm->p_dst2)
				    /*||( symb == 'S' && dst2 <= mprm->s_dst2) */
				    ){

			      ligand_add(&(self->ligand),&(self->prox), resindx, offset, 'P', rule, sqrt(dst2));
			      //                        map->add_site(this, &self->ligand, ligand.size-1);
			      //printf("Added PHP\n");
			      //molecule_metal_add_site(map, rna, site, self->proxres[resindx]);
			}
		  }
	    }
      }
}
void site_comp_protein(struct site* self, char rule,
	    struct current_params* mprm){
      struct atom* currmetal = self->metal;

      //mprm->print();
      double dst2;
      for(int resindx=0; resindx<self->prox.size; ++resindx){
	    // if(prox.restype[resindx] != 'N') continue;


	    if(self->prox.restype[resindx] != 'P') continue; // Not Protein

	    //int start = rna->resbeg[self->proxres[resindx]];
	    //int end   = rna->resend[self->proxres[resindx]];
	    struct residue* resi = prox_get_residue(&(self->prox),resindx); //  ->residue+resindx;
	    for(int offset=0; offset<resi->size; ++offset){
		  struct atom* curatom = resi->atom+offset;
		  //currmetal->print();
		  //curatom->print();
		  dst2 = distsqr(currmetal->center, curatom->center)+0.00001;
		  //printf("dst=%lf\n", dst2);

		  char symb = curatom->loc[0];
		  //if(len == 2){ //Then it is Nucleobase for Cor dataset
		  if((symb == 'O' && dst2 <= mprm->o_dst2)
			      ||( symb == 'N' && dst2 <= mprm->n_dst2)
			      ||( symb == 'C' && dst2 <= mprm->c_dst2)
			      ||( symb == 'S' && dst2 <= mprm->s_dst2)){

			ligand_add(&(self->ligand),&(self->prox), resindx, offset, 'P', rule, sqrt(dst2));
			//                        map->add_site(this, &self->ligand, ligand.size-1);
			//printf("Added NUC\n");
		  }
		  //}
	    }
      }
}

void site_comp_water(struct site* self, char rule, struct current_params *mprm) {
      struct atom* currmetal = self->metal;
      double dst2;
      for(int resindx=0; resindx<self->prox.size; ++resindx){
	    if(self->prox.restype[resindx] != 'W') continue; // Not Water acid
	    struct residue* resi = prox_get_residue(&(self->prox),resindx); //  ->residue+resindx;
	    for(int offset=0; offset<resi->size; ++offset){
		  struct atom* curatom = resi->atom+offset;
		  dst2 = distsqr(currmetal->center, curatom->center)+0.00001;
		  if(strcmp(curatom->resname, "HOH") == 0 && dst2 <= mprm->hoh_dst2){
			ligand_add(&(self->ligand),&(self->prox), resindx, offset, 'W', rule, sqrt(dst2));
		  }
	    }
      }
}

void site_comp_water_mediated(struct site* self, char rule, struct current_params *mprm) {
      double metdist2=0.0, wdist2=0.0, includist2=0.0;
      includist2 = mprm->hb_dst + mprm->o_dst2 + 1.0;
      includist2 *= includist2;
      double hbdst2 = (mprm->hb_dst +0.0001) * (mprm->hb_dst + 0.0001);

      
      
      for(int i=0; i<self->ligand.size; ++i){
	    if(self->ligand.restype[i] != 'W') continue;
	    struct residue* wres = self->ligand.mol[i]->residue + self->ligand.resindx[i];

	    struct atom* water = wres->atom + self->ligand.offset[i];
	    for(int j=0; j<self->prox.size; ++j){
		  struct residue* res = self->prox.mol[j]->residue+ self->prox.resindx[j];
		  for(int k=0; k<res->size; ++k){
			struct atom* curratm = res->atom + k;
			metdist2 = distsqr(self->metal->center, curratm->center);
			
			
			wdist2 = distsqr(water->center, curratm->center);
			int flag = 1;
			if(rule != '0'){
			      if(curratm->loc[0] == 'N' || curratm->loc[0] == 'O')
				    flag = 1;
			      else
				    flag = 0;
			}
			double angle = angle3d(self->metal->center, water->center, curratm->center);
			if(metdist2 <= includist2 &&
				    metdist2 > hbdst2 &&
				    wdist2 < hbdst2   &&
				    angle >= -520.0   &&
				    flag == 1
			  ){
			      char loc= get_res_loc(curratm, self->prox.restype[j]);
			      mediated_add(&(self->wmed), &(self->prox), j, k, loc, rule, i, wdist2, angle);
			}
		  }
	    }
      }
}
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
	    struct parameters* prm){

      //site_init(site);

      struct current_params mprm;
      //printf("struct atomic %d\n", self->atomic_no);
      currprm_adjst(&mprm, prm, self->atomic_no,rule, 0.0001);


      //metal_params_adjust(&mprm, prm, self->atomic_no, 0.5);



      site_fill_proximity(self, rna, 'N');
      site_comp_nucleic(self, rule, &mprm);
      site_fill_proximity(self, pro, 'P');
      site_comp_protein(self, rule, &mprm);
      site_fill_proximity(self, hoh, 'W');
      site_comp_water(self, rule, &mprm);

      site_comp_water_mediated(self, rule, &mprm);

      //site_addproxatm(site, hoh, 'W');
      //site_addhoh(site, hoh, rule, &mprm);

}

void site_fprint_chelate(struct site* self, FILE* fp, int* flag, int detail_flag)
{
      int size = self->ligand.size;
      struct residue* res1;
      struct residue* res2;
      struct atom* atom1;
      struct atom* atom2;
      int numcount = 0;
      int allowed = 0;
      for(int i=0; i<size; ++i){
	    char restype = self->ligand.restype[i];
	    char bindloc = self->ligand.loc[i];
	    res1 = self->ligand.mol[i]->residue + self->ligand.resindx[i];
	    atom1 = res1->atom + self->ligand.offset[i];
	    if(allowed ==0){
		  if(detail_flag == 0 && !(restype == 'N' && (bindloc == 'N' ||(bindloc == 'S' && strcmp(atom1->loc, "O2*") == 0) ))) continue;
		  if(detail_flag == 1 && restype != 'N' ) continue;
	    }
	    allowed = 1;

	    int counter = 0;
	    char bindloc2;
	    int c=0;
	    for(int j=i+1; j<size; ++j){
		  res2 = self->ligand.mol[j]->residue + self->ligand.resindx[j];
		  atom2 = res2->atom + self->ligand.offset[j];
		  if(atom1->resid == atom2->resid && 
			      strcmp(atom1->chain, atom2->chain) ==0 && 
			      strcmp(atom1->ins, atom2->ins) == 0){
			c=1;
			bindloc2 = self->ligand.loc[j];
			if(*flag == 0){
			      *flag = 1;
			      fprintf(fp, "           Metal Detail   |   ligand Detail    | positions and type |\n");
			      fprintf(fp, "       _____________________________________________________________\n");
			      fprintf(fp, "         resid  chn  met  |   resid  chn res   |  locations     type\n");
			      fprintf(fp, "       _____________________________________________________________\n");
			}
			if(counter == 0){
			      counter ++;
			      fprintf(fp, "\n");
			      fprintf(fp, "CHLT    %6d  %3s  %3s  |  %6d  %3s %3s   |  %3s", 
					  self->metal->resid,
					  self->metal->chain,
					  self->metal->resname,
					  atom1->resid,
					  atom1->chain,
					  atom1->resname,
					  atom1->loc);
			}
			fprintf(fp, "  %3s", 
				    atom2->loc);
			double d1 = self->ligand.dist[i];
			double d2 = self->ligand.dist[j];
			fprintf(fp, "  %8.3lf  %8.3lf", d1, d2);
		  }
	    }
	    numcount = numcount + c;
	    if(counter != 0){
		  if(counter == 1){
			char name[20];
			if(restype == 'N' && (bindloc =='N' && bindloc2 == 'N')){
			      strcpy(name,"alpha");
			}else if(restype == 'N' && (bindloc =='N' && bindloc2 == 'S')){
			      strcpy(name,"beta");
			}else if(restype == 'N' && (bindloc =='S' && bindloc2 == 'N')){
			      strcpy(name,"beta");
			}else if(restype == 'N' && (bindloc =='N' && bindloc2 == 'P')){
			      strcpy(name,"beta");
			}else if(restype == 'N' && (bindloc =='P' && bindloc2 == 'N')){
			      strcpy(name,"beta");
			}else if(restype == 'N' && (bindloc =='S' && bindloc2 == 'S')){
			      strcpy(name,"gamma");
			}else if(restype == 'N' && (bindloc =='P' && bindloc2 == 'P')){
			      strcpy(name,"gamma");
			}else if(restype == 'N' && (bindloc =='S' && bindloc2 == 'P')){
			      strcpy(name,"gamma");
			}else if(restype == 'N' && (bindloc =='P' && bindloc2 == 'S')){
			      strcpy(name,"gamma");
			}else if(restype == 'P' ){
			      strcpy(name,"P-alpha");
			}else{
			      fprintf(stderr, "Error in function %s()... Invalid binding location found %c %c %c.\n", __func__, restype, bindloc, bindloc2);
			      exit(EXIT_FAILURE);
			}
			fprintf(fp, "   %c-%s, ", restype, name);
			fprintf(fp, " bidentate");
		  }else if(counter == 2){
			fprintf(fp, " tridentate");
		  }else if(counter == 3){
			fprintf(fp, " tetradentate");
		  }else if(counter == 4){
			fprintf(fp, " pentadentate");
		  }else{
			fprintf(stderr, "Error in function %s()... too many dencity encountered.\n", __func__);
			exit(EXIT_FAILURE);
		  }
		  //		  fprintf(fp, "\n");
	    }
      }

      if(numcount >1 ){
	    fprintf(fp, "   *");
      }
      for(int k=1; k<numcount; ++k){
	    fprintf(fp, "*");
      }

}


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
	    struct sysparams* syspar) {
      int nsites = mol_atom_count(met);
      struct atom** met_atoms = (struct atom**) malloc (nsites * sizeof(struct atom*)); //new struct atom*[nsites];
      int count; // should be same as nsites;
      //    met->get_all_atoms(met_atoms, &count);

      
      mol_get_all_atoms(met, met_atoms, &count);
      struct site* sites = (struct site*) malloc ( nsites * sizeof(struct site));
      //struct site* sites = new struct site[nsites];//(struct Site*) malloc(met->size * sizeof(struct Site));
      //int nsites = 10*met->size;
      for(int i=0; i<nsites; ++i){
	    site_init(sites+i);
      }
      //    Map rnamap = Map(rna);
      //    Map promap = Map(pro);
      //molecule_metal_init(&rnamap, rna);
      for(int i=0; i<nsites; ++i){


	    struct atom* metal = met_atoms[i];
	    site_add_metal(sites+i,metal);
	    //site_init(sites+i, met->atom+met->resbeg[i], get_atomic_no(metal->loc));
	    site_populate(sites+i,rna, pro, hoh, met, rule, prm); //site_populate(sites+i, rna, pro,

	    // hoh,
	    // met, rule,
	    //  prm, &rnamap);

	    //printf("Total atoms in site %s %d\n",sites[i].metal->resname, sites[i].ligand.size);

      }


      for(int i=0; i<nsites; ++i){
	    if(i == 0){
		  fprintf(runpar->summaryfp, "\n\n\n\n            +---------------------- S U M M A R Y  R E P O R T "
			      "---------------------+\n\n\n\n");
		  fprintf(runpar->summaryfp, "      |  Metal Detail  |   Inner Coordination     |     Water Mediated "
			      "         |\n");

		  fprintf(runpar->summaryfp, "      ---------------------------------------------------------------------------\n");
		  fprintf(runpar->summaryfp, "       resid  chn  mtl  cord   NUC  PRO  H2O  MET    cnt     NUC  PRO  H2O  MET \n");
		  fprintf(runpar->summaryfp, "      ---------------------------------------------------------------------------\n");
	    
	    
/*	    fprintf(runpar->summaryfp, "\n\n\n\n            +---------------------- S U M M A R Y  R E P O R T "
			      "---------------------+\n\n\n\n");
		  fprintf(runpar->summaryfp, "      |  Metal Detail  |   Inner Coordination     "
			      "         |\n");

		  fprintf(runpar->summaryfp, "      ---------------------------------------------\n");
		  fprintf(runpar->summaryfp, "       resid  chn  mtl  cord   NUC  PRO  H2O  MET \n");
		  fprintf(runpar->summaryfp, "      ---------------------------------------------\n");*/
	   
	    }
	    site_fprint_summary(sites+i,runpar->summaryfp);
      }
      fprintf(runpar->summaryfp, "#\n");
      //fprintf(runpar->metalfp, "\n\n\n\n        +--------------------- D E T A I L  R E P O R T "
	//	  "---------------------+\n\n\n\n");

      int* bparray = (int*) malloc (nsites * sizeof(int));
      for(int k=0; k<nsites; ++k){

	    bparray[k] = 0;
      }
      for(int i=0; i<nsites; ++i){
	    int bpflag = 0;
	    site_bp_flag(sites+i,rnabp,runpar->detailflag, &bpflag, runpar->allbaseflag);
	    if(bpflag == 1){
		  bparray[i] = 1;
	    }
      }
      int flag=0;
      /*for(int i=0; i<nsites; ++i){
	    int bpflag = 0;
	    site_fprint_basepair(sites+i,runpar->metalfp, rnabp, &flag, runpar->detailflag, &bpflag, runpar->allbaseflag);
	    if(bpflag == 1){
		  bparray[i] = 1;
	    }
      }
      fprintf(runpar->metalfp,"#\n");


      flag=0;
      for(int i=0; i<nsites; ++i){
	    site_fprint_chelate(sites+i, runpar->metalfp, &flag, runpar->detailflag);
      }
      fprintf(runpar->metalfp,"\n#\n");



      flag=0;
      for(int i=0; i<nsites; ++i){
	    if(runpar->detailflag == 0 && bparray[i] == 0) continue; 
	    site_fprint_basepair_motifs(sites+i,runpar->metalfp, rnabp, sec, &flag, runpar->detailflag);
      }
      fprintf(runpar->metalfp,"#\n");
      for(int i=0; i<nsites; ++i){
	    if(runpar->detailflag == 0 && bparray[i] == 0) continue; 
	    if(runpar->detailflag == 1 && bparray[i] == 0) continue; 
	    site_fprint_inligand(sites+i,runpar);
	    site_fprint_angle(sites+i,runpar->metdetailfp);
	    //        sites[i].fprint_wmed(fp);
      }
*/

      char pmlfp_file_name[512];
      file_name_join(pmlfp_file_name, file_path, file_name, ".pml");

      FILE	*pmlfp;										

      pmlfp	= fopen( pmlfp_file_name, "w" );
      if ( pmlfp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			pmlfp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }

      fprintf(pmlfp, "load %s.cif\n", file_name);
      if(strcmp(syspar->chainparam, "-ML") == 0){
          fprintf(pmlfp, "select polyall, chain %s\n", syspar->chainvalparam);
      }else{
          fprintf(pmlfp, "select polyall, polymer\n");
      }

      fprintf(pmlfp, "select waterall, solvent\n");

      fprintf(pmlfp, "hide everything,%s\n", file_name);
      fprintf(pmlfp, "set cartoon_ring_mode, 0, polyall\n");
      fprintf(pmlfp, "set cartoon_ladder_mode, 0, polyall\n");
      fprintf(pmlfp, "show cartoon, polyall\n");
      fprintf(pmlfp, "set cartoon_nucleic_acid_color, red\n");

      //    fprintf(pmlfp, "select protein, polymer.protein\n");
	  //  fprintf(pmlfp, "select nucleic, polymer.nucleic\n");
	  //  fprintf(pmlfp, "select water, solvent\n");
	  //  fprintf(pmlfp, "show cartoon, %s\n", file_name);
      for(int i=0; i<nsites; ++i){

	    
//	    if(runpar->detailflag == 0 && bparray[i] == 0) continue; 
//	    if(runpar->detailflag == 1 && bparray[i] == 0) continue; 
//	    site_fprint_pml(sites+i, rnabp, pmlfp);
	    site_fprint_pml_wmed(sites+i, rnabp, pmlfp);

      }

      if( fclose(pmlfp) == EOF ) {			
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			pmlfp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }  



      
      for(int k=0; k<nsites; ++k){

	    bparray[k] = 0;
      }

      flag=0;
      for(int i=0; i<nsites; ++i){
	    int bpflag = 0;
	    site_fprint_wmed_basepair(sites+i,runpar->hohfp, rnabp, &flag, runpar->detailflag, &bpflag, runpar->allbaseflag);
	    if(bpflag == 1){
		  bparray[i] = 1;
	    }
      }

      if(syspar->nuparm_bp == 'T'){
	    flag=0;
	    for(int i=0; i<nsites; ++i){
		  int bpflag = 0;
		  site_fprint_wmed_bporient(sites+i, runpar->hohfp, rnabp, nuparm, &flag, runpar->detailflag, &bpflag, runpar->allbaseflag);
		  if(bpflag == 1){
			bparray[i] = 1;
		  }
	    }
      }
      flag=0;
      for(int i=0; i<nsites; ++i){

	    if(runpar->detailflag == 0 && bparray[i] == 0) continue; 
	    if(runpar->detailflag == 1 && bparray[i] == 0) continue; 
	    site_fprint_wmed_basepair_motifs(sites+i,runpar->hohfp, rnabp, sec, &flag, runpar->detailflag);
      }
      fprintf(runpar->hohfp,"#\n");
      //fprintf(fp,"#\n");

      //rnabp->fprint_bp(fp);
      for(int i=0; i<nsites; ++i){
	    if(runpar->detailflag == 0 && bparray[i] == 0) continue; 
	    if(runpar->detailflag == 1 && bparray[i] == 0) continue; 
	    //        sites[i].fprint_inligand(runpar);
	    //        sites[i].fprint_angle(fp);
	    site_fprint_wmed(sites+i,runpar->hohdetailfp);
      }
      //    fprintf(fp,"TER\n");

      //gen_verna(rna,rnabp, &rnamap, sec, sites, nsites);




      //    char nuc_file[512];
      //
      //    file_name_join(nuc_file, file_path,file_name, ".met");


      //    rnamap.gen_verna(nuc_file);


      free(bparray);
      bparray = NULL;

      free(sites);
      free(met_atoms);
}




