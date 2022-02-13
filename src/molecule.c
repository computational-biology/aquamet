//
// Created by parthajit on 10/7/20.
//

#include "molecule.h"




void mol_init(struct molecule* self)
{
    self->size = 0;
    self->residue = NULL;
}
int mol_atom_count(struct molecule* self){
      int count = 0;
      for(int i=0; i<self->size; ++i){
	    count = count + self->residue[i].size;
      }
      return count;
}

void mol_get_all_atoms(struct molecule* self, struct atom* atmarray[], int* atmcount){
      int count =0;
      for(int i=0; i<self->size; ++i){
	    for(int j=0; j<self->residue[i].size; ++j){
		  atmarray[count] = self->residue[i].atom+j;
		  count ++;
	    }
      }
      *atmcount = count;
}

void mol_polulate(struct molecule* self, struct atom* atoms, int size)
{

      if(size ==0) return;
      int resid = atoms[0].resid;
      char chain[5];
      char ins[5];
      strcpy(ins, atoms[0].ins);
      strcpy(chain, atoms[0].chain);
     
      
      int resindex = 0;
      for(int i=0; i<size; ++i){
	    if(atoms[i].resid == resid && strcmp(atoms[i].chain, chain) == 0 && strcmp(atoms[i].ins, ins) == 0){
		  resi_add(self->residue + resindex, atoms + i);
	    }else{
		  resid = atoms[i].resid;
		  strcpy(chain, atoms[i].chain);
		  strcpy(ins, atoms[i].ins);

		  resindex ++;
		  resi_add(self->residue + resindex, atoms + i);
	    }
      }
      self->size = resindex + 1;
}

void mol_scan_cif(struct molecule* self, char *ciffile,  int (*pf)(char *)) {
      //residue = new Residue[size];
      //self->residue = new  Residue[size];
      //self->max = size;
      self ->size = 0;


      char tokary[100][100];
      FILE* fp = fopen(ciffile, "r");
      assert(fp != NULL);
      char line[1024];
      int grouploc=-1;  // group_PDB
      int type_symbloc = -1;

      int idloc = -1;   // id
      int locloc = -1;  // label_atom_id
      int resnameloc = -1;  // label_comp_id
      int chainloc = -1; // auth_asym_id
      int residloc = -1; // auth_seq_id
      int modelloc = -1; // pdbx_PDB_model_num
      int occuloc  = -1; //occupancy
      int xloc = -1;
      int yloc = -1;
      int zloc = -1;
      int count = -1;
      while(fgets(line, 1024, fp) != NULL){
	    if(strncmp(line, "_atom_site.", 11) != 0 && count == -1) continue;
	    if(strncmp(line, "_atom_site.", 11) != 0 && count != -1) break;
	    count ++;
	    if(strncmp(line,"_atom_site.group_PDB", 20) == 0 && line[20] != '_' ){
		  grouploc = count;
	    }else if(strncmp(line,"_atom_site.id", 13) == 0 && line[13] != '_' ){
		  idloc = count;
	    }else if(strncmp(line,"_atom_site.type_symbol", 22) == 0 && line[22] != '_' ){
		  type_symbloc = count;
	    }else if(strncmp(line,"_atom_site.label_atom_id", 24) == 0 && line[24] != '_' ){
		  locloc = count;
	    }else if(strncmp(line,"_atom_site.label_comp_id", 24) == 0 && line[24] != '_' ){
		  resnameloc = count;
	    }else if(strncmp(line,"_atom_site.auth_asym_id", 23) == 0 && line[23] != '_' ){
		  chainloc = count;
	    }else if(strncmp(line,"_atom_site.auth_seq_id", 22) == 0 && line[22] != '_' ){
		  residloc = count;
	    }else if(strncmp(line,"_atom_site.Cartn_x", 18) == 0 && line[18] != '_' ){
		  xloc = count;
	    }else if(strncmp(line,"_atom_site.Cartn_y", 18) == 0 && line[18] != '_' ){
		  yloc = count;
	    }else if(strncmp(line,"_atom_site.Cartn_z", 18) == 0 && line[18] != '_' ){
		  zloc = count;
	    }else if(strncmp(line,"_atom_site.occupancy", 20) == 0 && line[20] != '_' ){
		  occuloc = count;
	    }else if(strncmp(line,"_atom_site.pdbx_PDB_model_num", 29) == 0 && line[29] != '_' ){
		  modelloc = count;
	    }else{
		  ;
	    }

      }
      char sep[10] = "\t \"\n";
      char* token;
      int prevresid = -999;
      char prevchain[10] = "-";
      int resindx = -1;
      int modelflag=0;
      int modelnum = -1;
      while(line[0] != '#'){
	    int indx=0;
	    token = strtok(line, sep);
	    while(token != NULL){
		  //	         printf("%s, ",token);
		  strcpy(tokary[indx], token);
		  indx++;
		  token = strtok(NULL, sep);
	    }

	    if(tokary[type_symbloc][0] == 'H'){ // We dont consider H
		  if(fgets(line, 1024, fp) == NULL) break;
		  continue;
	    }
	    if(pf(tokary[resnameloc]) < 0 ){
		  if(fgets(line, 1024, fp) == NULL) break;
		  continue;
	    }

	    int resid = atol(tokary[residloc]);
	    char chn[10];
	    strcpy(chn, tokary[chainloc]);
	    if(modelflag==0){
		  modelnum = atoi(tokary[modelloc]);
		  modelflag = 1;
	    }
	    if(modelnum != atoi(tokary[modelloc])){
		  fclose(fp);
		  return;
	    }
	    if(prevresid != resid || strcmp(prevchain, chn) != 0){
		  resindx ++;
		  self->size ++;
		  assert(self->size < self->max );
		  prevresid = resid;
		  strcpy(prevchain, chn);
	    }

	    struct atom atm;
	    char atmhetatm[20];
	    strcpy(atmhetatm, tokary[grouploc]);

	    atm.id = atol(tokary[idloc]);
	    atm.resid = atol(tokary[residloc]);
	    strcpy(atm.loc, tokary[locloc]);
	    strcpy(atm.resname, tokary[resnameloc]);
	    strcpy(atm.chain, tokary[chainloc]);
	    if(strncmp(tokary[grouploc], "HETATM", 5) == 0)
		  atm.type = 'H';
	    else
		  atm.type = 'A';

	    atm.center.x = atof(tokary[xloc]);
	    atm.center.y = atof(tokary[yloc]);
	    atm.center.z = atof(tokary[zloc]);
	    atm.occu = atof(tokary[occuloc]);

	    resi_add(self->residue + resindx, &atm);

	    if(fgets(line, 1024, fp) == NULL) break;
	    //printf("%s",line);
      }
      fclose(fp);
      //if(self->size>0)
      //printf("Moltest1=%s %lf\n", self->residue[0].atom[0].resname, residue[0].atom[0].x);
      //else{
      //   printf("No Mol found\n");
      //}

}

void mol_scan_rna(struct molecule* self, char *corfile, struct rnabp *rnabp) {
      FILE* fp = fopen(corfile, "r");


      if(fp == NULL){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s(). (File: %s, Line %d)... Cannot open <accn>_rna.pdb\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }

      self->size = rnabp->nres;
      self->hetatm = rnabp->hetatm;
      self->atm = rnabp->atm;

      self->residue = (struct residue*) malloc (self->size * sizeof(struct residue));
      for(int i=0; i<self->size; ++i){
	    resi_init(self->residue + i);
      }
      //cout<<"SIZE="<<size<<endl;
      char line[1024];

      char sep[] = "\t \n";
      char *token;
      int prev_resid = -999;
      char chain[5];
      strcpy(chain,"--NA");
      char ins[5];
      strcpy(ins, "--NA");
      int i= -1;
      char tmptoken[2];
      struct atom atom ;//= new Atom();
      while(fgets(line,1024, fp) != NULL){
	    if(line[0] == '#') continue;
	    int resid;
	    if(strncmp(line,"END", 3)==0){
		  break;
	    }
	    token = strtok(line, sep); /* skip "ATOM"*/

	    token = strtok(NULL, sep);
	    atom.id = atol(token);
	    if(line[16] != ' '){
		  tmptoken[0] = line[16];
		  tmptoken[1] = '\0';
		  line[16] = ' ';
		  token = strtok(NULL, sep);
		  strcpy(atom.loc, token);
		  strcat(atom.loc, tmptoken);
		  // printf("Found %s in %s\n", atom->loc, corfile);
	    }else{

		  token = strtok(NULL, sep);
		  strcpy(atom.loc, token);
	    }
	    if(atom.loc[0] == 'H') continue; // We dont consider Hydrogen

	    token = strtok(NULL, sep);
	    strcpy(atom.resname, token);

	    token = strtok(NULL, sep);
	    resid = atol(token);
	    atom.resid = rnabp->bp[resid-1].cifid;

	    char token1[20];
	    memset(token1, '\0', 20);
	    strncpy(token1, line+30,8);
	    atom.center.x = atof(token1);

	    memset(token1, '\0', 20);
	    strncpy(token1, line+38,8);
	    atom.center.y = atof(token1);



	    memset(token1, '\0', 20);
	    strncpy(token1, line+46,8);
	    atom.center.z = atof(token1);
	    memset(token1, '\0', 10);
	    strncpy(token1, line+54,6);
	    atom.occu = atof(token1);
	    strcpy(atom.chain, rnabp->bp[resid-1].chain); // cor file does not contain chain info
	    strcpy(atom.ins, rnabp->bp[resid-1].ins); // cor file does not contain chain info
	    atom.type ='A'; // normal atom. not hetatm
	    if(prev_resid != resid) {
		  prev_resid = resid;
		  i++;
	    }
	    resi_add(self->residue + i, &atom);
      }
      assert(i == self->size-1);
      //    int last= residue[i].size-1;
      //    int last1= residue[i-1].size-1;
      //printf("z=%ld ",    self->residue[i-1].atom[last1].resid);
      //printf("z=%ld ",    self->residue[i].atom[0].resid);
      // printf(" z=%ld\n", self->residue[i].atom[last].resid);
      fclose(fp);
}


void mol_free(struct molecule* self) {
      if(self->residue != NULL){
            free(self->residue);
      }
      // printf("Free invoked for Molecule %ld\n", max);

}

void mol_copy(struct molecule* dst, struct molecule* src) {
      dst->size = src->size;
      
      dst->max = src->max;
      if(src->size == 0) return;
      dst->atm = src->atm;
      dst->hetatm = src->hetatm;
      dst->residue = (struct residue*) malloc (dst->size * sizeof(struct residue));
      for(int i=0; i<dst->size; ++i){
	    resi_copy(dst->residue + i, src->residue+i);
      }
}

void mol_create(struct molecule* self, int maxsize, int atm, int hetatm) {
      //Molecule::Molecule(int maxsize, bool atm, bool hetatm) {
      self->max = maxsize;
      self->size = 0;
      if(maxsize == 0) return;
      //    self->residue = new Residue[self->max];
      self->residue = (struct residue*) malloc (self->max * sizeof(struct residue));
      self->atm = atm;
      self->hetatm = hetatm;
}

void mol_reset(struct molecule* self, int atm, int hetatm) {
      for(int i=0; i<self->max; ++i){
	    self->residue[i].size = 0;
      }
      self->size = 0;
      self->atm = atm;
      self->hetatm = hetatm;
}



