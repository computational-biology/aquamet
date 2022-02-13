//
// Created by parthajit on 13/7/20.
//


#include "struct.h"

void structure_init(struct structure* self, char *dat_file, int size) {
    self->num_res = size;
    if(self->num_res <= 0) return;
    self->secseq = (char*) malloc ((self->num_res+1) * sizeof(char));
    //self->dbn = (char*) malloc ((self->num_res+200+1) * sizeof(char)); // 200 extra to accomodate & (cnain break)
    //self->priseq = (char*) malloc ((self->num_res+200+1) * sizeof(char)); // 200 extra to accomodate & (cnain break)
    self->num_chain = 0;
    FILE* fpdat = fopen(dat_file, "r");
    assert( fpdat != NULL);
    char line[1024];
    char sep[]="\t :\n";
    char* token;
    long int index = 0;
    while(fgets(line, 1024, fpdat) != NULL){
        if(line[0] == '>'){
            token = strtok(line, sep);
            token = strtok(NULL, sep);//second one is the chain info
            strcpy(self->chain[self->num_chain], token);
            self->num_chain++;
        }else{
            token = strtok(line, sep);
            int len = strlen(token);
            assert(len>0);
            strcpy(self->secseq+index, token);
            index = index + len;
        }
    }
    //fprintf(stdout, "size= %ld: %s\n",this->num_res, this->secseq);
    assert(self->num_res == index);

    fclose(fpdat);
   /* FILE* fpdbn = fopen(dbn_file,"r");
    assert(fpdbn != NULL);

    char* longline = (char*) malloc((this->num_res+200) * sizeof(char));

    fgets(longline, this->num_res+200, fpdbn);

    fgets(longline, this->num_res+200, fpdbn);
    token = strtok(longline,"\n");
    strcpy(this->priseq, token);

    fgets(longline, this->num_res+200, fpdbn);
    token = strtok(longline,"\n");
    strcpy(this->dbn, token);*/

    //fprintf(stdout, "%s\n",secst->dbn);
    //fclose(fpdbn);
}


void structure_free(struct structure* self) {
    if(self->num_res > 0){
        free(self->secseq);
        free(self->dbn);
        free(self->priseq);
    }
    //printf("Sec free called\n");
}
