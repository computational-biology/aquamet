//
// Created by parthajit on 16/7/20.
//

#ifndef CPPMET_TREE_H
#define CPPMET_TREE_H
#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#define SNIPPET_SIZE 10
struct counting_tree_t{
    char text[SNIPPET_SIZE];
    int txtlen;
    long freq;
    struct counting_tree_t* left;
    struct counting_tree_t* right;
};


struct counting_tree_t* counting_tree_getnode(char* line);


void counting_tree_insert(struct counting_tree_t* root, char* line);
long counting_tree_node_count(struct counting_tree_t* t);

void counting_tree_fprintf(FILE* fp, struct counting_tree_t* t, int* count, char newline);
void counting_tree_free(struct counting_tree_t* t);


int counting_tree_search(struct counting_tree_t* t, char* word);
#endif //CPPMET_TREE_H
