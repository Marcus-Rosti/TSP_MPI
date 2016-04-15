#ifndef LIST_H
#define LIST_H

typedef struct {
    struct Node * next;
    int * path;
    int path_size;
    int value;
} Node;

typedef struct {
    Node * head;
    int size;
} List;

void add_node(List l, const int value);

void delete_node(List l, const int value);

#endif
