#ifndef LIST_H
#define LIST_H

typedef struct {
    struct Node * next;
    int value;
} Node;

typedef struct {
    Node * head;
} List;

void add_node(List l, const int value);

void delete_node(List l, const int value);

#endif
