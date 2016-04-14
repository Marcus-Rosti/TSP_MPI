#include <stdlib.h>
#include <list.h>

void add_node(List l, const int value)
{
    Node * iter = l.head;
    while (iter != NULL) iter = (Node *) iter->next;
    iter = malloc(sizeof(int) + sizeof(Node *));
    iter->value = value;
}

void delete_node(List l, const int value)
{
    Node * iter = l.head;
    while ( ((Node *)(iter -> next)) -> value  != value) iter = (Node *) iter->next;
    Node *new_next = (Node *) iter->next;
    free(iter->next);
    iter->next = new_next;
}