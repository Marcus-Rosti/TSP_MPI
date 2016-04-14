typedef struct {
  struct Node * next;
  int value;
} Node;

typedef struct {
  Node * head; 
} List;

void add_node(List l, const int value)
{
    Node * iter = l.head;
    while (iter != NULL) iter = iter->next;
    iter = malloc(sizeof(int) + sizeof(Node *));
    iter->value = value;
}

void delete_node(List l, const int value)
{
    Node * iter = l.head;
    while (iter->next->value != value) iter = iter->next;
    temp * next_next = iter->next;
    free(iter->next);
    iter->next = new_next
}