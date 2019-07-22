#include "linked.h"

// int main()
// {
//   List l = Init_Empty_List();
//   printf("\nCreated Empty List:\n\n" );
//   print_list(l);
//
//   printf("\nAppended Three and Prepended Four\n\n" );
//   l = append("Three",l);
//   l = prepend("Four",l);
//
//   print_list(l);
//
//   printf("\nAppended 'Two' and 'Blast Off' and Prepended 'Six'\n\n" );
//   l = append("Two",l);
//   l = append("Blast Off",l);
//   l = prepend("Six",l);
//
//   print_list(l);
//
//   printf("\nInserted 'Five' after 'Six' and 'One' before 'Blast Off'\n\n" );
//
//   l = insert_after("Six", "Five", l);
//   l = insert_before("Blast Off","One", l);
//
//   print_list(l);
//
//   l = reverse(l);
//   printf("\nNow printing reversed\n\n");
//   print_list(l);
//
//
//   printf("\nThis should be the original list:\n\n");
//   l = reverse(l);
//   print_list(l);
//
//   free_list(l);
//
//   return 0;
// }
double* findListLineSetLabel(List l, int n , int newLabel)
{
  Cell* temp = l.head;
  while (temp != NULL) {
    if (temp->label == n) {
      temp->label = newLabel;
      return temp->line;
    }
    temp = temp->next;
  }
  fprintf(stderr, "n (%d) not found\n",n );
  exit(1);
}

double* findListLine(List l, int n)
{
  Cell* temp = l.head;
  while (temp != NULL) {
    if (temp->label == n) {
      return temp->line;
    }
    temp = temp->next;
  }
  fprintf(stderr, "n (%d) not found\n",n );
  exit(1);
}


List Init_Empty_List()
/* Empty List Created*/
{
  // List declared and memory allocated:
  List l;
  l.head = malloc(sizeof(Cell));
  l.tail = malloc(sizeof(Cell));

  // Cells arranged:
  l.tail->next = NULL;
  l.head->next = l.tail;
  l.head->prev = NULL;
  l.tail->prev = l.head;
  l.head->label =  l.tail->label = -1;
  l.head->line = l.tail->line = NULL;

  return l;
}


List append(struct denseData *ds, List l, int n)
/* Function to append a linked list with a new string. */
{
  printf("n is %d\n",n );

  // Check if l.head, l.tail have words:
  if (l.head->line == NULL)
  {
    l.head->label = n;
    l.head->line = malloc(sizeof(double)*ds->nInstances);
    for (int i = 0; i < ds->nInstances; i++) {
      l.head->line[i] = 0.0;
      for (int j = 0; j < ds->nFeatures; j++) {
        l.head->line[i] += ds->data[i][j]*ds->data[n][j];
      }
      l.head->line[i]*= ds->y[i]*ds->y[n];
    }
    return l;
  }
  if (l.tail->line == NULL)
  {
    l.tail->label = n;
    l.tail->line = malloc(sizeof(double)*ds->nInstances);
    for (int i = 0; i < ds->nInstances; i++) {
      l.tail->line[i] = 0.0;
      for (int j = 0; j < ds->nFeatures; j++) {
        l.tail->line[i] += ds->data[i][j]*ds->data[n][j];
      }
      l.tail->line[i]*= ds->y[i]*ds->y[n];
    }
    return l;
  }

  /* If l.head & l.tail are filled in
  Create a new cell and arrange it after tail: */

  Cell *new = malloc(sizeof(Cell));
  new->line = malloc(sizeof(double)*ds->nInstances);
  for (int i = 0; i < ds->nInstances; i++) {
    new->line[i] = 0.0;
    for (int j = 0; j < ds->nFeatures; j++) {
      new->line[i] += ds->data[i][j]*ds->data[n][j];
    }
    new->line[i]*= ds->y[i]*ds->y[n];
  }

  new->label = n;
  new->prev = l.tail;
  l.tail->next = new;
  new->next = NULL;
  l.tail = new;

  return l;
}

void print_list(List l)
/* Function to print a linked list (each entry on new line)*/
{
  //l.head->prev is used as temp - ends up back to being NULL:

  l.head->prev = l.head;
  while(l.head->prev != NULL)
  {
    if (l.head->prev->line == NULL)
    /* Ensures nothing prints for empty list and null not printed
    for list of one entry */
    {
      l.head->prev = NULL;
      return;
    }
    printf("%lf\n",l.head->prev->line[0] );
    l.head->prev = l.head->prev->next;
  }
}

void free_list(List l)
/* Function to free memory from a linked list*/
{
  //l.head->prev used as tmp again
  l.head->prev = l.head->next->next;
  while (l.head->prev != NULL)
  {
    free(l.head->prev->prev->line);
    free(l.head->prev->prev);
    l.head->prev = l.head->prev->next;
  }

  //l.head and l.tail not freed in while loop:
  free(l.tail->line);
  free(l.tail);
  free(l.head->line);
  free(l.head);

  return;
}

List delete(int find, List l)
/* Function to insert a new_word after a given word in linked list.*/
{
  // Checks case if list is of length 1.
  if (l.tail == NULL)
  {
    if (l.head->label == find)
    {
      free(l.head->line);
      free(l.head);
      l.head = NULL;
      return l;
    }
    else
    {
      printf("Word not found.\n");
      exit(1);
    }
  }

  l.head->prev = l.head;
  while(l.head->prev->label != find)
  {
    if (l.head->prev == NULL)
    {
      printf("Word not found.\n");
      exit(191);
    }
    l.head->prev = l.head->prev->next;
  }
  if (l.head->prev == l.head) {
    l.head = l.head->next;
    free(l.head->prev->line);
    free(l.head->prev);
    l.head->prev = NULL;
  }
  else if (l.head->prev == l.tail){
    l.tail = l.tail->prev;
    l.tail->next = NULL;
    free(l.head->prev->line);
    free(l.head->prev);
    l.head->prev = NULL;
  }
  else{
    l.head->prev->prev->next = l.head->prev->next;
    l.head->prev->next->prev = l.head->prev->prev;
    free(l.head->prev->line);
    free(l.head->prev);
    l.head->prev = NULL;
  }
  return l;
}
