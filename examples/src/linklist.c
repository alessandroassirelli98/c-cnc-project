/* example of a linked list implementation */
/* Lesson of 201118                        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 1. Define the structs representing the objects:

// A structire representing the elements to be stored in the
// linked list
typedef struct {
  char *name;
  void *next; // if NULL, this is the last element
  void *prev; // if NULL, this is the first element
} element_t;

// A structure representing the list itself
typedef struct {
  element_t *first, *last;
} list_t;

// 2. Define the functions for dealing with the list

// Create a new list, with one element newly created with "name"
list_t *list_new(char *name) {
  list_t *l = malloc(sizeof(list_t)); // memory space for the list
  element_t *e = malloc(sizeof(element_t)); // memory for its first element
  e->name = malloc(strlen(name)+1); // memory for the name string in the element
  strncpy(e->name, name, strlen(name)); // copy the parameter into e->name
  e->prev = NULL; // these shall already be NULL, but different compilers 
  e->next = NULL; // can behave differently, so better safe than sorry
  l->first = e; // e is the only element, so it is both the first
  l->last =e;   // and the last
  return l; // return the new list (as a pointer)
}

// apped an existing element e to the list
void list_append_element(list_t *list, element_t *e) {
  list->last->next = e; // point the next field of the last element to e
  e->prev = list->last; // previous item of e is the current last element
  list->last = e; // last element of the list is e
  e->next = NULL; // just to be sure
}

// create and append a new element with "name"
void list_append(list_t *list, char *name) {
  // create the element
  element_t *e = malloc(sizeof(element_t));
  e->name = malloc(strlen(name) + 1);
  strncpy(e->name, name, strlen(name));
  // append it to the list
  list_append_element(list, e);
}

// create and insert a new element after the first element named after
void list_insert(list_t *list, char *name, char *after) {
  element_t *new = malloc(sizeof(element_t)); // the new element to be created
  element_t *e; // this is gonna be the element named "after" 
  // initialize new
  new->name = malloc(strlen(name) + 1);
  strncpy(new->name, name, strlen(name));
  // loop and search for the first element with name=after
  e = list->first; // start looping form the beginning of the list
  do {
    if (strcmp(e->name, after) == 0) {
      new->prev = e;
      new->next = e->next;
      ((element_t *)e->next)->prev = new; // the cast is needed for e->next is of type void*
      e->next = new;
      break; // stop searching after the first item is found
    }
  } while ((e = e->next));
}

// free memory
void list_free(list_t *list) {
  element_t *e, *next;
  // first, free list contents
  e = list->first; // start looping form the beginning of the list
  do {
    next = e->next; // temporary storage, for we are going to free(e)
    free(e->name);
    free(e);
  } while ((e = next));
  // finally, free the list itself
  free(list);
}


int main() {
  char names[4][10] = {"one", "two", "three", "four"};
  element_t *e;
  size_t i;
  // create a list:
  list_t *list = list_new("zero");

  // populate the list with the entries of names array
  for (i = 0; i < 4; i++) {
    // create element and append it:
    e = malloc(sizeof(element_t));
    e->name = malloc(strlen(names[i]) + 1);
    strncpy(e->name, names[i], strlen(names[i]));
    list_append_element(list, e);
    // or, more easily:
    // list_append(list, names[i]);
  }

  // insert a new element in between "two" and "three"
  list_insert(list, "two.five", "two");

  // loop from first to last element:
  e = list->first;
  i = 0;
  do {
    printf("element %ld: %s\n", i++, e->name);
  } while ((e = e->next));

  // reverse loop:
  e = list->last;
  i = 0;
  do {
    printf("element %ld: %s\n", i++, e->name);
  } while ((e = e->prev));

  // free memory:
  list_free(list);

}