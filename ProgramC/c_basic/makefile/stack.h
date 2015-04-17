/* stack.h */
#ifndef STACK_H
#define STACK_H

#include "main.h" /* provides definition for item_t */

extern void push(item_t);
extern item_t pop(void);
extern int is_empty(void);

#endif

