/* stack.c */
#include "stack.h"

static item_t stack[512];
static int top = 0;

void push(item_t p)
{
	stack[top++] = p;
}

item_t pop(void)
{
	return stack[--top];
}

int is_empty(void)
{
	return top == 0;
}

