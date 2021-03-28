#ifndef _FT_API_INTERNAL_H
#define _FT_API_INTERNAL_H
typedef struct
{
    int (*func)(POINTER, double*, double*);
    void *params;
}GENERIC_VELOCITY_PARAMS;
#endif
