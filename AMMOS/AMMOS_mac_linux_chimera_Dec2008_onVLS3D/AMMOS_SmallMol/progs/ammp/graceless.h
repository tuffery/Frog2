
#ifdef GRACELESS
/* dec malloc controls */
#include <sys/types.h>
unsigned long __noshrink = 0; /* force release with free */
size_t __minshrink = 16384;
double __minshrinkfactor = 0.001;
size_t __mingrow = 49152;
double __mingrowfactor = 0.1;
unsigned long __madvisor = 1;
unsigned long __small_buff = 0;
int __fast_free_max = 13;
#endif

