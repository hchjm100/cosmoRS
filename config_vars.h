#ifndef _CONFIG_VARS_H_
#define _CONFIG_VARS_H_
#include <stdint.h>

// todo
#define string(a,b)  extern char *  a;
#define real(a,b)    extern double  a;
#define real3(a,b)   extern double  a[3];
#define integer(a,b) extern int64_t a;

#include "config.template.h"
#include "config_part_mass_ZXY.h"

#undef string
#undef integer
#undef real
#undef real3

#endif /* _CONFIG_VARS_H_ */
