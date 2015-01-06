#include "m_pd.h"
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_rng.h>

static t_class * mt19937_class;

typedef struct _mt19937 {
  t_object x_obj;
  long seed;
  gsl_rng * r;
} t_mt19937;


void mt19937_free(t_mt19937 * x){
  char info[80];
  sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
  post(info);
  gsl_rng_free(x->r);
}

// constructor
void * mt19937_new(t_symbol * s, int argc, t_atom * argv){
  t_mt19937 * x = (t_mt19937 *) pd_new(mt19937_class);

  x->seed = (argc==0)? rand()*LONG_MAX : (long)atom_getint(argv);
  x->r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(x->r,x->seed);

  char info[80];
  sprintf(info,"%s generator with seed %d has been initiated.\n", gsl_rng_name(x->r), x->seed);
  post(info);

  return (void *) x;
}

void mt19937_setup(void){
  mt19937_class = class_new(gensym("mt19937"),
                            (t_newmethod) mt19937_new,
                            (t_method) mt19937_free,
                            sizeof(t_mt19937),
                            CLASS_DEFAULT,
                            A_GIMME, 0);
  class_sethelpsymbol(mt19937_class, gensym("mt19937-help"));

}
