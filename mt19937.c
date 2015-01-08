#include "m_pd.h"
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * mt19937_class;

typedef struct _mt19937 {
  t_object x_obj;
  gsl_rng * r;
} t_mt19937;




// poisson sampling
void mt19937_poisson(t_mt19937 * x, t_floatarg mu){
  unsigned int k = gsl_ran_poisson(x->r, (double)mu);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}

// normal sampling
void mt19937_gaussian(t_mt19937 * x, t_floatarg sigma, t_floatarg mu){
  double s = (sigma==0)? 1.0 : (double)sigma;
  double k = gsl_ran_gaussian(x->r, s);
  outlet_float(x->x_obj.ob_outlet, mu+(float)k);
}
// exponential
void mt19937_exponential(t_mt19937 * x, t_floatarg lambda){
  double l = (lambda==0)? 1.0 : (double)lambda;
  double k = gsl_ran_exponential(x->r, l);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}
// gamma
void mt19937_gamma(t_mt19937 * x, t_floatarg a, t_floatarg b){
  double kappa = (a==0)? 2.0 : (double)a;
  double theta = (b==0)? 2.0 : (double)b;
  double k = gsl_ran_gamma(x->r, kappa, theta);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}





void mt19937_seed(t_mt19937 * x, t_floatarg s){
  gsl_rng_set(x->r,(long)s);
}

// destroyer
void mt19937_free(t_mt19937 * x){
  char info[80];
  sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
  post(info);
  gsl_rng_free(x->r);
}


// constructor
void * mt19937_new(t_symbol * s, int argc, t_atom * argv){
  t_mt19937 * x = (t_mt19937 *) pd_new(mt19937_class);


  long seed = (argc==0)? (long)(rand()*LONG_MAX) : (long)atom_getint(argv);

  x->r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(x->r,seed);

  char info[80];
  sprintf(info,"%s generator with seed %d has been initiated.\n", gsl_rng_name(x->r), (int)seed);
  post(info);

  outlet_new(&x->x_obj, &s_float);
  return (void *) x;
}

void mt19937_setup(void){
  mt19937_class = class_new(gensym("mt19937"),
                            (t_newmethod) mt19937_new,
                            (t_method) mt19937_free,
                            sizeof(t_mt19937),
                            CLASS_DEFAULT,
                            A_GIMME, 0);
  // poisson
  class_addmethod(mt19937_class,
                  (t_method)mt19937_poisson, gensym("poisson"),
                  A_DEFFLOAT, 0);
  // gaussian
  class_addmethod(mt19937_class,
                  (t_method)mt19937_gaussian, gensym("gaussian"),
                  A_DEFFLOAT, A_DEFFLOAT, 0);

  // exponential
  class_addmethod(mt19937_class,
                  (t_method)mt19937_exponential, gensym("exponential"),
                  A_DEFFLOAT, 0);
  // gamma
  class_addmethod(mt19937_class,
                  (t_method)mt19937_gamma, gensym("gamma"),
                  A_DEFFLOAT, A_DEFFLOAT, 0);

  // seed
  class_addmethod(mt19937_class,
                  (t_method)mt19937_seed, gensym("seed"),
                  A_DEFFLOAT, 0);

  // help file
  class_sethelpsymbol(mt19937_class, gensym("mt19937-help"));

}
