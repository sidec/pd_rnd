#include "m_pd.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * rnd_class;

typedef struct _rnd {
  t_object x_obj;
  gsl_rng * r;
} t_rnd;




// poisson sampling
void rnd_poisson(t_rnd * x, t_floatarg mu){
  unsigned int k = gsl_ran_poisson(x->r, (double)mu);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}

// normal sampling
void rnd_gaussian(t_rnd * x, t_floatarg sigma, t_floatarg mu){
  double s = (sigma==0)? 1.0 : (double)sigma;
  double k = gsl_ran_gaussian(x->r, s);
  outlet_float(x->x_obj.ob_outlet, mu+(float)k);
}
// exponential
void rnd_exponential(t_rnd * x, t_floatarg lambda){
  double l = (lambda==0)? 1.0 : (double)lambda;
  double k = gsl_ran_exponential(x->r, l);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}
// gamma
void rnd_gamma(t_rnd * x, t_floatarg a, t_floatarg b){
  double kappa = (a==0)? 2.0 : (double)a;
  double theta = (b==0)? 2.0 : (double)b;
  double k = gsl_ran_gamma(x->r, kappa, theta);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}
// uniform
void rnd_uniform(t_rnd * x){
  double k = gsl_rng_uniform(x->r);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}
void rnd_uniform_int(t_rnd * x, t_floatarg n){
  unsigned long int k = gsl_rng_uniform_int(x->r, (unsigned long int)n);
  outlet_float(x->x_obj.ob_outlet, (float)k);
}



void rnd_seed(t_rnd * x, t_floatarg s){
  if(s==0){
    gsl_rng_set(x->r,(long)(rand()*INT_MAX));
  } else {
    gsl_rng_set(x->r,(long)s);
  }
}

// destroyer
void rnd_free(t_rnd * x){
  char info[80];
  sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
  post(info);
  gsl_rng_free(x->r);
}


// constructor
void * rnd_new(t_symbol * s, int argc, t_atom * argv){
  t_rnd * x = (t_rnd *) pd_new(rnd_class);

  long seed = (argc==2)? (long)atom_getint(argv+1):(long)(rand()*INT_MAX);

  t_symbol * algorithm = atom_getsymbol(argv);

  if (strcmp(algorithm->s_name, "gfsr4")==0){
    x->r = gsl_rng_alloc(gsl_rng_gfsr4);
  } else if (strcmp(algorithm->s_name,"ranlxs")==0){
    x->r = gsl_rng_alloc(gsl_rng_ranlxs2);
  } else if (strcmp(algorithm->s_name,"taus")==0){
    x->r = gsl_rng_alloc(gsl_rng_taus2);
  }  else {
    x->r = gsl_rng_alloc(gsl_rng_mt19937);
  }

  gsl_rng_set(x->r,seed);

  char info[80];
  sprintf(info,"%s generator with seed %d has been initiated.\n", gsl_rng_name(x->r), (int)seed);
  post(info);

  outlet_new(&x->x_obj, &s_float);
  return (void *) x;
}

void rnd_setup(void){
  rnd_class = class_new(gensym("rnd"),
                            (t_newmethod) rnd_new,
                            (t_method) rnd_free,
                            sizeof(t_rnd),
                            CLASS_DEFAULT,
                            A_GIMME, 0);
  // poisson
  class_addmethod(rnd_class,
                  (t_method)rnd_poisson, gensym("poisson"),
                  A_DEFFLOAT, 0);
  // gaussian
  class_addmethod(rnd_class,
                  (t_method)rnd_gaussian, gensym("gaussian"),
                  A_DEFFLOAT, A_DEFFLOAT, 0);

  // exponential
  class_addmethod(rnd_class,
                  (t_method)rnd_exponential, gensym("exponential"),
                  A_DEFFLOAT, 0);
  // gamma
  class_addmethod(rnd_class,
                  (t_method)rnd_gamma, gensym("gamma"),
                  A_DEFFLOAT, A_DEFFLOAT, 0);
  // uniform
  class_addmethod(rnd_class,
                  (t_method)rnd_uniform, gensym("uniform"),0);
  // uniform int
  class_addmethod(rnd_class,
                  (t_method)rnd_uniform_int, gensym("uniform_int"),
                  A_DEFFLOAT, 0);

  // seed
  class_addmethod(rnd_class,
                  (t_method)rnd_seed, gensym("seed"),
                  A_DEFFLOAT, 0);

  // help file
  class_sethelpsymbol(rnd_class, gensym("rnd-help"));

}
