#include "m_pd.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


static t_class * gaussian_class;

typedef struct _gaussian {
  t_object x_obj;
  // t_sample freq;


  /* double sigma; */
  t_float  f;
  t_sample freq;
  double sigma;

  t_float count;

  t_sample rnd_x;

  gsl_rng * r;

  t_inlet * sigma_in;
  t_outlet * x_out;
} t_gaussian;


t_int *gaussian_perform(t_int *w){
  t_gaussian *x =  (t_gaussian *)(w[1]);
  t_sample *s_freq = (t_sample *)(w[2]);
  t_sample *s_sigma = (t_sample *)(w[3]);
  t_sample *x_out = (t_sample *)(w[4]);
  int n = (int)(w[5]);


  int i;
  for(i = 0; i < n; i++){
    float freq = fabs(s_freq[i]);
    x->freq = (freq==0)?  x->freq : freq;
    if(x->count >= (sys_getsr()/x->freq)){
      x->count = 0;
      double sigma = fabs(s_sigma[i]);
      x->sigma = (sigma==0) ? x->sigma : sigma;
      x->rnd_x = (t_sample)(gsl_ran_gaussian(x->r,x->sigma));
    } else {
      x->count++;
    }
    // out
    x_out[i] = x->rnd_x;

  }

  return(w+6);
}


// dsp
void gaussian_dsp(t_gaussian *x, t_signal **sp){
  dsp_add(gaussian_perform, 5, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[2]->s_n);

}




// constructor
void * gaussian_new(t_symbol * s, int argc, t_atom * argv){
  t_gaussian * x = (t_gaussian *) pd_new(gaussian_class);
  // random generator
  x->r = gsl_rng_alloc(gsl_rng_taus2);
  // seed
  gsl_rng_set(x->r,0);

  // default sigma, normal
  double sigma = 1;

  //
  t_float f = 440;
  switch(argc){
  default:
  case 2:
    sigma  = atom_getfloat(argv+1);
  case 1:
    f = atom_getfloat(argv);
    f = (f==0)? 440 : f;
    break;
  }

  x->sigma = sigma;
  x->freq = f;
  x->count = 0;

  //sigma inlet
  x->sigma_in = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  //outlet
  x->x_out=outlet_new(&x->x_obj, &s_signal);

  return (void *) x;
}
// destroyer
void gaussian_free(t_gaussian  * x){
  //char info[80];
  //sprintf(info,"%s generator has been freed.", gsl_rng_name(x->r));
  //post(info);
  gsl_rng_free(x->r);
  inlet_free(x->sigma_in);
  outlet_free(x->x_out);
}

void gaussian_tilde_setup(void){
  gaussian_class = class_new(gensym("gaussian~"),
                            (t_newmethod) gaussian_new,
                           (t_method) gaussian_free,
                            sizeof(t_gaussian),
                            CLASS_DEFAULT, A_GIMME, 0);
  class_addmethod(gaussian_class,
                  (t_method)gaussian_dsp,
                  gensym("dsp"),0);
  // freq inlet
  CLASS_MAINSIGNALIN(gaussian_class, t_gaussian, f);

  class_sethelpsymbol(gaussian_class, gensym("gaussian~-help"));

}
