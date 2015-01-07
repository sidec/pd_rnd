#include "m_pd.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static t_class * pdf_class;

typedef struct _pdf {
  t_object x_obj;
  t_int len;
  t_outlet * l_out;
} t_pdf;


void pdf_gaussian(t_pdf * x){
  int i;
  double a;
  double step = 6.0/x->len;
  double j = -3.0;
  for(i = 0; i < x->len; i++){
    a = gsl_ran_gaussian_pdf (j, 1.0);
    t_atom  ai[2];
    SETFLOAT(ai, (float)a);
    SETFLOAT(ai+1, (float)i);
    outlet_list(x->l_out, gensym("list"), 2, ai );
    j += step;
  }
}

void pdf_gamma(t_pdf * x, t_floatarg a, t_floatarg b){
  double kappa = (a==0)? 2.0 : (double)a;
  double theta = (b==0)? 2.0 : (double)b;
  double g;
  double step = 20.0/x->len;
  double j = 0.0;
  int i;
  for(i = 0; i < x->len; i++){
    g = gsl_ran_gamma_pdf (j, kappa,theta);
    t_atom ai[2];
    SETFLOAT(ai, (float)g);
    SETFLOAT(ai+1, (float)i);
    outlet_list(x->l_out, gensym("list"), 2, ai );
    j += step;
  }
}



// constructor
void * pdf_new(t_floatarg len){
  t_pdf * x = (t_pdf *) pd_new(pdf_class);
  x->len = len;
  x->l_out = outlet_new(&x->x_obj, &s_list);
  return (void *) x;
}

void pdf_setup(void){
  pdf_class = class_new(gensym("pdf"),
                        (t_newmethod) pdf_new,
                        0,
                        sizeof(t_pdf),
                        CLASS_DEFAULT,
                        A_DEFFLOAT, 0);
  class_addmethod(pdf_class,(t_method)pdf_gaussian,gensym("gaussian"),0);
  class_addmethod(pdf_class,
                  (t_method)pdf_gamma, gensym("gamma"),
                  A_DEFFLOAT, A_DEFFLOAT, 0);

    // help file
  class_sethelpsymbol(pdf_class, gensym("pdf-help"));

}
