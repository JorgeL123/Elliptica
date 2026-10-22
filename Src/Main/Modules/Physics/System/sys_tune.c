/*
// Alireza Rashti
// November 2020
*/


/* collection of functions to deal with tuning of the physical system */

#include "sys_tune.h"

/* tuning P_ADMs to be zero. */
int sys_tune_ADM_momenta(Physics_T *const phys)
{
  FUNC_TIC
  
  const char *const par = sysGets("P_ADM_control_method");
  char *adjust[3] = {0};
  double adm_p[3] = {0};
  double adm_j[3] = {0};
  
  printf(Pretty0"method        = %s\n",par);
  printf(Pretty0"update weight = %e\n",
        sysGetd("P_ADM_control_update_weight"));
  
  parse_adjust_parameter(par,adjust);
  
  void (*P_ADM_control_0)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[0]);
  
  void (*P_ADM_control_1)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[1]);
                              
  void (*P_ADM_control_2)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[2]);
  
  
  observe(phys,"ADM(P)",Gets("Observe_ADM_P"),adm_p);
  observe(phys,"ADM(J)",Gets("Observe_ADM_J"),adm_j);
  
  printf(Pretty0"Current P_ADM = (%e,%e,%e)\n",adm_p[0],adm_p[1],adm_p[2]);
  printf(Pretty0"Current J_ADM = (%e,%e,%e)\n",adm_j[0],adm_j[1],adm_j[2]);
  
  sysSetd("Px_ADM",adm_p[0]);
  sysSetd("Py_ADM",adm_p[1]);
  sysSetd("Pz_ADM",adm_p[2]);
  
  sysSetd("Jx_ADM",adm_j[0]);
  sysSetd("Jy_ADM",adm_j[1]);
  sysSetd("Jz_ADM",adm_j[2]);
  
  if (P_ADM_control_0)
    P_ADM_control_0(phys);
  
  if (P_ADM_control_1)
    P_ADM_control_1(phys);
    
  if (P_ADM_control_2)
    P_ADM_control_2(phys);
  
  Jz_ADM_control(phys);
  E_ADM_control(phys);
    
  Free(adjust[0]);
  Free(adjust[1]);
  Free(adjust[2]);
  
  FUNC_TOC
  
  return EXIT_SUCCESS;
}

/* getting adjustment str, returns the relevant function. */
static fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust)
{
  fAdjustment_t *f = 0;
  
  if (!adjust)
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"none"))
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"x_CM"))
  {
    f = Py_ADM_is0_by_x_CM;
  }
  else if (strcmp_i(adjust,"y_CM"))
  {
    f = Px_ADM_is0_by_y_CM;
  }
  else
    Error0(NO_OPTION);
  
  return f;
}

/* parsing adjust parameter value and fill components consequently */
static void parse_adjust_parameter(const char *const par,char *adjust[3])
{
  if (!strstr_i(par,"adjust(") && !strstr_i(par,"none"))
    Errors("Syntax error for '%s'.\n",par);
  
  /* if it is none */  
  if (strcmp_i(par,"none"))
  {
    adjust[0] = 0;
    adjust[1] = 0;
    adjust[2] = 0;
    
    return;
  }
  
  /* parse if not none */
  char *str = dup_s(par);
  char *save;
  char *main_str = sub_s(str,'(',')',&save);
  char *tok  = tok_s(main_str,',',&save);
  
  adjust[0] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[1] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[2] = dup_s(tok);
  
  free(str);
}

/* find y_CM by demanding Px_ADM = 0 */
static void Px_ADM_is0_by_y_CM(Physics_T *const phys)
{
  double dy_CM = 0,y_CM_new;
  const double px    = sysGetd("Px_ADM");
  const double py    = sysGetd("Py_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double P_thr = sysGetd("P_ADM_control_threshold");
  const double Omega = sysGetd("angular_velocity");
  const double a     = infall_rate(phys);
  const double y_CM0 = sysGetd("y_CM");
  const double admM  = sysGetd("ADM_mass");
  const double P_err = (a == 0.) ? fabs(px) : sqrt(px*px+py*py);
  
  printf(Pretty0"adjusting Px_ADM by y_CM ...\n");
  
  /* changing center of mass */
  dy_CM    = (a*py-Omega*px)/(admM*(a*a+Omega*Omega));
  y_CM_new = y_CM0+dy_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(P_err,dP) && LSS(P_err,P_thr))
  {
    sysSetd("y_CM",y_CM_new);
  }
}
  
/* find x_CM by demanding Py_ADM = 0 */
static void Py_ADM_is0_by_x_CM(Physics_T *const phys)
{
  double  dx_CM = 0,x_CM_new;
  const double px    = sysGetd("Px_ADM");
  const double py    = sysGetd("Py_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double P_thr = sysGetd("P_ADM_control_threshold");
  const double Omega = sysGetd("angular_velocity");
  const double a     = infall_rate(phys);
  const double x_CM0 = sysGetd("x_CM");
  const double admM  = sysGetd("ADM_mass");
  const double P_err = (a == 0.) ? fabs(py) : sqrt(px*px+py*py);
  
  printf(Pretty0"adjusting Py_ADM by x_CM ...\n");
  
  /* changing center of mass */
  dx_CM    = (a*px+Omega*py)/(admM*(a*a+Omega*Omega));
  x_CM_new = x_CM0+dx_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(P_err,dP) && LSS(P_err,P_thr))
  {
    sysSetd("x_CM",x_CM_new);
  }
}

static double infall_rate(Physics_T *const phys)
{
  double Vr = 0., D = 0.;
  
  sprintf(phys->par,"%s_%s",phys->ssys,"infall_velocity");
  if (PgetdEZ(phys->par) != DBL_MAX) Vr = PgetdEZ(phys->par);
  sprintf(phys->par,"%s_%s",phys->ssys,"separation");
  if (PgetdEZ(phys->par) != DBL_MAX) D = PgetdEZ(phys->par);
  
  return (D > 0.) ? Vr/D : 0.;
}



/* double parameter of the system with a default value if it is not set */
static double sys_par_d(Physics_T *const phys,const char *const name,
                        const double def)
{
  sprintf(phys->par,"%s_%s",phys->ssys,name);
  const double v = PgetdEZ(phys->par);
  
  return (v == DBL_MAX) ? def : v;
}

/* the J_ADM and E_ADM controls are meant for orbits set by hand.
// for quasi-circular data Elliptica sets Omega itself via the force
// balance equation. */
static void assert_orbit_fixed_by_hand(const char *const control)
{
  const char *const fb_pars[] = {"NS_force_balance_equation",
                                 "NS1_force_balance_equation",
                                 "NS2_force_balance_equation",0};
  Uint i;
  
  for (i = 0; fb_pars[i]; ++i)
  {
    const char *const fb = PgetsEZ(fb_pars[i]);
    if (fb && !strcmp_i(fb,"none"))
    {
      char msg[1000];
      sprintf(msg,"%s needs the orbit set by hand, but %s = %s.\n"
             "These controls are for scattering data; for quasi-circular data\n"
             "the force balance sets Omega. Set %s = none or the control to none.\n",
             control,fb_pars[i],fb,fb_pars[i]);
      Errors("%s",msg);
    }
  }
}

/* ->: 1 if a J/E control should update this iteration, otherwise 0
// after printing the reason (weight 0, within tolerance, above threshold). */
static int control_should_update(const char *const control,const double W,
                                 const double rel_err,const double tol,
                                 const double thr)
{
  if (W == 0.)
  {
    printf(Pretty0"%s control: no update, weight = 0.\n",control);
    return 0;
  }
  if (!(fabs(rel_err) > tol))
  {
    printf(Pretty0"%s control: no update, |rel. error| = %.6e <= "
                  "tolerance = %.6e.\n",control,fabs(rel_err),tol);
    return 0;
  }
  if (!(fabs(rel_err) < thr))
  {
    printf(Pretty0"%s control: no update, |rel. error| = %.6e >= "
                  "threshold = %.6e.\n",control,fabs(rel_err),thr);
    return 0;
  }

  return 1;
}

/* drive Jz_ADM to the target "Jz_ADM_target" by adjusting the angular
// velocity (for runs where Omega is not fixed by force balance. 
// J_ADM ~ mu D^2 Omega, hence:
// Omega_new = Omega*(1 + W*(Jz_target/Jz - 1)),
// with the relative change limited to "J_ADM_control_max_step".
// parameters:
// J_ADM_control_method       : none (default) | adjust(angular_velocity)
// Jz_ADM_target              : target Jz_ADM
// J_ADM_control_update_weight: W (default 0, can be iterative)
// J_ADM_control_tolerance    : no update if |Jz/Jz_target-1| < tol (default 1E-5)
// J_ADM_control_threshold    : no update if |Jz/Jz_target-1| > thr (default 0.5)
// J_ADM_control_max_step     : max |dOmega/Omega| per iteration (default 0.01) */
static void Jz_ADM_control(Physics_T *const phys)
{
  sprintf(phys->par,"%s_%s",phys->ssys,"J_ADM_control_method");
  const char *const method = PgetsEZ(phys->par);
  
  if (!method || strcmp_i(method,"none"))
    return;
  
  if (!strcmp_i(method,"adjust(angular_velocity)"))
    Errors("Unknown option '%s' for J_ADM_control_method.\n",method);
  
  assert_orbit_fixed_by_hand("J_ADM_control_method");
  
  const double Jt    = sys_par_d(phys,"Jz_ADM_target",DBL_MAX);
  const double W     = sys_par_d(phys,"J_ADM_control_update_weight",0.);
  const double tol   = sys_par_d(phys,"J_ADM_control_tolerance",1E-5);
  const double thr   = sys_par_d(phys,"J_ADM_control_threshold",0.5);
  const double dmax  = sys_par_d(phys,"J_ADM_control_max_step",0.01);
  const double Jz    = sysGetd("Jz_ADM");
  const double Omega = sysGetd("angular_velocity");
  
  if (Jt == DBL_MAX)
    Error0("J_ADM_control_method is set but Jz_ADM_target is not.\n");
  if (Jt == 0.)
    Error0("Jz_ADM_target = 0: the relative error is undefined.\n");

  if (Jz == 0. || Omega == 0.)
  {
    printf(Pretty0"J_ADM control: SKIPPED, Jz = %.10e, Omega = %.10e "
                  "(Omega cannot be rescaled).\n",Jz,Omega);
    return;
  }

  const double rel_err = Jz/Jt-1.;

  printf(Pretty0"J_ADM control: Jz = %.10e, target = %.10e, rel. error = %.6e, "
                "weight = %e\n",Jz,Jt,rel_err,W);

  if (!control_should_update("J_ADM",W,rel_err,tol,thr))
    return;

  double step = W*(Jt/Jz-1.);
  if (step >  dmax) step =  dmax;
  if (step < -dmax) step = -dmax;

  const double Omega_new = Omega*(1.+step);
  printf(Pretty0"adjusting Jz_ADM by angular velocity: %.10e -> %.10e\n",
         Omega,Omega_new);
  
  sysSetd("angular_velocity",Omega_new);
}

/* drive E_ADM to the target "E_ADM_target" by adjusting the infall 
// (radial) velocity v_r (for orbits set by hand).
// E_ADM is measured with "E_ADM_control_observe" (default S_inf,default;
// the default S+V,conformal is biased in non-quasi-circular runs). 
// Newtonian energy, E ~ M + mu v^2/2 - mu M/D, gives
// dE/d(v_r^2) = mu/2, hence:
// v_r^2 <- v_r^2 + 2 W (E_target - E)/mu, keeping the sign of v_r
// (infall if v_r = 0), with mu = E q/(1+q)^2 from "mass_ratio" q and 
// |dv_r| limited to max_step*max(|v_r|,1E-3).
// parameters:
// E_ADM_control_method       : none (default) | adjust(infall_velocity)
// E_ADM_target               : target E_ADM
// E_ADM_control_observe      : ADM(M) method (default S_inf,default)
// E_ADM_control_update_weight: W (default 0, can be iterative)
// E_ADM_control_tolerance    : no update if |E/E_target-1| < tol (default 1E-6)
// E_ADM_control_threshold    : no update if |E/E_target-1| > thr (default 1E-2)
// E_ADM_control_max_step     : max |dv_r|/|v_r| per iteration (default 0.02) */
static void E_ADM_control(Physics_T *const phys)
{
  sprintf(phys->par,"%s_%s",phys->ssys,"E_ADM_control_method");
  const char *const method = PgetsEZ(phys->par);
  
  if (!method || strcmp_i(method,"none"))
    return;
  
  if (!strcmp_i(method,"adjust(infall_velocity)"))
    Errors("Unknown option '%s' for E_ADM_control_method.\n",method);
  
  assert_orbit_fixed_by_hand("E_ADM_control_method");
  
  const double Et   = sys_par_d(phys,"E_ADM_target",DBL_MAX);
  const double W    = sys_par_d(phys,"E_ADM_control_update_weight",0.);
  const double tol  = sys_par_d(phys,"E_ADM_control_tolerance",1E-6);
  const double thr  = sys_par_d(phys,"E_ADM_control_threshold",1E-2);
  const double dmax = sys_par_d(phys,"E_ADM_control_max_step",0.02);
  const double q    = sys_par_d(phys,"mass_ratio",DBL_MAX);
  const double vr   = sysGetd("infall_velocity");
  char obs_method[400] = {'\0'};
  double E = 0.;
  
  if (Et == DBL_MAX)
    Error0("E_ADM_control_method is set but E_ADM_target is not.\n");
  if (Et == 0.)
    Error0("E_ADM_target = 0: the relative error is undefined.\n");
  if (q == DBL_MAX || !(q > 0.))
    Error0("E_ADM_control needs the mass ratio (<sys>_mass_ratio).\n");
  
  sprintf(phys->par,"%s_%s",phys->ssys,"E_ADM_control_observe");
  sprintf(obs_method,"%s",
          PgetsEZ(phys->par) ? PgetsEZ(phys->par) : "S_inf,default");
  
  observe(phys,"ADM(M)",obs_method,&E);
  
  const double rel_err = E/Et-1.;
  const double mu      = E*q/((1.+q)*(1.+q));
  
  printf(Pretty0"E_ADM control: E(%s) = %.10e, target = %.10e, rel. error = %.6e, "
                "weight = %e, v_r = %.10e\n",obs_method,E,Et,rel_err,W,vr);

  if (!(E > 0.))
  {
    printf(Pretty0"E_ADM control: SKIPPED, E = %.10e is not positive.\n",E);
    return;
  }

  if (!control_should_update("E_ADM",W,rel_err,tol,thr))
    return;
  
  double vr2 = vr*vr + 2.*W*(Et-E)/mu;
  if (vr2 < 0.) vr2 = 0.;
  
  double vr_new = (vr > 0. ? 1. : -1.)*sqrt(vr2);
  const double dv_max = dmax*(fabs(vr) > 1E-3 ? fabs(vr) : 1E-3);
  if (vr_new - vr >  dv_max) vr_new = vr + dv_max;
  if (vr_new - vr < -dv_max) vr_new = vr - dv_max;
  
  printf(Pretty0"adjusting E_ADM by infall velocity: %.10e -> %.10e\n",
         vr,vr_new);
  
  sysSetd("infall_velocity",vr_new);
}
