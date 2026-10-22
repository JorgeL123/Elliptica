#include "sys_header.h"


typedef void fAdjustment_t (Physics_T *const phys);


int sys_tune_ADM_momenta(Physics_T *const phys);
static fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust);
static void parse_adjust_parameter(const char *const par,char *adjust[3]);
static void Px_ADM_is0_by_y_CM(Physics_T *const phys);
static void Py_ADM_is0_by_x_CM(Physics_T *const phys);
static double infall_rate(Physics_T *const phys);
static void Jz_ADM_control(Physics_T *const phys);
static void E_ADM_control(Physics_T *const phys);
static void assert_orbit_fixed_by_hand(const char *const control);
static int control_should_update(const char *const control,const double W,
                                 const double rel_err,const double tol,
                                 const double thr);
static double sys_par_d(Physics_T *const phys,const char *const name,
                        const double def);




