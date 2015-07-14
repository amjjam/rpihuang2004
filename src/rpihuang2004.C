#include "../include/rpihuang2004.H"

/*=============================================================================
  RPIHUANG2004() - constructor
  ============================================================================*/
RPIHUANG2004::RPIHUANG2004(){
  A=4833;
  B=3.64;
  C=0.2;
  D=0.03;
  gamma=-0.14;
  alpha=1.25;
}


/*=============================================================================
  ~RPIHUANG2004() - destructor
  ============================================================================*/
RPIHUANG2004::~RPIHUANG2004(){

}


/*=============================================================================
  double getDensity(double L, double lambda) - get the density, in
  cm^-3 as a function of L-shell and magnetic latitude of a point on
  the field line.
  
  double L - L-shell in Earth radii
  double lambda - magnetic latitude in degrees
  ============================================================================*/
double RPIHUANG2004::getDensity(double L, double lambda){
  double lambdainv=acos(1.0/sqrt(L))/M_PI*180;
  return getDensity(L,lambda,lambdainv);
}


/*=============================================================================
  double getDensity(double L, double lambda, double lambdainv) - get
  the densit, in cm^-3 as a function of L-shell and magnetic latitude
  of a point on the field line and invariant latitude of the field
  line. If the field line is dipolar then simply use the call form
  that does not include invariant latitude as it is then computed
  using a dipole model.

  double L - L-shell in Earth radii
  double lambda - magnetic latitude of point on field line in degrees
  double labdainv - invariant latitude of field line in degrees
  ============================================================================*/
double RPIHUANG2004::getDensity(double L, double lambda, double lambdainv){
  double N0=A*(B/L-1);
  double beta=C+D*L;
  
  
  return N0*(1+gamma*lambda/lambdainv)
    *pow(1.0/cos(M_PI/2*alpha*lambda/lambdainv),beta);
}
