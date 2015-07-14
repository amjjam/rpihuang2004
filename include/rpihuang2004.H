/******************************************************************************
 * This is the class RPIHUANG2004 which implements the Huang et al. 2004      *
 * RPI data based electron number density model of plasmaspheric plasma.      *
 ******************************************************************************/

#include <math.h>

class RPIHUANG2004{
public:
  RPIHUANG2004();
  ~RPIHUANG2004();
  double getDensity(double L, double lambda);
  double getDensity(double L, double lambda, double lambdainv);
private:
  double A,B,C,D,gamma,alpha;
};
