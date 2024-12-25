#include "magic.h"
#include "const.h"

class cal{

  private:
    int TotalNumber;
    void output(double *x);

  public:
    void Input();
    void Cal_Matrix( double **K, double *b, double *a, int num_nodes );
    void CG_Method( double **K, double *b, double *a, int num_nodes );

    void setTotalNumber(int num){
      TotalNumber = num;
    }
    int getTotalNumber(){
      return TotalNumber;
    }

};