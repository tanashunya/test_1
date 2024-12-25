#include "cal.h"

void Cal( double **K, double *b, double *a, int num_nodes){
  // cout << "Cal\n";

  // cout << "K[0][0] = " << K[0][0] << endl;
  // cout << K[1][1] << endl;
  // cout << b[0] << endl;
  // cout << a[0] << endl;
  
  cal Cal;

  // cout << num_nodes << endl;
  
  Cal.setTotalNumber(num_nodes);

  Cal.Cal_Matrix(K,b,a,num_nodes);
  
}