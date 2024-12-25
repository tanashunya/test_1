
#include "cal.h"

void cal::Cal_Matrix( double **K, double *b, double *a, int num_nodes){

  // double **K = new double*[num_of_matrixelement];
  // for(int i=0;i<num_of_matrixelement;i++){
  //   K[i] = new double[num_of_matrixelement];
  // }
  // cout << "Cal_Matrix\n";
  // cout << K[0][0] << endl;
  // cout << K[1][1] << endl;
  // cout << b[0] << endl;
  // cout << a[0] << endl;









  // double Ktmp[2][2]={{1.0,-2.0},
  //                 {-2.0,6.0}};

  // double btmp[2]={-5.0,16.0};

  // for(int i=0;i<num_of_matrixelement;i++){
  //   for(int j=0;j<num_of_matrixelement;j++){
  //     K[i][j] = Ktmp[i][j];
  //   }
  // }
  // double *b = new double[num_of_matrixelement];
  // for(int i=0;num_of_matrixelement;i++){
  //  b[i] = btmp[i];
  // }



  // 解ベクトル
  // double *a = new double[num_of_matrixelement];
  // for(int i=0;i<num_of_matrixelement;i++){
  //   a[i] = 0.0;
  // }

  // CG法実行
  CG_Method( K , b , a, num_nodes);

  // 解ベクトルの表示
  // cout << "---- delta B ----" << endl;
  // for(int i=0;i<num_of_matrixelement;i++){
  //   cout << a[i] << endl;
  // }

  // メモリ開放
  // delete [] a;
  // delete [] b;

  // for(int i=0;i<num_of_matrixelement;i++){
  //   delete [] K[i];
  // }
  // delete[] K;
}


void cal::CG_Method(double **K, double *b, double *a, int num_nodes){

  // cout << "CG_Method\n";
  // cout << K[0][0] << endl;
  // cout << K[1][1] << endl;
  // cout << b[0] << endl;
  // cout << a[0] << endl;

  // cout << "CG_Method2\n";
  // cout << num_of_matrixelement << endl;
  // double *r = new double[num_of_matrixelement];
  // double *p = new double[num_of_matrixelement];
  // double *Kp = new double[num_of_matrixelement];
  // double *Ka = new double[num_of_matrixelement];

  double r[num_nodes];
  double p[num_nodes];
  double Kp[num_nodes];
  double Ka[num_nodes];


  // cout << "CG_Method3\n";

  for( int i=0; i<num_nodes;i++){
    // cout << i << endl;
    r[i] = 0.0;
    p[i] = 0.0;
    Kp[i] = 0.0;
    Ka[i] = 0.0;
    // cout << "Ka[" << i << "] = " << Ka[i] << endl;
  }

  // cout << "CG_Method\n";

  for( int i=0;i<num_nodes;i++){
    for( int j=0;j<num_nodes;j++){
      Ka[i] += K[i][j] * a[j];
    }
    r[i] = b[i] - Ka[i];
    p[i] = r[i];
  }

  double alpha, ke;
  int loop;
  for(loop=0;loop<10000;loop++){
    double pKp = 0;
    double rp = 0;
    double rKp = 0;


    for(int i=0;i<num_nodes;i++){
      Kp[i] = 0;
      rp += r[i] * p[i];

      for(int j=0;j<num_nodes;j++){
        Kp[i] += K[i][j] * p[j];
      }
      pKp += p[i] * Kp[i];
    }
  
    alpha = rp/pKp;

    // cout << "rp    = " << rp << endl;
    // cout << "pKp   = " << pKp << endl;
    // cout << "alpha = " << alpha << endl;

    double beta;
    for (int i=0;i<num_nodes;i++){
      a[i] = a[i] + alpha * p[i];
      r[i] = r[i] - alpha * Kp[i];
      rKp += r[i] * Kp[i];
    }
    beta = -rKp/pKp;

    for(int i=0;i<num_nodes;i++){
      p[i] = r[i] + beta * p[i];
    }

    double normR = 0.0;
    double normB = 0.0;
    for(int i=0;i<num_nodes;i++){
      normR = normR + r[i] * r[i];
      normB = normB + b[i] * b[i];
    }

    if(sqrt(normR) < sqrt(normB) * 10e-10) break;

    // cout << sqrt(normR) << endl;
    ke = loop;
  }
  // cout << ke <<endl;

  // delete[] r;
  // delete[] p;
  // delete[] Kp;
  // delete[] Ka;
  cout << "CG loop = " << loop+1 << endl;
}






