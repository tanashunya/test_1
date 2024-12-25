/* **********************************************************************
// 非対称行列クラス  GenMatrix 宣言部
// n×mの非対称マトリックス用 + OpenMP
// Ver. 3.10.0 2010.9.27  K. Watanabe
********************************************************************** */

#include "magic.h"

#ifndef ___GENMATRIX
#define ___GENMATRIX

#define ___DBGPRINTG 1  /*デバッグデータ出力  1:する 0:しない*/
#define ___DEFAULT_BW 10 /*初期バンド幅*/

#define ___OPENMP_SUPPORT 1 /* OpenMPの使用 1:する 0:しない */

#if ___OPENMP_SUPPORT
#include <omp.h>
#endif

using namespace std;

/*テンプレートが目障りなので*/
typedef complex<double> Complex;

class GenMatrix{
	
	friend class SymGenMethod;
	
	private:
		/* member */
//
		/* method */
		void newRegion(int n, int m);       /*領域確保  */
		void delRegion();                   /*領域削除  */
		static void copy(const GenMatrix* x, GenMatrix* y); /* x -> yへのコピー*/

	protected:
		double** value;   /* 値格納配列の各行の先頭アドレスポインタ*/
		int**    column;  /* 列番号格納配列の各行の先頭アドレスポインタ*/
		int* szColumn;    /* 各行の確保されている領域サイズ*/
		int* numNZ;       /* 各行の非零要素の数*/
		
		/*numNZ > szColumn となると自動的に領域を増やしてszColumnを増やす*/
		
		int  numRow;      /* 行の数   */
		int  numColumn;   /* 列の数   */
//		int  numMaxBW;    /* 非零要素数の最大値（最大バンド幅）*/

		int sortFlag;     /* 列データをソートしたかどうかのフラグ (1:ソート済み)*/
		
		int searchNumMaxBW () const; /*最大バンド幅を調べる*/
		
	public:
		GenMatrix();                     /*デフォルトコンストラクタ */
		GenMatrix(int n);                /*コンストラクタ(n×n行列) */
		GenMatrix(int n, int m);         /*コンストラクタ(n×m行列) */
		GenMatrix(GenMatrix& x);         /*コピーコンストラクタ     */
		~GenMatrix(void);                /*デストラクタ             */

		double set(int i, int j, double val); /* データセットA(i,j)  = val */
		double add(int i, int j, double val); /* データ追加  A(i,j) += val */
		double get(int i, int j);             /* データ取り出し*/
		
		void sortColumn();             /*列データのソート*/
		void removeZero(double eps);   /*絶対値がepsより小さい要素の除外*/
		void quickAllZero();           /*ゼロ行列にする（メモリ解放はしない）*/
		void clearZeroRow(int i);      /*i行目の成分をゼロにする*/
		
		int getNumRow()   const {return numRow;};       /*行の数を返却*/
		int getNumColumn()const {return numColumn;};    /*列の数を返却*/
		int getNumszColumn(int i)const {return szColumn[i];};    /*列の使用サイズを返却*/
//		int getNumMaxBW() const {return numMaxBW;};     /*最大バンド幅を返却*/

		int getNumTotalEntry() const; /*格納データ数を返却（CRS形式変換時の便利）*/
		
		/*AにBを加える*/
		void set(const GenMatrix& B, int offsetRow, int offsetColumn);

		/* ************************************************************************* */
		/* ***以下行列演算Method *** */
		
		/* 転置行列の作成*/
		static void transpose(const GenMatrix& A, GenMatrix& B);

		/* [行列]と{ベクトル}の積*/
		static void productAX(const GenMatrix& A, const double* x, double* ax );

		/* [行列]と{ベクトル}の積(複素ベクトル用)*/
		static void productAX(const GenMatrix& A, const Complex* x, Complex* ax );

		/* [行列]の転置行列と{ベクトル}の積*/
		static void productTranAX(const GenMatrix& A, const double* x, double* ax );

		/* [行列]の転置行列と{ベクトル}の積(複素ベクトル用)*/
		static void productTranAX(const GenMatrix& A, const Complex* x, Complex* ax );

		/* [行列]と[行列]の積*/
		static void productGG(GenMatrix& A, GenMatrix& B, GenMatrix& AB );

		/* [行列]と[行列の転置行列]の積*/
		static void productGGT(GenMatrix& A, GenMatrix& B, GenMatrix& AB );
		
		/* [行列]の転置行列と[行列]の積(Bの方が非零要素が少ない場合)*/
		static void productGtG(GenMatrix& A, GenMatrix& B, GenMatrix& AB );
		
		/* [行列]の転置行列と[行列]の積(Bの方が非零要素が少ない場合&&add value to AB)*/
		static void productGtG_(GenMatrix& A, GenMatrix& B, GenMatrix& AB );

		/* 逆行列を求める（正方行列）*/
		static int inv(GenMatrix& A);

		/* 対角スケーリング（正方行列）*/
		static void diagnalScaling(GenMatrix& A);
		static void diagnalScaling(GenMatrix& A, const int an);// for CEFC2010

		/* convert to 1-dim. array for CLAPACK & MKL */
		static void makeuppmat(const GenMatrix& A, double *up);

		/* convert to CSR format for Intel MKL */
		void convCSR(double* values, int* columns, int* rowIndex) const;

		/* 非圧縮の配列に格納*/
		void convFullMatrix(double* fullmat) const;


		/* 対角項を抽出し、1次元配列に格納(ただし、正方行列のみ)*/
		/* 正方行列では無い場合は例外OverIdxErrを投げる*/
		static void getDiagVector(const GenMatrix& A, double* diag );

		/*  指定したi行目を全てゼロにする*/
		void clearRow(int i);

		/*  有限要素法の境界条件設定用method(高速版)*/
		void setBoundaryCondition(double* b,       vector<int>& bclist, double  val);
		void setBoundaryCondition(double* b, const vector<int>& bclist, double* val);/*同上*/

		static int find_index(vector<int>& bclist, int val);/*検索関数*/

		/* Imcomplete Cholesky decomposition of A */
		static void icdcmp(const GenMatrix& A, GenMatrix& LD, double ga);

		/* solve [L][D][Lt]{x} = {b} for ICCG solver */
		static void icSolv(const GenMatrix& LD, const double* b, double* x, double* y);

		/* ICCG Method */
		static int iccgSolv(const GenMatrix& A, const GenMatrix& L, const double* b,
                      double* x, double* r, double eps, int numItr);

		/* Bi-CG STAB Solver with ILU precondition*/
		/* 係数行列が正方行列では無い場合は例外OverIdxErrを投げる*/
		/* 戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）*/
		static int biCGStabSolv(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, double ga, bool pivotFlag =false);

		/* Bi-CG STAB Solver with Diagnal precondition*/
		static int biCGStabSolvWithDiagP(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr);

		/* Bi-CG STAB Solver with GS precondition*/
		static int biCGStabSolvWithGS(GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr);

		/* GPBiCS Solver with ILU preconditioned*/
		static int GPBiCGSolv(const GenMatrix& A, const double* b, double* x,
                          double eps, int maxItr, double ga, bool pivotFlag=false);

		/* GPBiCS Solver with Diagnal preconditioned*/
		static int GPBiCGSolvWithDiagP(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr);
		

		/* Imcomplete LU decomposition of A */
		/* 不完全LU分解（Bi-CG Stab法で使用）*/
		static void iLUcmp(const GenMatrix& A, GenMatrix& LU, int* swapData, double ga,
		                   bool pivotFlag);

		/*不完全ILU分解された係数行列で連立方程式を解く（Bi-CG Stab法用）*/
		static void iLUSolv(const GenMatrix& LU, const double* b, double* x);
		static void iLUSolv(const GenMatrix& LU, const double* b, double* x,
		                    int* swapData);

		/* Jacobi smoother for multigrid */
		static double JacobiSmoother(const GenMatrix& A, const double* b, double* x,
                              double* r, int maxItr);

		/* Gauss-Seidel Solver */
		/* 係数行列が正方行列では無い場合は例外OverIdxErrを投げる*/
		/* 戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）*/
		static int gsSolv(const GenMatrix& A, const double* b, double* x,
                       double eps, int maxItr, bool increStopFlag =false, int incCount = 1);

		/* Gauss-Seidel type precondition */
		static void gsPrecondition(GenMatrix& A, const double* p, double* pt, double* q);

		/* IDR(s) method */
		static int IDR_sSolv(const GenMatrix& A, const double* b, double* x,
                        const int s, double eps, int maxItr);

		/* IDR(s) smoother */
		static double IDR_sSmoother(const GenMatrix& A, const double* b, double* x,
                          double* r, const int s, int maxItr);

		static int AdaptiveIDR_sJacobi(const GenMatrix& A, const double* b, double* x,
                        const int s, double eps, int maxItr);


		/* IDR-based SOR method */
		static int iSORSolv(const GenMatrix& A, const double* b, double* x,
                         double eps, int maxItr, double omega, bool increStopFlag = false);

		/* IDR-based SOR method with best X */
		static int iSORSolvBest(const GenMatrix& A, const double* b, double* x,
                         double eps, int maxItr, double omega);

		/* IDR-based adaptive Jacobi method */
		static int adaptiveIJacobi(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr);

		/* IDR-based adaptive Jacobi smoother for multigrid */
		static double adaptiveIJacobiSmoother(const GenMatrix& A, const double* b, double* x,
                              double* r, int maxItr);

		/* IDR-based relaxed Jacobi method (do not work correctly)*/
		static int relaxedIJacobi(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, double omega);

		/* restriction [R][A][R]t */
		static void restriction(const GenMatrix& A, int* R, int Rsize, GenMatrix& AS);

		/* restriction [P1]t[A][P2] */
		static void restriction(const GenMatrix& A, int* P1, int* P2,
                            int r1size, int r2size, GenMatrix& AS);

		/* restriction [R1][A][R2]t */
		static void restriction(const GenMatrix& aM, int* rM1, int r1Size,
                                                     int* rM2, int r2Size, GenMatrix& asM);

		static void restriction(int* rM, int rSize,  double* x,  double* xr);

		static void restriction(int* rM, int rSize, int* x, int* rx);
		
		static void prolongation_add(int* rM, int rRow,           double* x,  double* px);

		static void prolongation(int* rM, int rRow, int rColumn,  double* x,  double* px);

		/* matrix(compress form)-matrix multiplication */
		static void multipleRtA(const GenMatrix& A, const int rn, const int rm,
                        double* Rvalue, int* Rcolumn, GenMatrix& RtA);

		/* matrix-matrix(compress form) multiplication */
		static void multipleAR(const GenMatrix& A, const int rn, const int rm,
                        double* Rvalue, int* Rcolumn, GenMatrix& AR);

		/* modified Gram Schmidt */
		void modGramSchmidt();
		static void modGramSchmidt(double** A, const int n, const int m);

		static int GMRES_mSolvWithDiagP(GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, const int m);


		/* ************************************************************************* */
		/* *******以下ファイル入出力Method *********** */
		/* ファイルに出力(Ver1.0) Ver2.0推奨 */
		static int writeFile_V1(const char* filename, const GenMatrix& A);

		/* ファイルに出力(ver2.0) */
		static int writeFile(const char* filename, const GenMatrix& A);

		/* ファイルから取りこみ(ver1.0) */
		static int readFile_V1(const char* filename, GenMatrix& A);

		/* ファイルから取りこみ(Ver2.0) ver1.0にも対応 */
		static int readFile(const char* filename, GenMatrix& A);

		/* 1次元(double)配列をファイルに出力 */
		static int vectWriteFile(const char* filename, const double* x, int n);
		
		/* 1次元(double)配列をファイルから取りこみ */
		static int vectReadFile(const char* filename, double** x, int* n);

		
		/* ***例外送出用クラス*** */
		/* メモリエラー*/
		class MemErr{
			private:
				GenMatrix* ident;
			public:
				MemErr(){}
				MemErr(GenMatrix* p){}
		};
		/* 添字範囲エラー */
		class OverIdxErr{
			private:
				GenMatrix* ident;
			public:
				OverIdxErr(){}
				OverIdxErr(GenMatrix* p){}
		};
		/* I/Oエラー */
		class IOErr{
			public:
				IOErr(){}
		};

		/* ***演算子の多重定義*** */
		GenMatrix& operator=(const GenMatrix& x){/*代入演算子の多重定義*/
			if (&x != this){
				copy(&x, this);
			}
			return *this;
		}

		/* ************************************************************************* */
		/* *******テスト用などの特殊関数 *********** */
		static int gsSolv_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax);

		static int adaptiveIJacobi_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax);

		static int IDR_sSolv_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, int s, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax);


}; /* *** 宣言部終了 *** */

#endif

