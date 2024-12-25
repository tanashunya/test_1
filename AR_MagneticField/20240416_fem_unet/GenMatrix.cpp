/* **********************************************************************
// 非対称行列クラス  GenMatrix 本体
// n×mの非対称マトリックス用
// Ver. 3.10.0 2010.9.27  K. Watanabe
 ********************************************************************** */

#include "GenMatrix.h"
using namespace std;


/* ************ private and protected Method  *************** */
/* ********************************************************** */

/* ********* *領域確保* ********* */
/* ****************************** */
void GenMatrix::newRegion(int n, int nbw){ 
	int i;
	try{
		value    = new double* [n];
		column   = new    int* [n];
		szColumn = new    int  [n];
		numNZ    = new    int  [n];
		
		for (i = 0; i < n; i++){ 
			value [i] = new double[nbw]; /*各行の領域を確保*/
			column[i] = new    int[nbw]; /*各行の領域を確保*/
		}
	}
	catch(bad_alloc) {
		throw MemErr(this);
	}
}

/* ********* *領域削除* ********* */
/* ****************************** */
void GenMatrix::delRegion(){
	int i;
	for (i = 0; i < numRow; i++){
		delete[]  value[i];
		delete[] column[i];
	}
	delete[] value;
	delete[] column;
	delete[] szColumn;
	delete[] numNZ;
}

/* ********* *データコピー* ********* */
/* ****************************** */
void GenMatrix::copy(const GenMatrix* x, GenMatrix* y){ /* x -> yへのコピー*/
	int i, mu, n;
	try{
		(*y).delRegion();     /*コピー先の領域を一旦削除*/
		n = (*y).numRow = (*x).numRow;
		(*y).numColumn  = (*x).numColumn;
//		(*y).numMaxBW   = (*x).numMaxBW;
		
		(*y).sortFlag = (*x).sortFlag;
		
		/*領域確保*/
		(*y).value    = new double* [n];
		(*y).column   = new    int* [n];
		(*y).szColumn = new    int  [n];
		(*y).numNZ    = new    int  [n];
		
		/* **データコピー** */
		for (i = 0; i < n; i++){
			(*y).value [i] = new double[(*x).szColumn[i]];
			(*y).column[i] = new    int[(*x).szColumn[i]];
			
			for(mu = 0; mu < (*x).numNZ[i]; mu++){
				(*y).value [i][mu] = (*x).value [i][mu];
				(*y).column[i][mu] = (*x).column[i][mu];
			}
			(*y).szColumn[i] = (*x).szColumn[i];
			(*y).numNZ[i]    = (*x).numNZ[i];
		}
	}
	catch(bad_alloc) {
		throw MemErr(y);
	}
}

/* ************* constructor and destructor ***************** */
/* ********************************************************** */

/* ***デフォルトコンストラクタ*** */
/* ****************************** */
GenMatrix::GenMatrix(){
	int i, nbw;
	
	numRow    = 1;
	numColumn = 1;
//	numMaxBW  = 0;
	sortFlag  = 0;

	nbw = ___DEFAULT_BW;

	try{
		newRegion(numRow, nbw); /*領域確保*/
		for (i = 0; i < numRow; i++){ 
			numNZ[i]    = 0; /*ゼロに初期化*/
			szColumn[i] = nbw;
		}
	}
	catch(MemErr) {
		throw MemErr(this);
		
	}
}

/* ***コンストラクタ(引数2個)*** */
/*   n×mの非正方行列生成        */
/* ***************************** */
GenMatrix::GenMatrix(int n, int m){

	int i, nbw;
	
	numRow    = n;
	numColumn = m;
//	numMaxBW  = 0;
	sortFlag  = 0;

	if (m < ___DEFAULT_BW){
		nbw = m;
	}
	else{
		nbw = ___DEFAULT_BW;
	}

	try{
		newRegion(n, nbw); /*領域確保*/
		for (i = 0; i < n; i++){ 
			numNZ[i]  = 0; /*ゼロに初期化*/
			szColumn[i] = nbw;
		}
	}
	catch(MemErr) {
		throw MemErr(this);
		
	}
}

/* ***コンストラクタ(引数1個)*** */
/*   n×nの正方行列生成          */
/* ***************************** */
GenMatrix::GenMatrix(int n){

	int i, nbw;
	
	numRow    = n;
	numColumn = n;
//	numMaxBW  = 0;
	sortFlag  = 0;

	if (n < ___DEFAULT_BW){
		nbw = n;
	}
	else{
		nbw = ___DEFAULT_BW;
	}

	try{
		newRegion(n, nbw); /*領域確保*/
		for (i = 0; i < n; i++){ 
			numNZ[i]  = 0; /*ゼロに初期化*/
			szColumn[i] = nbw;
		}
	}
	catch(MemErr) {
		throw MemErr(this);
		
	}
}


/* ***コピーコンストラクタ*** */
/* ************************** */

GenMatrix::GenMatrix(GenMatrix& x){ 

	numRow    = 1;
	numColumn = 1;
//	numMaxBW = 0;
	sortFlag = 0;

	try{
		newRegion(1, 1);  /*領域仮確保(下のcopyで改めて領域確保される)*/
		copy(&x, this);   /*データコピー*/
	}
	catch(MemErr) {
		throw MemErr(this);
		
	}
}


/* ***デストラクタ*** */
/* ****************** */
GenMatrix::~GenMatrix(void){ 
	delRegion();
}

/* ******************* public Method  *********************** */
/* ********************************************************** */

/* *** *データセット* *** */
/* * A[i,j] = val を行う。既存のA[i,j]のデータはvalに上書きされる */
/* * 処理速度最優先のため、エラーチェックは最低限しか行わない     */
/* * i,jが添え字範囲を超えている場合は例外OverIdxErrを投げる      */
/* * 戻り値はvalである                                            */
/* * Ver.1.11 K. Watanabe 2001 10.7                               */
/* ************************************************************* */
double GenMatrix::set(int i, int j, double val){
	
	int mu, nu, numnz;
	double* tmpVal;
	int* tmpCom;
	
	if (i >= numRow || j >= numColumn){
		throw OverIdxErr(this);
	}
	
	/* a[i,j]が既にあるかどうかチェック あれば上書き */
	numnz = numNZ[i];
	for (mu =0; mu < numnz; mu++){
		if (column[i][mu] == j){
			value[i][mu] = val;
			return val;
		}
	}
	/* a[i,j]を新規登録 */
	/*確保容量チェック*/
	if ( mu >= szColumn[i]){ /*おーばー*/
		tmpVal = new double[szColumn[i] + ___DEFAULT_BW];
		tmpCom = new    int[szColumn[i] + ___DEFAULT_BW];
		for (nu = 0; nu < numNZ[i]; nu++){
			tmpVal[nu] =  value[i][nu];
			tmpCom[nu] = column[i][nu];
		}
		delete[]  value[i];
		delete[] column[i];
		
		 value[i] = tmpVal;
		column[i] = tmpCom;
		szColumn[i] += ___DEFAULT_BW;

	}
	column[i][mu] = j;
	 value[i][mu] = val;
	numNZ[i]++;
//	if(numMaxBW < numNZ[i]){
//		numMaxBW = numNZ[i];
//	}
	sortFlag = 0;
	
	return val;
}


/* *** *データセット* *** */
/* * A[i,j] <= B[i+offsetRow, j+offsetColumn] を行う。            */
/* * 既存のA[i,j]のデータは上書きされる                           */
/* * i,jが添え字範囲を超えている場合は例外OverIdxErrを投げる      */
/* * 戻り値はvalである                                            */
/* * Ver.1.10 K. Watanabe 2010.5.3                                */
/* ************************************************************** */
void GenMatrix::set(const GenMatrix& B, int offsetRow, int offsetColumn){
//
//               offsetColumn
//           ____:________
//          |             |
//offsetRow ..   _______  |
//          |    |      | |
//          |    |   B  | |
//          |    |______| |
//          |_____________|
//	
	try{
	const int nb = B.numRow;
	const int mb = B.numColumn;
	
	if (numRow < (nb + offsetRow) || numColumn < (mb + offsetColumn)){
		throw OverIdxErr(this);
	}
	
	for(int ib = 0; ib < nb; ib++){
		int ia = ib + offsetRow;
		int numNZA =   numNZ[ia];
		int numNZB = B.numNZ[ib];

		/* remove A entry in B region */
		int numNZAtmp = 0;
		for(int nu = 0; nu < numNZA; nu++){
			int j = column[ia][nu];
			if(j < offsetColumn || j >= offsetColumn + mb){
				column[ia][numNZAtmp] = j;
				value [ia][numNZAtmp] = value[ia][nu];
				numNZAtmp++;
			}
		}
		numNZA = numNZAtmp;
			
		int numNZ_NEW = numNZA + numNZB;

		/* set B entry into A */
		/*check size */
		if(numNZ_NEW >= szColumn[ia]){
			double* tmpVal = new double[numNZ_NEW];
			int*    tmpCom = new    int[numNZ_NEW];
			for(int nu = 0; nu < numNZA; nu++){
				tmpVal[nu] =  value[ia][nu];
				tmpCom[nu] = column[ia][nu];
			}
			delete[]  value[ia];
			delete[] column[ia];
			
			 value[ia] = tmpVal;
			column[ia] = tmpCom;
			szColumn[ia] = numNZ_NEW;
		}
		for(int nu = 0; nu < numNZB; nu++){
			 value[ia][numNZA + nu] = B.value [ib][nu];
			column[ia][numNZA + nu] = B.column[ib][nu] + offsetColumn;
		}
		numNZ[ia] = numNZ_NEW;
	}
	sortFlag = 0;

	return;
	}
	catch(...){
		cout << "error at set" << endl;
		throw;
	}
}

/* *** *データ追加* *** */
/* * A[i,j] = A[i,j] + val を行う。                              */
/* * 既存のA[i,j]が無ければ新規登録される                        */
/* * 処理速度最優先のため、エラーチェックは最低限しか行わない    */
/* * i,jが添え字範囲を超えている場合は例外OverIdxErrを投げる      */
/* * 戻り値はvalである                                           */
/* * Ver.1.11 K. Watanabe 2001 10.7                              */
/* ************************************************************* */
double GenMatrix::add(int i, int j, double val){
	int mu, nu, numnz;
	double* tmpVal;
	int* tmpCom;
	
	if (i >= numRow || j >= numColumn){
		throw OverIdxErr(this);
	}
	
	/* a[i,j]が既にあるかどうかチェック あれば値を追加 */
	numnz = numNZ[i];
	for (mu =0; mu < numnz; mu++){
		if (column[i][mu] == j){
			value[i][mu] += val;
			return val;
		}
	}
	/* a[i,j]を新規登録 */
	/*確保容量チェック*/
	if ( mu >= szColumn[i]){ /*おーばー*/
		tmpVal = new double[szColumn[i] + ___DEFAULT_BW];
		tmpCom = new    int[szColumn[i] + ___DEFAULT_BW];
		for (nu = 0; nu < numNZ[i]; nu++){
			tmpVal[nu] =  value[i][nu];
			tmpCom[nu] = column[i][nu];
		}
		delete[]  value[i];
		delete[] column[i];
		
		value[i]  = tmpVal;
		column[i] = tmpCom;
		szColumn[i] += ___DEFAULT_BW;
		
	}
	column[i][mu] = j;
	value[i][mu] = val;
	numNZ[i]++;
//	if(numMaxBW < numNZ[i]){
//		numMaxBW = numNZ[i];
//	}
	sortFlag = 0;
	
	return val;
}


/* *** *データ取りだし* *** */
/* * A[i,j] を返却する。                                         */
/* * 処理速度最優先のため、エラーチェックは最低限しか行わない    */
/* * ver.1.1 K. Watanabe 2001 3.1                                */
/* ************************************************************* */
double GenMatrix::get(int i, int j){
	int mu;

	/* a[i,j]が既にあるかどうかチェック */
	for (mu =0; mu < numNZ[i]; mu++){
		if (column[i][mu] == j){
			return value[i][mu];
		}
	}
	/* 登録されていないので零を返却 */
	return 0.0;
}

/* * 列データのソート   */
/* * column[i] が昇順になるように並べ替える                      */
/* * 処理速度最優先のため、エラーチェックは最低限しか行わない    */
/* * ver.1.0 K. Watanabe 2001 3.1                                */
/* ************************************************************* */
void GenMatrix::sortColumn(){
	int i, mu, k;
	int    idxtmp;
	double valtmp;
	
	for (i = 0; i < numRow; i++){
		for (mu = 0; mu < numNZ[i]; mu++){
			for (k = mu + 1; k < numNZ[i]; k++){
				if (column[i][mu] > column[i][k]){
					idxtmp = column[i][mu];
					valtmp = value[i][mu];
					column[i][mu] = column[i][k];
					 value[i][mu] = value[i][k];
					column[i][k] = idxtmp;
					 value[i][k] = valtmp;
				}
			}
		}
	}
	sortFlag = 1; /* ソート済みであることを示すフラグを立てる */
	return;
}

/* * 絶対値がepsより小さい要素の除外     */
/* * |value[i][j]| <= eps である要素を除外                      */
/* * 例えば A.set(i,j,3.0) -> A.add(i,j,-3.0) をするとA(i,j)は   */
/* * ゼロにもかかわらず登録が残ったままになりメモリを消費する    */
/* * このように計算過程でゼロになった要素を除外する              */
/* * なお、この処理には時間がかかるので頻繁に呼び出さないように. */
/* * 処理速度最優先のため、エラーチェックは最低限しか行わない    */
/* * ワーク領域としてsizeof(double+int)*numMaxBW のメモリ使う    */
/* * この領域の確保に失敗した場合は例外MemErrを投げる            */
/* * ver.1.1 K. Watanabe 2001 9.26                               */
/* ************************************************************* */
void GenMatrix::removeZero(double eps){
	int i, mu;
	int cont;
	
	try{
		int numMaxBW = searchNumMaxBW();
		double* tmpVal = new double[numMaxBW];
		int*    tmpIdx = new int[numMaxBW];
		
		for (i = 0; i < numRow; i++){
			cont = 0;
			for (mu = 0; mu < numNZ[i]; mu++){
				if(fabs(value[i][mu]) > eps){
					tmpVal[cont] = value[i][mu];
					tmpIdx[cont] = column[i][mu];
					cont++;
				}
			}
			/* write back */
			for(mu = 0; mu < cont; mu++){
				 value[i][mu] = tmpVal[mu];
				column[i][mu] = tmpIdx[mu];
			}
			numNZ[i] = cont;
		}
		delete[] tmpVal;
		delete[] tmpIdx;
		
	}
	catch(bad_alloc){
		throw MemErr(this);
	}
}


/* * 零行列にする     */
/* * メモリの解放をせずに行列の全成分をゼロにする                */
/* * 具体的には numNZ[i] をゼロにするだけ                        */
/* * ver.1.0 K. Watanabe 2003 12.9                               */
/* ************************************************************* */
void GenMatrix::quickAllZero(){
	int i;
	
	for (i = 0; i < numRow; i++){
		numNZ[i] = 0;
	}
//	numMaxBW = 0;
}

/* * メモリの解放をせずにi行目の成分をゼロにする                 */
/* * 具体的には numNZ[i] をゼロにするだけ                        */
/* * ver.1.0 K. Watanabe 2008.5.3                                */
/* ************************************************************* */
void GenMatrix::clearZeroRow(int i){
	if(i < 0 || i >= numRow){
		throw OverIdxErr(this);
	}
	
	numNZ[i] = 0;
}

/* * 格納データ数を返却（CSR形式変換時の便利)                    */
/* * ver.1.0 K. Watanabe 2010.4.9                                */
/* ************************************************************* */
int GenMatrix::getNumTotalEntry() const{
	int size = 0;
	for(int i = 0; i < numRow; i++){
		size += numNZ[i];
	}
	return size;
}


/* * CSR格納データを作成（for Intel MKL zero-based indexing)     */
/* * ver.1.0 K. Watanabe 2010.4.13                               */
/* ************************************************************* */
void GenMatrix::convCSR(double* values, int* columns, int* rowIndex) const{
	
	try{
		int pos = 0;
		for(int i = 0; i < numRow; i++){
			int nNZ = numNZ[i];
			double*  valueP =  value[i];
			int*    columnP = column[i];

			rowIndex[i] = pos;
			for(int nu = 0; nu < nNZ; nu++){
				 values[pos] =  valueP[nu];
				columns[pos] = columnP[nu];
				pos++;
			}
		}
		rowIndex[numRow] = pos;// equal to the number of non-zeros
	
		return;
	}
	catch(...){
		cout << "error inconvCSR" << endl;
		throw;
	}
}

/* * 非圧縮データを作成                                          */
/* * ver.1.0 K. Watanabe 2010.4.13                               */
/* ************************************************************* */
void GenMatrix::convFullMatrix(double* fullmat) const{
	const int size = numRow * numColumn;
	for(int i = 0; i < size; i++){
		fullmat[i] = 0.0;
	}
	
	for(int i = 0; i < numRow; i++){
		for(int nu = 0; nu < numNZ[i]; nu++){
			fullmat[i*numColumn + column[i][nu]] = value[i][nu];
		}
	}
}


/* * 転置行列の作成                                              */
/* * [B] = [A]t                                                  */
/* * ver.1.0 K. Watanabe 2006.4.1                                */
/* ************************************************************* */
//     ==== input ====
//    [A]...... Matrix(nXm)
//
//    ==== output ====
//    [B].......Matrix (mXn)

void GenMatrix::transpose(const GenMatrix& A, GenMatrix& B ){
	try{

	int ia, ja, jap;
	int numNZa;
	int n, m;
	double val;
	
	/*check dimension*/
	n = A.getNumRow();
	m = A.getNumColumn();
	if( B.getNumRow() != m ||  B.getNumColumn() != n ){
		throw OverIdxErr();
	}

	/* [AB]の初期化 */
	B.quickAllZero();

	/*計算*/
	for(ia = 0; ia < n; ia++){
		numNZa = A.numNZ[ia];
		for(jap = 0; jap < numNZa; jap++){
			/*Aのia行ja列を取り出し*/
			ja  = A.column[ia][jap];
			val = A.value [ia][jap];

			B.set(ja, ia, val);
		}
	}
	}
	catch(...){
		cout << "error in transpose" << endl;
		throw;
	}
}



/* * [行列A]{ベクトルx}の計算   support OpenMP                  */
/* * {ax} = [A]{x}                                               */
/* * 処理速度最優先のため、エラーチェックはしない                */
/* * 当然{ax}はA.numRow 分の領域が確保されている必要がある       */
/* * ver.3.0 K. Watanabe 2009.4.17                               */
/* ************************************************************* */
void GenMatrix::productAX(const GenMatrix& A, const double* x, double* ax ){

	const int numRow = A.numRow;
//	for (i = 0; i < numRow; i++){
//		ax[i] = 0.0;
//	}

#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
	for(int i = 0; i < numRow; i++){
		int      numNZ = A.numNZ [i];
		int*   columnP = A.column[i];
		double* valueP = A.value [i];
		double a = 0.0;
		for(int mu = 0; mu < numNZ; mu++){
//			ax[i] += A.value[i][mu] * x[A.column[i][mu]];/*遅い*/
//			ax[i] += valueP[mu] * x[columnP[mu]];/*8%高速@ICC*/
			a += valueP[mu] * x[columnP[mu]];/*48%高速！＠ICC*/
		}
		ax[i] = a;
	}
}

/* * [行列A]{ベクトルx}の計算(複素ベクトル用)   */
/* * {ax} = [A]{x}                                               */
/* * 処理速度最優先のため、エラーチェックはしない                */
/* * 当然{ax}はA.numRow 分の領域が確保されている必要がある       */
/* * ver.3.0 K. Watanabe 2009.4.17                               */
/* ************************************************************* */

void GenMatrix::productAX(const GenMatrix& A, const Complex* x, Complex* ax ){

	const int numRow = A.numRow;
//	for (i = 0; i < numRow; i++){
//		ax[i] = 0.0;
//	}
#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
	for(int i = 0; i < numRow; i++){
		int      numNZ = A.numNZ [i];
		int*   columnP = A.column[i];
		double* valueP = A.value [i];
		Complex a = 0.0;
		for(int mu = 0; mu < numNZ; mu++){
//			ax[i] += A.value[i][mu] * x[A.column[i][mu]];
			a += valueP[mu] * x[columnP[mu]];
		}
		ax[i] = a;
	}
}

/* * [行列A]の転置行列と{ベクトルx}の積の計算   */
/* * {ax} = [A]T{x}                                              */
/* * 処理速度最優先のため、エラーチェックはしない                */
/* * 当然{ax}はA.numColumn 分の領域が確保されている必要がある    */
/* * ver.1.1 K. Watanabe 2001 7/19                               */
/* ************************************************************* */

void GenMatrix::productTranAX(const GenMatrix& A, const double* x, double* ax ){
	
	/* {ax}の初期化 */
	const int numRow    = A.numRow;
	const int numColumn = A.numColumn;
	for(int i = 0; i < numColumn; i++){
		ax[i] = 0.0;
	}

	for(int i = 0; i < numRow; i++){
		int numNZ = A.numNZ[i];
		for(int mu = 0; mu < numNZ; mu++){
			ax[A.column[i][mu]] += A.value[i][mu] * x[i];
		}
	}
}

/* * [行列A]の転置行列と{ベクトルx}の積の計算 (複素ベクトル用)  */
/* * {ax} = [A]T{x}                                              */
/* * 処理速度最優先のため、エラーチェックはしない                */
/* * 当然{ax}はA.numColumn 分の領域が確保されている必要がある    */
/* * ver.1.0 K. Watanabe 2001 7/19                               */
/* ************************************************************* */

void GenMatrix::productTranAX(const GenMatrix& A, const Complex* x, Complex* ax ){
	
	/* {ax}の初期化 */
	const int numRow    = A.numRow;
	const int numColumn = A.numColumn;
	for(int i = 0; i < numColumn; i++){
		ax[i] = 0.0;
	}
	
	for(int i = 0; i < numRow; i++){
		int numNZ = A.numNZ[i];
		for(int mu = 0; mu < numNZ; mu++){
			ax[A.column[i][mu]] += A.value[i][mu] * x[i];
		}
	}
}


/* * [行列]と[行列]の積                                          */
/* * [AB] = [A][B]                                               */
/* * ver.1.0 K. Watanabe 2004.5.31                               */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(nXm)
//     [B]...... Matrix(mXo)
//
//    ==== output ====
//    [AB].......Matrix (nXo)

void GenMatrix::productGG( GenMatrix& A, GenMatrix& B, GenMatrix& AB ){
	try{

	int ib, ja, jb;
	int numNZa, numNZb;
	double val;
	
	/*check dimension*/
	const int n = A.getNumRow();
	const int m = A.getNumColumn();
	const int o = B.getNumColumn();
	if( B.getNumRow() != m || 
	   AB.getNumRow() != n || AB.getNumColumn() != o){
		throw OverIdxErr();
	}

	/* [AB]の初期化 */
	AB.quickAllZero();

	/*計算*/
	for(int ia = 0; ia < n; ia++){
		numNZa = A.numNZ[ia];
		for(int jap = 0; jap < numNZa; jap++){
			/*Aのia行ja列を取り出し*/
			ja  = A.column[ia][jap];
			val = A.value [ia][jap];
			/*Bのja行を取り出し*/
			ib = ja;
			numNZb = B.numNZ[ib];
			for(int jbp = 0; jbp < numNZb; jbp++){
				jb = B.column[ib][jbp];
				AB.add(ia, jb, val*B.value[ib][jbp]);
			}
		}
	}
	}
	catch(...){
		cout << "error in productGG" << endl;
		throw;
	}
}


/* * [行列]と[行列の転置行列]の積                                */
/* * [AB] = [A][B]t                                              */
/* * ver.1.0 K. Watanabe 2003 12.12                               */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(nXm)
//     [B]...... Matrix(nXm)
//
//    ==== output ====
//    [AB].......Matrix (nXn)

void GenMatrix::productGGT( GenMatrix& A, GenMatrix& B, GenMatrix& AB ){
	int ia, ib, ja, jb, jpa, jpb;
	int numNZa, numNZb;

	try{
	int n, m;
	
	/*check dimension*/
	n = A.getNumRow();
	m = A.getNumColumn();
	if( B.getNumRow() != n ||  B.getNumColumn() != m ||
	   AB.getNumRow() != n || AB.getNumColumn() != n){
		throw OverIdxErr();
	}

	if(A.sortFlag == 0){
		A.sortColumn();
	}

	if(B.sortFlag == 0){
		B.sortColumn();
	}

	/* [AB]の初期化 */
	AB.quickAllZero();

	/*計算*/
	double sum;
	for (ia = 0; ia < n; ia++){
		for (ib = 0; ib <n ; ib++){
			jpa = 0;
			jpb = 0;
			numNZa = A.numNZ[ia];
			numNZb = B.numNZ[ib];
			sum = 0.0;
			while(jpa < numNZa && jpb < numNZb){
				ja = A.column[ia][jpa];
				jb = B.column[ib][jpb];
				if(ja == jb){
					sum += A.value[ia][jpa] * B.value[ib][jpb];
					jpa++;
					jpb++;
				}
				else if(ja < jb){
					jpa++;
				}
				else{
					jpb++;
				}
			}
			if (sum != 0.0){
				AB.add(ia, ib, sum);
			}
		}
	}
	}
	catch(...){
		cout << "error in productRA" << endl;
		throw;
	}
}

/* * [行列の転置行列]と[行列]の積(Bの方が非零要素が少ない場合)             */
/* * [AtB] = [A]t[B]                                               */
/* * ver.1.0 Y. Sato 2012.9.29                                   */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(mXn)
//     [B]...... Matrix(mXo)
//
//    ==== output ====
//    [AB].......Matrix (nXo)

void GenMatrix::productGtG( GenMatrix& A, GenMatrix& B, GenMatrix& AB ){
	try{

	int ib, ja, jb;
	int numNZa, numNZb;
	double val;
	
	/*check dimension*/
	const int n = A.getNumColumn();
	const int m = A.getNumRow();
	const int o = B.getNumColumn();
	if( B.getNumRow() != m || 
	   AB.getNumRow() != n || AB.getNumColumn() != o){
		throw OverIdxErr();
	}

	/* [AB]の初期化 */
	AB.quickAllZero();

	/*計算*/
	for(int ib = 0; ib < m; ib++){
		numNZb = B.numNZ[ib];
		for(int jbp = 0; jbp < numNZb; jbp++){
			/*Bのib行jb列を取り出し*/
			jb  = B.column[ib][jbp];
			val = B.value [ib][jbp];
			/*Aのjb行を取り出し*/
			int ia = ib;
			numNZa = A.numNZ[ia];
			for(int jap = 0; jap < numNZa; jap++){
				ja = A.column[ia][jap];
				AB.add( ja, jb , val*A.value[ia][jap]);
			}
		}
	}
	}
	catch(...){
		cout << "error in productGG" << endl;
		throw;
	}
}

/* * [行列の転置行列]と[行列]の積(Bの方が非零要素が少ない場合)             */
/* * [AtB] = [A]t[B]                                               */
/* * ver.1.0 Y. Sato 2014.9.12                                   */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(mXn)
//     [B]...... Matrix(mXo)
//
//    ==== output ====
//    [AB].......Matrix (nXo)

void GenMatrix::productGtG_( GenMatrix& A, GenMatrix& B, GenMatrix& AB ){
	try{

	int ib, ja, jb;
	int numNZa, numNZb;
	double val;
	
	/*check dimension*/
	const int n = A.getNumColumn();
	const int m = A.getNumRow();
	const int o = B.getNumColumn();
	if( B.getNumRow() != m || 
	   AB.getNumRow() != n || AB.getNumColumn() != o){
		throw OverIdxErr();
	}

	/* [AB]の初期化 */
//	AB.quickAllZero();

	/*計算*/
	for(int ib = 0; ib < m; ib++){
		numNZb = B.numNZ[ib];
		for(int jbp = 0; jbp < numNZb; jbp++){
			/*Bのib行jb列を取り出し*/
			jb  = B.column[ib][jbp];
			val = B.value [ib][jbp];
			/*Aのjb行を取り出し*/
			int ia = ib;
			numNZa = A.numNZ[ia];
			for(int jap = 0; jap < numNZa; jap++){
				ja = A.column[ia][jap];
				AB.add( ja, jb , val*A.value[ia][jap]);
			}
		}
	}
	}
	catch(...){
		cout << "error in productGG" << endl;
		throw;
	}
}


/* * 行列の逆行列を計算 Gaussian elimination Method              */
/* * 受け取った行列の逆行列に変換して返却                        */
/* * 行列サイズは数十以下が実用的,また密行列を前提にしている     */
/* * ver.1.03 K. Watanabe 2009.5.7                               */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(nXn)
//
//    ==== output ====
//    [A].......Matrix (nXn)
//   return value : 対角項が0となった行
//                  正常終了なら-1

int GenMatrix::inv( GenMatrix& A){
	try{
	
	/*check dimension*/
	const int n = A.getNumRow();
	const int m = A.getNumColumn();
	if( n != m ){
		throw OverIdxErr();
	}


	/*work領域確保と右側に単位行列作成*/
	const int mm = 2*m;
	double** w = new double*[n];
	for(int i = 0; i < n; i++){
		w[i] = new double[mm];
		for(int j = 0; j < mm; j++){
			w[i][j] = 0.0;
		}
		w[i][m+i] = 1.0;
	}
	
	/*右側に変換前のデータをコピー*/
	for(int i = 0; i < n; i++){
		int nz = A.numNZ[i];
		for(int jj = 0; jj < nz; jj++){
			w[i][A.column[i][jj]] = A.value[i][jj];
		}
	}
	
	/*掃きだし計算*/
	for(int i = 0; i < n; i++){
		double ad = w[i][i];/*対角項*/
		if(fabs(ad) < 1.0e-200){
			cout << "diag is very small value in " << i << "-th Row :" << ad << endl;
			return i;
		}

		for(int j = i+1; j < mm; j++){
			w[i][j] /= ad;
		}
		for(int ii =0; ii < n; ii++){
			if(ii != i){
				double mp = w[ii][i];
				for(int j = i+1; j < mm; j++){
					w[ii][j] -= mp*w[i][j];
				}
			}
		}
		
	}
	/*結果の書き戻し*/
	A.quickAllZero();
	for(int i = 0; i < n; i++){
		int sz = A.szColumn[i];
		if(sz < m){/*サイズチェック*/
			delete[] A.value [i];
			delete[] A.column[i];
			int newsz = m;
			A.column[i] = new    int[newsz];
			A.value [i] = new double[newsz];
			A.szColumn[i] = newsz;
		}
		for(int j = 0; j < m; j++){
			A.column[i][j] = j;
			A.value [i][j] = w[i][m+j];
		}
		A.numNZ[i] = m;
	}
	
	
	/*work領域解放*/
	for(int i = 0; i < n; i++){
		delete[] w[i];
	}
	delete[] w;
	return -1;/*正常終了*/
	
	}
	catch(...){
		cout << "error in inv()" << endl;
		throw;
	}
}



/* * 対角スケーリング　行列の固有値分布の調査用                  */
/* * ver.1.00 K. Watanabe 2010.3.5                               */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(nXn)
//
//    ==== output ====
//    [A].......Matrix (nXn)

void GenMatrix::diagnalScaling( GenMatrix& A){
	try{
	
	/*check dimension*/
	const int n = A.getNumRow();
	const int m = A.getNumColumn();
	if( n != m ){
		throw OverIdxErr();
	}

	/*get diagnal part of L*/
	double* diag = new double[n];
	getDiagVector(A, diag);

	
	/*scaling*/
	for(int i = 0; i < n; i++){
		const int    numNZ   = A.numNZ [i];
		const int*   columnP = A.column[i];
		      double* valueP = A.value [i];
		for(int nu = 0; nu < numNZ; nu++){
			int j = columnP[nu];
			valueP[nu] /= sqrt(diag[i]*diag[j]);
		}
	}

	delete[] diag;
	return;/*正常終了*/
	
	}
	catch(...){
		cout << "error in diagnalScaling()" << endl;
		throw;
	}
}



/* * 対角スケーリング　行列の固有値分布の調査用                  */
/* *an行までとそれ以降を分けてスケーリング for deflated ICCG     */
/* * ver.1.00 K. Watanabe 2010.5.26                               */
/* ************************************************************* */
//     ==== input ====
//     [A]...... Matrix(nXn)
//
//    ==== output ====
//    [A].......Matrix (nXn)

void GenMatrix::diagnalScaling( GenMatrix& A, const int an){
	try{
	
	/*check dimension*/
	const int n = A.getNumRow();
	const int m = A.getNumColumn();
	if( n != m || an > n){
		throw OverIdxErr();
	}

	/*get diagnal part of L*/
	double* diag = new double[n];
	getDiagVector(A, diag);

	
	/*scaling*/
	for(int i = 0; i < an; i++){
		const int    numNZ   = A.numNZ [i];
		const int*   columnP = A.column[i];
		double*       valueP = A.value [i];
		for(int nu = 0; nu < numNZ; nu++){
			int j = columnP[nu];
			if(j < an){
				valueP[nu] /= sqrt(diag[i]*diag[j]);
			}
			else{
				valueP[nu] /= sqrt(1.0*diag[j]);
			}
		}
	}
	for(int i = an; i < n; i++){
		const int    numNZ   = A.numNZ [i];
		const int*   columnP = A.column[i];
		double*       valueP = A.value [i];
		for(int nu = 0; nu < numNZ; nu++){
			int j = columnP[nu];
			if(j < an){
				valueP[nu] /= sqrt(diag[i]*1.0);
			}
			else{
				valueP[nu] /= sqrt(diag[i]*diag[j]);
			}
		}
	}

	delete[] diag;
	return;/*正常終了*/
	
	}
	catch(...){
		cout << "error in diagnalScaling()" << endl;
		throw;
	}
}



/* **************************************************************** */
/* 行列の上半分を埋めた一次元配列の作成 for CLAPACK & MKL           */
/*  対称正方行列のみ                                                */
/*   | 0,1,3,6,...                */
/*   |   2,4,7,...                */
/*   |     5,8,...                */
/*   |       9,...                */
/*   |       ,...                 */

/* input  : A       */
/* output : up                          */
/*     ver. 1.0 2010.3.5  K. Watanabe   */
/* **************************************************************** */
void GenMatrix::makeuppmat(const GenMatrix& A, double *up){
	try{
	
	const int n   = A.getNumRow();
	const int nup = n*(n+1)/2;
	
	/* ***Initialize *** */
	for(int i = 0; i < nup; i++){
		up[i] = 0.0;
	}
	
	for(int i = 0; i < n; i++){
		const int    numNZ   = A.numNZ [i];
		const int*   columnP = A.column[i];
		   double*    valueP = A.value [i];

		for(int mu = 0; mu < numNZ; mu++){
			int j = columnP[mu];
			if(j >= i){
				int idx = i +  (j + 1)*j/2;
				up[idx] = valueP[mu];
			}
		}
		
	}
	return;
	}
	catch (...){
		cout << "unknown error in makeuppmat" << endl;
	}
}




/* * 対角項を抽出し、1次元配列に格納(ただし、正方行列のみ)       */
/* * 正方行列で無い場合は例外OverIdxErrを投げる                  */
/* * 格納する配列の領域はあらかじめ確保するように                */
/* * Ver. 1.00 K. Watanabe 2001 9.25                             */
/* ************************************************************* */
void GenMatrix::getDiagVector(const GenMatrix& A, double* diag ){
	int i, mu;
	int numNZ;
	
	/*正方行列かどうかのチェック*/
	int numRow = A.numRow;
	if(numRow != A.numColumn){
		throw OverIdxErr();
	}
	
	/* 対角項の抽出*/
	for (i = 0; i < numRow; i++){
		numNZ = A.numNZ[i];
		diag[i] = 0.0;
		for (mu =0; mu < numNZ; mu++){
			if (A.column[i][mu] == i){
				diag[i] = A.value[i][mu];
			}
		}
	}
}

/* * 最大バンド幅（各行の非ゼロ要素数の最大値）を取得して        */
/* * 返却する                                                    */
/* * Ver. 1.00 K. Watanabe 2001 9.25                             */
/* ************************************************************* */
int GenMatrix::searchNumMaxBW() const {
	int i;
	int maxNZ;
	
	maxNZ = 0;
	for (i = 0; i < numRow; i++){
		if(numNZ[i] > maxNZ){
			maxNZ = numNZ[i];
		}
	}
	
	return maxNZ;
}

/* ***************************************************************** */
void GenMatrix::clearRow(int i){
/* ********************************************************************
//  指定したi行目を全てゼロにする
//
//
//     ver. 1.00 2002 10.2   K. Watanabe
********************************************************************** */
	 numNZ[i] = 0;   /*i行目を項全てを零にする*/

//	/*最大バンド幅を再調査*/
//	numMaxBW = searchNumMaxBW();

}


/* ***************************************************************** */
void GenMatrix::setBoundaryCondition(double* b, vector<int>& bclist, double val){
/* ********************************************************************
//  有限要素法の境界条件設定用method(高速版)
//  境界条件を設定したい行番号リストを配列(vector<int> bclist)に格納して、
//  配列bclistにリストされた行iに対して
//   ・i行i列の非対角項をゼロにする
//   ・i行の対角項を1.0にする
//   ・右辺ベクトルbのi行目をvalにする（b[i] = val）
//   ・i行目以外の行j (j=0,1,...,numRow-1)に対して b[j] -= this[j][i]*val
//  を実行する（this[j][i]はこのmethodの対象の行列のj行i列目）
//
// このメソッドを呼び出すとbclistが昇順にソートされることに注意
//
//
//     ver. 1.02 2003 12.11   K. Watanabe
********************************************************************** */
	int i, k, mu, numnz, bcsize;

	/*検索速度向上のために配列bclistをソーティング*/
	sort(bclist.begin(), bclist.end());

	/*i列目の非対角項を右辺ベクトルに移項*/
	for(i = 0; i < numRow; i++){
		numnz = numNZ[i];
		for(mu = 0; mu < numnz; mu++){
			if(binary_search(bclist.begin(), bclist.end(),column[i][mu])){
				b[i] -= value[i][mu]*val;
				value[i][mu] = 0.0;
			}
		}
	}
	bcsize = (int) bclist.size();
	for (k = 0; k < bcsize; k++){
		i = bclist[k];
		     b[i] = val;
		 numNZ[i] = 1;   /*i行目の対角項のみを残して非対角項を削除*/
		 value[i][0] = 1.0; /*対角項を1.0*/
		column[i][0] = i;
	}

	/*最大バンド幅を再調査*/
//	numMaxBW = searchNumMaxBW();

}

/* ***************************************************************** */
int GenMatrix::find_index(vector<int>& bclist, int val){
/* ********************************************************************
//  ソート済み配列bclistの中からvalと一致する要素を見つけ，その配列番号を
//  返す。もし万が一見つからなかった場合は-1を返す。
//  一致する要素が複数ある場合は，最初に見つかった要素の番号を返す
// 主に境界条件設定に用いる。
// C++標準の検索関数では戻り値の型がiteratorなのでとっても不便
//
//     ver. 1.00 2007.8.7   K. Watanabe
********************************************************************** */
	int begin, end, half;
	begin = 0;
	end = bclist.size();
	while(1){
		if(end - begin < 1){
			return -1;
		}
		else{
			half = begin + (end - begin)/2;
			if(bclist[half] == val){
				return half;
			}
			else if(bclist[half] > val){
				end = half;
			}
			else{
				begin = half+1;
			}
		}
	}
}


/* ***************************************************************** */
void GenMatrix::setBoundaryCondition(double* b, const vector<int>& bclistOrg, double* bcValOrg){
/* ********************************************************************
//  有限要素法の境界条件設定用method(高速版)
//  境界条件を設定したい行番号リストを配列(vector<int> bclist)に格納する。
// さらに，設定したい値を配列bcValに入れる。

//  配列bclistにリストされた行iに対して
//   ・i行i列の非対角項をゼロにする
//   ・i行の対角項を1.0にする
//   ・右辺ベクトルbのi行目をbcVal[i]にする（b[i] = bcVal[i]）
//   ・i行目以外の行j (j=0,1,...,numRow-1)に対して b[j] -= this[j][i]*bcVal[i]
//  を実行する（this[j][i]はこのmethodの対象の行列のj行i列目）
//
//     ver. 1.00 2009.2.2   K. Watanabe
********************************************************************** */

	/*行番号リストをコピー*/
	vector<int> bclist;
	const int n = bclistOrg.size();
	bclist.resize(n);
	for(int i = 0; i < n; i++){
		bclist[i] = bclistOrg[i];
	}
	
	/*検索速度向上のために配列bclistをソーティング*/
	sort(bclist.begin(), bclist.end());

	/*bcValもソートする*/
	double* bcVal = new double[n];
	int index;
	for(int i = 0; i < n; i++){
		index = find_index(bclist, bclistOrg[i]);
		bcVal[index] = bcValOrg[i];
	}

	/*j列目の非対角項を右辺ベクトルに移項*/
	for(int i = 0; i < numRow; i++){
		int numnz = numNZ[i];
		for(int mu = 0; mu < numnz; mu++){
			if((index = find_index(bclist, column[i][mu])) != -1){
				b[i] -= value[i][mu]*bcVal[index];
				value[i][mu] = 0.0;
			}
		}
	}

	const int nbc = (int) bclist.size();
	for(int k = 0; k < nbc; k++){
		int i = bclist[k];
		     b[i] = bcVal[k];
		 numNZ[i] = 1;   /*i行目の対角項のみを残して非対角項を削除*/
		 value[i][0] = 1.0; /*対角項を1.0*/
		column[i][0] = i;
	}

	delete[] bcVal;
}

/* *********************************************************** */
void GenMatrix::icdcmp(const GenMatrix& A, GenMatrix& LD, double ga){
//
//     Imcomplete Cholesky decomposition of A
//           into LDLt in ICCG method  for symetory matrix
//           (Lt = transpose of L)
//
//     Relation between L(K,K) and D(K) is
//           L(K,K)*D(K) = 1  (ICCG version)
//
/*     ver. 1.01 2006.12.23   K. Watanabe                      */
/* *********************************************************** */


	int n, k, mu, i;
	int mui, muk;
	double t;
	int* columnkP;
	int* columniP;
	
	int numRow = A.numRow;
	n = numRow;

	try{

		double* diag = new double[n];//diagnal part of A

		/*check */
		if(numRow != A.numColumn){
			cout << "GenMatrix::icdcmp() Matrix must be symmetry matrix" << endl;
			throw;
		}

		/* copy A -> LD with sorting*/
		LD = A;
		LD.sortColumn();

		/*get diagnal part*/
		getDiagVector(LD, diag);
				
		for (k = 0; k < n; k++){
			columnkP = LD.column[k];
			for (mu = 0; (i = columnkP[mu]) < k; mu++){
//				i = LD.column[k][mu];
				columniP = LD.column[i];

				mui = 0;
				muk = 0;
				while ((columniP[mui] < i )&&(columnkP[muk] < k)){
					if(columniP[mui] > columnkP[muk]){
						muk++;
					}
					else if(columniP[mui] < columnkP[muk]){
						mui++;
					}
					else{
						LD.value[k][mu] -= diag[columnkP[muk]] 
						                 * LD.value[i][mui] * LD.value[k][muk];
						mui++;
						muk++;
					}
				}
			}
			t = ga * diag[k];
			for (mu = 0; (i = columnkP[mu]) < k; mu++){
				t = t - diag[columnkP[mu]] * LD.value[k][mu] * LD.value[k][mu];
			}
			diag[k] = 1.0/t;

			/*write back diagnal part*/
			if(i != k){ cout << "error ! " <<i<<" : "<<k<< endl; exit(1);}
			LD.value[k][mu] = diag[k];
			
			/*discard upper parts*/
			LD.numNZ[k] = mu+1;
			
		}
		
		delete[] diag;
	}
	catch(...){
		cout << "error in icdcmp()" << endl;
		throw;
	}
}


/* *********************************************************** */
void GenMatrix::icSolv(const GenMatrix& LD,
                        const double* b, double* x, double* y){
//
//              Solve system of linear equations
//                  (L*D*Lt)*X = B
//             in ICCG method (Lt = transpose of L)
//             LDLt is incomplete Cholesky decomposition
//             arranged by ICDCMP
//
//     ver. 1.03 2009.4.17  K. Watanabe
/* *********************************************************** */
// input  LD : incomplete Cholesky decomposited matrix
//        b  : RHS vector
// in/out x  : solution vector


	double* LDpv;
	int*    LDpi;
	
	const int n = LD.numRow;
	try{
//		double* y = new double[n];/*関数の外で用意するよりも速い*/
/*     ---- forward substitution ----*/
		/* solve [L]{y} = {b} */
		for(int k = 0; k < n; k++){
			double t = b[k];
			int numNZ = LD.numNZ[k];
			LDpv  = LD.value [k];
			LDpi  = LD.column[k];
			for(int mu = 0; mu < numNZ-1; mu++){// numNZ-1 is the diagnal !!
//				t -= LD.value[k][mu] * y[LD.column[k][mu]];
				t -= LDpv[mu] * y[LDpi[mu]];
			}
//			y[k] = diag[k] * t;
//			y[k] = LD.value[k][numNZ-1] * t;
			y[k] = LDpv[numNZ-1] * t;
		}
	
		/* init x[] */
		for(int k = 0; k < n; k++){
			x[k] = 0.0;
		}

/*     ---- backward substitution ----*/
		/* solve [D][L]{x} = {y} */
		/* ** d(i,i)*L(i,i) = 1 ** */
		for(int k = n-1; k >= 0; k--){
			int numNZ = LD.numNZ[k];
			LDpv = LD.value[k];
			LDpi = LD.column[k];
//			x[k] = y[k] - diag[k] * x[k];
//			x[k] = y[k] - LD.value[k][numNZ-1] * x[k];
			x[k] = y[k] - LDpv[numNZ-1] * x[k];
			for(int nu = 0; nu < numNZ-1; nu++){
//				x[LD.column[k][nu]] += LD.value[k][nu] * x[k];
				x[LDpi[nu]] += LDpv[nu] * x[k];
			}
		}
		
//		delete[] y;
	}
	catch(...){
		throw ;
	}

}

/* Conjugate gradient method for ICCG*/
int GenMatrix::iccgSolv(const GenMatrix& A, const GenMatrix& L, const double* b,
                      double* x, double* r, double eps, int maxItr){
/*     Solution of system of linear equation                     */
/*             A*X = B (A:symmetric)                             */
/*      Conjugate gradient method for ICCG                       */
/*                                                               */
/*     ver. 1.11 2009.8.13  K. Watanabe                          */
/* ************************************************************* */
//    ==== input ====
//     A......Matrix
//     L......Matrix [L][D][Lt] = A

//     b[i]...right hand side vector(i=0,1,2,...n-1)
//      eps...absolute tolerance for residual vector
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== output ====
//    r[i]...residual vector (i=0,1,2,...n-1)
//
//    ==== return value ====
//        converged : Number of iteration
//    not converged : -1

	int convFlag  = 0;
	double  eps2;
	double r2sum;
	double rur0, rur1, pap, alpha, beta;

	double* p;
	double* ap;
	double* ru;
	double* work;
	const int n = A.numRow;
	
	eps2 = eps * eps;
	
	int numItr = 0;
	try{

	/*get diagnal part of L*/
//	double* diag = new double[n];
//	getDiagVector(L, diag);


#if ___DBGPRINTG
//	fstream fp("result_ICCG.csv", ios::out);
//	if( fp.fail() == true ){
//		return -1;
//	}
//	double firstr2sum;
	
#endif

	ap   = new double[n];
	 p   = new double[n];
	ru   = new double[n];
	work = new double[n];
	
	/***Cal. {AP} = [A]{x0}***/
	productAX(A, x, ap);
	

	/***Cal. residual vector {r0} = {b} - [A]{X0} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] = b[i] - ap[i];
		r2sum += r[i]*r[i];
	}
#if ___DBGPRINTG
//		firstr2sum = r2sum;
#endif
	/*** check convergence**/
	if (r2sum <= eps2){
#if ___DBGPRINTG
//		cout << "---- converged ------" << endl;
//		cout << "Number of iteration = " << 0 << endl;
#endif
		/***free memory***/
//		delete[] diag;

		delete[]  p;
		delete[] ap;
		delete[] ru;
		delete[] work;
	
		/***  termination ***/
		return 0;
	}

	/***solve [L][D][Lt]{ru} = {r}***/
//	icSolv(L, r, ru);
	icSolv(L, r, ru, work);

	rur0 = 0.0;
	for(int i = 0; i < n; i++){
		p[i] = ru[i];
		rur0 += r[i] * ru[i];
	}
	
	/**** iteration loop ****/
	for(int k = 0; k < maxItr; k++){
		
		
		/* *** {ap} = [A]{p} *** */
		productAX(A, p, ap);
		pap = 0.0;
		for(int i = 0; i < n; i++){
			pap += p[i] * ap[i];
		}

		alpha = rur0/pap;
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			x[i] += alpha *  p[i];
			r[i] -= alpha * ap[i];
			r2sum += r[i]*r[i];
		}
		
//cout << k << " , " << sqrt(r2sum/firstr2sum) << endl;

#if ___DBGPRINTG
//		cout << k <<  ":norm{r} = " << sqrt(r2sum)  << endl;
//		if(k%10==0){
//			fp << k << " , " << sqrt(r2sum/firstr2sum) << endl;
//		}
#endif
		/*** check convergence**/
		if (r2sum <= eps2){
#if ___DBGPRINTG
			//cout <<"---- converged ------" << endl;
			//cout <<"Number of iteration = " << k+1 << endl;
			
#endif
			numItr = k+1;
			convFlag  = 1;
			
			break;
		}
		
		/***solve [L][D][Lt]{ru} = {r}***/
//		icSolv(L, r, ru);
		icSolv(L, r, ru, work);

		rur1 = 0.0;
		for(int i = 0; i < n; i++){
			rur1 += r[i] * ru[i];
		}
		beta = rur1/rur0;
		rur0 = rur1;
		
		for(int i = 0; i < n; i++){
			p[i] = ru[i] + beta * p[i];
		}
	}

	/*** not convergence**/
	if(convFlag == 0){
		numItr = -1;
#if ___DBGPRINTG
		
		cout << "--- not converged ---" << endl;
	
#endif
	}

	/***free memory***/
//	delete[] diag;
	
	delete[]  p;
	delete[] ap;
	delete[] ru;
	delete[] work;

	}
	catch(...){
		throw ;
	}

	/***  termination ***/
	return numItr;
	
}


/* *********************************************************************
//    Bi-CG STAB Solver
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned BiConjugate Gradient Stabilized method
//
//     Ver. 1.3 2008.7.24  K. Watanabe
********************************************************************** */
int GenMatrix::biCGStabSolv(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, double ga, bool pivotFlag){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     ga........acc. ratio for ILU decomposition
//     pivotFlag.. true = pivot serection in ILU decomposition
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{
	double r2sum;

	int convFlag  = 0;
	double alpha = 0.0;
	double beta;
	double omega = 1.0;/* non zero value */
	double rou0  = 1.0;/* non zero value */
	double rou, ts, tt, r0v;
	
	const int n = A.numRow; /*係数行列の行数を取得*/

	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*ILU分解の際のピボット選択による行入れ替え情報配列用意*/
	int* swapData = new int[n];
	cout << "start:ILU decomp." << endl;

	/*不完全LU分解行列生成*/
	GenMatrix LU;
	iLUcmp(A, LU, swapData, ga, pivotFlag);
	cout << "finish:ILU decomp." << endl;
	
	double*  y  = new double[n];
	double*  r  = new double[n];
	double*  r0 = new double[n];
	
	double*  p  = new double[n];
	double*  pt = new double[n];
	double*  s  = new double[n];
	double*  st = new double[n];
	double*  t  = new double[n];
	double*  v  = new double[n];

	/*初期残差ベクトルの計算*/
	/*** {AP} = [A]{x}***/
	productAX(A, x, y);

	/*** residual vector {r0} = {b} - [A]{X} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] =  b[i] - y[i];
		p[i] = r0[i] = r[i];
		r2sum += r[i]*r[i];
	}


	/*** check convergence**/
	if(r2sum < eps2) {/*収束判定*/
#if ___DBGPRINTG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
		/***  termination ***/
		convFlag  = 1;
	}

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for (itr = 1; itr < maxItr; itr++){
//		cout << itr << "-th iteration" << endl;


		/* rou = {r0}.{r} */
		rou = 0.0;
		for(int i =0; i < n; i++){
			rou += r0[i]*r[i];
		}
		if(rou == 0.0){
			cout << "Error? BiCGStab method fail !" << endl;
			break;
		}

		/*comp. beta*/
		beta = (rou/rou0) * (alpha/omega);

		for(int i =0; i < n; i++){
			p[i] = r[i] + beta * (p[i] - omega*v[i]);
		}

		/* [M]{pt} == {p} を解く（M：前処理行列）*/
		if(pivotFlag){
			iLUSolv(LU, p, pt, swapData);/* M : 不完全LU分解 */
		}
		else{
			iLUSolv(LU, p, pt);
		}
//for(int i =0; i < n; i++){
//	pt[i] = p[i];
//}
		
		/***{v} = [A]{pt}***/
		productAX(A, pt, v);

		/* alpha = rou/ {r0}.{v} */
		r0v = 0.0;
		for(int i =0; i < n; i++){
			r0v += r0[i]*v[i];
		}
		alpha = rou/r0v;
//		cout << "alpha =" << alpha << endl;

		/* {s} = {r} - alpha*{v} */
		for(int i =0; i < n; i++){
			s[i] = r[i] - alpha * v[i];
		}

		/* [M]{st} == {s} を解く（M：前処理行列）*/
		if(pivotFlag){
			iLUSolv(LU, s, st, swapData);/* M : 不完全LU分解 */
		}
		else{
			iLUSolv(LU, s, st);
		}
//for(int i =0; i < n; i++){
//	st[i] = s[i];
//}
		
		/* {t} = [A]{st}*/
		productAX(A, st, t);
		
		/* omega = {t}.{s}/{t}.{t} */
		ts = 0.0;
		tt = 0.0;
		for(int i =0; i < n; i++){
			ts += t[i]*s[i];
			tt += t[i]*t[i];
		}
		omega = ts/tt;
//		cout << "omega =" << omega << endl;
		
		/*解の更新*/
		r2sum = 0.0;
		for(int i =0; i < n; i++){
			x[i] += alpha * pt[i] + omega * st[i];
			r[i] = s[i] - omega * t[i];
			r2sum += r[i]*r[i];
		}
//		cout << "norm{r} = " << sqrt(r2sum) << endl;
		
#if ___DBGPRINTG
		cout << "norm{r} = " << sqrt(r2sum) << endl;
#endif


		if( r2sum < eps2){ /*収束判定*/
			itr++;
			convFlag = 1;
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;

			break;
		}

		rou0 = rou;

		if (omega == 0.0){
			cout << "Error: omega==0" << endl;
			itr = -1; 
			break;
		}
	}/**************************** end loop ***********************/
	
	/*** not convergence**/
	if (convFlag == 0){
		itr = -1;
		cout << "--- not converged ---" << endl;
	}
	else{
		cout << "norm_cg{r} = " << sqrt(r2sum) << endl;
		/*真の残差を計算*/
		productAX(A, x, r);
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r[i] = b[i] - r[i];
			r2sum += r[i]*r[i];
		}
		cout << "norm_true{r} = " << sqrt(r2sum) << endl;
	}

	/***free memory***/
	delete[]  y;
	delete[]  r;
	delete[]  r0;
	delete[]  p;
	delete[]  pt;
	delete[]  s;
	delete[]  st;
	delete[]  t;
	delete[]  v;

	return itr;
	}
	catch(...){
		throw IOErr();
	}
}


/* *********************************************************************
//    Bi-CG STAB Solver
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned BiConjugate Gradient Stabilized method
//       Preconditon Matrix : M = diag[A];
//
//     Ver. 1.31 2010.8.26  K. Watanabe
********************************************************************** */
int GenMatrix::biCGStabSolvWithDiagP(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{
	double r2sum;

	int convFlag  = 0;
	double alpha = 0.0;
	double beta;
	double omega = 1.0;/* non zero value */
	double rou0  = 1.0;/* non zero value */
	double rou, ts, tt, r0v;
	
	const int n = A.numRow; /*係数行列の行数を取得*/

	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*対角情報配列用意*/
	double* diag = new double[n];

	/*対角成分取り出し*/
	getDiagVector(A, diag);
	
	double*  y  = new double[n];
	double*  r  = new double[n];
	double*  r0 = new double[n];
	
	double*  p  = new double[n];
	double*  pt = new double[n];
	double*  s  = new double[n];
	double*  st = new double[n];
	double*  t  = new double[n];
	double*  v  = new double[n];

	/*初期残差ベクトルの計算*/
	/*** {AP} = [A]{x}***/
	productAX(A, x, y);

	/*** residual vector {r0} = {b} - [A]{X} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] =  b[i] - y[i];
		p[i] = r0[i] = r[i];
		r2sum += r[i]*r[i];
	}
	cout << "norm|r|= " << sqrt(r2sum) << endl;

	/*** check convergence**/
	if(r2sum < eps2) {/*収束判定*/
#if ___DBGPRINTG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
		/***  termination ***/
		convFlag  = 1;
	}

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for(itr = 1; itr < maxItr; itr++){
//		cout << itr << "-th iteration" << endl;


		/* rou = {r0}.{r} */
		rou = 0.0;
		for(int i =0; i < n; i++){
			rou += r0[i]*r[i];
		}
		if(rou == 0.0){
			cout << "itr =" << itr << endl;
			cout << "Error? rou==0 : BiCGStab method fail !" << endl;
			break;
		}

		/*comp. beta*/
		beta = (rou/rou0) * (alpha/omega);

		for(int i =0; i < n; i++){
			p[i] = r[i] + beta * (p[i] - omega*v[i]);
		}

		/* [M]{pt} == {p} を解く（M：前処理行列）*/
		for(int i =0; i < n; i++){
			pt[i] = p[i]/diag[i];
		}
		
		/***{v} = [A]{pt}***/
		productAX(A, pt, v);

		/* alpha = rou/ {r0}.{v} */
		r0v = 0.0;
		for(int i =0; i < n; i++){
			r0v += r0[i]*v[i];
		}
		alpha = rou/r0v;
//		cout << "alpha =" << alpha << endl;

		/* {s} = {r} - alpha*{v} */
		for(int i =0; i < n; i++){
			s[i] = r[i] - alpha * v[i];
		}

		/* [M]{st} == {s} を解く（M：前処理行列）*/
		for(int i =0; i < n; i++){
			st[i] = s[i]/diag[i];
		}
		
		/* {t} = [A]{st}*/
		productAX(A, st, t);
		
		/* omega = {t}.{s}/{t}.{t} */
		ts = 0.0;
		tt = 0.0;
		for(int i =0; i < n; i++){
			ts += t[i]*s[i];
			tt += t[i]*t[i];
		}
		omega = ts/tt;
//		cout << "omega =" << omega << endl;
		
		/*解の更新*/
		r2sum = 0.0;
		for(int i =0; i < n; i++){
			x[i] += alpha * pt[i] + omega * st[i];
			r[i] = s[i] - omega * t[i];
			r2sum += r[i]*r[i];
		}
//		cout << "norm{r} = " << sqrt(r2sum) << endl;
		
#if ___DBGPRINTG
		cout << "norm{r} = " << sqrt(r2sum) << endl;
#endif


		if( r2sum < eps2){ /*収束判定*/
			itr++;
			convFlag = 1;
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;

			break;
		}

		rou0 = rou;

		if (omega == 0.0){
			cout << "Error: omega==0" << endl;
			itr = -1; 
			break;
		}
	}/**************************** end loop ***********************/
	
	/*** not convergence**/
	if (convFlag == 0){
		itr = -1;
		cout << "--- not converged ---" << endl;
	}
	else{
		cout << "norm_cg{r} = " << sqrt(r2sum) << endl;
		/*真の残差を計算*/
		productAX(A, x, r);
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r[i] = b[i] - r[i];
			r2sum += r[i]*r[i];
		}
		cout << "norm_true{r} = " << sqrt(r2sum) << endl;
	}

	/***free memory***/
	delete[] diag;

	delete[]  y;
	delete[]  r;
	delete[]  r0;
	delete[]  p;
	delete[]  pt;
	delete[]  s;
	delete[]  st;
	delete[]  t;
	delete[]  v;

	return itr;
	}
	catch(...){
		throw IOErr();
	}
}


/* *********************************************************************
//    Bi-CG STAB Solver with GS type preconditioner
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned BiConjugate Gradient Stabilized method
//
//     Ver. 1.00 2010.8.26  K. Watanabe
********************************************************************** */
int GenMatrix::biCGStabSolvWithGS(GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{
	double r2sum;

	int convFlag  = 0;
	double alpha = 0.0;
	double beta;
	double omega = 1.0;/* non zero value */
	double rou0  = 1.0;/* non zero value */
	double rou, ts, tt, r0v;
	
	const int n = A.numRow; /*係数行列の行数を取得*/

	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
//	/*対角情報配列用意*/
//	double* diag = new double[n];

//	/*対角成分取り出し*/
//	getDiagVector(A, diag);
	
	double*  y  = new double[n];
	double*  r  = new double[n];
	double*  r0 = new double[n];
	
	double*  p  = new double[n];
	double*  pt = new double[n];
	double*  pq = new double[n];
	double*  s  = new double[n];
	double*  st = new double[n];
	double*  sq = new double[n];
	double*  t  = new double[n];
	double*  v  = new double[n];

	/*初期残差ベクトルの計算*/
	/*** {AP} = [A]{x}***/
	productAX(A, x, y);

	/*** residual vector {r0} = {b} - [A]{X} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] =  b[i] - y[i];
		p[i] = r0[i] = r[i];
		r2sum += r[i]*r[i];
	}
	cout << "norm|r|= " << sqrt(r2sum) << endl;

	/*** check convergence**/
	if(r2sum < eps2) {/*収束判定*/
#if ___DBGPRINTG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
		/***  termination ***/
		convFlag  = 1;
	}

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for(itr = 1; itr < maxItr; itr++){
//		cout << itr << "-th iteration" << endl;


		/* rou = {r0}.{r} */
		rou = 0.0;
		for(int i =0; i < n; i++){
			rou += r0[i]*r[i];
		}
		if(rou == 0.0){
			cout << "itr =" << itr << endl;
			cout << "Error? rou==0 : BiCGStab method fail !" << endl;
			break;
		}

		/*comp. beta*/
		beta = (rou/rou0) * (alpha/omega);

		for(int i =0; i < n; i++){
			p[i] = r[i] + beta * (p[i] - omega*v[i]);
		}

		/* {pt} <= [U]inv[L+D]{p}（GS処理）*/
		GenMatrix::gsPrecondition(A, p, pt, pq);

		for(int i =0; i < n; i++){
			v[i] = p[i] - pq[i];
		}
		
		/* alpha = rou/ {r0}.{v} */
		r0v = 0.0;
		for(int i =0; i < n; i++){
			r0v += r0[i]*v[i];
		}
		alpha = rou/r0v;
//		cout << "alpha =" << alpha << endl;

		/* {s} = {r} - alpha*{v} */
		for(int i =0; i < n; i++){
			s[i] = r[i] - alpha * v[i];
		}

		/* {st} <= [U]inv[L+D]{s}（GS処理）*/
		GenMatrix::gsPrecondition(A, s, st, sq);

		for(int i =0; i < n; i++){
			t[i] = s[i] - sq[i];
		}
		
		/* omega = {t}.{s}/{t}.{t} */
		ts = 0.0;
		tt = 0.0;
		for(int i =0; i < n; i++){
			ts += t[i]*s[i];
			tt += t[i]*t[i];
		}
		omega = ts/tt;
//		cout << "omega =" << omega << endl;
		
		/*解の更新*/
		r2sum = 0.0;
		for(int i =0; i < n; i++){
			x[i] += alpha * pt[i] + omega * st[i];
			r[i] = s[i] - omega * t[i];
			r2sum += r[i]*r[i];
		}
//		cout << "norm{r} = " << sqrt(r2sum) << endl;
		
#if ___DBGPRINTG
		cout << "norm{r} = " << sqrt(r2sum) << endl;
#endif


		if( r2sum < eps2){ /*収束判定*/
			itr++;
			convFlag = 1;
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;

			break;
		}

		rou0 = rou;

		if (omega == 0.0){
			cout << "Error: omega==0" << endl;
			itr = -1; 
			break;
		}
	}/**************************** end loop ***********************/
	
	/*** not convergence**/
	if (convFlag == 0){
		itr = -1;
		cout << "--- not converged ---" << endl;
	}
	else{
		cout << "norm_cg{r} = " << sqrt(r2sum) << endl;
		/*真の残差を計算*/
		productAX(A, x, r);
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r[i] = b[i] - r[i];
			r2sum += r[i]*r[i];
		}
		cout << "norm_true{r} = " << sqrt(r2sum) << endl;
	}

	/***free memory***/
	delete[]  y;
	delete[]  r;
	delete[]  r0;
	delete[]  p;
	delete[]  pt;
	delete[]  pq;
	delete[]  s;
	delete[]  st;
	delete[]  sq;
	delete[]  t;
	delete[]  v;

	return itr;
	}
	catch(...){
		throw IOErr();
	}
}

/* *********************************************************************
//    GPBiCG Solver with ILU decomposition
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned GPBiCG method
//       Preconditon Matrix : M = diag[A];
//
//     Ver. 1.00 2010.1.30  K. Watanabe
********************************************************************** */
int GenMatrix::GPBiCGSolv(const GenMatrix& A, const double* b, double* x,
                          double eps, int maxItr, double ga, bool pivotFlag){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     ga........acc. ratio for ILU decomposition
//     pivotFlag.. true = pivot serection in ILU decomposition
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{

	bool convFlag = false;
	const int n = A.numRow; /*係数行列の行数を取得*/
	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*ILU分解の際のピボット選択による行入れ替え情報配列用意*/
	int* swapData = new int[n];
	cout << "start:ILU decomp." << endl;

	/*不完全LU分解行列生成*/
	GenMatrix LU;
	iLUcmp(A, LU, swapData, ga, pivotFlag);
	cout << "finish:ILU decomp." << endl;

	double alpha = 0.0;
	double beta  = 0.0;
	double eta   = 1.0;
	double zeta, tau;
	double omega = 1.0;/* non zero value */
	double rou0  = 1.0;/* non zero value */
	double rou,   r0v,   r2sum;

	double mu  [5];
	double zeta_eta[2];

	double*  work = new double[n];//working vector

	double*  y  = new double[n];
	double*  r  = new double[n];
	double*  rh = new double[n];
	double*  r00= new double[n];
	double*  r0 = new double[n];
	
	double*  p  = new double[n];
	double*  t  = new double[n];
	double*  th = new double[n];
	double*  th0= new double[n];
	double*  u  = new double[n];
	double*  ap = new double[n];
	double*  aph= new double[n];
	double*  at = new double[n];
	double*  ath= new double[n];
	double*  w  = new double[n];
	double*  z  = new double[n];

	/***Cal. {AP} = [A]{x0}***/
	GenMatrix::productAX(A, x, work);
	
	/***Cal. residual vector {r0} = {b} - [A]{X0} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] = b[i] - work[i];
		r2sum += r[i]*r[i];
	}
	double firstr2sum = r2sum;
	cout << "first norm[r] = " << sqrt(r2sum) << endl;
	/*** check convergence**/
	if(r2sum <= eps2){
#if ___DEBUGPRINT_GPBICG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
	
		/***  termination ***/
		convFlag  = true;
	}

	/*** save initial residual ***/
	for(int i = 0; i < n; i++){
		p[i] = r0[i] = r[i];
	}
	
	/*** init***/
	for(int i = 0; i < n; i++){
		t[i] = th0[i] = w[i] = 0.0;
	}

	/***solve [M]{rh} == {r}***/
//	GenMatrix::icSolv(L, r, rh, work);
//	for(int i =0; i < n; i++){
//		rh[i] = r[i]/diag[i];
//	}
	if(pivotFlag){
		iLUSolv(LU, r, rh, swapData);/* M : 不完全LU分解 */
	}
	else{
		iLUSolv(LU, r, rh);
	}

	/* ***シャドー残差r0の選択 *** */
	
	/*r0 = r*/
	for(int i = 0; i < n; i++){
		r0[i] = r[i];
	}

//	/*シャドー残差を定数ベクトルにしてみる　効果があるといいなぁ*/
//	for(int i = 0; i < n; i++){
//		r0[i] = 1.0;
//	}

//	/*シャドー残差を[0,1]乱数ベクトルにしてみる　効果があるといいなぁ*/
//	for(int i = 0; i < n; i++){
//		r0[i] = rand()/(double)RAND_MAX;
//	}



#if CR_METHOD
	
	/*** {r00] = [A]{r0} ***/
	GenMatrix::productAX(A, r0, r00);
	delete[] r0;
	r0 = r00;
	
#endif
	
	/***rou = {r0}{r}***/
//	rou0 = rou = r2sum;
	rou = 0.0;
	for(int i = 0; i < n; i++){
		rou += r0[i]*r[i];
	}
	
	rou0 = rou;

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for(itr = 1; itr < maxItr; itr++){
//#if ___DEBUGPRINT_GPBICG
//		if(myid == 0){ cout << itr << "-th iteration" << endl;}
//#endif
		
		/*p の更新**/
		if(itr != 1){ // if itr==1 then p = r
			for(int i =0; i < n; i++){
				p[i] = rh[i] + beta * (p[i] - u[i]);
			}
		}

		/***compute {ap} = [A]{p}***/
		GenMatrix::productAX(A, p, ap);

		/***solve [M]{aph} == {ap}***/
//		GenMatrix::icSolv(LS, ap, aph, work);
//		for(int i =0; i < n; i++){
//			aph[i] = ap[i]/diag[i];
//		}
		if(pivotFlag){
			iLUSolv(LU, ap, aph, swapData);/* M : 不完全LU分解 */
		}
		else{
			iLUSolv(LU, ap, aph);
		}

		/***compute alpha ***/
		r0v = 0.0;
		for(int i = 0; i < n; i++){
			r0v += r0[i] * ap[i];
		}
		alpha = rou/r0v;
//		cout << "alpha = " << alpha << endl;
		
		/* {y} = {t} - {r} -alpha*({w} -{v}) */
		for(int i = 0; i < n; i++){
			y[i] = t[i] -r[i] - alpha*(w[i] - ap[i]);
		}

		/* {t} = {r} - alpha*{v} */
		for(int i = 0; i < n; i++){
			t[i] = r[i] - alpha*ap[i];
		}
		/* {th} = {rh} - alpha*{v} */
		for(int i = 0; i < n; i++){
			th[i] = rh[i] - alpha*aph[i];
		}
		
		/* *** {ath} = [A]{th} *** */
		GenMatrix::productAX(A, th, ath);

		/*eta and zeta の更新*/
//		mu1_s = mu2_s = mu3_s = mu4_s = mu5_s = 0.0;
		for(int i = 0; i < 5; i++){
			mu[i] = 0.0;
		}
		
		for(int i =0; i < n; i++){
			mu[0] +=   y[i] *   y[i];
			mu[1] += ath[i] *   t[i];
			mu[2] +=   y[i] *   t[i];
			mu[3] += ath[i] *   y[i];
			mu[4] += ath[i] * ath[i];
		}
//		cout << " mu0= " << mu[0] << " mu1= " << mu[1] << " mu2= " << mu[2] 
//		     << " mu3= " << mu[3] << " mu4= " << mu[4] << endl;

		if(itr ==1){
			zeta = 1.0;// mu[1]/mu[4];// !!!!????
			eta  = 0.0;
		}
		else{
			tau  =  mu[4]*mu[0] - mu[3]*mu[3];
			zeta = (mu[0]*mu[1] - mu[2]*mu[3])/tau;
			eta  = (mu[4]*mu[2] - mu[3]*mu[1])/tau;
		}

//#if ___DEBUGPRINT_GPBICG
//		cout << " tau= " << tau << endl;
//		cout << "zeta= " << zeta << endl;
//		cout << "eta = " << eta  << endl;
//#endif

		/***{u} = zeta{ap} + eta({t} - {r} + beta{u}) ***/
		for(int i =0; i < n; i++){
			if(itr == 1){
				u[i] = zeta*aph[i];
			}
			else{
				u[i] = zeta*aph[i] + eta*(th0[i] - rh[i] + beta*u[i]);
			}
//cout << "aph[" << i << "]=" << aph[i] << endl;
//cout << "u[" << i << "]=" << u[i] << endl;

		}
		
		/***{z} = zeta{r} + eta*{z} - alpha{u} ***/
		for(int i =0; i < n; i++){
			z[i] = zeta*rh[i] + eta*z[i] - alpha*u[i];
		}

		/*解の更新*/
		r2sum = 0.0;
		for(int i =0; i < n; i++){
			x[i] += alpha * p[i] + z[i];
			r[i] = t[i] - eta*y[i] - zeta*ath[i];
//cout << "rh[" << i << "]=" << rh[i] << endl;
//cout << "u[" << i << "]=" << u[i] << endl;
			r2sum += r[i]*r[i];
		}

#if ___DBGPRINTG
		cout << itr << " :normalized norm{r} = " << sqrt(r2sum/firstr2sum)  << endl;
//		if(itr % INTVL_PRINT == 0){
//			fp << itr << " , " << sqrt(r2sum/firstr2sum) << endl;
//		}
#endif

		/*** check convergence**/
		if (r2sum <= eps2){
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;
			convFlag = true;
			break;
		}

		/***solve [M]{rh} == {r}***/
//		GenMatrix::icSolv(L, r, rh, work);
//		for(int i =0; i < n; i++){
//			rh[i] = r[i]/diag[i];
//		}
		if(pivotFlag){
			iLUSolv(LU, r, rh, swapData);/* M : 不完全LU分解 */
		}
		else{
			iLUSolv(LU, r, rh);
		}

		/*** comp. {r0}{r}****/
		rou = 0.0;
		for(int i = 0; i < n; i++){
			rou += r0[i] * r[i];
		}
		
		if(rou == 0.0){
			cout << "GPBiCG method fail!" << endl;
			convFlag = false;
			break;
		}

		beta = (rou/rou0) * (alpha/zeta);
//#if ___DEBUGPRINT_GPBICG
//			cout << "rou  = " << rou  << endl;
//			cout << "rou0 = " << rou0 << endl;
//			cout << "beta = " << beta << endl;
//#endif

		rou0 = rou;

		/***{w} = {at} + beta{ap} ***/
		for(int i =0; i < n; i++){
			w[i] = ath[i] + beta*ap[i];
		}
		
		for(int i =0; i < n; i++){
			th0[i] = th[i];
		}
		

	}/****************★★ end : iteration loop ★★********************/

	/*** not convergence**/
	if(convFlag == false){
		itr = -1;
		cout << "--- not converged ---" << endl;
		cout << "Number of iteration = " << itr << endl;
	}
	
	/*** check real residual {r} = {bs} - [A]{Xs} ***/
	GenMatrix::productAX(A, x, y);

	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] = b[i] - y[i];
		r2sum += r[i]*r[i];
	}
	cout << "**real norm{r}       = " << sqrt(r2sum) << endl;
	cout << "**normalized norm{r} = " << sqrt(r2sum/firstr2sum)  << endl;

	delete[] swapData;
	delete[] work;
	delete[] y;
	delete[] r;
	delete[] rh;
	delete[] r00;
	delete[] r0;
	delete[] p;
	delete[] t;
	delete[] th;
	delete[] th0;
	delete[] u;
	delete[] ap;
	delete[] aph;
	delete[] at;
	delete[] ath;
	delete[] w;
	delete[] z;

	return itr;
	}
	catch(...){
		throw IOErr();
	}
}




/* *********************************************************************
//    GPBiCG Solver
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned GPBiCG method
//       Preconditon Matrix : M = diag[A];
//
//     Ver. 1.01 2009.11.17  K. Watanabe
********************************************************************** */
int GenMatrix::GPBiCGSolvWithDiagP(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{

	bool convFlag = false;
	const int n = A.numRow; /*係数行列の行数を取得*/
	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*対角情報配列用意*/
	double* diag = new double[n];

	/*対角成分取り出し*/
	getDiagVector(A, diag);

	double alpha = 0.0;
	double beta  = 0.0;
	double eta   = 1.0;
	double zeta, tau;
	double omega = 1.0;/* non zero value */
	double rou0  = 1.0;/* non zero value */
	double rou,   r0v,   r2sum;

	double mu  [5];
	double zeta_eta[2];

	double*  work = new double[n];//working vector

	double*  y  = new double[n];
	double*  r  = new double[n];
	double*  rh = new double[n];
	double*  r00= new double[n];
	double*  r0 = new double[n];
	
	double*  p  = new double[n];
	double*  t  = new double[n];
	double*  th = new double[n];
	double*  th0= new double[n];
	double*  u  = new double[n];
	double*  ap = new double[n];
	double*  aph= new double[n];
	double*  at = new double[n];
	double*  ath= new double[n];
	double*  w  = new double[n];
	double*  z  = new double[n];

	/***Cal. {AP} = [A]{x0}***/
	GenMatrix::productAX(A, x, work);
	
	/***Cal. residual vector {r0} = {b} - [A]{X0} ***/
	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] = b[i] - work[i];
		r2sum += r[i]*r[i];
	}
	double firstr2sum = r2sum;
	cout << "first norm[r] = " << sqrt(r2sum) << endl;
	/*** check convergence**/
	if(r2sum <= eps2){
#if ___DEBUGPRINT_GPBICG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
	
		/***  termination ***/
		convFlag  = true;
	}

	/*** save initial residual ***/
	for(int i = 0; i < n; i++){
		p[i] = r0[i] = r[i];
	}
	
	/*** init***/
	for(int i = 0; i < n; i++){
		t[i] = th0[i] = w[i] = 0.0;
	}

	/***solve [M]{rh} == {r}***/
//	GenMatrix::icSolv(L, r, rh, work);
	for(int i =0; i < n; i++){
		rh[i] = r[i]/diag[i];
	}

	/* ***シャドー残差r0の選択 *** */
	
	/*r0 = r*/
	for(int i = 0; i < n; i++){
		r0[i] = r[i];
	}

//	/*シャドー残差を定数ベクトルにしてみる　効果があるといいなぁ*/
//	for(int i = 0; i < n; i++){
//		r0[i] = 1.0;
//	}

//	/*シャドー残差を[0,1]乱数ベクトルにしてみる　効果があるといいなぁ*/
//	for(int i = 0; i < n; i++){
//		r0[i] = rand()/(double)RAND_MAX;
//	}



#if CR_METHOD
	
	/*** {r00] = [A]{r0} ***/
	GenMatrix::productAX(A, r0, r00);
	delete[] r0;
	r0 = r00;
	
#endif
	
	/***rou = {r0}{r}***/
//	rou0 = rou = r2sum;
	rou = 0.0;
	for(int i = 0; i < n; i++){
		rou += r0[i]*r[i];
	}
	
	rou0 = rou;

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for(itr = 1; itr < maxItr; itr++){
//#if ___DEBUGPRINT_GPBICG
//		if(myid == 0){ cout << itr << "-th iteration" << endl;}
//#endif
		
		/*p の更新**/
		if(itr != 1){ // if itr==1 then p = r
			for(int i =0; i < n; i++){
				p[i] = rh[i] + beta * (p[i] - u[i]);
			}
		}

		/***compute {ap} = [A]{p}***/
		GenMatrix::productAX(A, p, ap);

		/***solve [M]{aph} == {ap}***/
//		GenMatrix::icSolv(LS, ap, aph, work);
		for(int i =0; i < n; i++){
			aph[i] = ap[i]/diag[i];
		}

		/***compute alpha ***/
		r0v = 0.0;
		for(int i = 0; i < n; i++){
			r0v += r0[i] * ap[i];
		}
		alpha = rou/r0v;
//		cout << "alpha = " << alpha << endl;
		
		/* {y} = {t} - {r} -alpha*({w} -{v}) */
		for(int i = 0; i < n; i++){
			y[i] = t[i] -r[i] - alpha*(w[i] - ap[i]);
		}

		/* {t} = {r} - alpha*{v} */
		for(int i = 0; i < n; i++){
			t[i] = r[i] - alpha*ap[i];
		}
		/* {th} = {rh} - alpha*{v} */
		for(int i = 0; i < n; i++){
			th[i] = rh[i] - alpha*aph[i];
		}
		
		/* *** {ath} = [A]{th} *** */
		GenMatrix::productAX(A, th, ath);

		/*eta and zeta の更新*/
//		mu1_s = mu2_s = mu3_s = mu4_s = mu5_s = 0.0;
		for(int i = 0; i < 5; i++){
			mu[i] = 0.0;
		}
		
		for(int i =0; i < n; i++){
			mu[0] +=   y[i] *   y[i];
			mu[1] += ath[i] *   t[i];
			mu[2] +=   y[i] *   t[i];
			mu[3] += ath[i] *   y[i];
			mu[4] += ath[i] * ath[i];
		}
//		cout << " mu0= " << mu[0] << " mu1= " << mu[1] << " mu2= " << mu[2] 
//		     << " mu3= " << mu[3] << " mu4= " << mu[4] << endl;

		if(itr ==1){
			zeta = 1.0;// mu[1]/mu[4];// !!!!????
			eta  = 0.0;
		}
		else{
			tau  =  mu[4]*mu[0] - mu[3]*mu[3];
			zeta = (mu[0]*mu[1] - mu[2]*mu[3])/tau;
			eta  = (mu[4]*mu[2] - mu[3]*mu[1])/tau;
		}

//#if ___DEBUGPRINT_GPBICG
//		cout << " tau= " << tau << endl;
//		cout << "zeta= " << zeta << endl;
//		cout << "eta = " << eta  << endl;
//#endif

		/***{u} = zeta{ap} + eta({t} - {r} + beta{u}) ***/
		for(int i =0; i < n; i++){
			if(itr == 1){
				u[i] = zeta*aph[i];
			}
			else{
				u[i] = zeta*aph[i] + eta*(th0[i] - rh[i] + beta*u[i]);
			}
//cout << "aph[" << i << "]=" << aph[i] << endl;
//cout << "u[" << i << "]=" << u[i] << endl;

		}
		
		/***{z} = zeta{r} + eta*{z} - alpha{u} ***/
		for(int i =0; i < n; i++){
			z[i] = zeta*rh[i] + eta*z[i] - alpha*u[i];
		}

		/*解の更新*/
		r2sum = 0.0;
		for(int i =0; i < n; i++){
			x[i] += alpha * p[i] + z[i];
			r[i] = t[i] - eta*y[i] - zeta*ath[i];
//cout << "rh[" << i << "]=" << rh[i] << endl;
//cout << "u[" << i << "]=" << u[i] << endl;
			r2sum += r[i]*r[i];
		}

#if ___DBGPRINTG
		cout << itr << " :normalized norm{r} = " << sqrt(r2sum/firstr2sum)  << endl;
//		if(itr % INTVL_PRINT == 0){
//			fp << itr << " , " << sqrt(r2sum/firstr2sum) << endl;
//		}
#endif

		/*** check convergence**/
		if (r2sum <= eps2){
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;
			convFlag = true;
			break;
		}

		/***solve [M]{rh} == {r}***/
//		GenMatrix::icSolv(L, r, rh, work);
		for(int i =0; i < n; i++){
			rh[i] = r[i]/diag[i];
		}

		/*** comp. {r0}{r}****/
		rou = 0.0;
		for(int i = 0; i < n; i++){
			rou += r0[i] * r[i];
		}
		
		if(rou == 0.0){
			cout << "GPBiCG method fail!" << endl;
			convFlag = false;
			break;
		}

		beta = (rou/rou0) * (alpha/zeta);
//#if ___DEBUGPRINT_GPBICG
//			cout << "rou  = " << rou  << endl;
//			cout << "rou0 = " << rou0 << endl;
//			cout << "beta = " << beta << endl;
//#endif

		rou0 = rou;

		/***{w} = {at} + beta{ap} ***/
		for(int i =0; i < n; i++){
			w[i] = ath[i] + beta*ap[i];
		}
		
		for(int i =0; i < n; i++){
			th0[i] = th[i];
		}
		

	}/****************★★ end : iteration loop ★★********************/

	/*** not convergence**/
	if(convFlag == false){
		itr = -1;
		cout << "--- not converged ---" << endl;
		cout << "Number of iteration = " << itr << endl;
	}
	
	/*** check real residual {r} = {bs} - [A]{Xs} ***/
	GenMatrix::productAX(A, x, y);

	r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		r[i] = b[i] - y[i];
		r2sum += r[i]*r[i];
	}
	cout << "**real norm{r}       = " << sqrt(r2sum) << endl;
	cout << "**normalized norm{r} = " << sqrt(r2sum/firstr2sum)  << endl;

	delete[] diag;
	delete[] work;
	delete[] y;
	delete[] r;
	delete[] rh;
	delete[] r00;
	delete[] r0;
	delete[] p;
	delete[] t;
	delete[] th;
	delete[] th0;
	delete[] u;
	delete[] ap;
	delete[] aph;
	delete[] at;
	delete[] ath;
	delete[] w;
	delete[] z;

	return itr;
	}
	catch(...){
		throw IOErr();
	}
}


/* *********************************************************** */
void GenMatrix::iLUcmp(const GenMatrix& A, GenMatrix& LU, int* swapData, double ga,
                       bool pivotFlag){
//
//     Imcomplete LU decomposition of A
//  双共役勾配法の前処理用不完全LU分解
//
//     Ver. 1.4 2008.7.23  K. Watanabe
/* *********************************************************** */
/*
//  input
//      A  : 分解前の行列
//      ga : 加速係数
//  output
//    LU: 分解後の行列（ピボット選択による行の入れ替えに注意）
//    swapData[i] : i行目が元の行列の何行目だったかを記憶する配列
//                  あらかじめ領域を確保しておくこと。pivotFlag==trueのときのみ使用
*/

	
	try{
		const int n = A.numRow;

		/*とりあえず分解後の行列LUにAをコピー*/
		LU = A;
		
		/*LUの列データのソート*/
		LU.sortColumn();
		
		/*対角項がゼロになるのを防ぐために対角項を少し大きくする*/
		for(int i = 0; i < n; i++){
			int numNZ = LU.numNZ[i];
			for(int mu = 0; mu < numNZ; mu++){
				if(LU.column[i][mu] == i){
					LU.value[i][mu] *= ga;
				}
			}
		}
		
		
		/*各行の未処理項の先頭（最も行番号が若い）を指すポインタ*/
		int**    fpColumn = new    int*[n];
		double** fpValue  = new double*[n];
		for(int i = 0; i < n; i++){
			fpColumn[i] = LU.column[i]; /* =&(LU.column[i][0]) */
			fpValue [i] = LU.value [i]; /* =&(LU.value [i][0]) */
		}
		
		/*分解開始*/
		for(int i = 0; i < n; i++){
//			cout << "now in" << i << "-th row" << endl;
			/*for debug */
//			if(*fpColumn[i] != i){
//				cout << "unknown Error!!" << endl;
//				cout << " i = " << i << " fpColumn[i]=" << *fpColumn[i];
//				cout << " swapData[i]=" << swapData[i];
//				cout << " numNZ[i]=" << LU.numNZ[i] << endl;
//				break;
//			}
			
			/*ピボット選択*/
			if(pivotFlag){
				double maxDiag    = 0.0;
				int maxDiagRow = -1;
				for(int ii = i; ii < n; ii++){
					if(*fpColumn[ii] == i){/*非零要素発見*/
						double mm = fabs(*fpValue[ii]);
						if(mm > maxDiag){
							maxDiag = mm;
							maxDiagRow = ii;
						}
					}
				}
				
				if(maxDiagRow == -1 || maxDiag == 0.0){
					cout << "diagnal[" << i << "] = 0!!" << endl;
					exit(1);
				}
				
				/*ピボット行入れ替え*/
				int*    intptmp = LU.column[i];
				double* doubleptmp = LU.value [i];
				LU.column[i]  = LU.column[maxDiagRow];
				LU.value [i]  = LU.value [maxDiagRow];
				LU.column[maxDiagRow] =    intptmp;
				LU.value [maxDiagRow] = doubleptmp;
				
				   intptmp = fpColumn[i];
				doubleptmp = fpValue [i];
				fpColumn[i] = fpColumn[maxDiagRow];
				fpValue [i] = fpValue [maxDiagRow];
				fpColumn[maxDiagRow] = intptmp;
				fpValue [maxDiagRow] = doubleptmp;
				
				int inttmp = LU.numNZ[i];
				LU.numNZ[i] = LU.numNZ[maxDiagRow];
				LU.numNZ[maxDiagRow] = inttmp;

				inttmp = LU.szColumn[i];
				LU.szColumn[i] = LU.szColumn[maxDiagRow];
				LU.szColumn[maxDiagRow] = inttmp;

				swapData[i] = maxDiagRow;
			}
			else{
				swapData[i] = i;
			}
			
			
			/* [L] および[U]の作成開始 */
			for(int ii = i+1; ii < n; ii++){
//				cout << "now in" << ii << "-th row" << endl;
				if(*fpColumn[ii] == i){/*非零要素発見*/
					double m = *fpValue[ii] / *fpValue[i];
					*fpValue[ii] = m; /*[L]作成*/
					
					fpColumn[ii]++;
					fpValue [ii]++;
					
					/*[U]作成*/
					int numNZi  = LU.numNZ[i];
					int numNZii = LU.numNZ[ii];
					int mui, muii;
					for(mui = 0; (mui < numNZi)&&(LU.column[i][mui] <= i); mui++){
						/*カウンタmuiがi列より大きくなるまでカウントアップ*/
					}
					for(muii = 0; (muii < numNZii)&&(LU.column[ii][muii] <= i); muii++){
						/*カウンタmuiiがi列より大きくなるまでカウントアップ*/
					}
					while((mui < numNZi) && (muii < numNZii)){
//						cout << "mui=" << mui << " muii=" << muii << endl;
						int ik  = LU.column[i ][mui];
						int iik = LU.column[ii][muii];
						if(ik == iik){
							LU.value[ii][muii] -= m*LU.value[i][mui];
							mui++;
							muii++;
						}
						else if(ik > iik){
							muii++;
						}
						else{ /* if(ik < iik) */
							mui++;
						}
					}
				}
			}
		}

		delete[] fpColumn;
		delete[] fpValue;
	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		exit(1);
	}
}

/* ******************************************************************* */
void GenMatrix::iLUSolv(const GenMatrix& LU, const double* b, double* x){
//              Solve system of linear equations
//                  [L*U]*X = B
//             LU is incomplete LU decomposition
//   Bi-CG法用の不完全LU分解された前処理行列に対して連立方程式を解く
//   不完全LU分解はiLCcmpを使用すること（columnデータがソートされている
//    必要がある） （ピボット選択なし）
//   ピボット選択をしていない場合はこっちの方が高速
//     ver. 1.2 2008.7.23  K. Watanabe
/* ******************************************************************* */

	
	try{
		const int n = LU.numRow;
		double* y       = new double[n];
		int*    diagIdx = new    int[n];/*対角項の位置を示すポインタ*/
		
/*     ---- forward substitution ----*/
		/* solve [L]{y} = {b} */
		for(int i = 0; i < n; i++){
			double t = b[i];
			int*   columnP = LU.column[i];
			double* valueP =  LU.value[i];
			int j, mu;
			for(mu = 0; (j=columnP[mu]) < i; mu++){
				t -= valueP[mu] * y[j];
			} /* がj=iでforループ脱出*/

			diagIdx[i] = mu; /*対角項の場所を記憶*/
			y[i] = t; /*[L]の対角項＝1*/
		}
	
		/* init x[] */
		for(int i = 0; i < n; i++){
			x[i] = 0.0;
		}

/*     ---- backward substitution ----*/
		/* solve [U]{xswap} = {y} */
		for(int i = n-1; i >= 0; i--){
			int     numNZ  =  LU.numNZ[i];
			double* valueP =  LU.value[i];
			double       t =  y[i];
			for(int mu = diagIdx[i]+1; mu < numNZ; mu++){
				t -= valueP[mu] * x[i];
			}
			x[i] = t / valueP[diagIdx[i]];
		}
		
		delete[] y;
		delete[] diagIdx;
	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(...){
		cout << "unknown Error" << endl;
		exit(1);
	}
}

/* ******************************************************************* */
void GenMatrix::iLUSolv(const GenMatrix& LU, const double* b, double* x
                        ,int* swapData){
//              Solve system of linear equations
//                  [L*U]*X = B
//             LU is incomplete LU decomposition
//   Bi-CG法用の不完全LU分解された前処理行列に対して連立方程式を解く
//   不完全LU分解はiLCcmpを使用すること（columnデータがソートされている
//    必要がある） （ピボット選択あり）
//     ver. 1.2 2008.7.23  K. Watanabe
/* ******************************************************************* */

	
	try{
		const int n = LU.numRow;
		double* y       = new double[n];
		double* xswap   = new double[n];
		int*    diagIdx = new    int[n];/*対角項の位置を示すポインタ*/
		
/*     ---- forward substitution ----*/
		/* solve [L]{y} = {b} */
		for(int i = 0; i < n; i++){
			double t = b[swapData[i]];
			int j, mu;
			for(mu = 0; (j=LU.column[i][mu]) < i; mu++){
				t -= LU.value[i][mu] * y[j];
			} /* がj=iでforループ脱出*/

			diagIdx[i] = mu; /*対角項の場所を記憶*/
			y[i] = t; /*[L]の対角項＝1*/
		}
	
		/* init xswap[] */
		for(int i = 0; i < n; i++){
			xswap[i] = 0.0;
		}

/*     ---- backward substitution ----*/
		/* solve [U]{xswap} = {y} */
		for(int i = n-1; i >= 0; i--){
			int numNZ = LU.numNZ[i];
			double t = y[i];
			for(int mu = diagIdx[i]+1; mu < numNZ; mu++){
				t -= LU.value[i][mu] * xswap[i];
			}
			xswap[i] = t / LU.value[i][diagIdx[i]];
		}
		
		/* 解の並べ替え*/
		for(int i = 0; i < n; i++){
			x[swapData[i]] = xswap[i];
		}
		
		
		delete[] y;
		delete[] xswap;
		delete[] diagIdx;
	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(...){
		cout << "unknown Error" << endl;
		exit(1);
	}
}

double GenMatrix::JacobiSmoother(const GenMatrix& A, const double* b, double* x,
                              double* r, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   Jacobi Method で連立方程式を解く for multigrid smoother
//     Ver. 1.00 2010.3.1  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of smoothing
//
//    ==== input ====
//      r[i].... residual vector
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は残差ノルム
*/
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		normR += r[i]*r[i];
	}
	normR = sqrt(normR);
	double normR0 = normR;
	cout << "first norm[r] = " << normR0 << endl;
	double minNormR = normR0;
	
	int k;
	for(k = 0; k < maxItr; k++){
	
		for(int i = 0; i < n; i++){
//			x[i] = (r[i] + x[i]/invDiag[i])*invDiag[i];
			x[i] += r[i]*invDiag[i];
		}

		double normR = 0.0;
		productAX(A, x, r);
		for(int i = 0; i < n; i++){
			r[i] = b[i] - r[i];
			normR += r[i]*r[i];
		}
		
		normR = sqrt(normR);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
	
	}

	/***free memory***/
	delete[] invDiag;

	
	/***  termination ***/
	return normR;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

/* ******************************************************************* */
int GenMatrix::gsSolv(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, bool increStopFlag, int incCount){
//              Solve system of linear equations
//                  [A]*X = B
//   ガウスザイデル法で連立方程式を解く
//     Ver. 1.02 2004.10.26  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     increStopFlag... 残差が増加したときに反復を停止するかどうか(default = false)
//     incCount...残差がこの回数増加したら停止する(default = 1)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double eps2 = eps * eps;
	double* r   = new double[n];
	int numItr = -1;
	double r2sum_old;
	int numInc = 0;
	int k;
	for(k = 0; k < maxItr; k++){
	
		for(int i = 0; i < n; i++){
			r[i] = b[i];
		}
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				r[i] -= ap[mu] * x[j];
			}
			x[i] += r[i]*invDiag[i]; /*x[i] = (r[i] + a.diag[i]*x[i])/a.diag[i] */
			r2sum += r[i]*r[i];
		}
#if ___DBGPRINTG
		double normR0;
		if(k == 0){
			normR0 = sqrt(r2sum);
			cout << "first norm[r] = " << normR0 << endl;
		}
		else{
			cout << k << " :Normalized|r| = " << sqrt(r2sum)/normR0 << endl;
		}
#endif

		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << sqrt(r2sum) << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		if(increStopFlag && k != 0 && r2sum > r2sum_old){
			numInc++;
			if(numInc == incCount){
				cout << "norm of residual increase!! stop iteration" << endl;
				cout << "norm{r} = " << sqrt(r2sum) << endl;
				cout << "Number of iteration = " << k+1 << endl;
				numItr = k+1;
				convFlag = 1;
				break;
			}
		}
		r2sum_old = r2sum;
	
	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}
	
	/***free memory***/
	delete[] r;
	delete[] invDiag;
	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}


/* ******************************************************************* */
void GenMatrix::gsPrecondition(GenMatrix& A, const double* p,
                               double* pt, double* q){
//  Gauss Seidel type preconditioner
//       [U] inv([L]+[D]){p} -> {q}
//    step1. solve ([L]+[D]){pt} == {p}
//    step2. [U]{pt} -> {q}
//   
//     Ver. 1.00 2010.8.25  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     p[i]
//
//    ==== output ====
//     pt[i]
//      q[i]
//
//
*/
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	/* sorting */
	if(A.sortFlag != 1){
		A.sortColumn();
	}

	int* columnIdx = new int[n];
	for(int i = 0; i < n; i++){
		columnIdx[i] = -1;
		q[i] = 0.0;
	}
	
	/* **** step1. solve ([L]+[D]){pt} == {p} **** */
	for(int i = 0; i < n; i++){
		pt[i] = p[i];

		int numNZ  = A.numNZ[i];
		double* ap = A.value[i];
		int*    cp = A.column[i];
		for(int mu = 0; mu < numNZ; mu++){
			int j = cp[mu];
			if(j == i){
				columnIdx[i] = mu;
				break;
			}
			pt[i] -= ap[mu] * pt[j];
		}
		if(columnIdx[i] == -1){
			cout << "error: diagonal part is missing" << endl;
		}
		pt[i] /= ap[columnIdx[i]];
		columnIdx[i]++;
	}
	
	/* **** step2. [U]{pt} -> {q} **** */
	for(int i = 0; i < n; i++){

		int numNZ  = A.numNZ[i];
		double* ap = A.value[i];
		int*    cp = A.column[i];
		for(int mu = columnIdx[i]; mu < numNZ; mu++){
			int j = cp[mu];
			q[i] += ap[mu] * pt[j];
		}
	}
	
	delete[] columnIdx;
	
	return;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(...){
		cout << "unknown Error" << endl;
		throw;
	}
}



/* ******************************************************************* */
int GenMatrix::IDR_sSolv(const GenMatrix& A, const double* b, double* x,
                        const int s, double eps, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR Method で連立方程式を解く
//     Ver. 1.00 2010.5.24  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	if(s < 1){
		cout << "parameter s should be s >= 1 " << endl;
		return -1;
	}

	
	/*係数行列から対角項を抽出*/
//	double* invDiag = new double[n];
//	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
//	for(int i = 0; i < n ; i++){
//		if(invDiag[i] == 0.0){
//			cout << "error! diagnal [" << i << " ] is zero" << endl;
//		}
//		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
//	}
	
	double eps2 = eps * eps;

	double* r   = new double[n];
	double* v   = new double[n];
	double* t   = new double[n];
	double* q   = new double[n];
	double* e   = new double[n];

	double** E = new double*[s];
	double** P = new double*[s];
	double** Q = new double*[s];
	for(int i = 0; i < s; i++){
		E[i] = new double[n];
		P[i] = new double[n];
		Q[i] = new double[n];
	}


	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int j = 0; j < n; j++){
		r[j] = b[j] - r[j];
		normR0 += r[j]*r[j];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;

	/* make P matrix */
//	GenMatrix PG(s, n);

	for(int j = 0; j < n; j++){
		P[0][j] = r[j];
//		P[0][j] = (double)rand();
//		PG.set(0, j, r[j]);
//		PG.set(0, j, (double)rand());
	}
	for(int i = 1; i < s; i++){
		for(int j = 0; j < n; j++){
//			P[i][j] = 0.0;
			P[i][j] = (double)rand();
		}
//		P[i][i] = 1.0;
//		PG.set(i, i, 1.0);
	}
//	PG.modGramSchmidt();
	GenMatrix::modGramSchmidt(P, s, n);

//	GenMatrix::writeFile("P.d", PG);
//	for(int i = 0; i < s; i++){
//		for(int j = 0; j < n; j++){
//			P[i][j] = PG.get(i, j);
//		}
//	}

/*
	int* num = new int[s];
	for(int i = 0; i < s; i++){
		num[i] = n/s;
		if(i < (n%s)){
			num[i]++;
		}
	}
	for(int i = 0; i < s; i++){
		for(int j = 0; j < n; j++){
			P[i][j] = 0.0;
		}
	}
	for(int j = 0; j < n; j++){
		PG.set(j%s, j, 1.0/num[j%s]);
		P[j%s][j] = 1.0/num[j%s];
	}
	delete[] num;
*/
	/*initial loop */
	double omega = 0.0;
	for(int k = 0; k < s; k++){
		GenMatrix::productAX(A, r, v);
		
		double vr = 0.0;
		double vv = 0.0;
		for(int j = 0; j < n; j++){
			vr += v[j]*r[j];
			vv += v[j]*v[j];
		}
		omega = vr/vv;
		
		for(int j = 0; j < n; j++){
			double qj =  omega*r[j];
			double ej = -omega*v[j];
			r[j] += ej;
			x[j] += qj;
			E[s-1-k][j] = ej;
			Q[s-1-k][j] = qj;
		}
		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif

	}
	
	


	double* c   = new double[s];
	double* ptr = new double[s];
//	double* rs  = new double[s];

	int k;
	for(k = s; k < maxItr; k++){
		
		/*make PtE */
		GenMatrix PtE(s, s);
		for(int i = 0; i < s; i++){
			for(int j = 0; j < s; j++){
				double pte = 0.0;
				for(int k = 0; k < n; k++){
					pte += P[i][k] * E[j][k];
				}
				PtE.set(i, j, pte);
			}
		}
//	GenMatrix::writeFile("PtE.d", PtE);
//	exit(1);
		
		GenMatrix PtEinv;
		PtEinv = PtE;
		GenMatrix::inv(PtEinv);
//	GenMatrix::writeFile("PtEinv.d", PtEinv);

		
		/*solve PtE{c} = Pt{r} */
		double ptr_norm = 0.0;
		for(int i = 0; i < s; i++){
			ptr[i] = 0.0;
			for(int j = 0; j < n; j++){
				ptr[i] += P[i][j]*r[j];
			}
			ptr_norm += ptr[i]*ptr[i];
			c[i] = 0.0;
		}
//		ptr_norm = sqrt(ptr_norm);
//		int itrs = GenMatrix::iSORSolv(PtE, ptr, c, eps*ptr_norm, maxItr, 0.1);
//		cout << "number of itr for c = " << itrs << endl;
		GenMatrix::productAX(PtEinv, ptr, c);
			
		/* v = r - E{c} */
		for(int j = 0; j < n; j++){
			v[j] = r[j];
		}
		for(int i = 0; i < s; i++){
			for(int j = 0; j < n; j++){
				v[j] -= E[i][j]*c[i];
			}
		}
		
		if(k%(s+1) == s){
			GenMatrix::productAX(A, v, t);
			
			double tv = 0.0;
			double tt = 0.0;
			for(int j = 0; j < n; j++){
				tv += t[j]*v[j];
				tt += t[j]*t[j];
			}
			omega = tv/tt;
			
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
				e[j] = -omega*t[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
					e[j] -= E[i][j]*c[i];
				}
			}
		}
		else{
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
				}
			}
			GenMatrix::productAX(A, q, e);
			for(int j = 0; j < n; j++){
				e[j] = -e[j];
			}
		}
		
		/*re-construct E={e(k), e(k-1), ...,e(k+1-s)} and Q*/
//for(int i = 0; i < s; i++){
//	cout << " E[" << i << "]=" << E[i];
//}
//cout << endl;

		double* e0p = E[s-1];
		double* q0p = Q[s-1];
		for(int i = s-1; i > 0; i--){
			E[i] = E[i-1];
			Q[i] = Q[i-1];
		}
		E[0] = e0p;
		Q[0] = q0p;
		
		for(int j = 0; j < n; j++){
			r[j] += e[j];
			x[j] += q[j];
			E[0][j] = e[j];
			Q[0][j] = q[j];
		}
//for(int i = 0; i < s; i++){
//	cout << " E[" << i << "]=" << E[i];
//}
//cout << endl;
//exit(1);		
		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
	delete[] r;
	delete[] v;
	delete[] t;
	delete[] q;
	delete[] e;
	for(int i = 0; i < s; i++){
		delete[] E[i];
		delete[] P[i];
		delete[] Q[i];
	}
	
	delete[] E;
	delete[] P;
	delete[] Q;

	delete[] c;
	delete[] ptr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}


/* ******************************************************************* */
double  GenMatrix::IDR_sSmoother(const GenMatrix& A, const double* b, double* x,
                          double* r, const int s, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR Method で連立方程式を解く Multigrid smoother
//     Ver. 1.00 2010.5.24  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は残差ノルム
*/
	int convFlag = 0;
	double r2sum, normR;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	if(s < 1){
		cout << "parameter s should be s >= 1 " << endl;
		return -1;
	}

	
	/*係数行列から対角項を抽出*/
//	double* invDiag = new double[n];
//	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
//	for(int i = 0; i < n ; i++){
//		if(invDiag[i] == 0.0){
//			cout << "error! diagnal [" << i << " ] is zero" << endl;
//		}
//		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
//	}
	
//	double eps2 = eps * eps;

//	double* r   = new double[n];
	double* v   = new double[n];
	double* t   = new double[n];
	double* q   = new double[n];
	double* e   = new double[n];

	double** E = new double*[s];
	double** P = new double*[s];
	double** Q = new double*[s];
	for(int i = 0; i < s; i++){
		E[i] = new double[n];
		P[i] = new double[n];
		Q[i] = new double[n];
	}


	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int j = 0; j < n; j++){
		r[j] = b[j] - r[j];
		normR0 += r[j]*r[j];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;

	/* make P matrix */
//	GenMatrix PG(s, n);

	for(int j = 0; j < n; j++){
		P[0][j] = r[j];
//		P[0][j] = (double)rand();
//		PG.set(0, j, r[j]);
//		PG.set(0, j, (double)rand());
	}
	for(int i = 1; i < s; i++){
		for(int j = 0; j < n; j++){
//			P[i][j] = 0.0;
			P[i][j] = (double)rand();
		}
//		P[i][i] = 1.0;
//		PG.set(i, i, 1.0);
	}
//	PG.modGramSchmidt();
	GenMatrix::modGramSchmidt(P, s, n);

//	GenMatrix::writeFile("P.d", PG);
//	for(int i = 0; i < s; i++){
//		for(int j = 0; j < n; j++){
//			P[i][j] = PG.get(i, j);
//		}
//	}

	
	/*initial loop */
	double omega = 0.0;
	for(int k = 0; k < s; k++){
		GenMatrix::productAX(A, r, v);
		
		double vr = 0.0;
		double vv = 0.0;
		for(int j = 0; j < n; j++){
			vr += v[j]*r[j];
			vv += v[j]*v[j];
		}
		omega = vr/vv;
		
		for(int j = 0; j < n; j++){
			double qj =  omega*r[j];
			double ej = -omega*v[j];
			r[j] += ej;
			x[j] += qj;
			E[s-1-k][j] = ej;
			Q[s-1-k][j] = qj;
		}
	}
	
	


	double* c   = new double[s];
	double* ptr = new double[s];
//	double* rs  = new double[s];

	int k;
	for(k = s; k < maxItr; k++){
		
		/*make PtE */
		GenMatrix PtE(s, s);
		for(int i = 0; i < s; i++){
			for(int j = 0; j < s; j++){
				double pte = 0.0;
				for(int k = 0; k < n; k++){
					pte += P[i][k] * E[j][k];
				}
				PtE.set(i, j, pte);
			}
		}
//	GenMatrix::writeFile("PtE.d", PtE);
//	exit(1);
		
		GenMatrix PtEinv;
		PtEinv = PtE;
		GenMatrix::inv(PtEinv);
//	GenMatrix::writeFile("PtEinv.d", PtEinv);

		
		/*solve PtE{c} = Pt{r} */
		double ptr_norm = 0.0;
		for(int i = 0; i < s; i++){
			ptr[i] = 0.0;
			for(int j = 0; j < n; j++){
				ptr[i] += P[i][j]*r[j];
			}
			ptr_norm += ptr[i]*ptr[i];
			c[i] = 0.0;
		}
//		ptr_norm = sqrt(ptr_norm);
//		double eps = 1.0e-6;
//		int itrs = GenMatrix::iSORSolv(PtE, ptr, c, eps*ptr_norm, maxItr, 0.1);
//		cout << "number of itr for c = " << itrs << endl;
		GenMatrix::productAX(PtEinv, ptr, c);
			
		/* v = r - E{c} */
		for(int j = 0; j < n; j++){
			v[j] = r[j];
		}
		for(int i = 0; i < s; i++){
			for(int j = 0; j < n; j++){
				v[j] -= E[i][j]*c[i];
			}
		}
		
		if(k%(s+1) == s){
			GenMatrix::productAX(A, v, t);
			
			double tv = 0.0;
			double tt = 0.0;
			for(int j = 0; j < n; j++){
				tv += t[j]*v[j];
				tt += t[j]*t[j];
			}
			omega = tv/tt;
			
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
				e[j] = -omega*t[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
					e[j] -= E[i][j]*c[i];
				}
			}
		}
		else{
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
				}
			}
			GenMatrix::productAX(A, q, e);
			for(int j = 0; j < n; j++){
				e[j] = -e[j];
			}
		}
		
		/*re-construct E={e(k), e(k-1), ...,e(k+1-s)} and Q*/
		double* e0p = E[s-1];
		double* q0p = Q[s-1];
		for(int i = s-1; i > 0; i--){
			E[i] = E[i-1];
			Q[i] = Q[i-1];
		}
		E[0] = e0p;
		Q[0] = q0p;
		
		for(int j = 0; j < n; j++){
			r[j] += e[j];
			x[j] += q[j];
			E[0][j] = e[j];
			Q[0][j] = q[j];
		}
		
		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}
		normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
//		/*収束判定*/
//		if(r2sum <= eps2){
//			cout << "---- converged ------" << endl;
//			cout << "norm{r} = " << normR << endl;
//			cout << "Number of iteration = " << k+1 << endl;
//			numItr = k+1;
//			convFlag = 1;
//			break;
//		}
	}

	/*** not convergence**/
//	if (convFlag == 0){
//		cout << "--- not converged ---" << endl;
//		cout << "norm{r} = " << sqrt(r2sum) << endl;
//		cout << "Number of iteration = " << k << endl;
//	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
//	delete[] r;
	delete[] v;
	delete[] t;
	delete[] q;
	delete[] e;
	for(int i = 0; i < s; i++){
		delete[] E[i];
		delete[] P[i];
		delete[] Q[i];
	}
	
	delete[] E;
	delete[] P;
	delete[] Q;

	delete[] c;
	delete[] ptr;

	
	/***  termination ***/
	return normR;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}



/* ******************************************************************* */
int GenMatrix::AdaptiveIDR_sJacobi(const GenMatrix& A, const double* b, double* x,
                        const int s, double eps, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR(s) based adaptive Jacobi Method で連立方程式を解く
//     Ver. 1.00 2010.5.24  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	if(s < 1){
		cout << "parameter s should be s >= 1 " << endl;
		return -1;
	}

	double* diag = new double[n];
	getDiagVector(A, diag);

	/*処理速度向上のため、対角項の逆数も保存*/
	double* invDiag = new double[n];
	for(int i = 0; i < n ; i++){
		if(diag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/diag[i];
	}

	const double eps2 = eps * eps;

	double* r   = new double[n];
	double* t   = new double[n];
	double* at  = new double[n];
	double* sv  = new double[n];
	double* as  = new double[n];
	double* ec  = new double[n];
	double* qc  = new double[n];

	double* dx  = new double[n];
	double* dr  = new double[n];

	double** E = new double*[s];
	double** P = new double*[s];
	double** Q = new double*[s];
	for(int i = 0; i < s; i++){
		E[i] = new double[n];
		P[i] = new double[n];
		Q[i] = new double[n];
	}


	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int j = 0; j < n; j++){
		r[j] = b[j] - r[j];
		normR0 += r[j]*r[j];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;


	for(int j = 0; j < n; j++){
		P[0][j] = r[j];
//		P[0][j] = (double)rand();
//		PG.set(0, j, r[j]);
//		PG.set(0, j, (double)rand());
	}
	for(int i = 1; i < s; i++){
		for(int j = 0; j < n; j++){
//			P[i][j] = 0.0;
			P[i][j] = (double)rand();
		}
//		P[i][i] = 1.0;
//		PG.set(i, i, 1.0);
	}
//	PG.modGramSchmidt();
	GenMatrix::modGramSchmidt(P, s, n);

	/*initial loop */
	double gamma = 0.0;
	double omega = 0.0;
	for(int k = 0; k < s; k++){
		
		if(k >=1){
			for(int j = 0; j < n; j++){
				t[j] = r[j] - gamma*dr[j];
			}
		}
		else{
			for(int j = 0; j < n; j++){
				t[j] = r[j];
			}
		}

		GenMatrix::productAX(A, t, at);
		
		double omega1 = 0.0;
		double omega2 = 0.0;
		for(int j = 0; j < n; j++){
			double invDat = invDiag[j] * at[j];
			omega1 += invDat * t[j];
			omega2 += invDat * invDat;
		}
		omega = omega1/omega2;

		for(int j = 0; j < n; j++){
			sv[j] = omega * invDiag[j] * t[j];
		}

		GenMatrix::productAX(A, sv, as);
		
		double pr = 0.0;
		double pe = 0.0;
		for(int j = 0; j < n; j++){
			dx[j] =  sv[j] - gamma*dx[j];
			dr[j] = -as[j] + sv[j]*diag[j]/omega -r[j];

			x[j] += dx[j];
			r[j] += dr[j];

			Q[k][j] = dx[j];
			E[k][j] = dr[j];

			pr += P[0][j]* r[j];
			pe += P[0][j]*dr[j];
		}
		gamma = pr/pe;

		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
	}

	double* c   = new double[s];
	double* ptr = new double[s];
//	double* rs  = new double[s];

	int k;
	for(k = s; k < maxItr; k++){
		int ks = k%s;// 0,1,2, ... s-1, 0, 1, 2, ...
		
		/*make PtE */
		GenMatrix PtE(s, s);
		for(int i = 0; i < s; i++){
			for(int j = 0; j < s; j++){
				double pte = 0.0;
				for(int k = 0; k < n; k++){
					pte += P[i][k] * E[j][k];
				}
				PtE.set(i, j, pte);
			}
		}
//	GenMatrix::writeFile("PtE.d", PtE);
//	exit(1);
		
		GenMatrix PtEinv;
		PtEinv = PtE;
		GenMatrix::inv(PtEinv);
//	GenMatrix::writeFile("PtEinv.d", PtEinv);

		
		/*solve PtE{c} = Pt{r} */
		double ptr_norm = 0.0;
		for(int i = 0; i < s; i++){
			ptr[i] = 0.0;
			for(int j = 0; j < n; j++){
				ptr[i] += P[i][j]*r[j];
			}
//			ptr_norm += ptr[i]*ptr[i];
//			c[i] = 0.0;
		}
//		ptr_norm = sqrt(ptr_norm);
//		int itrs = GenMatrix::iSORSolv(PtE, ptr, c, eps*ptr_norm, maxItr, 0.1);
//		cout << "number of itr for c = " << itrs << endl;
		GenMatrix::productAX(PtEinv, ptr, c);
			
		for(int j = 0; j < n; j++){
			ec[j] = 0.0;
			qc[j] = 0.0;
		}
		for(int i = 0; i < s; i++){
			for(int j = 0; j < n; j++){
				ec[j] += E[i][j]*c[i];
				qc[j] += Q[i][j]*c[i];
			}
		}

		for(int j = 0; j < n; j++){
			t[j] = r[j] - gamma*ec[j];
		}

		GenMatrix::productAX(A, t, at);

		double omega1 = 0.0;
		double omega2 = 0.0;
		for(int j = 0; j < n; j++){
			double invDat = invDiag[j] * at[j];
			omega1 += invDat * t[j];
			omega2 += invDat * invDat;
		}
		omega = omega1/omega2;

		for(int j = 0; j < n; j++){
			sv[j] = omega * invDiag[j] * (r[j] - ec[j]);
		}

		GenMatrix::productAX(A, sv, as);
		
		double pr = 0.0;
		double pe = 0.0;
		for(int j = 0; j < n; j++){
			double qj =   sv[j] - qc[j];
			double ej = -as[j] + sv[j]*diag[j]/omega -r[j];

			x[j] += qj;
			r[j] += ej;

			Q[ks][j] = qj;
			E[ks][j] = ej;

			pr += P[0][j]*r[j];
			pe += P[0][j]*ej;
		}
		gamma = pr/pe;

		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/

	delete[] diag;
	delete[] invDiag;
	
	delete[] r;
	delete[] t;
	delete[] at;
	delete[] sv;
	delete[] as;
	delete[] ec;
	delete[] qc;
	
	delete[] dx;
	delete[] dr;
	
	for(int i = 0; i < s; i++){
		delete[] E[i];
		delete[] P[i];
		delete[] Q[i];
	}
	
	delete[] E;
	delete[] P;
	delete[] Q;

	delete[] c;
	delete[] ptr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}


/* ******************************************************************* */
int GenMatrix::iSORSolv(const GenMatrix& A, const double* b, double* x,
                        double eps, int maxItr, double omega, bool increStopFlag){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based SOR Method で連立方程式を解く
//     Ver. 1.00 2010.1.30  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     omega... 加速係数
//     increStopFlag... 残差が増加したときに反復を停止するかどうか
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	double r2sum_old;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double eps2 = eps * eps;
	double* r   = new double[n];
	double* p   = new double[n];
	double* s   = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];
	int numItr = -1;
	double r2sum_old;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		normR0 += r[i]*r[i];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;
	double minNormR = normR0;
	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double gamma = 0.0;

	int k;
	for(k = 0; k < maxItr; k++){
	
		/* cal {s} */
		for(int i = 0; i < n; i++){
			s[i] = r[i] + gamma*dr[i];
		}
		for(int i = 0; i < n; i++){
			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				if(j < i){
					s[i] -= ap[mu] * s[j];
				}
			}
			s[i] *= omega * invDiag[i];
		}
		
		/* cal dx and dr */
		for(int i = 0; i < n; i++){
			dx[i] = s[i] + gamma * dx[i];
			dr[i] = -r[i];

			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				if(j > i){
					dr[i] -= ap[mu] * s[j];
				}
			}
			dr[i] -= (1.0 - 1.0/omega)/invDiag[i] * s[i];
		}
		
		/* cal update r and x */
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum += r[i]*r[i];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		if(increStopFlag && k != 0 && r2sum > r2sum_old){
			cout << "norm of residual increase!! stop iteration" << endl;
			cout << "norm{r} = " << sqrt(r2sum) << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		r2sum_old = r2sum;
		
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
			pdr += p[i] * dr[i];
		}
		gamma = -pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
	delete[] invDiag;
	delete[] r;
	delete[] s;
	delete[] dx;
	delete[] dr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

/* ******************************************************************* */
int GenMatrix::iSORSolvBest(const GenMatrix& A, const double* b, double* x,
                        double eps, int maxItr, double omega){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based SOR Method で連立方程式を解く 最小残差解を返す
//     Ver. 1.00 2010.5.22  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     omega... 加速係数
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	double r2sum_old;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	double* bestX = new double[n];
	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double eps2 = eps * eps;
	double* r   = new double[n];
	double* p   = new double[n];
	double* s   = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];
	int numItr = -1;
	double r2sum_old;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		bestX[i] = x[i];
		r[i] = b[i] - r[i];
		normR0 += r[i]*r[i];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;
	double minNormR = normR0;
	double bestNormR  = normR0;
	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double gamma = 0.0;

	int k;
	for(k = 0; k < maxItr; k++){
	
		/* cal {s} */
		for(int i = 0; i < n; i++){
			s[i] = r[i] + gamma*dr[i];
		}
		for(int i = 0; i < n; i++){
			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				if(j < i){
					s[i] -= ap[mu] * s[j];
				}
			}
			s[i] *= omega * invDiag[i];
		}
		
		/* cal dx and dr */
		for(int i = 0; i < n; i++){
			dx[i] = s[i] + gamma * dx[i];
			dr[i] = -r[i];

			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				if(j > i){
					dr[i] -= ap[mu] * s[j];
				}
			}
			dr[i] -= (1.0 - 1.0/omega)/invDiag[i] * s[i];
		}
		
		/* cal update r and x */
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum += r[i]*r[i];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		if(normR < bestNormR){
			bestNormR = normR;
			for(int i = 0; i < n; i++){
				bestX[i] = x[i];
			}
		}
			
		
		
		r2sum_old = r2sum;
		
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
			pdr += p[i] * dr[i];
		}
		gamma = -pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "best of norm{r} = " << bestNormR << endl;
		cout << "Normalized|r| = " << bestNormR/normR0 << endl;
		for(int i = 0; i < n; i++){
			x[i] = bestX[i];
		}

	}

#if ___DBGPRINTG
	/* check true residual */
//	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
	delete[] invDiag;
	delete[] r;
	delete[] s;
	delete[] dx;
	delete[] dr;
	delete[] bestX;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

/* ******************************************************************* */
int GenMatrix::adaptiveIJacobi(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based adaptive I-Jacobi Method で連立方程式を解く
//     Ver. 1.00 2010.1.30  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double eps2 = eps * eps;
	double* r   = new double[n];
	double* p   = new double[n];
	double* s   = new double[n];
	double* t   = new double[n];
	double* as  = new double[n];
	double* at  = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];
	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		normR0 += r[i]*r[i];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;
	double minNormR = normR0;
	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double beta = 0.0;
	

	int k;
	for(k = 0; k < maxItr; k++){
	
		/* cal {s} */
		for(int i = 0; i < n; i++){
			t[i] = r[i] - beta*dr[i];
		}
		GenMatrix::productAX(A, t, at);
		double omega1 = 0.0;
		double omega2 = 0.0;
		for(int i = 0; i < n; i++){
			double invDat = invDiag[i] * at[i];
			omega1 += invDat * t[i];
			omega2 += invDat * invDat;
		}
		double omega = omega1 / omega2;
#if ___DBGPRINTG
		cout <<"   omega = " << omega << endl;
#endif
		if(omega > 2.0){omega = 2.0;}
		if(omega < 0.1){omega = 0.1;}
		
		for(int i = 0; i < n; i++){
			s[i] = omega * invDiag[i] * t[i];
		}
		
		/* cal dx and dr */
		GenMatrix::productAX(A, s, as);
		for(int i = 0; i < n; i++){
			dx[i] = s[i] - beta * dx[i];
			dr[i] = -as[i] + s[i]/invDiag[i]/omega - r[i];
		}
		
		/* cal update r and x */
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum += r[i]*r[i];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
			pdr += p[i] * dr[i];
		}
		beta = pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
	delete[] invDiag;
	delete[] r;
	delete[] s;
	delete[] t;
	delete[] as;
	delete[] at;
	delete[] dx;
	delete[] dr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

/* ******************************************************************* */
double GenMatrix::adaptiveIJacobiSmoother(const GenMatrix& A, const double* b, double* x,
                              double* r, int maxItr){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based adaptive I-Jacobi Method で連立方程式を解く for multigrid smoother
//     Ver. 1.00 2010.1.31  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of smoothing
//
//    ==== input ====
//      r[i].... residual vector
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は残差ノルム
*/
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
//	double* x   = new double[n];

	double* p   = new double[n];
	double* s   = new double[n];
	double* t   = new double[n];
	double* as  = new double[n];
	double* at  = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];
	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
//		x[i] = xBest[i];
		r[i] = b[i] - r[i];
		normR += r[i]*r[i];
	}
	normR = sqrt(normR);
	double normR0 = normR;
	cout << "first norm[r] = " << normR0 << endl;
	double minNormR = normR0;
//	double normRbest = normR0;
//	int bestI = 0;
	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double beta = 0.0;

	int k;
	for(k = 0; k < maxItr; k++){
	
		/* cal {s} */
#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
		for(int i = 0; i < n; i++){
			t[i] = r[i] - beta*dr[i];
		}

		GenMatrix::productAX(A, t, at);
		double omega1 = 0.0;
		double omega2 = 0.0;
#if ___OPENMP_SUPPORT
#pragma omp parallel for reduction(+:omega1)
#endif
		for(int i = 0; i < n; i++){
			double invDat = invDiag[i] * at[i];
			omega1 += invDat * t[i];
		}
#if ___OPENMP_SUPPORT
#pragma omp parallel for reduction(+:omega2)
#endif
		for(int i = 0; i < n; i++){
			double invDat = invDiag[i] * at[i];
			omega2 += invDat * invDat;
		}
		double omega = omega1 / omega2;
#if ___DBGPRINTG
//		cout <<"   omega = " << omega << endl;
#endif
		if(omega > 2.0){omega = 2.0;}
		if(omega < 0.1){omega = 0.1;}
		
#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
		for(int i = 0; i < n; i++){
			s[i] = omega * invDiag[i] * t[i];
		}
		
		/* cal dx and dr */
		GenMatrix::productAX(A, s, as);

#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
		for(int i = 0; i < n; i++){
			dx[i] = s[i] - beta * dx[i];
			dr[i] = -as[i] + s[i]/invDiag[i]/omega - r[i];
		}
		
		/* cal update r and x */
#if ___OPENMP_SUPPORT
#pragma omp parallel for
#endif
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		double r2sum = 0.0;
#if ___OPENMP_SUPPORT
#pragma omp parallel for reduction(+:r2sum)
#endif
		for(int i = 0; i < n; i++){
			r2sum += r[i]*r[i];
		}
		normR = sqrt(r2sum);

#if ___DBGPRINTG
//		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
		cout << k << " norm|r| = " << normR << endl;
#endif
//		if(normRbest > normR){
//			bestI = k;
//			normRbest = normR;
//			for(int i = 0; i < n; i++){
//				xBest[i] = x[i];
//			}
//		}
		
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
#if ___OPENMP_SUPPORT
#pragma omp parallel  for reduction(+:pr)
#endif
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
		}

#if ___OPENMP_SUPPORT
#pragma omp parallel  for reduction(+:pdr)
#endif
		for(int i = 0; i < n; i++){
			pdr += p[i] * dr[i];
		}
		beta = pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}


#if ___DBGPRINTG
	/* check true residual */
//	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
//	productAX(A, x, r);
//	r2sum = 0.0;
//	for(int i = 0; i < n; i++){
//		r[i] = b[i] - r[i];
//		r2sum += r[i]*r[i];
//	}
//	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
//	delete[] x;

	delete[] invDiag;
	delete[] s;
	delete[] t;
	delete[] as;
	delete[] at;
	delete[] dx;
	delete[] dr;

	
	/***  termination ***/
	return normR;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

/* ******************************************************************* */
int GenMatrix::relaxedIJacobi(const GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, double omega){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based relaxed I-Jacobi Method で連立方程式を解く
//     Ver. 1.00 2010.1.30  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	double r2sum;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double eps2 = eps * eps;
	double* r   = new double[n];
	double* p   = new double[n];
	double* s   = new double[n];
	double* as  = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];
	int numItr = -1;

	/*cal r0 = b - Ax0 */
	double normR0 = 0.0;
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		normR0 += r[i]*r[i];
	}
	normR0 = sqrt(normR0);
	cout << "first norm[r] = " << normR0 << endl;
	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double beta = 0.0;
	

	int k;
	for(k = 0; k < maxItr; k++){
	
		/* cal {s} */
		for(int i = 0; i < n; i++){
			s[i] = omega * invDiag[i] * (r[i] - beta * dr[i]);
		}
		
		/* cal dx and dr */
		GenMatrix::productAX(A, s, as);
		for(int i = 0; i < n; i++){
			dx[i] = s[i] - beta * dx[i];
			dr[i] = -as[i] - r[i];
		}

		for(int i = 0; i < n; i++){
			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				if(j > i){
					dr[i] += ap[mu] * s[j] / omega;
				}
			}
		}
		
		/* cal update r and x */
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum += r[i]*r[i];
		}
		double normR = sqrt(r2sum);

#if ___DBGPRINTG
		cout << k << " :Normalized|r| = " << normR/normR0 << endl;
#endif
		
		/*収束判定*/
		if(r2sum <= eps2){
			cout << "---- converged ------" << endl;
			cout << "norm{r} = " << normR << endl;
			cout << "Number of iteration = " << k+1 << endl;
			numItr = k+1;
			convFlag = 1;
			break;
		}
		
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
			pdr += p[i] * dr[i];
		}
		beta = pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}

	/*** not convergence**/
	if (convFlag == 0){
		cout << "--- not converged ---" << endl;
		cout << "norm{r} = " << sqrt(r2sum) << endl;
		cout << "Number of iteration = " << k << endl;
	}

#if ___DBGPRINTG
	/* check true residual */
	cout << "     **norm|r| = " << sqrt(r2sum) << endl;
	productAX(A, x, r);
	r2sum = 0.0;
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
		r2sum += r[i]*r[i];
	}
	cout << " true{norm|r|} = " << sqrt(r2sum) << endl;
#endif
	
	/***free memory***/
	delete[] invDiag;
	delete[] r;
	delete[] s;
	delete[] as;
	delete[] dx;
	delete[] dr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}

void GenMatrix::restriction(const GenMatrix& aM, int* rM, int rSize, GenMatrix& asM){
/*     restriction oparation                                     */
/*             [R][A][R]t  => AS                                 */
/*                                                               */
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
/* ************************************************************* */
//    ==== input ====
//     A......Matrix[nxn]
//     R......int[>Rsize]: i-th row, R[i]-th column in Restrition Mat.
//
//      0  1  2  3
//    ----------- 
//   0| 1         |
//   1|       1   |    => [R] = {0,2,3}
//   2|          1|    => [P] = [R]t = {0,-1,1,2}
//
//    ==== output ====
//    AS:  dim(AS) = Rsize;

	try{

	const int n =  aM.getNumRow();
	const int m = asM.getNumRow();
	
	if(m != rSize){
		cout <<"error dimesion of matrix is different in restriction" << endl;
		return;
	}
	
	/* make [P] = [R]t */
	int* pM = new int[n];
	for(int i = 0; i < n; i++){
		pM[i] = -1;
	}
	
	for(int i = 0; i < m; i++){
		pM[rM[i]] = i;
	}
	
	/*cal. [R][A][R]t*/
	for(int i = 0; i < n; i++){
		int ii = pM[i];
		if(ii != -1){
			int numNZ = aM.numNZ[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = aM.column[i][mu];
				int jj = pM[j];
				if(jj != -1){
					asM.set(ii, jj, aM.value[i][mu]);
				}
			}
		}
	}
	
	delete[] pM;
	return;
	}
	catch(...){
		cout << "error in restriction" << endl;
		throw;
	}
	/***  termination ***/

}

void GenMatrix::restriction(const GenMatrix& aM, int* pM1, int* pM2,
                            int r1size, int r2size, GenMatrix& asM){
/*     restriction oparation                                     */
/*             [P1]t[A][P2]  => AS                               */
/*                                                               */
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
/* ************************************************************* */
//    ==== input ====
//     A......Matrix[nxn]
//     P1 = [R1]t, P2 = [R2]t
//      0  1  2  3
//    ----------- 
//   0| 1         |
//   1|       1   |    => [R] = {0,2,3}
//   2|          1|    => [P] = [R]t = {0,-1,1,2}
//
//    ==== output ====
//    AS:  Matrix[R1size x R2size]

	try{

	const int an = aM.getNumRow();
	
	
	/*cal. [P1]t[A][P2]*/
	for(int i = 0; i < an; i++){
		int ii = pM1[i];
		if(ii != -1){
			int numNZ = aM.numNZ[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = aM.column[i][mu];
				int jj = pM2[j];
				if(jj != -1){
					asM.set(ii, jj, aM.value[i][mu]);
				}
			}
		}
	}
	
	return;
	}
	catch(...){
		cout << "error in restriction" << endl;
		throw;
	}
	/***  termination ***/

}



void GenMatrix::restriction(const GenMatrix& aM, int* rM1, int r1Size,
                                                 int* rM2, int r2Size, GenMatrix& asM){
/*     restriction oparation                                     */
/*             [R][A][R]t  => AS                                 */
/*                                                               */
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
/* ************************************************************* */
//    ==== input ====
//     A......Matrix[nxn]
//     R......int[>Rsize]: i-th row, R[i]-th column in Restrition Mat.
//
//      0  1  2  3
//    ----------- 
//   0| 1         |
//   1|       1   |    => [R] = {0,2,3}
//   2|          1|    => [P] = [R]t = {0,-1,1,2}
//
//    ==== output ====
//    AS:  dim(AS) = Rsize;

	try{

	const int n  =  aM.getNumRow();
	const int m  = asM.getNumRow();
	const int mm = asM.getNumColumn();
	
	if((m != r1Size) || (mm != r2Size)){
		cout <<"error dimesion of matrix is different in restriction" << endl;
		return;
	}
	
	/* make [P] = [R]t */
	int* pM1 = new int[n];
	int* pM2 = new int[n];
	for(int i = 0; i < n; i++){
		pM1[i] = -1;
		pM2[i] = -1;
	}
	
	for(int i = 0; i < m; i++){
		pM1[rM1[i]] = i;
	}
	for(int i = 0; i < mm; i++){
		pM2[rM2[i]] = i;
	}
	
	/*cal. [R][A][R]t*/
	for(int i = 0; i < n; i++){
		int ii = pM1[i];
		if(ii != -1){
			int numNZ = aM.numNZ[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = aM.column[i][mu];
				int jj = pM2[j];
				if(jj != -1){
					asM.set(ii, jj, aM.value[i][mu]);
				}
			}
		}
	}
	
	delete[] pM1, pM2;
	return;
	}
	catch(...){
		cout << "error in restriction" << endl;
		throw;
	}
	/***  termination ***/

}

// *********************************************************
// Restriction  for Complex vector
//  [R]{x} => {rx}
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
// *********************************************************
void GenMatrix::restriction(int* rM, int rSize, double* x, double* rx){
//
// input:
//   rM....[R]
//   rSize ... # of row in [R]
//   x    .....vector{x}
//  output:
//   xr .......=[R]{x}
//
	for(int i = 0; i < rSize; i++){
		rx[i] = x[rM[i]];
	}
}

// *********************************************************
// Restriction for int vector
//  [R]{x} => {rx}
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
// *********************************************************
void GenMatrix::restriction(int* rM, int rSize, int* x, int* rx){
//
// input:
//   rM....[R]
//   rSize ... # of row in [R]
//   x    .....vector{x}
//  output:
//   xr .......=[R]{x}
//
	for(int i = 0; i < rSize; i++){
		rx[i] = x[rM[i]];
	}
}

// *********************************************************
// Prologation without init. 
//  [R]t{x} +=> {rx}
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
// *********************************************************
void GenMatrix::prolongation_add(int* rM, int rRow, double* x, double* px){
//
// input:
//   rM....[M]
//   rRow ...... #of row    in [R]
//   x    .....vector{x}
//  output:
//   xr .......=[R]t{x}
//
	
	for(int i = 0; i < rRow; i++){
		px[rM[i]] += x[i];
	}
}

// *********************************************************
// Prologation
//  [R]t{x} => {rx}
/*     ver. 1.01 2009.11.18  K. Watanabe                          */
// *********************************************************
// input:
//   rM....[M]
//   rRow ...... #of row    in [R]
//   rColumn ... #of column in [R]
//   x    .....vector{x}
//  output:
//   xr .......=[R]t{x}
//
void GenMatrix::prolongation(int* rM, int rRow, int rColumn, double* x, double* px){

	/*init*/
	for(int i = 0; i < rColumn; i++){
		px[i] = 0.0;
	}
	
	for(int i = 0; i < rRow; i++){
		px[rM[i]] = x[i];
	}
}

// *********************************************************
// matrix(compress form)-matrix multiplication
//  [R]t[A] => [RA]
/*     ver. 1.0 2010.4.30  K. Watanabe                  */
// *********************************************************
// input:
//   A....matrix[A] n X *
//   Rvalue....non-zero entry of [R]  m X n
//   RColumn ... column index
//  output:
//   rA .......=[R]t[A]
//
//(e.g.  Rt =|1 0 3 0 5|  => Rvalue = {1,3,5,   Rcolumn = {0,2,4}
//           |2 0 4 0 6|               2,4,6}   rn = 2, rm = 3 (not 5!!)
//
void GenMatrix::multipleRtA(const GenMatrix& A, const int rn, const int rm,
                        double* Rvalue, int* Rcolumn, GenMatrix& RtA){

	const int an = A.getNumRow();
	const int am = A.getNumColumn();


	if(rn != RtA.getNumRow() || am != RtA.getNumColumn() ){
		cout << "error dimension mismatch" << endl;
		return;
	}

	RtA.quickAllZero();
	

	for(int rtRow = 0; rtRow < rn; rtRow++){
		for(int j = 0; j < rm; j++){
			double rValue  = Rvalue[rtRow*rm + j];
			int    rColumn = Rcolumn[j];
			int aRow = rColumn;
			double* pAv = A.value [aRow];
			int*    pAc = A.column[aRow];
			int numNZ   = A.numNZ [aRow];
			for(int nu = 0; nu < numNZ; nu++){
				int    aColumn = pAc[nu];
				double aValue  = pAv[nu];
				RtA.add(rtRow, aColumn, rValue*aValue);
			}
		}
	}
}


// *********************************************************
// matrix-matrix(compress form) multiplication
//  [A][R] => [AR]
/*     ver. 1.0 2010.4.30  K. Watanabe                  */
// *********************************************************
// input:
//   A....matrix[A] an X am
//   Rvalue....non-zero entry of [R]  rm X rn
//   RColumn ... column index
//  output:
//   AR .......=[A][R]  an X rn
//
//(e.g.  Rt =|1 0 3 0 5|  => Rvalue = {1,3,5,   Rcolumn = {0,2,4}
//           |2 0 4 0 6|               2,4,6}   rn = 2, rm = 3 (not 5!!)
//
void GenMatrix::multipleAR(const GenMatrix& A, const int rn, const int rm,
                        double* Rvalue, int* Rcolumn, GenMatrix& AR){

	try{
	const int an = A.getNumRow();
	const int am = A.getNumColumn();


	if(an != AR.getNumRow() || rn != AR.getNumColumn() ){
		cout << "error dimension mismatch" << endl;
		return;
	}

	AR.quickAllZero();
	/*ARを密行列と仮定する*/
	AR.delRegion();
	AR.newRegion(an, an);
	for(int i = 0; i < an; i++){
		AR.numNZ[i] = an;
		AR.szColumn[i] = an;
		for(int j = 0; j < an; j++){
			AR.column[i][j] = j;
			AR.value [i][j] = 0.0;
		}
	}
	
	
	double* rValue = new double[am];

	for(int rtRow = 0; rtRow < rn; rtRow++){
		/* make un-compress i-th row of [R]t */
		for(int jj = 0; jj < am; jj++){
			rValue[jj] = 0.0;
		}
		for(int j = 0; j < rm; j++){
			rValue[Rcolumn[j]] = Rvalue[rtRow*rm + j];
		}
		
		for(int aRow = 0; aRow < an; aRow++){
			double* pAv = A.value [aRow];
			int*    pAc = A.column[aRow];
			int numNZ   = A.numNZ [aRow];
			for(int nu = 0; nu < numNZ; nu++){
				int    aColumn = pAc[nu];
				double aValue  = pAv[nu];
//				AR.add(aRow, rtRow, rValue[aColumn]*aValue);
				AR.value[aRow][rtRow] += rValue[aColumn]*aValue;
			}
		}
	}
	delete[] rValue;
	}
	catch(...){
		cout << "error at multipleAR" << endl;
		throw;
	}
}


/* ***以下ファイル入出力Method *** */
/* ******************************* */

/* ファイルに出力 */
/* ********************************************************** */
int GenMatrix::writeFile_V1(const char* filename, const GenMatrix& A){
/* ************************************************************ */
//    ファイルへデータの書き出し(Ver.1形式用)
//    テキスト形式ですが、独自フォーマットです
//    容量を食うので使用しないで下さい（ver2.0推奨）
//    szColumnは出力しない
//
/*     ver. 1.20 2002 9.25  K. Watanabe                           */
/* ************************************************************* */
	try{
	int i, mu;
	
	/* ファイルオープン*/
	fstream fp(filename, ios::out);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* 書式設定 */
	fp << setiosflags(ios::scientific) /* 指数表示 */
	   << setprecision(15);            /* 有効数字16桁(小数点以下15桁) */

	
	/* ヘッダ情報出力 */
	fp << "##GenMatrix writeFile Method Ver 1.0" <<endl;
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}
	
	/* サイズ情報出力 */
	fp << "numRow= "    << A.numRow
	   << " numColumn= "<< A.numColumn
	   << " numMaxBW= " << A.searchNumMaxBW()
	   << " sortFlag= " << A.sortFlag
	   << endl;
	
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}

	/* データ出力(diag, numNZ) */
	for (i = 0; i < A.numRow; i++){
		fp << "DATA1"
		   << " i= "     << i
		   << " numNZ= " << A.numNZ[i]
		   << endl;
		if( fp.fail() == true ){
			cout << "can not write" << endl;
			throw IOErr();
		}
	}

	/* データ出力(index, value) */
	for (i = 0; i < A.numRow; i++){
		for (mu = 0; mu < A.numNZ[i]; mu++){
			fp << "DATA2"
			   << " i= "  << i
			   << " mu= " << mu
			   << " j= "  << A.column[i][mu]
			   << " val= "<< A.value[i][mu]
			   << endl;
			if( fp.fail() == true ){
				cout << "can not write" << endl;
				throw IOErr();
			}
		}
	}

	fp.close();
	return 0;
	}
	catch(...){
		throw IOErr();
	}
}

/* ********************************************************** */
int GenMatrix::writeFile(const char* filename, const GenMatrix& A){
/* ************************************************************ */
//    ファイルへデータの書き出し(ver2.0)
//    テキスト形式ですが、独自フォーマットです
//    szColumnは出力しない
//
/*     ver. 2.00 2002 10.3  K. Watanabe                           */
/* ************************************************************* */
	try{
	int i, mu;
	
	/* ファイルオープン*/
	fstream fp(filename, ios::out);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* 書式設定 */
	fp << setiosflags(ios::scientific) /* 指数表示 */
	   << setprecision(15);            /* 有効数字16桁(小数点以下15桁) */

	
	/* ヘッダ情報出力 */
	fp << "##GenMatrix writeFile Method Ver 2.0" <<endl;
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}
	
	/* サイズ情報出力 */
	fp << "numRow= "    << A.numRow
	   << " numColumn= "<< A.numColumn
	   << " numMaxBW= " << A.searchNumMaxBW()
	   << " sortFlag= " << A.sortFlag
	   << endl;
	
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}

	/* データ出力(diag, numNZ) */
	for (i = 0; i < A.numRow; i++){
		fp << i
		   << " " << A.numNZ[i]
		   << endl;
		if( fp.fail() == true ){
			cout << "can not write" << endl;
			throw IOErr();
		}
	}

	/* データ出力(index, value) */
	for (i = 0; i < A.numRow; i++){
		fp << i << endl;
		for (mu = 0; mu < A.numNZ[i]; mu++){
			fp << mu
			   << " " << A.column[i][mu]
			   << " " << A.value[i][mu]
			   << endl;
			if( fp.fail() == true ){
				cout << "can not write" << endl;
				throw IOErr();
			}
		}
	}

	fp.close();
	return 0;
	}
	catch(...){
		throw IOErr();
	}
}

/* ファイルから取りこみ */
/* ********************************************************** */
int GenMatrix::readFile_V1(const char* filename, GenMatrix& A){
/* ************************************************************ */
//    ファイルからデータの取りだし(Ver.1形式用)
//    writeFIleメソッドにより生成されたファイルからデータを
//    読みこむ
//   旧形式（ver1.0）で書かれたデータを取りこむためのもの
//
/*     ver. 1.21 2004.10.26   K. Watanabe                          */
/* ************************************************************* */
	int i, mu;
	int iDummy, muDummy;

	string len, strDummy1, strDummy2, strDummy3, strDummy4, strDummy5;
	const string VERSION="##GenMatrix writeFile Method Ver 1.0";
	
	try{
	
	/* ファイルオープンチェック*/
	fstream fp(filename, ios::in);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* ヘッダ情報取りこみと確認 */
	getline(fp, len);
	
	if( fp.fail() == true || len != VERSION ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	
	/* delete region A */
	A.delRegion();

	/* サイズ情報取得 */
	fp >> strDummy1 >> A.numRow
	   >> strDummy2 >> A.numColumn
	   >> strDummy3 >> iDummy
	   >> strDummy4 >> A.sortFlag;
	
	if( fp.fail() == true ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	

	A.newRegion(A.numRow,1);


	/* データ取りこみ(diag, numNZ) */
	for (i = 0; i < A.numRow; i++){
		fp >> strDummy1
		   >> strDummy2 >> iDummy
		   >> strDummy3 >> A.numNZ[i];
		
		if( fp.fail() == true || i != iDummy){
			cout << "format incompatible" << endl;
			throw IOErr();
		}
	}
	
	/*各行の領域確保*/
	int numNZ;
	for (i = 0; i < A.numRow; i++){
		numNZ = A.numNZ[i];
		delete[] A.value[i];
		delete[] A.column[i];
		A.value[i]  = new double[numNZ];
		A.column[i] = new    int[numNZ];
		A.szColumn[i] = numNZ;
	}

	/* データ取りこみ(index, value) */
	for (i = 0; i < A.numRow; i++){
		for (mu = 0; mu < A.numNZ[i]; mu++){
			fp >> strDummy1
			   >> strDummy2 >> iDummy
			   >> strDummy3 >> muDummy
			   >> strDummy4 >> A.column[i][mu]
			   >> strDummy5 >> A.value[i][mu];
			   
			if( fp.fail() == true || i != iDummy || mu != muDummy){
				cout << "format incompatible" << endl;
				throw IOErr();
			}
		}
	}
	fp.close();

	}
	catch(MemErr x){
		throw x;
	}
	catch(OverIdxErr x){
		throw x;
	}

	return 0;
}


/* ********************************************************** */
int GenMatrix::readFile(const char* filename, GenMatrix& A){
/* ************************************************************ */
//    ファイルからデータの取りだし(Ver2.0)
//    writeFIleメソッドにより生成されたファイルからデータを
//    読みこむ
//   旧形式（ver1.0）で書かれたデータの場合はreadFile_V1を呼び出す
//
/*     Ver. 2.11 2004.10.26   K. Watanabe                          */
/* ************************************************************* */
	int i, mu;
	int iDummy, muDummy;

	string len, strDummy1, strDummy2, strDummy3, strDummy4;
	const string VERSION_V1="##GenMatrix writeFile Method Ver 1.0";
	const string VERSION_V2="##GenMatrix writeFile Method Ver 2.0";
	
	try{
	
	/* ファイルオープンチェック*/
	fstream fp(filename, ios::in);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* ヘッダ情報取りこみと確認 */
	getline(fp, len);

	if( fp.fail() == true || (len != VERSION_V1 && len != VERSION_V2) ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	if(len == VERSION_V1){
		fp.close();
		GenMatrix::readFile_V1(filename, A);
		return 0;
	}
	
	
	/* delete region A */
	A.delRegion();

	/* サイズ情報取得 */
	fp >> strDummy1 >> A.numRow
	   >> strDummy2 >> A.numColumn
	   >> strDummy3 >> iDummy
	   >> strDummy4 >> A.sortFlag;
	
	if( fp.fail() == true ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	

	A.newRegion(A.numRow,1);


	/* データ取りこみ( numNZ) */
	for (i = 0; i < A.numRow; i++){
		fp >> iDummy
		   >> A.numNZ[i];
		
		if( fp.fail() == true || i != iDummy){
			cout << "format incompatible" << endl;
			throw IOErr();
		}
	}
	
	/*各行の領域確保*/
	int numNZ;
	for (i = 0; i < A.numRow; i++){
		numNZ = A.numNZ[i];
		delete[] A.value[i];
		delete[] A.column[i];
		A.value [i] = new double[numNZ];
		A.column[i] = new    int[numNZ];
		A.szColumn[i] = numNZ;
	}

	/* データ取りこみ(index, value) */
	for (i = 0; i < A.numRow; i++){
		fp >> iDummy;
		if( fp.fail() == true || i != iDummy){
			cout << "format incompatible" << endl;
			throw IOErr();
		}

		numNZ = A.numNZ[i];
		for (mu = 0; mu < numNZ; mu++){
			fp >> muDummy
			   >> A.column[i][mu]
			   >> A.value[i][mu];
			   
			if( fp.fail() == true || mu != muDummy){
				cout << "format incompatible" << endl;
				throw IOErr();
			}
		}
	}
	fp.close();

	}
	catch(MemErr x){
		throw x;
	}
	catch(OverIdxErr x){
		throw x;
	}

	return 0;
}



/* ベクトルデータ(1次元配列)をファイルに出力 */
/* ********************************************************** */
int GenMatrix::vectWriteFile(const char* filename, const double* x, int n){
/* ********************************************************** */
//    ファイルへデータの書き出し
//    テキスト形式ですが、独自フォーマットです
//
/*     ver. 2.01 2002 9.25  K. Watanabe                           */
/* ************************************************************* */
	try{
	int i;
	
	/* ファイルオープン*/
	fstream fp(filename, ios::out);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* 書式設定 */
	fp << setiosflags(ios::scientific) /* 指数表示 */
	   << setprecision(15);            /* 有効数字16桁(小数点以下15桁) */

	
	/* ヘッダ情報出力 */
	fp << "##vectWriteFile Method Ver 1.0" <<endl;
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}
	
	/* サイズ情報出力 */
	fp << "numDim= " << n
	   << endl;
	
	if( fp.fail() == true ){
		cout << "can not write" << endl;
		throw IOErr();
	}

	/* データ出力 */
	for (i = 0; i < n; i++){
		fp << " i= "   << setw(6) << i
		   << " val= " << x[i]
		   << endl;
		if( fp.fail() == true ){
			cout << "can not write" << endl;
			throw IOErr();
		}
	}

	fp.close();
	}
	catch(...){
		throw IOErr();
	}

	return 0;
}

/* ベクトルデータ(1次元配列)ファイルから取りこみ */
/* ********************************************************** */
int GenMatrix::vectReadFile(const char* filename, double** x, int* n){
/* ************************************************************ */
//    ファイルからデータの取りだし
//    writeFIleメソッドにより生成されたファイルからデータを
//    読みこむ
//
/*     ver. 2.01 2002 9.25  K. Watanabe                           */
/* ************************************************************* */
	try{
	int i, iDummy;
	string len, strDummy1, strDummy2;
	const string VERSION="##vectWriteFile Method Ver 1.0";
	
	
	/* ファイルオープンチェック*/
	fstream fp(filename, ios::in);
	if( fp.fail() == true ){
		cout << "can not open file" << endl;
		throw IOErr();
	}
	
	/* ヘッダ情報取りこみと確認 */
	getline(fp, len);
	
	if( fp.fail() == true || len != VERSION ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	
	/* サイズ情報取得 */
	fp >> strDummy1 >> *n;
	
	if( fp.fail() == true ){
		cout << "format incompatible" << endl;
		throw IOErr();
	}
	
	/* malloc region */
	*x = new double[*n];

	/* データ取りこみ(diag, numNZ) */
	for (i = 0; i < *n; i++){
		fp >> strDummy1 >> iDummy
		   >> strDummy2 >> (*x)[i];
		
		if( fp.fail() == true || i != iDummy){
			cout << "format incompatible" << endl;
			throw IOErr();
		}
	}

	fp.close();

	}
	catch(bad_alloc ){
		throw MemErr();
	}
	catch(...){
		throw IOErr();
	}

	return 0;
}


/* *********************************************************** */
void GenMatrix::modGramSchmidt(){
//
//  modified Gram-Schmidt method
//
/*     ver. 1.00 2010.5.23   K. Watanabe                      */
/* *********************************************************** */
/*
// 行ベクトルが正規直交するようにする（注意：計算速度の観点から
// 列ベクトルではなく，行ベクトルを直交化する）
// 0行目は単に正規化するだけ
*/


	try{

		const int n = numRow;
		const int m = numColumn;
		double* vi = new double[m];
		double* vj = new double[m];
		for(int i = 0; i < n; i++){
			/*expand vi*/
			for(int k = 0; k < m; k++){
				vi[k] = 0.0;
			}
			int*   columnPi = column[i];
			double* valuePi = value[i];
			int      numNZi = numNZ[i];
			for(int nu = 0; nu < numNZi; nu++){
				vi[columnPi[nu]] = valuePi[nu];
			}
			
			
			for(int j = 0; j < i; j++){
				/*expand vj*/
				for(int k = 0; k < m; k++){
					vj[k] = 0.0;
				}
				int*   columnPj = column[j];
				double* valuePj = value[j];
				int      numNZj = numNZ[j];
				for(int nu = 0; nu < numNZj; nu++){
					vj[columnPj[nu]] = valuePj[nu];
				}
				/*cal (vi, vj)*/
				double innerIJ = 0.0;
				for(int k = 0; k < m; k++){
					innerIJ += vi[k]*vj[k];
				}
				
				/* do modG-S */
				for(int k = 0; k < m; k++){
					vi[k] -= innerIJ*vj[k];
				}
			}
			
			/*normarize */
			double normVi = 0.0;
			for(int k = 0; k < m; k++){
				normVi += vi[k]*vi[k];
			}
			normVi = sqrt(normVi);
			if(normVi == 0.0){
				cout << "error norm{vi} == 0" << endl;
				return;
			}
			for(int k = 0; k < m; k++){
				vi[k] /= normVi;
			}
			
			/*count non-zero entry */
			int numNZnew = 0;
			for(int k = 0; k < m; k++){
				if(vi[k] != 0.0){
					numNZnew++;
				}
			}
			if(numNZnew > numNZi){
				delete[] column[i];
				delete[] value[i];
				column[i] = new    int[numNZnew];
				 value[i] = new double[numNZnew];
				 numNZ[i]    = numNZnew;
				 szColumn[i] = numNZnew;
			}
			int nu = 0;
			for(int k = 0; k < m; k++){
				if(vi[k] != 0.0){
					column[i][nu] = k;
					value[i][nu] = vi[k];
					nu++;
				}
			}
		}
		
		delete[] vi;
		delete[] vj;
	}
	catch(...){
		cout << "error in modGramSchmidt()" << endl;
		throw;
	}
}

/* *********************************************************** */
void GenMatrix::modGramSchmidt(double** A, const int n, const int m){
//
//  modified Gram-Schmidt method for two-dimensional array
//
/*     ver. 1.00 2010.5.25   K. Watanabe                      */
/* *********************************************************** */
/*
// 行ベクトルが正規直交するようにする（注意：計算速度の観点から
// 列ベクトルではなく，行ベクトルを直交化する）
// 0行目は単に正規化するだけ
*/


	try{

		double* vi;
		double* vj;
		for(int i = 0; i < n; i++){
			vi = A[i];
			
			for(int j = 0; j < i; j++){
				vj = A[j];
				
				/*cal (vi, vj)*/
				double innerIJ = 0.0;
				for(int k = 0; k < m; k++){
					innerIJ += vi[k]*vj[k];
				}
				
				/* do modG-S */
				for(int k = 0; k < m; k++){
					vi[k] -= innerIJ*vj[k];
				}
			}
			
			/*normarize */
			double normVi = 0.0;
			for(int k = 0; k < m; k++){
				normVi += vi[k]*vi[k];
			}
			normVi = sqrt(normVi);
			if(normVi == 0.0){
				cout << "error norm{vi} == 0" << endl;
				return;
			}
			for(int k = 0; k < m; k++){
				vi[k] /= normVi;
			}
			
		}
	}
	catch(...){
		cout << "error in modGramSchmidt()" << endl;
		throw;
	}
}


/* ******************************************************************* */
int GenMatrix::gsSolv_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax){
//   ガウスザイデル法で連立方程式を解く for CEFC2010
//     Ver. 1.00 201.5.27  K. Watanabe
/* ******************************************************************* */

/*
// 残差履歴と、低周波成分の残差履歴の調査
// ev は最小非零固有値に対応する固有ベクトル
*/
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double* r   = new double[n];

	/*cal r0 = b - Ax0 */
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
	}

	double r2sum     = 0.0;
	double evr0sum   = 0.0;
	double evrmaxsum = 0.0;
	for(int i = 0; i < n; i++){
		r2sum     += r[i]*r[i];
		evr0sum   += r[i]*ev0[i];
		evrmaxsum += r[i]*evmax[i];
	}

	double normR0     = sqrt(r2sum);
	double absEVR0    = fabs(evr0sum);
	double absEVRMAX0 = fabs(evrmaxsum);

	cout << "first norm[r]       = " << normR0 << endl;
	cout << "first |(r, ev0)  |  = " << absEVR0 << endl;
	cout << "first |(r, evmax)|  = " << absEVRMAX0 << endl;
	r_ev0[0]   += 1.0;
	r_evmax[0] += 1.0;


	for(int k = 1; k < maxItr; k++){
	
		for(int i = 0; i < n; i++){
			r[i] = b[i];
		}
		double r2sum = 0.0;
		double evr0sum = 0.0;
		double evrmaxsum = 0.0;
		for(int i = 0; i < n; i++){
			int numNZ = A.numNZ[i];
			double* ap = A.value[i];
			int*    cp = A.column[i];
			for(int mu = 0; mu < numNZ; mu++){
				int j = cp[mu];
				r[i] -= ap[mu] * x[j];
			}
			x[i] += r[i]*invDiag[i]; /*x[i] = (r[i] + a.diag[i]*x[i])/a.diag[i] */
			r2sum   += r[i]*r[i];
			evr0sum   += r[i]*ev0[i];
			evrmaxsum += r[i]*evmax[i];
		}
		cout << k << ", Normalized|r|= , " << sqrt(r2sum)/normR0 
		     << " , Normalized|(r, ev0)  |= , "<< fabs(evr0sum)/absEVR0
		     << " , Normalized|(r, evmax)|= , "<< fabs(evrmaxsum)/absEVRMAX0 <<  endl;
		r_ev0  [k] += fabs(evr0sum)/absEVR0;
		r_evmax[k] += fabs(evrmaxsum)/absEVRMAX0;
	}

	/***free memory***/
	delete[] r;
	delete[] invDiag;
	
	/***  termination ***/
	return 0;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}


/* ******************************************************************* */
int GenMatrix::adaptiveIJacobi_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR-based adaptive I-Jacobi Method で連立方程式を解く for multigrid smoother
//     Ver. 1.00 2010.1.31  K. Watanabe
/* ******************************************************************* */

/*
// 残差履歴と、低周波成分の残差履歴の調査
// ev は最小非零固有値に対応する固有ベクトル
*/

	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}

	
	/*係数行列から対角項を抽出*/
	double* invDiag = new double[n];
	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
	for(int i = 0; i < n ; i++){
		if(invDiag[i] == 0.0){
			cout << "error! diagnal [" << i << " ] is zero" << endl;
		}
		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
	}
	
	double* r   = new double[n];

	double* p   = new double[n];
	double* s   = new double[n];
	double* t   = new double[n];
	double* as  = new double[n];
	double* at  = new double[n];
	double* dx  = new double[n];
	double* dr  = new double[n];

	/*cal r0 = b - Ax0 */
	productAX(A, x, r);
	for(int i = 0; i < n; i++){
		r[i] = b[i] - r[i];
	}

	double r2sum     = 0.0;
	double evr0sum   = 0.0;
	double evrmaxsum = 0.0;
	for(int i = 0; i < n; i++){
		r2sum     += r[i]*r[i];
		evr0sum   += r[i]*ev0[i];
		evrmaxsum += r[i]*evmax[i];
	}

	double normR0     = sqrt(r2sum);
	double absEVR0    = fabs(evr0sum);
	double absEVRMAX0 = fabs(evrmaxsum);

	cout << "first norm[r]       = " << normR0 << endl;
	cout << "first |(r, ev0)  |  = " << absEVR0 << endl;
	cout << "first |(r, evmax)|  = " << absEVRMAX0 << endl;
	r_ev0[0]   += 1.0;
	r_evmax[0] += 1.0;

	
	/* set p vector (random vector) */
	for(int i = 0; i < n; i++){
		p[i] = r[i];
//		p[i] = 1.0;
	}

	double beta = 0.0;

	for(int k = 1; k < maxItr; k++){
	
		/* cal {s} */
		for(int i = 0; i < n; i++){
			t[i] = r[i] - beta*dr[i];
		}

		GenMatrix::productAX(A, t, at);
		double omega1 = 0.0;
		double omega2 = 0.0;
		for(int i = 0; i < n; i++){
			double invDat = invDiag[i] * at[i];
			omega1 += invDat * t[i];
		}
		for(int i = 0; i < n; i++){
			double invDat = invDiag[i] * at[i];
			omega2 += invDat * invDat;
		}
		double omega = omega1 / omega2;
		if(omega > 2.0){omega = 2.0;}
		if(omega < 0.1){omega = 0.1;}
		
		for(int i = 0; i < n; i++){
			s[i] = omega * invDiag[i] * t[i];
		}
		
		/* cal dx and dr */
		GenMatrix::productAX(A, s, as);

		for(int i = 0; i < n; i++){
			dx[i] = s[i] - beta * dx[i];
			dr[i] = -as[i] + s[i]/invDiag[i]/omega - r[i];
		}
		
		/* cal update r and x */
		for(int i = 0; i < n; i++){
			r[i] += dr[i];
			x[i] += dx[i];
		}
		
		double r2sum    = 0.0;
		double evr0sum   = 0.0;
		double evrmaxsum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum     += r[i]*r[i];
			evr0sum   += r[i]*ev0[i];
			evrmaxsum += r[i]*evmax[i];
		}
//		normR = sqrt(r2sum);
		cout << k << ", Normalized|r|= , " << sqrt(r2sum)/normR0 
		     << " , Normalized|(r, ev0)|= , "<< fabs(evr0sum)/absEVR0
		     << " , Normalized|(r, evmax)|= , "<< fabs(evrmaxsum)/absEVRMAX0 <<  endl;
		r_ev0  [k] += fabs(evr0sum)/absEVR0;
		r_evmax[k] += fabs(evrmaxsum)/absEVRMAX0;
	
		/* cal gamma type 1 */
		double pr  = 0.0;
		double pdr = 0.0;
		for(int i = 0; i < n; i++){
			pr  += p[i] *  r[i];
		}

		for(int i = 0; i < n; i++){
			pdr += p[i] * dr[i];
		}
		beta = pr/pdr;

		/* cal gamma type 2 */
//		if(normR < minNormR){
//			minNormR = normR;
//		}
//		gamma = minNormR;

	}


	
	/***free memory***/
	delete[] r;

	delete[] invDiag;
	delete[] s;
	delete[] t;
	delete[] as;
	delete[] at;
	delete[] dx;
	delete[] dr;

	
	/***  termination ***/
	return 0;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}



/* ******************************************************************* */
int GenMatrix::IDR_sSolv_tmp(const GenMatrix& A, const double* b, double* x,
                              int maxItr, int s, double* ev0, double* evmax,
                              double* r_ev0, double* r_evmax){
//              Solve system of linear equations
//                  [A]*X = B
//   IDR Method で連立方程式を解く
//     Ver. 1.00 2010.5.24  K. Watanabe
/* ******************************************************************* */

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	int convFlag = 0;
	
	try{
	const int n = A.numRow; /*係数行列の行数を取得*/

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	if(s < 1){
		cout << "parameter s should be s >= 1 " << endl;
		return -1;
	}

	
	/*係数行列から対角項を抽出*/
//	double* invDiag = new double[n];
//	getDiagVector(A, invDiag);

	/*処理速度向上のため、対角項の逆数を保存*/
//	for(int i = 0; i < n ; i++){
//		if(invDiag[i] == 0.0){
//			cout << "error! diagnal [" << i << " ] is zero" << endl;
//		}
//		invDiag[i] = 1.0/invDiag[i];
		
//		cout << invDiag[i] << endl;
//	}
	

	double* r   = new double[n];
	double* v   = new double[n];
	double* t   = new double[n];
	double* q   = new double[n];
	double* e   = new double[n];

	double** E = new double*[s];
	double** P = new double*[s];
	double** Q = new double*[s];
	for(int i = 0; i < s; i++){
		E[i] = new double[n];
		P[i] = new double[n];
		Q[i] = new double[n];
	}


	int numItr = -1;

	/*cal r0 = b - Ax0 */
	productAX(A, x, r);
	for(int j = 0; j < n; j++){
		r[j] = b[j] - r[j];
	}

	double r2sum     = 0.0;
	double evr0sum   = 0.0;
	double evrmaxsum = 0.0;
	for(int i = 0; i < n; i++){
		r2sum     += r[i]*r[i];
		evr0sum   += r[i]*ev0[i];
		evrmaxsum += r[i]*evmax[i];
	}

	double normR0     = sqrt(r2sum);
	double absEVR0    = fabs(evr0sum);
	double absEVRMAX0 = fabs(evrmaxsum);

	cout << "first norm[r]       = " << normR0 << endl;
	cout << "first |(r, ev0)  |  = " << absEVR0 << endl;
	cout << "first |(r, evmax)|  = " << absEVRMAX0 << endl;
	r_ev0[0]   += 1.0;
	r_evmax[0] += 1.0;


	/* make P matrix */
//	GenMatrix PG(s, n);

	for(int j = 0; j < n; j++){
		P[0][j] = r[j];
//		P[0][j] = (double)rand();
//		PG.set(0, j, r[j]);
//		PG.set(0, j, (double)rand());
	}
	for(int i = 1; i < s; i++){
		for(int j = 0; j < n; j++){
//			P[i][j] = 0.0;
			P[i][j] = (double)rand();
		}
//		P[i][i] = 1.0;
//		PG.set(i, i, 1.0);
	}
//	PG.modGramSchmidt();
	GenMatrix::modGramSchmidt(P, s, n);

	/*initial loop */
	double omega = 0.0;
	for(int k = 0; k < s; k++){
		GenMatrix::productAX(A, r, v);
		
		double vr = 0.0;
		double vv = 0.0;
		for(int j = 0; j < n; j++){
			vr += v[j]*r[j];
			vv += v[j]*v[j];
		}
		omega = vr/vv;
		
		for(int j = 0; j < n; j++){
			double qj =  omega*r[j];
			double ej = -omega*v[j];
			r[j] += ej;
			x[j] += qj;
			E[s-1-k][j] = ej;
			Q[s-1-k][j] = qj;
		}
		r2sum = 0.0;
		for(int j = 0; j < n; j++){
			r2sum += r[j]*r[j];
		}


	}
	
	


	double* c   = new double[s];
	double* ptr = new double[s];
//	double* rs  = new double[s];

	int k;
	for(k = s; k < maxItr; k++){
		
		/*make PtE */
		GenMatrix PtE(s, s);
		for(int i = 0; i < s; i++){
			for(int j = 0; j < s; j++){
				double pte = 0.0;
				for(int k = 0; k < n; k++){
					pte += P[i][k] * E[j][k];
				}
				PtE.set(i, j, pte);
			}
		}
//	GenMatrix::writeFile("PtE.d", PtE);
//	exit(1);
		
		GenMatrix PtEinv;
		PtEinv = PtE;
		GenMatrix::inv(PtEinv);
//	GenMatrix::writeFile("PtEinv.d", PtEinv);

		
		/*solve PtE{c} = Pt{r} */
		double ptr_norm = 0.0;
		for(int i = 0; i < s; i++){
			ptr[i] = 0.0;
			for(int j = 0; j < n; j++){
				ptr[i] += P[i][j]*r[j];
			}
			ptr_norm += ptr[i]*ptr[i];
			c[i] = 0.0;
		}
		GenMatrix::productAX(PtEinv, ptr, c);
			
		/* v = r - E{c} */
		for(int j = 0; j < n; j++){
			v[j] = r[j];
		}
		for(int i = 0; i < s; i++){
			for(int j = 0; j < n; j++){
				v[j] -= E[i][j]*c[i];
			}
		}
		
		if(k%(s+1) == s){
			GenMatrix::productAX(A, v, t);
			
			double tv = 0.0;
			double tt = 0.0;
			for(int j = 0; j < n; j++){
				tv += t[j]*v[j];
				tt += t[j]*t[j];
			}
			omega = tv/tt;
			
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
				e[j] = -omega*t[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
					e[j] -= E[i][j]*c[i];
				}
			}
		}
		else{
			for(int j = 0; j < n; j++){
				q[j] =  omega*v[j];
			}
			for(int i = 0; i < s; i++){
				for(int j = 0; j < n; j++){
					q[j] -= Q[i][j]*c[i];
				}
			}
			GenMatrix::productAX(A, q, e);
			for(int j = 0; j < n; j++){
				e[j] = -e[j];
			}
		}
		
		/*re-construct E={e(k), e(k-1), ...,e(k+1-s)} and Q*/

		double* e0p = E[s-1];
		double* q0p = Q[s-1];
		for(int i = s-1; i > 0; i--){
			E[i] = E[i-1];
			Q[i] = Q[i-1];
		}
		E[0] = e0p;
		Q[0] = q0p;
		
		for(int j = 0; j < n; j++){
			r[j] += e[j];
			x[j] += q[j];
			E[0][j] = e[j];
			Q[0][j] = q[j];
		}

		double normR = sqrt(r2sum);
		double r2sum    = 0.0;
		double evr0sum   = 0.0;
		double evrmaxsum = 0.0;
		for(int i = 0; i < n; i++){
			r2sum     += r[i]*r[i];
			evr0sum   += r[i]*ev0[i];
			evrmaxsum += r[i]*evmax[i];
		}
		cout << k << ", Normalized|r|= , " << sqrt(r2sum)/normR0 
		     << " , Normalized|(r, ev0)|= , "<< fabs(evr0sum)/absEVR0
		     << " , Normalized|(r, evmax)|= , "<< fabs(evrmaxsum)/absEVRMAX0 <<  endl;
		r_ev0  [k] += fabs(evr0sum)/absEVR0;
		r_evmax[k] += fabs(evrmaxsum)/absEVRMAX0;

		
	}


	/***free memory***/
	delete[] r;
	delete[] v;
	delete[] t;
	delete[] q;
	delete[] e;
	for(int i = 0; i < s; i++){
		delete[] E[i];
		delete[] P[i];
		delete[] Q[i];
	}
	
	delete[] E;
	delete[] P;
	delete[] Q;

	delete[] c;
	delete[] ptr;

	
	/***  termination ***/
	return numItr;

	}
	catch(bad_alloc){
		throw MemErr();
	}
	catch(OverIdxErr x){
		throw x;
	}
	catch(...){
		cout << "unknown Error" << endl;
		return -1;
	}
}



/* *********************************************************************
//    update solutions for GMRES 
//     Ver. 1.00 2010.9.13  K. Watanabe
********************************************************************** */
void Update(double* x, const int k, double** h, double* s, double** v,
       const int n, const int m){

	double* y = new double[m+1];
	for(int i = 0; i < m+1; i++){
		y[i] = s[i];
	}

	// Backsolve:  
	for(int i = k; i >= 0; i--) {
		y[i] /= h[i][i];
		for(int j = i - 1; j >= 0; j--){
			y[j] -= h[j][i] * y[i];
		}
	}

	for(int j = 0; j <= k; j++){
		for(int i = 0; i < n; i++){
			x[i] += v[j][i] * y[j];
		}
	}
	delete[] y;
}

/* *********************************************************************
//    GeneratePlaneRotation for GMRES 
//     Ver. 1.00 2010.9.13  K. Watanabe
********************************************************************** */
void GeneratePlaneRotation(double* dx, double* dy, double* cs, double *sn){
	if(*dy == 0.0) {
		*cs = 1.0;
		*sn = 0.0;
	} else if(fabs(*dy) > fabs(*dx)) {
		double temp = *dx / *dy;
		*sn = 1.0 / sqrt( 1.0 + temp*temp );
		*cs = temp * (*sn);
	} else{
		double temp = *dy / *dx;
		*cs = 1.0 / sqrt( 1.0 + temp*temp );
		*sn = temp * (*cs);
	}
}

/* *********************************************************************
//    ApplyPlaneRotation for GMRES 
//     Ver. 1.00 2010.9.13  K. Watanabe
********************************************************************** */
void ApplyPlaneRotation(double* dx, double* dy, double* cs, double *sn){
	double temp  =  (*cs) * (*dx) + (*sn) * (*dy);
	*dy = -(*sn) * (*dx) + (*cs) * (*dy);
	*dx = temp;
}



/* *********************************************************************
//    GMRES with restart Solver with GS type preconditioner
//     Solution of system of linear equation
//             A*X = B
//     by Preconditioned BiConjugate Gradient Stabilized method
//
//     Ver. 1.00 2010.8.26  K. Watanabe
********************************************************************** */
int GenMatrix::GMRES_mSolvWithDiagP(GenMatrix& A, const double* b, double* x,
                              double eps, int maxItr, const int m){

/*
//    ==== input ====
//     A.........Matrix
//     b[i]......right hand side vector(i=0,1,2,...n-1)
//     maxItr... Max(limit) number of iterations
//     eps.......tolerance(下記注意書き参照)
//     m ........restart
//
//    ==== in/output ====
//     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
//             output : solution
//
//    ==== return value ====
//     number of iterations
//
// ・係数行列が正方行列では無い場合は例外OverIdxErrを投げる
// ・戻り値は収束に要した反復回数（maxItrを超えた場合は-1を返す）
// ・収束判定は一般的な norm(r)/norm(b) < eps ではなくnorm(r) < eps
//   である（rは残差ベクトル）。
*/
	try{
	double r2sum, Mr2sum;

	bool convFlag  = false;
//	double alpha = 0.0;
//	double beta;
//	double omega = 1.0;/* non zero value */
//	double rou0  = 1.0;/* non zero value */
//	double rou, ts, tt, r0v;
	
	const int n = A.numRow; /*係数行列の行数を取得*/

	const double eps2 = eps * eps;

	/*正方行列かどうかのチェック*/
	if(n != A.numColumn){
		throw OverIdxErr();
	}
	
	/*対角情報配列用意*/
	double* invDiag = new double[n];

	/*対角成分取り出し*/
	getDiagVector(A, invDiag);
	for(int i = 0; i < n; i++){
		invDiag[i] = 1.0/invDiag[i];
	}
	
	double*   r  = new double[n];
	double*  Mr  = new double[n];
	double*   w  = new double[n];

	double** v = new double*[m+1];
	for(int i = 0; i < m+1; i++){
		v[i] = new double[n];
	}

	double**  H  = new double*[m+1];//   H[m+1][m];
	for(int j = 0; j < m+1; j++){
		H[j] = new double[m];
	}
	double*   s  = new double [m+1];
	double*  cs  = new double [m+1];
	double*  sn  = new double [m+1];


	/*初期残差ベクトルの計算*/
	/*** {AP} = [A]{x}***/
	productAX(A, x, r);

	/*** residual vector {r0} = {b} - [A]{X} ***/
	Mr2sum = 0.0;
	 r2sum = 0.0;
	for(int i = 0; i < n ; i++){
		 r[i] =  b[i] - r[i];
		Mr[i] = invDiag[i] * r[i];
		 r2sum +=  r[i] *  r[i];
		Mr2sum += Mr[i] * Mr[i];
	}
	cout << "norm|r|= " << sqrt(r2sum) << endl;
	double beta = sqrt(Mr2sum);

	/*** check convergence**/
	if(r2sum < eps2) {/*収束判定*/
#if ___DBGPRINTG
		cout << "---- converged ------" << endl;
		cout << "Number of iteration = " << 0 << endl;
#endif
		/***  termination ***/
		convFlag  = 1;
	}

	/***************★★★ iteration loop ★★★******************************/
	int itr;
	for(itr = 1; itr < maxItr; itr++){
//		cout << itr << "-th iteration" << endl;
		for(int i = 0; i < n; i++){
			v[0][i] = Mr[i] /beta;
		}

		for(int j = 0; j < m+1; j++){
			s[j] = 0.0;
		}
		s[0] = beta;
		for(int j = 0; j < m && itr <= maxItr; j++, itr++) {
			productAX(A, v[j], w);
			for(int i = 0; i < n; i++){
				w[i] = invDiag[i] * w[i];
			}
			
			for(int k = 0; k <= j; k++){
				H[k][j] = 0.0;
				for(int i = 0; i < n; i++){
					H[k][j] += w[i] * v[k][i];
				}
				for(int i = 0; i < n; i++){
					w[i] -= H[k][j] * v[k][i];
				}
			}
			double normW = 0.0;
			for(int i = 0; i < n; i++){
				normW += w[i] * w[i];
			}
			H[j+1][j] = sqrt(normW);
			for(int i = 0; i < n; i++){
				v[j+1][i] = w[i] / H[j+1][j];
			}

			for(int k = 0; k < j; k++){
				ApplyPlaneRotation(&H[k][j], &H[k+1][j], &cs[k], &sn[k]);
			}
      
			GeneratePlaneRotation(&H[j][j], &H[j+1][j], &cs[j], &sn[j]);
			   ApplyPlaneRotation(&H[j][j], &H[j+1][j], &cs[j], &sn[j]);
			   ApplyPlaneRotation(&s[j],    &s[j+1],    &cs[j], &sn[j]);
      
			if(fabs(s[j+1]) < eps){ // /*要再検討*/
				Update(x, j, H, s, v, n, m);
				convFlag = true;
				break;
			}
		}
		Update(x, m - 1, H, s, v, n, m);

		productAX(A, x, r);
		Mr2sum = 0.0;
		 r2sum = 0.0;
		for(int i = 0; i < n ; i++){
			 r[i] =  b[i] - r[i];
			Mr[i] = invDiag[i] * r[i];
			 r2sum +=  r[i] *  r[i];
			Mr2sum += Mr[i] * Mr[i];
		}
#if ___DBGPRINTG
		cout << "norm{r} = " << sqrt(r2sum) << endl;
#endif
		beta = sqrt(Mr2sum);
		if( r2sum < eps2){ /*収束判定*/
			convFlag = true;
			cout <<"---- converged ------" << endl;
			cout <<"Number of iteration = " << itr << endl;
			break;
		}
	}/**************************** end loop ***********************/
	
	/*** not convergence**/
	if (convFlag == false){
		itr = -1;
		cout << "--- not converged ---" << endl;
	}
	else{
		cout << "norm_GMRES{r} = " << sqrt(r2sum) << endl;
		/*真の残差を計算*/
		productAX(A, x, r);
		r2sum = 0.0;
		for(int i = 0; i < n; i++){
			r[i] = b[i] - r[i];
			r2sum += r[i]*r[i];
		}
		cout << "norm_true{r} = " << sqrt(r2sum) << endl;
	}

	delete[] invDiag;

	/***free memory***/
	delete[]   r;
	delete[]  Mr;
	delete[]   w;
	delete[]   s;
	delete[]  cs;
	delete[]  sn;

	for(int j = 0; j < m+1; j++){
		delete[] H[j];
		delete[] v[j];
	}
	delete[] H;
	delete[] v;


	return itr;
	}
	catch(...){
		throw IOErr();
	}
}
