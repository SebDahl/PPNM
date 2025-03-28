#include"matrix.h"
#include<string>
#include<algorithm>
#define SELF (*this)
#define FOR_V(i,v) for(int i=0;i<v.size();i++)
#define FOR_COLS(i,A) for(int i=0;i<A.sizecol();i++)
namespace pp{

vector& vector::operator+=(const vector& other) {
	FOR_V(i,SELF) data[i]+=other.data[i];
	return SELF; }

vector& vector::operator-=(const vector& other) {
	FOR_V(i,SELF) data[i]-=other.data[i];
	return SELF; }

vector& vector::operator*=(NUMBER x) {
	FOR_V(i,SELF) data[i]*=x;
	return SELF; }

vector& vector::operator/=(NUMBER x) {
	FOR_V(i,SELF) data[i]/=x;
	return SELF; }

void vector::print(std::string s,FILE* stream) const {
	fprintf(stream,"%s\n",s.c_str());
	for(int i=0;i<size();i++)fprintf(stream,"%9.4g ",(double)SELF[i]);
	fprintf(stream,"\n");
	}

vector operator/(const vector& v, NUMBER x){
	vector result; result.resize(v.size());
	FOR_V(i,result) result.data[i]=v.data[i]/x;
	return result; }

vector operator*(const vector& a, NUMBER x){
	vector result; result.resize(a.size());
	FOR_V(i,result) result.data[i]=a.data[i]*x;
	return result; }

vector operator*(NUMBER x,const vector& a){
	vector result; result.resize(a.size());
	FOR_V(i,result) result.data[i]=a.data[i]*x;
	return result; }

vector operator+(const vector& a, const vector& b){
	vector result; result.resize(a.size());
	FOR_V(i,result) result.data[i]=a.data[i]+b.data[i];
	return result; }

vector operator-(const vector& a, const vector& b){
	vector result; result.resize(a.size());
	FOR_V(i,result) result.data[i]=a.data[i]-b.data[i];
	return result; }

bool operator==(const vector& a, const vector& b) {
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++) {
		if (std::abs(a[i] - b[i]) > 1e-10) return false;
	}
	return true;
}

void matrix::resize(int n, int m){
	cols.resize(m);
	for(int i=0;i<m;++i)cols[i].resize(n);
	}

matrix matrix::transpose() const{
    matrix R; R.resize(sizecol(),sizerow());
    for(int j=0;j<R.sizecol();++j)
    for(int i=0;i<R.sizerow();++i)
        R[i,j]=SELF[j,i];
    return R;
    }

matrix matrix::identity(int n){
	matrix I; I.resize(n,n);
	for(int i=0;i<n;i++)I(i,i)=1;
	return I;
	}



matrix& matrix::operator+=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]+=other[i];
	return SELF; }

matrix& matrix::operator-=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]-=other[i];
	return SELF; }

matrix& matrix::operator*=(NUMBER x) {
	FOR_COLS(i,SELF) SELF[i]*=x;
	return SELF; }

matrix& matrix::operator/=(NUMBER x) {
	FOR_COLS(i,SELF) SELF[i]/=x;
	return SELF; }

matrix operator/(const matrix& A,NUMBER x){
	matrix R; R.resize(A.sizerow(),A.sizecol());
	FOR_COLS(i,R) R[i]=A[i]/x;
	return R; }

matrix operator*(const matrix& A,NUMBER x){
	matrix R; R.resize(A.sizerow(),A.sizecol());
	FOR_COLS(i,R) R[i]=A[i]*x;
	return R; }

matrix operator*(NUMBER x,const matrix& A){
	matrix R; R.resize(A.sizerow(),A.sizecol());
	FOR_COLS(i,R) R[i]=A[i]*x;
	return R; }

matrix operator+(const matrix& A, const matrix& B){
	matrix R; R.resize(A.sizerow(),A.sizecol());
	FOR_COLS(i,R) R[i]=A[i]+B[i];
	return R; }

matrix operator-(const matrix& A, const matrix& B){
	matrix R; R.resize(A.sizerow(),A.sizecol());
	FOR_COLS(i,R) R[i]=A[i]-B[i];
	return R; }

vector operator*(const matrix& M, const vector& v){
	vector r; r.resize(M.sizerow());
	for(int i=0;i<r.size();i++){
		NUMBER sum=0;
		for(int j=0;j<v.size();j++)sum+=M[i,j]*v[j];
		r[i]=sum;
		}
	return r;
	}

matrix operator*(const matrix& A, const matrix& B){
	matrix R; R.resize(A.sizerow(),B.sizecol());
	for(int k=0;k<A.sizecol();k++)
	for(int j=0;j<B.sizecol();j++)
		{
		for(int i=0;i<A.sizerow();i++)R[i,j]+=A[i,k]*B[k,j];
		}
	return R;
	}

bool operator==(const matrix& A, const matrix& B) {
	if (A.sizerow() != B.sizerow() || A.sizecol() != B.sizecol()) return false;
	for (int i = 0; i < A.sizerow(); i++) {
		for (int j = 0; j < A.sizecol(); j++) {
			if (std::abs(A(i, j) - B(i, j)) > 1e-10) return false;
		}
	}
	return true;
}
/*
vector matrix::get_col(int j){
	vector cj=SELF[j];
	return cj;
	}

void matrix::set_col(int j,vector& cj){
	SELF[i]=cj;
	}
*/

void matrix::print(std::string s,FILE* stream){
	fprintf(stream,"%s\n",s.c_str());
	for(int i=0;i<sizerow();i++){
		for(int j=0;j<sizecol();j++)fprintf(stream,"%9.4g ",(double)SELF[i,j]);
		fprintf(stream,"\n");
		}
	}




}//pp