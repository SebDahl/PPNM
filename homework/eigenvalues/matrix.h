#ifndef HAVE_MATRIX_H
#define HAVE_MATRIX_H
#ifdef LONG_DOUBLE
	#define NUMBER long double
#else
	#define NUMBER double
#endif
#include<string>
#include<vector>
namespace pp{
struct vector {
	std::vector<NUMBER> data;
	vector(int n) : data(n) {}
	vector()			=default;
	vector(const vector&)		=default;
	vector(vector&&)		=default;
	vector& operator=(const vector&)=default;
	vector& operator=(vector&&)	=default;
	int size() const {return data.size();}
	void resize(int n) {data.resize(n);}
	NUMBER& operator[](int i) {return data[i];}
	const NUMBER& operator[](int i) const {return data[i];}
	vector& operator+=(const vector&);
	vector& operator-=(const vector&);
	vector& operator*=(NUMBER);
	vector& operator/=(NUMBER);
	void print(std::string s="",FILE* stream=stdout) const;
	void append(const double a);

	static void write(const vector& v, const std::string& filename);
	static vector read(const std::string& filename);
};

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const vector&, NUMBER);
vector operator*(NUMBER, const vector&);
vector operator/(const vector&, NUMBER);
bool operator==(const vector&, const vector&);
bool operator!=(const vector&, const vector&);
bool approx_equal(const vector&, const vector&, double tol=1e-5);


struct matrix {
	std::vector<vector> cols;
	matrix()=default;
	matrix(int n,int m) : cols(m, vector(n)) {}
	matrix(const matrix& other)=default;
	matrix(matrix&& other)=default;
	matrix& operator=(const matrix& other)=default;
	matrix& operator=(matrix&& other)=default;
	int sizerow() const {return cols.empty() ? 0 : cols[0].size(); }
	int sizecol() const {return cols.size();}
	void resize(int n, int m);
	
	NUMBER get (int i, int j) const {return cols[j][i];}
	void set(int i, int j, NUMBER value){cols[j][i] = value;}
	NUMBER& operator()(int i, int j) { return cols[j][i]; }       // Non-const (modifies matrix)
	const NUMBER& operator()(int i, int j) const { return cols[j][i]; }  // Const version (read-only)
	

	NUMBER& operator[](int i, int j){return cols[j][i];}
	const NUMBER& operator[](int i, int j) const {return cols[j][i];}
	vector& operator[](int i){return cols[i];}
	const vector& operator[](int i) const {return cols[i];}
	vector get_col(int j);
	void set_col(int j,vector& cj);
	matrix transpose() const;
	static matrix identity(int n);
	static void write(const matrix& A, const std::string& filename);
	static matrix read(const std::string& filename);

	matrix& operator+=(const matrix&);
	matrix& operator-=(const matrix&);
	matrix& operator*=(const matrix&);
	matrix& operator*=(const NUMBER);
	matrix& operator/=(const NUMBER);
	matrix  operator^(int);

	//void print(const char* s="");
	void print(std::string s="",FILE* stream=stdout);
};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, NUMBER);
matrix operator*(NUMBER, const matrix&);
matrix operator/(const matrix&, NUMBER);
bool operator==(const matrix&, const matrix&);
bool operator!=(const matrix&, const matrix&);
bool approx_equal(const matrix&, const matrix&, double tol=1e-5);
vector operator*(const matrix&, const vector&);

}
#endif