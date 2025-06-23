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
	void append(const vector& a);
	double norm() const;
	vector copy() const;
	vector(std::initializer_list<NUMBER> init) : data(init) {}
	vector(const std::vector<NUMBER>& init) : data(init) {}
	vector(std::initializer_list<vector> init) {
		for (const auto& v : init)
			data.insert(data.end(), v.data.begin(), v.data.end());
	}
	double mean() const;
	double back() const;
	double std() const;
	int argmax() const;
	int argmin() const;
	double var() const;

	static void write(const vector& v, const std::string& filename);
	static vector loadtxt(const std::string& filename);
	// bool has_value();
};

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const vector&, NUMBER);
vector operator*(NUMBER, const vector&);
double operator*(const vector&, const vector&);
vector operator/(const vector&, NUMBER);
bool operator==(const vector&, const vector&);
bool operator!=(const vector&, const vector&);
bool approx_equal(const vector&, const vector&, double tol=1e-5);
bool has_value(const vector&, double tol=1e-5);

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
	static matrix loadtxt(const std::string& filename);
	pp::vector back() const;

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
bool has_value(const matrix&, double tol=1e-5);
}
#endif