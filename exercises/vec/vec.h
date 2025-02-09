#ifndef HAVE_VEC_H
#define HAVE_VEC_H
#include<iostream>
#include<string>
struct vec{
    double x,y,z;
    vec(double x,double y,double z); // parametrized constructor
    vec():vec(0,0,0){} // default constructor
    vec(const vec&)=default; // copy constructor
    vec(vec&&)=default; // move constructor
    ~vec()=default; // destructor
    vec& operator=(const vec&)=default; // copy assignment
    vec& operator=(vec&&)=default; // move assignment
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(double);
    vec& operator/=(double);
    double dot(const vec& v) const;
    vec cross(const vec& v) const;
    void set(double a,double b,double c){x=a;y=b;z=c;}
    void print(std::string s="") const; // for debugging
    friend std::ostream& operator<<(std::ostream&, const vec&);
};
vec operator-(const vec&);
vec operator-(const vec&, const vec&);
vec operator+(const vec&, const vec&);
vec operator*(const vec&, double);
vec operator*(double, const vec&);
vec operator/(const vec&, double);
bool approx(const vec&, const vec&, double acc=1e-6,double eps=1e-6);
#endif