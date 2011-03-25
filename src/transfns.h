#ifndef TRANSFNS_H
#define TRANSFNS_H



typedef double (*pt2trans)(double);


double trans_log2(double x);
double trans_loge(double x);
double trans_log10(double x);
double trans_sqrt(double x);
double trans_cuberoot(double x);

pt2trans transFunc(int code);


#endif
