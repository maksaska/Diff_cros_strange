#include <vector>
#include <iostream>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <ctime>
#include <chrono>
#include "TApplication.h"
#include "TRint.h"
#include "TText.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TSpline.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TPad.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TF1.h"
#include <getopt.h>
#include <fstream>

using namespace std;

extern double E0, Q2_min_, Q2_max_, W_min_, W_max_, cos_min, cos_max, phi_min, phi_max;
extern bool h_L, ratio_str, add_factor, W_sys;
extern int err_option;

extern double m_p, m_K, m_S, m_L;

extern vector<vector<double>> Data1;
extern vector<vector<double>> Data2;
extern vector<vector<double>> Data3;

extern vector<vector<double>> Data_Diff;
extern vector<vector<double>> Data_Sigma;

extern vector<vector<double>> Data_Q2cos;
extern vector<vector<double>> cosQ2_grid;

extern vector<vector<double>> W_syserr;

/*   ---   Interpolation.cpp   ---   */

void Reading(string Path, vector<vector<double>>& V);            /*   Reads data file   */
void input_check(int argc, char* argv[]);                        /*  input parameters handler   */
void Greet_message();
void Transfer(vector<vector<double>>& from_, vector<vector<double>>& to_); /*   Transfer data to appropriate format   */
void Print_data(vector<vector<double>>& V);                         /*   Prints vector<vector<double>>   */
double fRand(const double& fMin, const double& fMax);               /*  Gives random double from the interval   */
vector<double> cub_interp(vector<vector<double>>& V, const double x);
vector<double> cub_interp_err(vector<vector<double>>& V, const double x);
vector<double> interp_cub(vector<vector<double>>& V, const double& cos, const bool& statement);
vector<double> interp_cub_TT_LT(vector<vector<double>>& V, const double& cos_th);
vector<double> approx_cos_leg(vector<vector<double>>& V, const double& cos_th, const bool& statement, const bool& S_LT_check);
vector<double> approx_cos_leg_Sigma(vector<vector<double>>& V, const double& cos_th);
vector<double> extrapolate_higher(vector<vector<double>>& V, const double& cos_th, const bool& statement);
vector<double> Cubic_fit(vector<vector<double>>& V, const double& arg);
double eps(const double& W, const double& Q2);
double eps_beam(const double& W, const double& Q2, const double& E_beam);
bool isData1(const double& W, const double& Q2);
bool isData2(const double& W, const double& Q2);
bool isData3(const double& W, const double& Q2);
vector<double> linear_inter(const double& x1, const double& x2, const double& y1, const double& y2, const double& dy1, const double& dy2, const double& x);
vector<double> giveData1(const double& W, const double& Q2, const double& cos);
vector<double> giveData2(const double& W, const double& Q2, const double& cos);
vector<double> giveData3(const double& W, const double& Q2, const double& cos);
bool isW1(const double& W);
bool isW2(const double& W);
bool isW3(const double& W);
bool isWforQ2inter(const double& W, const double& Q2);
bool isSigma(const double& W);
vector<double> Photo_diff_fixed(const double& W, const double& cos_th);
vector<double> Photo_Sigma_fixed(const double& W, const double& cos_th);
vector<double> Photo_diff(const double& W, const double& cos_th);
vector<double> Photo_Sigma(const double& W, const double& cos_th);
vector<double> Photo_Sigma_TT(const double& W, const double& cos_th);
vector<double> lower_Q2_int(const double& W, const double& Q2, const double& cos_th);
vector<double> plane_inter(vector<vector<double>>& V, const double& x, const double& y);
bool triangular(vector<vector<double>>& V, const double& x, const double& y);
bool check_line(vector<vector<double>>& V);
vector<vector<double>> trig_points(const double& cos_th, const double& Q2);
vector<double> trig_points_plane(const double& cos_th, const double& Q2);
vector<double> lower_Q2_int_with_ratio(const double& W, const double& Q2, const double& cos_th);
vector<double> higher_Q2_int(const double& W, const double& Q2, const double& cos_th);
double get_sys_W_err_cos(const double& W, const double& cos_th);
double get_sys_W_err(const double& W, const double& cos_th);
vector<double> extrapolate_W(const double& W, const double& Q2, const double& cos_th);
vector<double> extrapolate_W_with_ratio(const double& W, const double& Q2, const double& cos_th);
double error_handler(vector<double>& V, const double& average);
vector<double> Point_diff_phi(const double& W, const double& Q2, const double& cos_th, const double& phi);

/*   ---   Functions you can use  ---   */

vector<double> Str_func_all(const double& W, const double& Q2, const double& cos_th); // str. functions St dSt Sl dSl Slt dSlt Stt dStt [nb/sr] in the point (W, Q2, cos)
vector<double> Point_diff(const double& W, const double& Q2, const double& cos, const double& phi); // cross section S0 dS0 [nb/sr] in the point (W, Q2, cos)
vector<double> Average_CS_stat(); // Average cross section stat. method [nb/sr]
vector<double> Average_CS(); // Average cross section [nb/sr]
vector<double> Average_CS_phi(); // integrted over phi + divided by 2PI cross section [nb/sr]
