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

extern double E0, Q2_min, Q2_max, W_min, W_max, cos_min, cos_max, W_p, Q2_p, cos_p, phi_p;
extern bool h_L, average;

extern double m_p, m_K, m_S, m_L; 

extern vector<vector<double>> Data;
extern vector<vector<double>> Data2;
extern vector<vector<double>> Data3;

/*   ---   Interpolation.cpp   ---   */

void Reading(string Path, vector<vector<double>>& V);
void input_check(int argc, char* argv[]);
void Transfer(vector<vector<double>>& from_, vector<vector<double>>& to_);
void Print_data(vector<vector<double>>& V);
double fRand(const double& fMin, const double& fMax);
vector<double> cub_interp(vector<vector<double>>& V, const double x);
vector<double> interp_cub(vector<vector<double>>& V, const double& cos, const bool& statement);
vector<double> approx_cos_leg(vector<vector<double>>& V, const double& cos, const bool& statement);
double eps(const double& W, const double& Q2);
bool isData(const double& W, const double& Q2);
bool isData2(const double& W, const double& Q2);
bool isData3(const double& W, const double& Q2);
vector<double> giveData(const double& W, const double& Q2, const double& cos);
vector<double> giveData2(const double& W, const double& Q2, const double& cos);
vector<double> giveData3(const double& W, const double& Q2, const double& cos);
vector<double> Str_func(const double& W, const double& Q2, const double& cos);
vector<double> Average_diff(const double& W, const double& Q2, const double& cos);
void Point_diff();
