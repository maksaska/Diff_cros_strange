#include "header.h"

using namespace std;

double E0(6.535), Q2_min_(0), Q2_max_(0), W_min_(0), W_max_(0), cos_min(0), cos_max(0), phi_min(0), phi_max(0);
bool h_L(true); double seed(time(0));

vector<vector<double>> Data1;
vector<vector<double>> Data2;
vector<vector<double>> Data3;

vector<vector<double>> Data_Diff;
vector<vector<double>> Data_Sigma;

double m_p(0.93827), m_K(0.498), m_S(1.1926), m_L(1.1157);

void Reading(string Path, vector<vector<double>>& V)
{
	string line; stringstream ss;
	ifstream File;
	double dub; vector<double> Numbers;

	File.open(Path,fstream::in | fstream::out | fstream::app);

	if (!File.is_open())
	{
		cout << "Can't open " << Path << " !" << endl;
	} else
	{
		while(!File.eof())
		{
			getline(File,line);
			ss << line;

			if(File.eof())
			{
				break;
			}

			while(ss>>dub)
			{
				Numbers.push_back(dub);
				while(ss.peek() == ' ')
            			ss.ignore();
			}

			V.push_back(Numbers);

			Numbers.clear();

			ss.clear();
		}
	}
	File.close();
}

void input_check(int argc, char* argv[])
{
	const char* short_options = "e:hz:x:c:v:b:n:u:i:o:p:m:l:r:"; int rez; int option_index;

	const struct option long_options[] = {
						{"beam_energy", required_argument, NULL, 'e'},
	        				{"KSigma0", no_argument, NULL, 'h'},
	        				{"Q2_min", required_argument, NULL, 'z'},
	        				{"Q2_max", required_argument, NULL, 'x'},
	        				{"W_min", required_argument, NULL, 'c'},
	        				{"W_max", required_argument, NULL, 'v'},
	        				{"cos_min", required_argument, NULL, 'b'},
	        				{"cos_max", required_argument, NULL, 'n'},
	        				{"phi_min", required_argument, NULL, 'm'},
	        				{"phi_max", required_argument, NULL, 'l'},
							{"SEED", required_argument, NULL, 'r'},
	        				{NULL, 0, NULL, 0}
								};
	while ((rez=getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
	{
		switch(rez)
		{
			case 'e': {
				E0 = atof(optarg);
				break;
			};
			case 'h': {
				h_L = false;
				break;
			};
			case 'z': {
				Q2_min_ = atof(optarg);
				break;
			};
			case 'x': {
				Q2_max_ = atof(optarg);
				break;
			};
			case 'c': {
				W_min_ = atof(optarg);
				break;
			};
			case 'v': {
				W_max_ = atof(optarg);
				break;
			};
			case 'b': {
				cos_min = atof(optarg);
				break;
			};
			case 'n': {
				cos_max = atof(optarg);
				break;
			};
			case 'm': {
				phi_min = M_PI*atof(optarg)/180;
				break;
			};
			case 'l': {
				phi_max = M_PI*atof(optarg)/180;
				break;
			};
			case 'r': {
				seed = atof(optarg);
				break;
			};
			case '?': default: {
				cerr << "Unkhown option" << endl;
				break;
			};
		};
	};

	double srand(seed);

    if(W_min_ < 1.5) W_min_ = 1.61;
	if(Q2_min_ < 0) Q2_min_ = 0;
	if(Q2_max_ < Q2_min_) Q2_max_ = Q2_min_ + 0.05;
	if(W_max_ < W_min_) W_max_ = W_min_ + 0.01;
	if(cos_max < cos_min) cos_max = cos_min + 0.01;
	if(cos_min < -1) cos_min = -1;
	if(cos_min > 1) cos_min = 1;
	if(cos_max < -1) cos_max = -1;
	if(cos_max > 1) cos_max = 1;

	vector<vector<double>> buff;

	if(h_L)
	{
		Reading("./Data/P1.dat", buff);
		Transfer(buff, Data1); buff.clear();
		Reading("./Data/P2.dat", buff);
		Transfer(buff, Data2); buff.clear();
		Reading("./Data/P3.dat", buff);
		Transfer(buff, Data3); buff.clear();
		Reading("./Data/Diff_L_Photo.dat", Data_Diff);
		Reading("./Data/Sigma_L_Photo.dat", Data_Sigma);
	}
	else
	{
		Reading("./Data/P4.dat", buff);
		Transfer(buff, Data1); buff.clear();
		Reading("./Data/P5.dat", buff);
		Transfer(buff, Data2); buff.clear();
		Reading("./Data/P6.dat", buff);
		Transfer(buff, Data3); buff.clear();
		Reading("./Data/Diff_S_Photo.dat", Data_Diff);
		Reading("./Data/Sigma_S_Photo.dat", Data_Sigma);
	}

	buff.clear();
}

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void Transfer(vector<vector<double>>& from_, vector<vector<double>>& to_)
{
	vector<double> buff;

	for(auto i:from_)
	{
		buff.push_back(i[0]);
		buff.push_back(i[1]);
		buff.push_back(i[3]);
		buff.push_back(i[4]/(1 + i[2]/5));
		buff.push_back(i[5]/(1 + i[2]/5));
		buff.push_back(i[4]/(5 + i[2]));
		buff.push_back(i[5]/(5 + i[2]));

		//buff.push_back(i[4]);
		//buff.push_back(i[5]);

		buff.push_back(i[6]);
		buff.push_back(i[7]);
		buff.push_back(i[8]);
		buff.push_back(i[9]);

		to_.push_back(buff); buff.clear();
	}
}

void Print_data(vector<vector<double>>& V)
{
	for(auto i:V)
	{
		for(auto j:i)
		{
			cout << j << "\t";
		}
		cout << endl;
	}
	cout << "\nNumber of rows: " << V.size() << endl;
}

vector<double> cub_interp(vector<vector<double>>& V, const double x)
{
	vector<double> result;

	double x1, x2, x3, y1, y2, y3, dy1, dy2, dy3;
	double a, da, b, db, c, dc;
	double delta;
	double f, df;

	x1 = V[0][0]; x2 = V[1][0]; x3 = V[2][0];
	y1 = V[0][1]; y2 = V[1][1]; y3 = V[2][1];
	dy1 = V[0][2] + y1; dy2 = V[1][2] + y2; dy3 = V[2][2] + y3;

	delta = (x2 - x3)*(x1 - x2)*(x1 - x3);

	a = ((x2 - x3)*y1 + (x3 - x1)*y2 + (x1 - x2)*y3)/delta;
	b = ((x3*x3 - x2*x2)*y1 + (x1*x1 - x3*x3)*y2 + (x2*x2 - x1*x1)*y3)/delta;
	c = (x2*x3*(x2 - x3)*y1 + x1*x3*(x3 - x1)*y2 + x1*x2*(x1 - x2)*y3)/delta;

	f = a*x*x + b*x + c;

	a = ((x2 - x3)*dy1 + (x3 - x1)*dy2 + (x1 - x2)*dy3)/delta;
	b = ((x3*x3 - x2*x2)*dy1 + (x1*x1 - x3*x3)*dy2 + (x2*x2 - x1*x1)*dy3)/delta;
	c = (x2*x3*(x2 - x3)*dy1 + x1*x3*(x3 - x1)*dy2 + x1*x2*(x1 - x2)*dy3)/delta;

	df = a*x*x + b*x + c - f;

	result.push_back(f);
	result.push_back(df);

	return result;
}

vector<double> cub_interp_err(vector<vector<double>>& V, const double x)
{
	vector<double> result;

	double x1, x2, x3, y1, y2, y3, dy1, dy2, dy3;
	double a, da, b, db, c, dc;
	double delta;
	double f, df;

	x1 = V[0][0]; x2 = V[1][0]; x3 = V[2][0];
	y1 = V[0][1]; y2 = V[1][1]; y3 = V[2][1];
	dy1 = V[0][2]; dy2 = V[1][2]; dy3 = V[2][2];

	delta = (x2 - x3)*(x1 - x2)*(x1 - x3);

	a = ((x2 - x3)*y1 + (x3 - x1)*y2 + (x1 - x2)*y3)/delta;
	b = ((x3*x3 - x2*x2)*y1 + (x1*x1 - x3*x3)*y2 + (x2*x2 - x1*x1)*y3)/delta;
	c = (x2*x3*(x2 - x3)*y1 + x1*x3*(x3 - x1)*y2 + x1*x2*(x1 - x2)*y3)/delta;

	f = a*x*x + b*x + c;

	da = sqrt(pow((x2 - x3)*dy1/delta, 2) + pow((x3 - x1)*dy2/delta, 2) + pow((x1 - x2)*dy3/delta, 2));
	db = sqrt(pow((x3*x3 - x2*x2)*dy1/delta, 2) + pow((x1*x1 - x3*x3)*dy2/delta, 2) + pow((x2*x2 - x1*x1)*dy3/delta, 2));
	dc = sqrt(pow(x2*x3*(x2 - x3)*dy1/delta, 2) + pow(x1*x3*(x3 - x1)*dy2/delta, 2) + pow(x1*x2*(x1 - x2)*dy3/delta, 2));

	df = sqrt(pow(x*x*da, 2) + pow(x*db, 2) + pow(dc, 2));

	result.push_back(f);
	result.push_back(df);

	return result;
}

vector<double> interp_cub(vector<vector<double>>& V, const double& cos_th, const bool& statement)
{
	vector<vector<double>> transit;
	vector<double> result, buff;
	double func, dfunc;
	int t = V.size() - 1;

	if(cos_th < V[1][0])
	{
		buff.push_back(V[0][0]); buff.push_back(V[0][1]); buff.push_back(V[0][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[1][0]); buff.push_back(V[1][1]); buff.push_back(V[1][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[2][0]); buff.push_back(V[2][1]); buff.push_back(V[2][2]);
		transit.push_back(buff); buff.clear();
	}

	if(cos_th >= V[t-1][0])
	{
		buff.push_back(V[t-2][0]); buff.push_back(V[t-2][1]); buff.push_back(V[t-2][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[t-1][0]); buff.push_back(V[t-1][1]); buff.push_back(V[t-1][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[t][0]); buff.push_back(V[t][1]); buff.push_back(V[t][2]);
		transit.push_back(buff); buff.clear();
	}

	if(cos_th >= V[1][0] and cos_th < V[t-1][0])
	{
		for(int i = 0; i < V.size(); i++)
		{
			if(cos_th >= V[i][0] and cos_th < V[i+1][0])
			{
				buff.push_back(V[i][0]); buff.push_back(V[i][1]); buff.push_back(V[i][2]);
				transit.push_back(buff); buff.clear();
				buff.push_back(V[i+1][0]); buff.push_back(V[i+1][1]); buff.push_back(V[i+1][2]);
				transit.push_back(buff); buff.clear();
				buff.push_back(V[i+2][0]); buff.push_back(V[i+2][1]); buff.push_back(V[i+2][2]);
				transit.push_back(buff); buff.clear();

				break;
			}
		}
	}

	buff = cub_interp(transit, cos_th);

	func = buff[0];
	dfunc = buff[1];

	if(func < 0 and statement){func = 0;}

	result.push_back(func);
	result.push_back(dfunc); transit.clear(); buff.clear();

	return result;
}

vector<double> approx_cos_leg(vector<vector<double>>& V, const double& cos_th, const bool& statement)
{
	vector<double> result;
	double func, dfunc;
	double arg = cos_th;
	vector<double> buff;

	double *X = new double[V.size()];
	double *Y = new double[V.size()];
	double *dY = new double[V.size()]; int count(0);

	for(auto i:V)
	{
			X[count] = i[0];
			Y[count] = i[1];
			dY[count] = i[2]; count++;
	}

	TGraphErrors* gr;
	gr = new TGraphErrors(V.size(), X, Y, NULL, dY);

	TF1 *ff;
ff = new TF1("ff", "[0] + [1]*x + [2]*0.5*(3*pow(x,2) - 1) + [3]*0.5*(5*pow(x,3) - 3*x) + [4]*(35*pow(x, 4) - 30*pow(x, 2) + 3)/8", -1, 1);
	gr->Fit(ff, "Q");

	double A = ff->GetParameter(0);
	double B = ff->GetParameter(1);
	double C = ff->GetParameter(2);
	double D = ff->GetParameter(3);
	double E = ff->GetParameter(4);

	double dA = ff->GetParError(0);
	double dB = ff->GetParError(1);
	double dC = ff->GetParError(2);
	double dD = ff->GetParError(3);
	double dE = ff->GetParError(4);

	delete [] X; delete [] Y; delete [] dY;
	ff->Clear(); delete ff;
	gr->Clear(); delete gr;

	func =  A + B*arg + C*0.5*(3*pow(arg,2) - 1) + D*0.5*(5*pow(arg,3) - 3*arg) + E*(35*pow(arg, 4) - 30*pow(arg, 2) + 3)/8;

	if(func < 0 and statement){func = 0;}
	buff = interp_cub(V, cos_th, statement);

	dfunc = buff[1];

	if(cos_th < V[0][0])
	{
		arg = V[0][0];
		//dfunc = sqrt(pow((func-buff[0])/2, 2) + dfunc*dfunc);
		func = A + B*arg + C*0.5*(3*pow(arg,2) - 1) + D*0.5*(5*pow(arg,3) - 3*arg) + E*(35*pow(arg, 4) - 30*pow(arg, 2) + 3)/8;
	}

	if(cos_th > V[V.size()-1][0])
	{
		arg = V[V.size()-1][0];
		//dfunc = sqrt(pow((func-buff[0])/2, 2) + dfunc*dfunc);
		func = A + B*arg + C*0.5*(3*pow(arg,2) - 1) + D*0.5*(5*pow(arg,3) - 3*arg) + E*(35*pow(arg, 4) - 30*pow(arg, 2) + 3)/8;
	}

	/*if(cos_th < V[0][0])
	{
		dfunc = V[0][2];
		func = V[0][1];
	}

	if(cos_th > V[V.size()-1][0])
	{
		dfunc = V[V.size()-1][2];
		func = V[V.size()-1][1];
	}*/

	result.push_back(func);
	result.push_back(dfunc); buff.clear();

	return result;
}

double eps(const double& W, const double& Q2)
{
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	return  1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));
}

double eps_beam(const double& W, const double& Q2, const double& E_beam)
{
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	return  1/(1 + 2*(nu*nu + Q2)/(4*(E_beam - nu)*E_beam - Q2));
}

bool isData1(const double& W, const double& Q2)
{
	if(h_L)
	{
		if(Q2 <= 1.0 and Q2 >= 0.65 and W <= 1.975 and W >= 1.65){ return true;}
		else{ return false;}
	} else
	{
		if(Q2 <= 1.0 and Q2 >= 0.65 and W <= 2.05 and W >= 1.725){ return true;}
		else{ return false;}
	}
}

bool isData2(const double& W, const double& Q2)
{
	if(h_L)
	{
	if((Q2 <= 2.55 and Q2 >= 1.0 and W <= 2.15 and W >= 1.65) or (Q2 <= 2.05 and Q2 >= 1.0 and W <= 2.35 and W >= 2.15)){ return true;}
	else{ return false;}
	} else
	{
	if((Q2 <= 2.55 and Q2 >= 1.0 and W <= 2.15 and W >= 1.75) or (Q2 <= 2.05 and Q2 >= 1.0 and W <= 2.35 and W >= 2.15)){ return true;}
	else{ return false;}
	}
}

bool isData3(const double& W, const double& Q2)
{
	if(h_L)
	{
		if(Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.575 and W >= 1.63){ return true;}
		else{ return false;}
	} else
	{
		if(Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.575 and W >= 1.695){ return true;}
		else{ return false;}
	}
}

bool isW1(const double& W)
{
	if(h_L)
	{
		if(W <= 1.975 and W >= 1.65){ return true;}
		else{ return false;}
	} else
	{
		if(W <= 2.05 and W >= 1.725){ return true;}
		else{ return false;}
	}
}

bool isW2(const double& W)
{
	if(h_L)
	{
	if(W >= 1.65 and W <= 2.35){ return true;}
	else{ return false;}
	} else
	{
	if(W >= 1.75 and W <= 2.35){ return true;}
	else{ return false;}
	}
}

bool isW3(const double& W)
{
	if(h_L)
	{
		if(W <= 2.575 and W >= 1.63){ return true;}
		else{ return false;}
	} else
	{
		if(W <= 2.575 and W >= 1.695){ return true;}
		else{ return false;}
	}
}

bool isSigma(const double& W)
{
	if(h_L)
	{
		if(W <= 2.18 and W >= 1.72){ return true;}
		else{ return false;}
	} else
	{
		if(W <= 2.17 and W >= 1.78){ return true;}
		else{ return false;}
	}
}

vector<double> linear_inter(const double& x1, const double& x2, const double& y1, const double& y2, const double& dy1, const double& dy2, const double& x)
{
	vector<double> result;

	if(x1 == x2)
	{
		result.push_back(y1);
		result.push_back(dy1);

		return result;
	}

	result.push_back((y2*(x - x1) + y1*(x2 - x))/(x2 - x1));
	result.push_back(sqrt(pow(dy2*(x - x1)/(x2 - x1), 2) + pow(dy1*(x2 - x)/(x2 - x1), 2)));

	return result;
}

vector<double> giveData1(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;
	double S, dS; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;

	Q2_min = 0.65;
	Q2_max = 1.0;

	W_min = double(floor((W - 0.025)*20))/20 + 0.025;
	W_max = W_min + 0.05;

	if(W < 1.725)
	{
		W_min = 1.65;
		W_max = 1.725;
	}

	if(W > 1.925)
	{
		W_min = 1.925;
		W_max = 1.975;
	}

	for(auto i:Data1)
	{
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_1.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_5.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_9.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_13.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_2.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_6.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_10.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_14.push_back(buff); buff.clear();
		}
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_3.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_7.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_11.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_15.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_4.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_8.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_12.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_16.push_back(buff); buff.clear();
		}
	}

	buff = approx_cos_leg(pack_1, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_2, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_3, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_4, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();

	buff = approx_cos_leg(pack_5, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_6, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_7, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_8, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();


	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();

	buff = approx_cos_leg(pack_9, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_10, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_11, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_12, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();

	buff = approx_cos_leg(pack_13, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_14, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_15, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_16, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();

	return result;
}

vector<double> giveData2(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;
	double S, dS; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;

	Q2_min = double(floor((Q2 - 0.05)*2))/2 + 0.05;
	Q2_max = Q2_min + 0.5;

	if(Q2 < 1.55)
	{
		Q2_min = 1.0;
		Q2_max = 1.55;
	}

	if(abs(Q2 - 2.55) < 1e-9)
	{
		Q2_min = 2.05;
		Q2_max = 2.55;
	}

	W_min = double(floor((W - 0.05)*10))/10 + 0.05;
	W_max = W_min + 0.1;

	if(abs(W - 2.15) < 1e-9)
	{
		W_min = 2.05;
		W_max = 2.15;
	}

	if(abs(W - 1.65) < 1e-9)
	{
		W_min = 1.65;
		W_max = 1.75;
	}

	if(abs(W - 1.75) < 1e-9)
	{
		W_min = 1.75;
		W_max = 1.85;
	}

	if(abs(W - 2.35) < 1e-9)
	{
		W_min = 2.25;
		W_max = 2.35;
	}

	for(auto i:Data2)
	{
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_1.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_5.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_9.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_13.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_2.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_6.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_10.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_14.push_back(buff); buff.clear();
		}
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_3.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_7.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_11.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_15.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_4.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_8.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_12.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_16.push_back(buff); buff.clear();
		}
	}

	buff = approx_cos_leg(pack_1, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_2, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_3, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_4, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();

	buff = approx_cos_leg(pack_5, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_6, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_7, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_8, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();


	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();

	buff = approx_cos_leg(pack_9, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_10, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_11, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_12, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();

	buff = approx_cos_leg(pack_13, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_14, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_15, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_16, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();

	return result;
}

vector<double> giveData3(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;
	double S, dS; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;

	W_min = double(floor((W - 0.025)*20))/20 + 0.025;
	W_max = W_min + 0.05;

	if(W < 1.675 and h_L)
	{
		W_min = 1.63;
		W_max = 1.675;
	}

	if(W < 1.725 and !h_L)
	{
		W_min = 1.695;
		W_max = 1.725;
	}

	if(abs(W - 2.575) < 1e-9)
	{
		W_min = 2.525;
		W_max = 2.575;
	}

	if(Q2 < 2.6)
	{
		Q2_min = 1.8;
		Q2_max = 2.6;
	}

	if(2.6 <= Q2)
	{
		Q2_min = 2.6;
		Q2_max = 3.45;
	}

	for(auto i:Data3)
	{
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_1.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_5.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_9.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_13.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_min - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_2.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_6.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_10.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_14.push_back(buff); buff.clear();
		}
		if(abs(W_min - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_3.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_7.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_11.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_15.push_back(buff); buff.clear();
		}
		if(abs(W_max - i[0]) < 1e-9 and abs(Q2_max - i[1]) < 1e-9)
		{
			buff.push_back(i[2]);
			buff.push_back(i[3]);
			buff.push_back(i[4]);
			pack_4.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[5]);
			buff.push_back(i[6]);
			pack_8.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[7]);
			buff.push_back(i[8]);
			pack_12.push_back(buff); buff.clear();

			buff.push_back(i[2]);
			buff.push_back(i[9]);
			buff.push_back(i[10]);
			pack_16.push_back(buff); buff.clear();
		}
	}

	buff = approx_cos_leg(pack_1, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_2, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_3, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_4, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();

	buff = approx_cos_leg(pack_5, cos_th, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_6, cos_th, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_7, cos_th, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_8, cos_th, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();


	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();

	buff = approx_cos_leg(pack_9, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_10, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_11, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_12, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();

	buff = approx_cos_leg(pack_13, cos_th, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_14, cos_th, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_15, cos_th, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_16, cos_th, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y1, y2, dy1, dy2, W);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();

	buff = linear_inter(W_min, W_max, y3, y4, dy3, dy4, W);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	buff = linear_inter(Q2_min, Q2_max, S1, S2, dS1, dS2, Q2);
	S = buff[0]; dS = buff[1]; buff.clear();

	result.push_back(S);
	result.push_back(dS);

	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();

	return result;
}

vector<double> Photo_diff_fixed(const double& W, const double& cos_th)
{
	vector<double> result, buff;
	vector<vector<double>> Block;

	for(int i = 0; i < Data_Diff.size(); i++)
	{
		if(abs(W - Data_Diff[i][0]) < 1e-8)
		{
			buff.push_back(Data_Diff[i][1]);
			buff.push_back(Data_Diff[i][2]);
			buff.push_back(Data_Diff[i][3]);
			Block.push_back(buff);
			buff.clear();
		}
	}

	result = approx_cos_leg(Block, cos_th, true);

	Block.clear(); buff.clear();

	return result;
}

vector<double> Photo_diff(const double& W, const double& cos_th)
{
	vector<double> result, buff;
	double S1, S2, dS1, dS2;

	double W_min = double(floor(W*100))/100 + 0.005;
	double W_max;

	if(W < W_min)
	{
		W_max = W_min;
		W_min = W_max - 0.01;
	}else
	{
		W_max = W_min + 0.01;
	}

	if(abs(W_min - 1.955) < 1e-8) W_min = 1.945;
	if(abs(W_max - 1.955) < 1e-8) W_max = 1.965;

	if(abs(W_min - 2.735) < 1e-8) W_min = 2.725;
	if(abs(W_max - 2.735) < 1e-8) W_max = 2.755;

	if(abs(W_min - 2.745) < 1e-8) W_min = 2.725;
	if(abs(W_max - 2.745) < 1e-8) W_max = 2.755;

	buff = Photo_diff_fixed(W_min, cos_th);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();
	buff = Photo_diff_fixed(W_max, cos_th);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	result.push_back((S2*(W - W_min) + S1*(W_max - W))/(W_max - W_min));
	result.push_back(sqrt(pow(dS2*(W - W_min)/(W_max - W_min), 2) + pow(dS1*(W_max - W)/(W_max - W_min), 2)));

	buff.clear();

	return result;
}

vector<double> Photo_Sigma_fixed(const double& W, const double& cos_th)
{
	vector<double> result, buff;
	vector<vector<double>> Block;

	for(int i = 0; i < Data_Sigma.size(); i++)
	{
		if(abs(W - Data_Sigma[i][0]) < 1e-8)
		{
			buff.push_back(Data_Sigma[i][1]);
			buff.push_back(Data_Sigma[i][2]);
			buff.push_back(Data_Sigma[i][3]);
			Block.push_back(buff);
			buff.clear();
		}
	}

	result = approx_cos_leg(Block, cos_th, true);

	Block.clear(); buff.clear();

	return result;
}

vector<double> Photo_Sigma(const double& W, const double& cos_th)
{
	vector<double> result, buff;
	double S1, S2, dS1, dS2;

	double W_min, W_max;

	if(h_L)
	{
		W_min = double(floor(W*50))/50;
		W_max = W_min + 0.02;
		if(W_max > 2.18)
		{
			W_max = 2.18;
			W_min = 2.16;
		}
	}else{
		W_min = double(floor((W-0.01)*25))/25 + 0.01;
		W_max  = W_min + 0.04;

		if(W >= 1.78 and W < 1.81)
		{
			W_min = 1.78;
			W_max = 1.81;
		}

		if(W > 2.13)
		{
			W_min = 2.13;
			W_max = 2.17;
		}
	}

	buff = Photo_Sigma_fixed(W_min, cos_th);
	S1 = buff[0]; dS1 = buff[1]; buff.clear();
	buff = Photo_Sigma_fixed(W_max, cos_th);
	S2 = buff[0]; dS2 = buff[1]; buff.clear();

	result.push_back((S2*(W - W_min) + S1*(W_max - W))/(W_max - W_min));
	result.push_back(sqrt(pow(dS2*(W - W_min)/(W_max - W_min), 2) + pow(dS1*(W_max - W)/(W_max - W_min), 2)));

	buff.clear();

	return result;
}

vector<double> Photo_Sigma_TT(const double& W, const double& cos_th)
{
	vector<double> result, buff1, buff2;

	buff1 = Photo_Sigma(W, cos_th);
	buff2 = Photo_diff(W, cos_th);

	result.push_back(-buff1[0]*buff2[0]);
	result.push_back(sqrt(pow(buff1[0]*buff2[1], 2) + pow(buff2[0]*buff1[1], 2)));
	//result.push_back(buff1[0]*buff2[0]*sqrt(pow(buff2[1]/buff2[0], 2) + pow(buff1[1]/buff1[0], 2)));

	buff1.clear(); buff2.clear();

	return result;
}

vector<double> lower_Q2_int(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;
	vector<double> buff1, buff2, buff, line;
	vector<vector<double>> Block;
	double Q2_1, Q2_2;

	buff = Photo_diff(W, cos_th);
	line.push_back(0);
	line.push_back(buff[0]);
	line.push_back(buff[1]);
	Block.push_back(line); line.clear(); buff.clear();

	if(isW1(W) and Q2 < 1.0)
	{
		Q2_1 = 0.65; Q2_2 = 1.0;
		buff1 = giveData1(W, Q2_1, cos_th);
		buff2 = giveData1(W, Q2_2, cos_th);
	}else{
		if(isW2(W))
		{
			Q2_1 = 1.0; Q2_2 = 1.55;
			buff1 = giveData2(W, Q2_1, cos_th);
			buff2 = giveData2(W, Q2_2, cos_th);
		}else{
			if(isW3(W))
			{
				Q2_1 = 1.8; Q2_2 = 2.6;
				buff1 = giveData3(W, Q2_1, cos_th);
				buff2 = giveData3(W, Q2_2, cos_th);
			}else
			{
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);

				return result;
			}
		}
	}

	line.push_back(Q2_1);
	line.push_back(buff1[0]);
	line.push_back(buff1[1]);
	Block.push_back(line); line.clear();


	line.push_back(Q2_2);
	line.push_back(buff2[0]);
	line.push_back(buff2[1]);
	Block.push_back(line); line.clear();

	result = cub_interp(Block, Q2); Block.clear();
	//result = cub_interp_err(Block, Q2); Block.clear();

	line.push_back(0);
	line.push_back(0);
	line.push_back(0);
	Block.push_back(line); line.clear();

	line.push_back(Q2_1);
	line.push_back(buff1[2]);
	line.push_back(buff1[3]);
	Block.push_back(line); line.clear();


	line.push_back(Q2_2);
	line.push_back(buff2[2]);
	line.push_back(buff2[3]);
	Block.push_back(line); line.clear();

	buff = cub_interp(Block, Q2); Block.clear();
	//buff = cub_interp_err(Block, Q2); Block.clear();
	result.push_back(buff[0]);
	result.push_back(buff[1]); buff.clear();

	line.push_back(0);
	line.push_back(0);
	line.push_back(0);
	Block.push_back(line); line.clear();

	line.push_back(Q2_1);
	line.push_back(buff1[4]);
	line.push_back(buff1[5]);
	Block.push_back(line); line.clear();


	line.push_back(Q2_2);
	line.push_back(buff2[4]);
	line.push_back(buff2[5]);
	Block.push_back(line); line.clear();

	buff = cub_interp(Block, Q2); Block.clear();
	//buff = cub_interp_err(Block, Q2); Block.clear();
	result.push_back(buff[0]);
	result.push_back(buff[1]); buff.clear();

	if(isSigma(W))
	{
		buff = Photo_Sigma_TT(W, cos_th);
		line.push_back(0);
		line.push_back(buff[0]);
		line.push_back(buff[1]);
		Block.push_back(line); line.clear(); buff.clear();
	}else
	{
		if(isW1(W))
		{
			buff = giveData1(W, (Q2_1 + Q2_2)/2, cos_th);
		}else{
			if(isW2(W))
			{
				buff = giveData2(W, (Q2_1 + Q2_2)/2, cos_th);
			}else{
				if(isW3(W))
				{
					buff = giveData3(W, (Q2_1 + Q2_2)/2, cos_th);
				}
			}
		}
		line.push_back((Q2_1 + Q2_2)/2);
		line.push_back(buff[6]);
		line.push_back(buff[7]);
		Block.push_back(line); line.clear(); buff.clear();
	}

	line.push_back(Q2_1);
	line.push_back(buff1[6]);
	line.push_back(buff1[7]);
	Block.push_back(line); line.clear();


	line.push_back(Q2_2);
	line.push_back(buff2[6]);
	line.push_back(buff2[7]);
	Block.push_back(line); line.clear();

	buff = cub_interp(Block, Q2); Block.clear();
	//buff = cub_interp_err(Block, Q2); Block.clear();
	result.push_back(buff[0]);
	result.push_back(buff[1]); buff.clear();

	buff1.clear(); buff2.clear();

	return result;
}

vector<double> extrapolate_higher(vector<vector<double>>& V, const double& Q2)
{
	vector<double> result;
	double func, dfunc;
	double arg = Q2;
	vector<double> buff;

	double *X = new double[V.size()];
	double *Y = new double[V.size()];
	double *dY = new double[V.size()]; int count(0);

	for(auto i:V)
	{
			X[count] = i[0];
			Y[count] = i[1];
			dY[count] = i[2]; count++;
	}

	TGraphErrors* gr;
	gr = new TGraphErrors(V.size(), X, Y, NULL, dY);

	TF1 *ff;
	//ff = new TF1("ff", "[0]*[0] + [1]*[1]/x + [2]*[2]/(x*x)", 1.8, 3.45);
	ff = new TF1("ff", "[0] + [1]/x + [2]/(x*x)", 1.8, 3.45);
	gr->Fit(ff, "Q");

	double A = ff->GetParameter(0);
	double B = ff->GetParameter(1);
	double C = ff->GetParameter(2);

	double dA = ff->GetParError(0);
	double dB = ff->GetParError(1);
	double dC = ff->GetParError(2);

	delete [] X; delete [] Y; delete [] dY;
	ff->Clear(); delete ff;
	gr->Clear(); delete gr;

	func =  A + B/arg + C/(arg*arg);
	//dfunc = sqrt(pow(dA, 2) + pow(dB/arg, 2) + pow(dC/(arg*arg), 2));

	buff = interp_cub(V, Q2, true);

	dfunc = buff[1];

	result.push_back(func);
	result.push_back(dfunc); buff.clear();

	return result;
}

vector<double> higher_Q2_int(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result, line;
	vector<double> buff1, buff2, buff3;
	vector<vector<double>> Block;

	buff1 = giveData3(W, 1.8, cos_th); //1.8
	buff2 = giveData3(W, 2.6, cos_th); // 2.6
	buff3 = giveData3(W, 3.45, cos_th); // 3.45

	for(int i = 0; i < 8; i += 2)
	{
		line.push_back(1.8);
		line.push_back(buff1[i]);
		line.push_back(buff1[i+1]);
		Block.push_back(line); line.clear();

		line.push_back(2.6);
		line.push_back(buff2[i]);
		line.push_back(buff2[i+1]);
		Block.push_back(line); line.clear();

		line.push_back(3.45);
		line.push_back(buff3[i]);
		line.push_back(buff3[i+1]);
		Block.push_back(line); line.clear();

		line = extrapolate_higher(Block, Q2);
		result.push_back(line[0]);
		result.push_back(line[1]); line.clear(); Block.clear();
	}

	buff1.clear(); buff2.clear(); buff3.clear();

	return result;
}

vector<double> Str_func(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result, buff;
	double St, dSt, Sl, dSl, Stt, dStt, Slt, dSlt;

	if(isData1(W, Q2))
	{
		buff = giveData1(W, Q2, cos_th);

		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();

		if(isData2(W, Q2))
		{
			buff = giveData2(W, Q2, cos_th);

			dSt = sqrt(dSt*dSt + buff[1]*buff[1] + pow(St - buff[0], 2))/2;
			St = (St + buff[0])/2;
			dSl = sqrt(dSl*dSl + buff[3]*buff[3] + pow(Sl - buff[2], 2))/2;
			Sl = (Sl + buff[2])/2;
			dSlt = sqrt(dSlt*dSlt + buff[5]*buff[5] + pow(Slt - buff[4], 2))/2;
			Slt = (Slt + buff[4])/2;
			dStt = sqrt(dStt*dStt + buff[7]*buff[7] + pow(Stt - buff[6], 2))/2;
			Stt = (Stt + buff[6])/2; buff.clear();
		}

		result.push_back(St);
		result.push_back(dSt);
		result.push_back(Sl);
		result.push_back(dSl);
		result.push_back(Slt);
		result.push_back(dSlt);
		result.push_back(Stt);
		result.push_back(dStt);

		return result;
	}

	if(isData2(W, Q2))
	{
		buff = giveData2(W, Q2, cos_th);

		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();

		/*if(isData3(W, Q2))
		{
			buff = giveData3(W, Q2, cos_th);

			dSt = sqrt(dSt*dSt + buff[1]*buff[1] + pow(St - buff[0], 2))/2;
			St = (St + buff[0])/2;
			dSl = sqrt(dSl*dSl + buff[3]*buff[3] + pow(Sl - buff[2], 2))/2;
			Sl = (Sl + buff[2])/2;
			dSlt = sqrt(dSlt*dSlt + buff[5]*buff[5] + pow(Slt - buff[4], 2))/2;
			Slt = (Slt + buff[4])/2;
			dStt = sqrt(dStt*dStt + buff[7]*buff[7] + pow(Stt - buff[6], 2))/2;
			Stt = (Stt + buff[6])/2; buff.clear();
		}*/

		result.push_back(St);
		result.push_back(dSt);
		result.push_back(Sl);
		result.push_back(dSl);
		result.push_back(Slt);
		result.push_back(dSlt);
		result.push_back(Stt);
		result.push_back(dStt);

		return result;
	}

	if(isData3(W, Q2))
	{
		buff = giveData3(W, Q2, cos_th);

		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();

		result.push_back(St);
		result.push_back(dSt);
		result.push_back(Sl);
		result.push_back(dSl);
		result.push_back(Slt);
		result.push_back(dSlt);
		result.push_back(Stt);
		result.push_back(dStt);

		return result;
	}

	if(Q2 < 3.45 and isW3(W))
	{
		buff = lower_Q2_int(W, Q2, cos_th);

		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();

		result.push_back(St);
		result.push_back(dSt);
		result.push_back(Sl);
		result.push_back(dSl);
		result.push_back(Slt);
		result.push_back(dSlt);
		result.push_back(Stt);
		result.push_back(dStt);

		return result;
	}

	if(Q2 > 3.45 and isW3(W))
	{
		buff = higher_Q2_int(W, Q2, cos_th);

		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();

		result.push_back(St);
		result.push_back(dSt);
		result.push_back(Sl);
		result.push_back(dSl);
		result.push_back(Slt);
		result.push_back(dSlt);
		result.push_back(Stt);
		result.push_back(dStt);

		return result;
	}

	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);

	return result;
}

vector<double> Cubic_fit(vector<vector<double>>& V, const double& W)
{
	vector<double> result;
	double func, dfunc;
	double arg = W;
	vector<double> buff;

	double *X = new double[V.size()];
	double *Y = new double[V.size()];
	double *dY = new double[V.size()]; int count(0);

	for(auto i:V)
	{
			X[count] = i[0];
			Y[count] = i[1];
			dY[count] = i[2]; count++;
	}

	TGraphErrors* gr;
	gr = new TGraphErrors(V.size(), X, Y, NULL, dY);

	double left, right;

	if(W < 2.0)
	{
		left = W;
		right = V[V.size()-1][0];
	}
	else
	{
		left = V[0][0];
		right = W;
	}

	TF1 *ff;
	ff = new TF1("ff", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", left, right);
	gr->Fit(ff, "Q");

	double A = ff->GetParameter(0);
	double B = ff->GetParameter(1);
	double C = ff->GetParameter(2);
	double D = ff->GetParameter(3);

	double dA = ff->GetParError(0);
	double dB = ff->GetParError(1);
	double dC = ff->GetParError(2);
	double dD = ff->GetParError(3);

	delete [] X; delete [] Y; delete [] dY;
	ff->Clear(); delete ff;
	gr->Clear(); delete gr;

	func =  A + B*arg + C*arg*arg + D*arg*arg*arg;
	//dfunc = sqrt(pow(dA, 2) + pow(dB/arg, 2) + pow(dC/(arg*arg), 2));

	buff = interp_cub(V, W, true);

	dfunc = buff[1];

	result.push_back(func);
	result.push_back(dfunc); buff.clear();

	return result;
}

vector<double> extrapolate_W(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;

	vector<double> buff1, buff2, buff3, buff4, buff;
	vector<vector<double>> Block;

	double W1, W2, W3, W4;

	if(W < 2.0)
	{
		if(h_L)
		{
			if(Q2 < 1.0)
			{
				W1 = 1.65; W2 = 1.725; W3 = 1.775; W4 = 1.825;
			}
			else if(Q2 >= 1.0 and Q2 < 1.8)
			{
				W1 = 1.65; W2 = 1.75; W3 = 1.85; W4 = 1.95;
			}
			else if(Q2 >= 1.8)
			{
				W1 = 1.63; W2 = 1.675; W3 = 1.725; W4 = 1.775;
			}
		}
		else
		{
			if(Q2 < 1.0)
			{
				W1 = 1.725; W2 = 1.775; W3 = 1.825; W4 = 1.875;
			}
			else if(Q2 >= 1.0 and Q2 < 1.8)
			{
				W1 = 1.75; W2 = 1.85; W3 = 1.95; W4 = 2.05;
			}
			else if(Q2 >= 1.8)
			{
				W1 = 1.695; W2 = 1.725; W3 = 1.775; W4 = 1.825;
			}

		}
	}
	else
	{
		//W1 = 2.425; W2 = 2.475; W3 = 2.525; W4 = 2.575;
		W1 = 2.025; W2 = 2.275; W3 = 2.425; W4 = 2.575;
	}

	buff1 = Str_func(W1, Q2, cos_th);
	buff2 = Str_func(W2, Q2, cos_th);
	buff3 = Str_func(W3, Q2, cos_th);
	buff4 = Str_func(W4, Q2, cos_th);

	for(int i = 0; i < 8; i += 2)
	{
		buff.push_back(W1);
		buff.push_back(buff1[i]);
		buff.push_back(buff1[i+1]);
		Block.push_back(buff); buff.clear();

		buff.push_back(W2);
		buff.push_back(buff2[i]);
		buff.push_back(buff2[i+1]);
		Block.push_back(buff); buff.clear();

		buff.push_back(W3);
		buff.push_back(buff3[i]);
		buff.push_back(buff3[i+1]);
		Block.push_back(buff); buff.clear();

		buff.push_back(W4);
		buff.push_back(buff4[i]);
		buff.push_back(buff4[i+1]);
		Block.push_back(buff); buff.clear();

		buff = Cubic_fit(Block, W);

		result.push_back(buff[0]);
		result.push_back(buff[1]);

		Block.clear();
	}

	buff1.clear(); buff2.clear(); buff3.clear(); buff4.clear(); Block.clear();
	buff.clear();

	return result;
}

vector<double> Str_func_all(const double& W, const double& Q2, const double& cos_th)
{
	vector<double> result;

	if(isW3(W)) result = Str_func(W, Q2, cos_th);
	else result = extrapolate_W(W, Q2, cos_th);

	return result;
}

vector<double> Point_diff(const double& W, const double& Q2, const double& cos_th, const double& phi)
{
	double f, df;
	vector<double> S, result;

	if(isW3(W)) S = Str_func(W, Q2, cos_th);
	else S = extrapolate_W(W, Q2, cos_th);

	f = S[0] + eps(W, Q2)*S[2] + eps(W, Q2)*S[6]*std::cos(2*phi) + sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[4]*std::cos(phi);

	if(S[0] < 0) f = eps(W, Q2)*S[2] + eps(W, Q2)*S[6]*std::cos(2*phi) + sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[4]*std::cos(phi);

	if(S[2] < 0) f = S[0] + eps(W, Q2)*S[6]*std::cos(2*phi) + sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[4]*std::cos(phi);

	if(S[0] < 0 and S[2] < 0) f = eps(W, Q2)*S[6]*std::cos(2*phi) + sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[4]*std::cos(phi);

	if(f < 0) f = 0;

	df = sqrt(S[1]*S[1] + pow(eps(W, Q2)*S[3], 2) + pow(eps(W, Q2)*S[7]*std::cos(2*phi), 2) + pow(sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[5]*std::cos(phi), 2));

	S.clear();

	result.push_back(f);
	result.push_back(df);

	return result;
}

vector<double> Point_diff_phi(const double& W, const double& Q2, const double& cos_th, const double& phi)
{
	double f, df;
	vector<double> S, result;

	if(isW3(W)) S = Str_func(W, Q2, cos_th);
	else S = extrapolate_W(W, Q2, cos_th);

	f = S[0] + eps(W, Q2)*S[2];

	if(f < 0) f = 0;

	df = sqrt(S[1]*S[1] + pow(eps(W, Q2)*S[3], 2));

	S.clear();

	result.push_back(f);
	result.push_back(df);

	return result;
}

double error_handler(vector<double>& V, const double& average) //stat. error
{
	double result(0);

	if(V.size() == 1) return 0;

	for(auto i:V)
	{
		result += pow(i - average, 2);
	}

	return sqrt(result/(V.size()*(V.size() - 1)));
	//return sqrt(result/(V.size() - 1));
}

vector<double> Average_CS_stat()
{
	auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	vector<double> result;

	int volume_W = int(20*(W_max_ - W_min_));
	int volume_Q2 = int(10*(Q2_max_ - Q2_min_));
	int volume_cos = int(10*(cos_max - cos_min));
 	int volume_phi = int(0.1*(phi_max - phi_min)*180/M_PI);

	double W, Q2, cos_th, phi;
	unsigned long long int j = 1; double cs(0), d_cs(0);
	vector<double> f, holder;

	unsigned long long int barier(1);

	if(W_max_ == W_min_) volume_W = 1;
	if(Q2_max_ == Q2_min_) volume_Q2 = 1;
	if(cos_max == cos_min) volume_cos = 1;
	if(phi_max == phi_min) volume_phi = 1;

	barier = volume_W*volume_Q2*volume_cos*volume_phi;

	if(barier == 0) barier++;

	while(j <= barier)
	{
		W = fRand(W_min_, W_max_);
		Q2 = fRand(Q2_min_, Q2_max_);
		cos_th = fRand(cos_min, cos_max);
		phi = fRand(phi_min, phi_max);

		f = Point_diff(W, Q2, cos_th, phi);

		holder.push_back(f[0]);
		cs = (cs*(j-1) + f[0])/j;
		//d_cs = sqrt(d_cs*d_cs*(j-1)*(j-1) + f[1]*f[1])/j;
		if(f[0] != 0) d_cs = (d_cs*(j-1) + f[1]/f[0])/j;

		j++; f.clear();

		finish = std::chrono::high_resolution_clock::now();
		elapsed = (barier - j)*(finish - start)/j; //totalPhysMem

		cout << std::fixed << std::setprecision(2) << "Progress: " << floor(10000*double(j)/double(barier))/100 << "%             Time remain: " << std::fixed << std::setprecision(0) <<   floor(elapsed.count()/3600) << " h ";
		cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
		cout << std::fixed << std::setprecision(1) << ceil(10*(elapsed.count() - 60*floor(elapsed.count()/60)))/10 << " s           \r" << flush;
	}
	cout << std::fixed << std::setprecision(2) << "                                                                                      ";
	//double error_stat = error_handler(holder, cs);

	result.push_back(cs);
	result.push_back(cs*d_cs);
	//result.push_back(sqrt(pow(d_cs, 2) + pow(error_stat, 2)));

	f.clear(); holder.clear();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h ";
	cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
	cout << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";

	return result;
}

vector<double> Average_CS()
{
	cout << "\nBeam energy E = " << E0 << " GeV" << endl;
	if(h_L) cout << "Channel: KL" << endl;
	else cout << "Channel: KS" << endl;
	if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Point" << endl;
		cout << "\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
		cout << "\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
		cout << "\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
		cout << "\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
		cout << "\tphi = " << phi_max*180/M_PI << " degree\n" << endl;
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] degree\n" << endl;
	}

	auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	vector<double> result;

	int volume_W = floor(20*(W_max_ - W_min_));
	int volume_Q2 = floor(10*(Q2_max_ - Q2_min_));
	int volume_cos = floor(10*(cos_max - cos_min));
 	int volume_phi = floor(0.1*(phi_max - phi_min)*180/M_PI);

	double W, Q2, cos_th, phi;
	double cs(0), d_cs(0);
	vector<double> f, holder;

	unsigned long long int barier(1);

	if(W_max_ == W_min_) volume_W = 1;
	if(Q2_max_ == Q2_min_) volume_Q2 = 1;
	if(cos_max == cos_min) volume_cos = 1;
	if(phi_max == phi_min) volume_phi = 1;

	barier = volume_W*volume_Q2*volume_cos*volume_phi;

	if(barier == 0) barier++;
	int count(1);

	for(int l = 1; l <= volume_phi; l++)
	{
		for(int k = 1; k <= volume_cos; k++)
		{
			for(int j = 1; j <= volume_Q2; j++)
			{
				for(int i = 1; i <= volume_W; i++)
				{
					W = W_min_ + (W_max_ - W_min_)*(2*i-1)/(2*volume_W);
					Q2 = Q2_min_ + (Q2_max_ - Q2_min_)*(2*j-1)/(2*volume_Q2);
					cos_th = cos_min + (cos_max - cos_min)*(2*k-1)/(2*volume_cos);
					phi = phi_min + (phi_max - phi_min)*(2*l-1)/(2*volume_phi);

					f = Point_diff(W, Q2, cos_th, phi);

					holder.push_back(f[0]);
					cs = (cs*(count-1) + f[0])/count;
					//d_cs = sqrt(d_cs*d_cs*(count-1)*(count-1) + f[1]*f[1])/count; f.clear();
					if(abs(f[0] - 0) > 1e-0) d_cs = (d_cs*(count-1) + f[1]/f[0])/count; f.clear();

					count++;

					finish = std::chrono::high_resolution_clock::now();
					elapsed = (barier - count)*(finish - start)/count; //totalPhysMem

					cout << std::fixed << std::setprecision(2) << "Progress: " << floor(10000*double(count)/double(barier))/100 << "%             Time remain: " << std::fixed << std::setprecision(0) <<   floor(elapsed.count()/3600) << " h ";
					cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
					cout << std::fixed << std::setprecision(1) << ceil(10*(elapsed.count() - 60*floor(elapsed.count()/60)))/10 << " s           \r" << flush;
				}
			}
		}
	}

	cout << std::fixed << std::setprecision(2) << "                                                                                      ";
	//double error_stat = error_handler(holder, cs);

	result.push_back(cs);
	result.push_back(cs*d_cs);
	//result.push_back(sqrt(pow(d_cs, 2) + pow(error_stat, 2)));

	f.clear(); holder.clear();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h ";
	cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
	cout << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";

	return result;
}


vector<double> Average_CS_phi()
{
	cout << "\nBeam energy E = " << E0 << " GeV" << endl;
	if(h_L) cout << "Channel: KL" << endl;
	else cout << "Channel: KS" << endl;
	if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Point" << endl;
		cout << "\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if((Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if((Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 = " << Q2_max_ << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
	} else if(!(Q2_max_ == Q2_min_) && (W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max_ << " GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos = ";
		cout << cos_max << endl;
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n";
	} else if(!(Q2_max_ == Q2_min_) && !(W_max_ == W_min_) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min_ << ", " << W_max_ << "] GeV\n\tQ2 in [" << Q2_min_ << ", " << Q2_max_ << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
	}

	auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	vector<double> result;

	int volume_W = ceil(20*(W_max_ - W_min_));
	int volume_Q2 = ceil(10*(Q2_max_ - Q2_min_));
	int volume_cos = ceil(10*(cos_max - cos_min));
 	int volume_phi = ceil(0.1*(phi_max - phi_min)*180/M_PI);

	double W, Q2, cos_th, phi;
	double cs(0), d_cs(0);
	vector<double> f, holder;

	unsigned long long int barier(1);

	if(W_max_ == W_min_) volume_W = 1;
	if(Q2_max_ == Q2_min_) volume_Q2 = 1;
	if(cos_max == cos_min) volume_cos = 1;
	volume_phi = 1;

	barier = volume_W*volume_Q2*volume_cos*volume_phi;

	if(barier == 0) barier++;
	int count(1);

	for(int k = 1; k <= volume_cos; k++)
	{
		for(int j = 1; j <= volume_Q2; j++)
		{
			for(int i = 1; i <= volume_W; i++)
			{
				W = W_min_ + (W_max_ - W_min_)*(2*i-1)/(2*volume_W);
				Q2 = Q2_min_ + (Q2_max_ - Q2_min_)*(2*j-1)/(2*volume_Q2);
				cos_th = cos_min + (cos_max - cos_min)*(2*k-1)/(2*volume_cos);
				phi = 0;

				f = Point_diff_phi(W, Q2, cos_th, phi);

				holder.push_back(f[0]);
				cs = (cs*(count-1) + f[0])/count;
				//d_cs = sqrt(d_cs*d_cs*(count-1)*(count-1) + f[1]*f[1])/count; f.clear();
				if(abs(f[0] - 0) > 1e-0) d_cs = (d_cs*(count-1) + f[1]/f[0])/count; f.clear();

				count++;

				finish = std::chrono::high_resolution_clock::now();
				elapsed = (barier - count)*(finish - start)/count; //totalPhysMem

				cout << std::fixed << std::setprecision(2) << "Progress: " << floor(10000*double(count)/double(barier))/100 << "%             Time remain: " << std::fixed << std::setprecision(0) <<   floor(elapsed.count()/3600) << " h ";
				cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
				cout << std::fixed << std::setprecision(1) << ceil(10*(elapsed.count() - 60*floor(elapsed.count()/60)))/10 << " s           \r" << flush;
			}
		}
	}


	cout << std::fixed << std::setprecision(2) << "                                                                                      ";
	//double error_stat = error_handler(holder, cs);

	result.push_back(cs);
	result.push_back(cs*d_cs);
	//result.push_back(sqrt(pow(d_cs, 2) + pow(error_stat, 2)));

	f.clear(); holder.clear();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h ";
	cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
	cout << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";

	return result;
}
