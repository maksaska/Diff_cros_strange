#include "header.h"

using namespace std;

double E0(2.567), Q2_min(0), Q2_max(0), W_min(0), W_max(0), cos_min(0), cos_max(0), phi_min(0), phi_max(0);
bool h_L(true);

int number(3);

vector<vector<double>> Data;
vector<vector<double>> Data2;
vector<vector<double>> Data3;

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
	const char* short_options = "e:hz:x:c:v:b:n:u:i:o:p:m:l:"; int rez; int option_index;
	
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
				Q2_min = atof(optarg);
				break;
			};
			case 'x': {
				Q2_max = atof(optarg);
				break;
			};
			case 'c': {
				W_min = atof(optarg);
				break;
			};
			case 'v': {
				W_max = atof(optarg);
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
			case '?': default: {
				cerr << "Unkhown option" << endl;
				break;
			};
		};
	};
	
	if(W_min == 0) W_min = 1.8;
	if(Q2_min == 0) Q2_min = 0.65;
	if(Q2_max < Q2_min) Q2_max = Q2_min + 0.05;
	if(W_max < W_min) W_max = W_min + 0.01;
	if(cos_max < cos_min) cos_max = cos_min + 0.01;
	if(cos_min < -1) cos_min = -1;
	if(cos_min > 1) cos_min = 1;
	if(cos_max < -1) cos_max = -1;
	if(cos_max > 1) cos_max = 1;
	if(phi_min < 0) phi_min = 0;
	if(phi_min > 2*M_PI) phi_min = 2*M_PI;
	if(phi_max < 0) phi_max = 0;
	if(phi_max > 2*M_PI) phi_max = 2*M_PI;
	
	cout << "\nBeam energy E = " << E0 << " GeV" << endl;
	if(h_L) cout << "Channel: KL" << endl;
	else cout << "Channel: KS" << endl;
	if((Q2_max == Q2_min) && (W_max == W_min) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Point" << endl;
		cout << "\tW = " << W_max << " GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " grad\n" << endl;
	} else if((Q2_max == Q2_min) && (W_max == W_min) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if((Q2_max == Q2_min) && (W_max == W_min) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n"; 
		cout << "\tphi = " << phi_max*180/M_PI << " grad\n" << endl;
	} else if((Q2_max == Q2_min) && !(W_max == W_min) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " grad\n" << endl; 
	} else if(!(Q2_max == Q2_min) && (W_max == W_min) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " grad\n" << endl; 
	} else if((Q2_max == Q2_min) && (W_max == W_min) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if((Q2_max == Q2_min) && !(W_max == W_min) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if(!(Q2_max == Q2_min) && (W_max == W_min) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if((Q2_max == Q2_min) && !(W_max == W_min) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n"; 
		cout << "\tphi = " << phi_max*180/M_PI << " grad\n" << endl;
	} else if(!(Q2_max == Q2_min) && (W_max == W_min) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n"; 
		cout << "\tphi = " << phi_max*180/M_PI << " grad\n" << endl;
	} else if(!(Q2_max == Q2_min) && !(W_max == W_min) && (cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi = " << phi_max*180/M_PI << " grad\n" << endl; 
	} else if((Q2_max == Q2_min) && !(W_max == W_min) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 = " << Q2_max << " GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if(!(Q2_max == Q2_min) && (W_max == W_min) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW = " << W_max << " GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t";
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if(!(Q2_max == Q2_min) && !(W_max == W_min) && (cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos = ";
		cout << cos_max << "\n\tphi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;	
	} else if(!(Q2_max == Q2_min) && !(W_max == W_min) && !(cos_min == cos_max) && (phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n"; 
		cout << "\tphi = " << phi_max*180/M_PI << " grad\n" << endl;
	} else if(!(Q2_max == Q2_min) && !(W_max == W_min) && !(cos_min == cos_max) && !(phi_min == phi_max))
	{
		cout << "\nDiff. cross section calculation option: Average" << endl;
		cout << "Chosen area:\n\tW in [" << W_min << ", " << W_max << "] GeV\n\tQ2 in [" << Q2_min << ", " << Q2_max << "] GeV2\n\tcos in [" << cos_min << ", " << cos_max << "]\n\t"; 
		cout << "phi in [" << phi_min*180/M_PI << ", " << phi_max*180/M_PI << "] grad\n" << endl;
	}	
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

vector<double> interp_cub(vector<vector<double>>& V, const double& cos, const bool& statement)
{
	vector<vector<double>> transit;
	vector<double> result, buff;
	double func, dfunc;
	int t = V.size() - 1;

	if(cos < V[1][0])
	{
		buff.push_back(V[0][0]); buff.push_back(V[0][1]); buff.push_back(V[0][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[1][0]); buff.push_back(V[1][1]); buff.push_back(V[1][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[2][0]); buff.push_back(V[2][1]); buff.push_back(V[2][2]);
		transit.push_back(buff); buff.clear();
	}

	if(cos >= V[t-1][0])
	{
		buff.push_back(V[t-2][0]); buff.push_back(V[t-2][1]); buff.push_back(V[t-2][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[t-1][0]); buff.push_back(V[t-1][1]); buff.push_back(V[t-1][2]);
		transit.push_back(buff); buff.clear();
		buff.push_back(V[t][0]); buff.push_back(V[t][1]); buff.push_back(V[t][2]);
		transit.push_back(buff); buff.clear();
	}
	
	if(cos >= V[1][0] and cos < V[t-1][0])
	{
		for(int i = 0; i < V.size(); i++)
		{
			if(cos >= V[i][0] and cos < V[i+1][0])
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
	
	buff = cub_interp(transit, cos);
	
	func = buff[0];
	dfunc = buff[1];
	
	if(func < 0 and statement){func = 0;}
	
	result.push_back(func);
	result.push_back(dfunc); transit.clear(); buff.clear();
	
	return result;
}

vector<double> approx_cos_leg(vector<vector<double>>& V, const double& cos, const bool& statement)
{
	vector<double> result;
	double func, dfunc;
	double arg; arg = cos;
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
	dfunc = sqrt(pow(dA, 2) + pow(arg*dB, 2) + pow(dC*0.5*(3*pow(arg,2) - 1), 2) + pow(dD*0.5*(5*pow(arg,3) - 3*arg), 2) + pow(dE*(35*pow(arg, 4) - 30*pow(arg, 2) + 3)/8, 2));
	
	if(func < 0 and statement){func = 0;}
	buff = interp_cub(V, cos, statement);
	dfunc = buff[1];
	if(cos < V[0][0])
	{
		if(number == 2)
		{
			func = buff[0];
		}
		
		if(number == 3)
		{
			dfunc = sqrt(pow((func-buff[0])/2, 2) + dfunc*dfunc);
			func = (func+buff[0])/2; 
		}				
	}
	
	
	result.push_back(func);
	result.push_back(dfunc); buff.clear();
	
	return result;
}

double eps(const double& W, const double& Q2)
{
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);	
	return  1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2)); 
}

bool isData(const double& W, const double& Q2)
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
	if((Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.175 and W >= 1.65) or (Q2 <= 2.6 and Q2 >= 1.8 and W <= 2.35 and W >= 2.175)){ return true;}
	else{ return false;}
	} else
	{
	if((Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.175 and W >= 1.695) or (Q2 <= 2.6 and Q2 >= 1.8 and W <= 2.35 and W >= 2.175)){ return true;}
	else{ return false;}
	}			
}

vector<double> giveData(const double& W, const double& Q2, const double& cos)
{
	vector<double> result;
	double Su, dSu; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;
	
	if(Q2 == 0.65)
	{
		Q2_min = 0.65;
		Q2_max = 0.65;
	} else if(Q2 == 1.0)
	{
		Q2_min = 1.0;
		Q2_max = 1.0;
	} else
	{
		Q2_min = 0.65;
		Q2_max = 1.0;
	}
	
	if(W < 1.725)
	{
		W_min = 1.65;
		W_max = 1.725;
	}
	
	if(1.725 < W and W < 1.775)
	{
		W_min = 1.725;
		W_max = 1.775;
	}
	
	if(1.775 < W and W < 1.825)
	{
		W_min = 1.775;
		W_max = 1.825;
	}
	
	if(1.825 < W and W < 1.875)
	{
		W_min = 1.825;
		W_max = 1.875;
	}
	
	if(1.875 < W and W < 1.925)
	{
		W_min = 1.875;
		W_max = 1.925;
	}
	
	if(1.925 < W and W < 1.975)
	{
		W_min = 1.925;
		W_max = 1.975;
	}
	
	if(1.975 < W)
	{
		W_min = 1.975;
		W_max = 2.05; 
	}
	
	if(W == 1.65 or W == 1.725 or W == 1.775 or W == 1.825 or W == 1.875 or W == 1.925 or W == 1.975 or W == 2.05)
	{
		W_max = W;
		W_min = W_max;
	}
	
	for(auto i:Data)
	{
		if(W_min == i[0] and Q2_min == i[1])
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
		if(W_max == i[0] and Q2_min == i[1])
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
		if(W_min == i[0] and Q2_max == i[1])
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
		if(W_max == i[0] and Q2_max == i[1])
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
	
	buff = approx_cos_leg(pack_1, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_2, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_3, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_4, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
		
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();

	buff = approx_cos_leg(pack_5, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_6, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_7, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_8, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_9, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_10, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_11, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_12, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_13, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_14, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_15, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_16, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();
	
	return result;
}

vector<double> giveData2(const double& W, const double& Q2, const double& cos)
{
	vector<double> result;
	double Su, dSu; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;
	
	if(W < 1.75)
	{
		W_min = 1.65;
		W_max = 1.75;
	}
	
	if(1.75 < W and W < 1.85)
	{
		W_min = 1.75;
		W_max = 1.85;
	}
	
	if(1.85 < W and W < 1.95)
	{
		W_min = 1.85;
		W_max = 1.95;
	}
	
	if(1.95 < W and W < 2.05)
	{
		W_min = 1.95;
		W_max = 2.05;
	}
	
	if(2.05 < W and W < 2.15)
	{
		W_min = 2.05;
		W_max = 2.15;
	}
	
	if(2.15 < W and W < 2.25)
	{
		W_min = 2.15;
		W_max = 2.25;
	}
	
	if(2.25 < W and W < 2.35)
	{
		W_min = 2.25;
		W_max = 2.35;
	}
	
	if(W == 1.65 or W == 1.75 or W == 1.85 or W == 1.95 or W == 2.05 or W == 2.15 or W == 2.25 or W == 2.35)
	{
		W_max = W;
		W_min = W_max;
	}
	
	if(Q2 < 1.55)
	{
		Q2_min = 1.0;
		Q2_max = 1.55;
	}
	
	if(1.55 < Q2 and Q2 < 2.05)
	{
		Q2_min = 1.55;
		Q2_max = 2.05;
	}
	
	if(2.05 < Q2 and Q2 < 2.55)
	{
		Q2_min = 2.05;
		Q2_max = 2.55;
	}
	
	if(Q2 == 1.0 or Q2 == 1.55 or Q2 == 2.05 or Q2 == 2.55)
	{
		Q2_max = Q2;
		Q2_min = Q2_max;
	}
	
	for(auto i:Data2)
	{
		if(W_min == i[0] and Q2_min == i[1])
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
		if(W_max == i[0] and Q2_min == i[1])
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
		if(W_min == i[0] and Q2_max == i[1])
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
		if(W_max == i[0] and Q2_max == i[1])
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
		
	buff = approx_cos_leg(pack_1, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_2, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_3, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_4, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_5, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_6, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();

	buff = approx_cos_leg(pack_7, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_8, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_9, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_10, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_11, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_12, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_13, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_14, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_15, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_16, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();
	
	return result;
}

vector<double> giveData3(const double& W, const double& Q2, const double& cos)
{
	vector<double> result;
	double Su, dSu; double W_min, W_max, Q2_min, Q2_max;
	double y1, y2, y3, y4, dy1, dy2, dy3, dy4, S1, S2, dS1, dS2;
	vector<vector<double>> pack_1, pack_2, pack_3, pack_4; //St
	vector<vector<double>> pack_5, pack_6, pack_7, pack_8; //Sl
	vector<vector<double>> pack_9, pack_10, pack_11, pack_12; //Slt
	vector<vector<double>> pack_13, pack_14, pack_15, pack_16; //Slt
	vector<double> buff;
	
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
	
	if(1.675 <= W and h_L)
	{
		W_min = floor(W*40); 
		if(int(W_min) % 2 == 0){W_min--;}
		W_min = W_min/40;
		
		W_max = ceil(W*40); 
		if(int(W_max) % 2 == 0){W_max++;}
		W_max = W_max/40;
	}
	
	if(1.725 <= W)
	{
		W_min = floor(W*40); 
		if(int(W_min) % 2 == 0){W_min--;}
		W_min = W_min/40;
		
		W_max = ceil(W*40); 
		if(int(W_max) % 2 == 0){W_max++;}
		W_max = W_max/40;
	}
	
	if(Q2 < 2.6)
	{
		Q2_min = 1.8;
		Q2_max = 2.6;
	}
	
	if(2.6 < Q2)
	{
		Q2_min = 2.6;
		Q2_max = 3.45;
	}
	
	if(Q2 == 1.8 or Q2 == 2.6 or Q2 == 3.45)
	{
		Q2_max = Q2;
		Q2_min = Q2_max;
	}
	
	for(auto i:Data3)
	{
		if(W_min == i[0] and Q2_min == i[1])
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
		if(W_max == i[0] and Q2_min == i[1])
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
		if(W_min == i[0] and Q2_max == i[1])
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
		if(W_max == i[0] and Q2_max == i[1])
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
	
	buff = approx_cos_leg(pack_1, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_2, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_3, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_4, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_1.clear(); pack_2.clear(); pack_3.clear(); pack_4.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_5, cos, true);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_6, cos, true);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_7, cos, true);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_8, cos, true);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_5.clear(); pack_6.clear(); pack_7.clear(); pack_8.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_9, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_10, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_11, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_12, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_9.clear(); pack_10.clear(); pack_11.clear(); pack_12.clear(); buff.clear();
	
	buff = approx_cos_leg(pack_13, cos, false);
	y1 = buff[0]; dy1 = buff[1]; buff.clear(); 
		
	buff = approx_cos_leg(pack_14, cos, false);
	y2 = buff[0]; dy2 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_15, cos, false);
	y3 = buff[0]; dy3 = buff[1]; buff.clear();
		
	buff = approx_cos_leg(pack_16, cos, false);
	y4 = buff[0]; dy4 = buff[1]; buff.clear();
	
	if(W_min == W_max)
	{
		if(Q2_min == Q2_max)
		{
			Su = y1;
			dSu = dy1;
		} else
		{
			Su = (y1 - y3)*Q2/(Q2_min - Q2_max) + (y3*Q2_min - y1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dy1, 2) + pow((Q2_min - Q2)*dy3, 2))/0.5;		
		}
	} else
	{
		S1 = (y1 - y2)*W/(W_min - W_max) + (y2*W_min - y1*W_max)/(W_min - W_max);
		dS1 = sqrt(pow((W - W_max)*dy1, 2) + pow((W_min - W)*dy2, 2))/abs(W_min - W_max);
		
		S2 = (y3 - y4)*W/(W_min - W_max) + (y4*W_min - y3*W_max)/(W_min - W_max);
		dS2 = sqrt(pow((W - W_max)*dy3, 2) + pow((W_min - W)*dy4, 2))/abs(W_min - W_max);
		
		if(Q2_min == Q2_max)
		{
			Su = S1;
			dSu = dS1;
		} else
		{
			Su = (S1 - S2)*Q2/(Q2_min - Q2_max) + (S2*Q2_min - S1*Q2_max)/(Q2_min - Q2_max);
			dSu = sqrt(pow((Q2 - Q2_max)*dS1, 2) + pow((Q2_min - Q2)*dS2, 2))/0.5;		
		} 	
	}	
	
	result.push_back(Su);
	result.push_back(dSu);
	
	pack_13.clear(); pack_14.clear(); pack_15.clear(); pack_16.clear(); buff.clear();
	
	return result;
}

vector<double> Str_func(const double& W, const double& Q2, const double& cos)
{
	vector<double> result, buff;
	double St, dSt, Sl, dSl, Stt, dStt, Slt, dSlt;
	
	if(isData(W, Q2))
	{
		buff = giveData(W, Q2, cos);
		
		St = buff[0];
		dSt = buff[1];
		Sl = buff[2];
		dSl = buff[3];	
		Slt = buff[4];
		dSlt = buff[5];
		Stt = buff[6];
		dStt = buff[7]; buff.clear();
		
		/*if(isData2(W, Q2))
		{
			buff = giveData2(W, Q2, cos);
			
			St = (St + buff[0])/2;
			dSt = sqrt(dSt*dSt + buff[1]*buff[1])/2;
			Sl = (Sl + buff[2])/2;
			dSl = sqrt(dSl*dSl + buff[3]*buff[3])/2;
			Slt = (Slt + buff[4])/2;
			dSlt = sqrt(dSlt*dSlt + buff[5]*buff[5])/2;
			Stt = (Stt + buff[6])/2;
			dStt = sqrt(dStt*dStt + buff[6]*buff[6])/2; buff.clear();
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
	
	if(isData2(W, Q2))
	{
		buff = giveData2(W, Q2, cos);
		
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
			buff = giveData3(W, Q2, cos);
			
			St = (St + buff[0])/2;
			dSt = sqrt(dSt*dSt + buff[1]*buff[1])/2;
			Sl = (Sl + buff[2])/2;
			dSl = sqrt(dSl*dSl + buff[3]*buff[3])/2;
			Slt = (Slt + buff[4])/2;
			dSlt = sqrt(dSlt*dSlt + buff[5]*buff[5])/2;
			Stt = (Stt + buff[6])/2;
			dStt = sqrt(dStt*dStt + buff[6]*buff[6])/2; buff.clear();
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
		buff = giveData3(W, Q2, cos);
		
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

vector<double> Point_diff(const double& W, const double& Q2, const double& cos, const double& phi)
{
	double f, df;
	vector<double> S, result;
	
	S = Str_func(W, Q2, cos);
	
	f = S[0] + eps(W, Q2)*S[2] + eps(W, Q2)*S[6]*std::cos(2*phi) + sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[4]*std::cos(phi);
	if(f < 0) f = 0;
	df = sqrt(S[1]*S[1] + pow(eps(W, Q2)*S[3], 2) + pow(eps(W, Q2)*S[7]*std::cos(2*phi), 2) + pow(sqrt(eps(W, Q2)*(eps(W, Q2) + 1))*S[5]*std::cos(phi), 2));

	S.clear();
	
	result.push_back(f);
	result.push_back(df);
		
	return result;
}

double error_handler(vector<double>& V, const double& average)
{
	double result(0);
	
	for(auto i:V)
	{
		result += pow(i - average, 2);
	}
	
	return sqrt(result/(V.size()*(V.size() - 1)));
}
