#include "header.h"

using namespace std;

int main(int argc, char* argv[])
{	
	auto start = std::chrono::high_resolution_clock::now();
	vector<vector<double>> buff;
	vector<double> f, holder;
	double W, Q2, cos, phi, increment(100), f_av, df_av, f_prev(0), df_prev(0), error_stat;
	int j(1);
		
	srand(time(0));
	
	input_check(argc, argv); 
	
	if((!(isData(W_min, Q2_min) or isData2(W_min, Q2_min) or isData3(W_min, Q2_min)) or !(isData(W_max, Q2_min) or isData2(W_max, Q2_min) or isData3(W_max, Q2_min)) or !(isData(W_max, Q2_max) or isData2(W_max, Q2_max) or isData3(W_max, Q2_max)) or !(isData(W_min, Q2_max) or isData2(W_min, Q2_max) or isData3(W_min, Q2_max))))
	{
		cout << "W or Q2 interval is out of the available range!" << endl;
		return 0;
	}
	
	if(h_L)
	{
		Reading("L1.dat", buff);
		Transfer(buff, Data); buff.clear();
		Reading("L2.dat", buff);
		Transfer(buff, Data2); buff.clear();
		Reading("L3.dat", buff);
		Transfer(buff, Data3); buff.clear();
	}
	else
	{
		Reading("L4.dat", buff);
		Transfer(buff, Data); buff.clear();
		Reading("L5.dat", buff);
		Transfer(buff, Data2); buff.clear();
		Reading("L6.dat", buff);
		Transfer(buff, Data3); buff.clear();
	}
	
	while(increment > 0.00000001)
	{
		W = fRand(W_min, W_max);
		Q2 = fRand(Q2_min, Q2_max);
		cos = fRand(cos_min, cos_max);
		phi = fRand(phi_min, phi_max);
			
		f = Point_diff(W, Q2, cos, phi);	
			
        	holder.push_back(f[0]);
		f_av = (f_prev*(j-1) + f[0])/j;
		df_av = sqrt(df_prev*df_prev*(j-1)*(j-1) + f[1]*f[1])/j;
			
		increment = abs(f_prev - f_av)/f_av;
			
		f_prev = f_av;
		df_prev = df_av; j++; 
			
		cout << "Convergence: " << floor((1 - increment)*10000000+0.1)/100000 << "%      \r" << flush;
	}	
	error_stat = error_handler(holder, f_av);	
	cout << "\ndS/dOmega_gamma_av = " << f_av << " +- " << sqrt(pow(df_av, 2) + pow(error_stat, 2)) << endl; 
	cout << "j = " << j << endl;
		
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h ";
	cout << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min ";
	cout << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";
	
	Data.clear(); Data2.clear(); Data3.clear(); f.clear(); buff.clear(); holder.clear();
	
	return 0;
}
