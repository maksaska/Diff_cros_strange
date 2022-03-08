#include "header.h"

using namespace std;

int main(int argc, char* argv[])
{
	input_check(argc, argv);

	/*  -----------------  You can write your code below  -----------------  */

	/*
		YOU CAN CHANGE:

		 ---> E0 			beam energy [GeV]
		 ---> Q2_min_ 		virtuality (left border) [GeV2]
		 ---> Q2_max_		virtuality (right border) [GeV2]
		 ---> W_min_ 		invariant mass (left border) [GeV]
		 ---> W_max_ 		invariant mass (right border) [GeV]
		 ---> cos_min      	theta in c.m. frame
		 ---> cos_max
		 ---> phi_min		phi in c.m. frame [rad]
		 ---> phi_max

		YOU CAN USE:

		(1) vector<double> calc = Str_func_all(W, Q2, cos); str. functions St dSt Sl dSl Slt dSlt Stt dStt [nb/sr] in the point (W, Q2, cos)

		Example: calc[0] - St -> dSigma/dOmega_T value
				 calc[1] - dSt -> error for dSigma/dOmega_T value

		(2) vector<double> calc = Point_diff(W, Q2, cos, phi); cross section S0 dS0 [nb/sr] in the point (W, Q2, cos)

		Example: calc[0] - S0 -> dSigma/dOmega_{gamma_virt} value
				 calc[1] - dS0 -> error for dSigma/dOmega_{gamma_virt} value

		(3) vector<double> calc = Average_CS();  Average cross main method [nb/sr]
		(4) vector<double> calc = Average_CS_stat();  Average cross section stat. method [nb/sr]
		(5) vector<double> calc = Average_CS_phi();  Average cross section [nb/sr] with phi in [-180, 180] degree (faster than Average_CS function)

		Note: dSigma/dOmega_{gamma_virt} = St + eps*Sl + eps*Stt*cos(2*phi) + sqrt(eps*(eps + 1))*Slt*cos(phi)
				 																		*/

	vector<double> calc = Average_CS(); //  final choice for average proceedure

	cout << "Average cs:\n\tdSigma/dOmega = " << calc[0] << " +- " << calc[1] << " [nb/sr]\n" << endl;

	/*  ----------------   You can write your code above   ------------------  */

	Data1.clear(); Data2.clear(); Data3.clear();
	Data_Diff.clear(); Data_Sigma.clear(); Data_Q2cos.clear();
	cosQ2_grid.clear();

	return 0;
}
