#include "header.h"

using namespace std;

int main(int argc, char* argv[])
{
	input_check(argc, argv);

	/*  -----------------  You can write your code below  -----------------  */

	/*
		You can change:
		 ---> E0 			beam energy [GeV]
		 ---> Q2_min_ 		virtuality (left border) [GeV2]
		 ---> Q2_max_		virtuality (right border) [GeV2]
		 ---> W_min_ 		invariant mass (left border) [GeV]
		 ---> W_max_ 		invariant mass (right border) [GeV]
		 ---> cos_min
		 ---> cos_max
		 ---> phi_min
		 ---> phi_max
		 												*/

	/*
		vector<double> calc = Str_func(W, Q2, cos); str. functions St dSt Sl dSl Slt dSlt Stt dStt [nb/sr] in the point (W, Q2, cos)
		Example: calc[0] - St -> dSigma/dOmega_T value
				 calc[1] - dSt -> error for dSigma/dOmega_T value

		vector<double> calc = Point_diff(W, Q2, cos, phi); cross section S0 dS0 [nb/sr] in the point (W, Q2, cos)
		Example: calc[0] - S0 -> dSigma/dOmega_{gamma_virt} value
				 calc[1] - dS0 -> error for dSigma/dOmega_{gamma_virt} value

		Note: dSigma/dOmega_{gamma_virt} = St + eps*Sl + eps*Stt*cos(2*phi) + sqrt(eps*(eps + 1))*Slt*cos(phi)
				 																		*/


	vector<double> calc = Average_CS(); // Average cross section [nb/sr]

	cout << "Average cs:\n\tdSigma/dOmega = " << calc[0] << " +- " << calc[1] << " [nb/sr]\n" << endl;

	/*  ----------------   You can write your code above   ------------------  */

	Data1.clear(); Data2.clear(); Data3.clear();
	Data_Diff.clear(); Data_Sigma.clear();

	return 0;
}
