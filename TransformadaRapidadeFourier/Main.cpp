#include "FFT.h"

#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

using namespace std;


int main() {


	/* Busco la se�al f(n) con N muestras*/
	ifstream f_n_txt;
	f_n_txt.open("xn_45_muestras.txt", ios::in);
	string line;


	if (f_n_txt.is_open())
	{
		//Aca tendria que asignarle con ALLOC/MALLOC el tama�o a los vectores f(n) y F(f)
		//ese tama�o cambia con la funcion checksize q es privada de la clase FastFourierTransform
		//Pero podria hacerla publica y obtener el tama;o y hacer con eso los vectores
		//pq estoy usando la funcion resize q tarda mucho

		vector<double> f;
		vector<complex<double>> F(size(f));

		while (getline(f_n_txt, line))
		{
			//aca le cargo los blaores a el vector f(n)
			f.push_back(stof(line));
		}
	}

	/* CALCULADORA DE TRANSFORMADA DE FOURIER DE f(n)*/

	FastFourierTransform SenalOne;
	//SenalOne.FFT_1(f, F);
	SenalOne.FFT_2(f, F);

	ofstream F_f_txt;
	int k;
	F_f_txt.open("Xfourier.txt", ios::out);
	if (F_f_txt.is_open())
	{
		for (k = 0; k <= F.size(); k++)
		{
			F_f_txt << F[k] << "\n";
			cout << F[k] << "\n";
		}
	}


}