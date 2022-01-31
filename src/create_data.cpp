#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

#include <random>
#include <chrono>


#include "../include/functions.h"


/*
compile with:

g++ test_2ndFF.cpp ../include/functions.h -std=c++17 `root-config --libs` `root-config --cflags` -lgsl -lgslcblas -o cdata

*/



#include <stdlib.h>
#include <gsl/gsl_math.h>



// Dumps data to text files for graphing in python
int main()
{ 
    int samples = int(1e4);
    double maxQ = 0.5;
    

    // LMD+V Form factor
    ofstream file1;
    file1.open("dataLMDV.txt");
    double q1;

    for(int n=0; n<samples; n++)
    {
        q1 = maxQ/samples * n;

        double val = form_factor(pow(q1,2), pow(q1,2));
        // double val = form_factor(pow(q1,2), 0);

        file1 << q1 << "," << val;
        file1 << std::endl;
    }
    file1.close();


    // TFF Form factor Q4 Expansion
    ofstream file2;
    file2.open("dataQ4.txt");
    q1 = 0;

    for(int n=0; n<samples; n++)
    {
        q1 = maxQ/samples * n;

        double val = form_factor_q4(pow(q1,2), pow(q1,2));
        // double val = form_factor2_q4(pow(q1,2), 0);

        file2 << q1 << "," << val;
        file2 << std::endl;
    }
    file2.close();


    // TFF Form factor Q6 Expansion
    ofstream file3;
    file3.open("dataQ6.txt");
    q1 = 0;

    for(int n=0; n<samples; n++)
    {
        q1 = maxQ/samples * n;

        double val = form_factor_q6(pow(q1,2), pow(q1,2));
        // double val = form_factor2_q6(pow(q1,2), 0);

        file3 << q1 << "," << val;
        file3 << std::endl;
    }
    file3.close();




    return 0;
}







