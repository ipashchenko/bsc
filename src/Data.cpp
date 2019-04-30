#include "Data.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

using namespace std;

// The static instance
Data Data::instance;

Data::Data() = default;


void Data::load(const char* filename)
{
    std::vector<int> _antennas;
    // Vectors to hold the data
    std::vector<double> _u;
    std::vector<double> _v;
    std::vector<double> _vis_real;
    std::vector<double> _vis_imag;
    std::vector<double> _sigma;

    // Open the file
    fstream fin(filename, ios::in);

    // Temporary variables
    // time, u, v, re, im, sigma, times_amp,  times_phase,
    double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
    // i, j, idx_amp_ant_i, idx_amp_ant_j, idx_phase_ant_i, idx_phase_ant_j
    int temp9, temp10, temp11, temp12, temp13, temp14;

    // Read until end of file
    while(fin>>temp1 && fin>>temp9 && fin>>temp10 && fin>>temp2 && fin>>temp3 && fin>>temp4 && fin>>temp5 && fin>>temp6
    && fin>>temp7 && fin>>temp11 && fin>>temp12 && fin>>temp8 && fin>>temp13 && fin>>temp14)
    {
        times.push_back(temp1);
        ant_i.push_back(temp9);
        ant_j.push_back(temp10);
        _u.push_back(temp2);
        _v.push_back(temp3);
        _vis_real.push_back(temp4);
        _vis_imag.push_back(temp5);
        _sigma.push_back(temp6);
        times_amp.push_back(temp7);
        idx_amp_ant_i.push_back(temp11);
        idx_amp_ant_j.push_back(temp12);
        times_phase.push_back(temp8);
        idx_phase_ant_i.push_back(temp13);
        idx_phase_ant_j.push_back(temp14);

    }

    // Close the file
    fin.close();

    // Vector of antenna numbers
    antennas.insert(antennas.end(), ant_i.begin(), ant_i.end());
    antennas.insert(antennas.end(), ant_j.begin(), ant_j.end());
    std::set<int> s( antennas.begin(), antennas.end() );
    antennas.assign( s.begin(), s.end() );

    // Generate the map between ant_i, ant_j and their position in antennas vector
    for (int i=0; i<antennas.size();i++) {
        antennas_map[antennas[i]] = i;
    }
    for (auto x : antennas_map) {
        std::cout << x.first << " -- " << x.second << std::endl;
    }

        // Copy the data to the valarrays
    u = valarray<double>(&_u[0], _u.size());
    v = valarray<double>(&_v[0], _v.size());
    vis_real = valarray<double>(&_vis_real[0], _vis_real.size());
    vis_imag = valarray<double>(&_vis_imag[0], _vis_imag.size());
    sigma = valarray<double>(&_sigma[0], _sigma.size());

}


