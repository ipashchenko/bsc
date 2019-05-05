#include <valarray>
#include <set>
#include <iostream>
#include <fstream>
#include <string.h>
#include "Gain.h"
#include "Gains.h"
#include "Data.h"
#include "DNestModel.h"
#include "SkyModel.h"
#include "DNest4.h"


void check_data();
using namespace DNest4;


//void check_random() {
//    std::valarray<double> r = make_normal_random(10);
//    for (int i=0; i<r.size(); i++) {
//        std::cout << r[i] << ", ";
//    }
//}


//void check_dnestmodel() {
//    Data::get_instance().load("/home/ilya/github/bsc/0716_gains.txt");
//    DNestModel model = DNestModel();
//    std::cout << "Original model : ";
//    model.print(std::cout);
//    DNestModel newmodel = model;
//    std::cout << "Copied model : ";
//    newmodel.print(std::cout);
//    double new_x = 100;
//    model.set_x_skymodel(new_x);
//    std::cout << "Original model after = : ";
//    model.print(std::cout);
//    std::cout << "Copied model after = : ";
//    newmodel.print(std::cout);
//    std::cout << "Done" << std::endl;
//}


//void check_skymodel() {
//    SkyModel sky_model = SkyModel();
//    int ncomp = 1;
//    for (int i=0; i<ncomp; i++) {
//        auto* comp = new CGComponent();
//        sky_model.add_component(comp);
//    }
//    sky_model.print(std::cout);
//    sky_model.set_x(1);
//    SkyModel new_sky_model = sky_model;
//    new_sky_model.print(std::cout);
//    sky_model.set_x(100);
//    sky_model.print(std::cout);
//    new_sky_model.print(std::cout);
//}


//void check_gain() {
//    // Gain has no pointers thus default copy ctor works fine
//    Data::get_instance().load("/home/ilya/github/bsc/0716_gains.txt");
//    Data data = Data::get_instance();
//    Gain gain = Gain(data.get_times_amp()[3], data.get_times_phase()[3]);
//    auto descr = gain.description();
//    std::cout << descr << std::endl;
//    gain.print_hp(std::cout);
//    Gain newgain = gain;
//    newgain.print_hp(std::cout);
//    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 2);
//    gain.set_hp_amp(hp_amplitudes);
//    gain.print_hp(std::cout);
//    newgain.print_hp(std::cout);
//}


//void check_gain_ptr() {
//    // Gain has no pointers thus default copy ctor works fine
//    Data::get_instance().load("/home/ilya/github/bsc/0716_gains.txt");
//    Data data = Data::get_instance();
//    Gain* gain = new Gain(data.get_times_amp()[3], data.get_times_phase()[3]);
//    gain->print_hp(std::cout);
//    Gain* newgain = new Gain(*gain);
//    newgain->print_hp(std::cout);
//    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 2);
//    gain->set_hp_amp(hp_amplitudes);
//    gain->print_hp(std::cout);
//    newgain->print_hp(std::cout);
//}


//void check_gains() {
//    Data::get_instance().load("/home/ilya/github/bsc/0716_gains.txt");
//    Gains* gains = new Gains(Data::get_instance());
//    auto descr = gains->description();
//    std::cout << descr << std::endl;
//    std::cout << "Gains before changing :" << std::endl;
//    gains->print_hp(std::cout);
//    auto* newgains = new Gains(*gains);
//    //*newgains = *gains;
//    std::cout << "Copied Gains before changing :" << std::endl;
//    newgains->print_hp(std::cout);
//    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 20);
//    gains->set_hp_amp(hp_amplitudes);
//    std::cout << "Gains after changing :" << std::endl;
//    gains->print_hp(std::cout);
//    std::cout << "Copied Gains after changing :" << std::endl;
//    newgains->print_hp(std::cout);
//
//}



//void check_description() {
//    //SkyModel sky_model = SkyModel();
//    //int ncomp = 1;
//    //for (int i=0; i<ncomp; i++) {
//    //    auto* comp = new CGComponent();
//    //    sky_model.add_component(comp);
//    //}
//    //std::string descr_skymodel = sky_model.description();
//    //std::cout << descr_skymodel << std::endl;
//    //
//    //Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    //Data data = Data::get_instance();
//    //Gain gain = Gain(data.get_times_amp(), data.get_times_phase());
//    //std::string descr_gain = gain.description();
//    //std::cout << descr_gain << std::endl;
//    //
//    //Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    //Gains* gains = new Gains(Data::get_instance());
//    //std::string descr_gains = gains->description();
//    //std::cout << descr_gains << std::endl;
//
//    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    DNestModel model = DNestModel();
//    std::string descr_dn = model.description();
//    std::cout << descr_dn << std::endl;
//    model.print(std::cout);
//}


//void check_data() {
//    Data::get_instance().load("/home/ilya/github/bsc/0716_gains.txt");
//    Data data = Data::get_instance();
//    auto times_amp = data.get_times_amp();
//    auto times_phases = data.get_times_phase();
//
//    std::cout << "Length time amplitudes 3 = " << times_amp[3].size() << std::endl;
//    for (auto x : times_amp[3]) {
//        std::cout << x << ", ";
//    }
//
//    std::cout << std::endl;
//    std::cout << "Length time amplitudes 5 = " << times_amp[5].size() << std::endl;
//    for (auto x : times_amp[5]) {
//        std::cout << x << ", ";
//    }
//
//    std::cout << std::endl;
//    std::cout << "Length time amplitudes 7 = " << times_amp[7].size() << std::endl;
//    for (auto x : times_amp[7]) {
//        std::cout << x << ", ";
//    }
//
//    std::cout << std::endl;
//    std::cout << "Length time amplitudes 8 = " << times_amp[8].size() << std::endl;
//    for (auto x : times_amp[8]) {
//        std::cout << x << ", ";
//    }
//
//    std::cout << "Length time phases 3 = " << times_phases[3].size() << std::endl;
//    for (auto x : times_phases[3]) {
//        std::cout << x << ", ";
//    }
//    std::cout << std::endl;
//
//    std::cout << "Length time phases 9 = " << times_phases[9].size() << std::endl;
//    for (auto x : times_phases[9]) {
//        std::cout << x << ", ";
//    }
//    std::cout << std::endl;
//}


//void check_read() {
//    // Open the file
//    std::fstream fin("/home/ilya/github/bsc/uvdata_gains.txt", std::ios::in);
//
//    // Temporary variables
//    // time, u, v, re, im, sigma, times_amp,  times_phase,
//    double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
//    // i, j, idx_amp_ant_i, idx_amp_ant_j, idx_phase_ant_i, idx_phase_ant_j
//    int temp9, temp10, temp11, temp12, temp13, temp14;
//
//    // Read until end of file
//    while(fin>>temp1 && fin>>temp9 && fin>>temp10 && fin>>temp2 && fin>>temp3 && fin>>temp4 && fin>>temp5 && fin>>temp6
//        && fin>>temp7 && fin>>temp11 && fin>>temp12 && fin>>temp8 && fin>>temp13 && fin>>temp14)
//    {
//
//        std::cout << "Adding to ant=" << temp9 << " amp time=" << temp7 << std::endl;
//        std::cout << "Adding to ant=" << temp10 << " amp time=" << temp7 << std::endl;
//
//
//    }
//
//    // Close the file
//    fin.close();
//}


//int main() {
//    //check_random();
//    //check_dnestmodel();
//    //check_skymodel();
//    //check_gain();
//    //check_gain_ptr();
//    //check_gains();
//    //check_description();
//    //check_data();
//    //check_read();
//    return 0;
//}


int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/bsc/uvdata_gains.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}