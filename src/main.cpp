#include <valarray>
#include <iostream>
#include "Gain.h"
#include "Gains.h"
#include "Data.h"
#include "DNestModel.h"
#include "SkyModel.h"
#include "DNest4.h"

using namespace DNest4;


void check_random() {
    std::valarray<double> r = make_normal_random(10);
    for (int i=0; i<r.size(); i++) {
        std::cout << r[i] << ", ";
    }
}


//void check_dnestmodel() {
//    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    DNestModel model = DNestModel();
//    std::cout << "Original model : ";
//    model.print(std::cout);
//    DNestModel newmodel = model;
//    std::cout << "New model : ";
//    newmodel.print(std::cout);
//    double new_x = 100;
//    model.set_x_skymodel(new_x);
//    std::cout << "Original model after = : ";
//    model.print(std::cout);
//    std::cout << "New model : ";
//    newmodel.print(std::cout);
//    std::cout << "Done" << std::endl;
//}


void check_skymodel() {
    SkyModel sky_model = SkyModel();
    int ncomp = 1;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model.add_component(comp);
    }
    sky_model.print(std::cout);
    sky_model.set_x(1);
    SkyModel new_sky_model = sky_model;
    new_sky_model.print(std::cout);
    sky_model.set_x(100);
    sky_model.print(std::cout);
    new_sky_model.print(std::cout);
}


//void check_gain() {
//    // Gain has no pointers thus default copy ctor works fine
//    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    Data data = Data::get_instance();
//    Gain gain = Gain(data.get_times_amp(), data.get_times_phase());
//    gain.print_hp(std::cout);
//    Gain newgain = gain;
//    newgain.print_hp(std::cout);
//    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 2);
//    gain.set_hp_amp(hp_amplitudes);
//    gain.print_hp(std::cout);
//    newgain.print_hp(std::cout);
//}


void check_gain_ptr() {
    // Gain has no pointers thus default copy ctor works fine
    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
    Data data = Data::get_instance();
    Gain* gain = new Gain(data.get_times_amp(), data.get_times_phase());
    gain->print_hp(std::cout);
    Gain* newgain = new Gain(*gain);
    newgain->print_hp(std::cout);
    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 2);
    gain->set_hp_amp(hp_amplitudes);
    gain->print_hp(std::cout);
    newgain->print_hp(std::cout);
}


void check_gains() {
    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
    Gains* gains = new Gains(Data::get_instance());
    std::cout << "Gains before changing :" << std::endl;
    gains->print_hp(std::cout);
    auto* newgains = new Gains(*gains);
    //*newgains = *gains;
    std::cout << "Copied Gains before changing :" << std::endl;
    newgains->print_hp(std::cout);
    std::valarray<double> hp_amplitudes = std::valarray<double>(100.0, 20);
    gains->set_hp_amp(hp_amplitudes);
    std::cout << "Gains after changing :" << std::endl;
    gains->print_hp(std::cout);
    std::cout << "Copied Gains after changing :" << std::endl;
    newgains->print_hp(std::cout);

}

int main() {
    //check_random();
    //check_dnestmodel();
    //check_skymodel();
    //check_gain();
    //check_gain_ptr();
    check_gains();
    return 0;
}


//int main(int argc, char** argv)
//{
//    Data::get_instance().load("/home/ilya/github/bsc/uv_data_gains.txt");
//    // set the sampler and run it!
//    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
//    sampler.run();
//
//    return 0;
//}