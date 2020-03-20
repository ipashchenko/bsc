#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"


using namespace DNest4;


int main(int argc, char** argv)
{
    //Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/RA/tests/BLLAC_STOKES_0_IF_0_amp120_phase60.txt");
    //Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/1502/1502_STOKES_0_IF_0_amp60_phase30.txt");
    Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/RA/tests/BLLAC_STOKES_0_IF_0_amp120_phase60_10ants_scan.txt");
    //Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/bllac/1IF/bllac_2015_02_20_STOKES_0_IF_0_aver30s_amp30_phase30.txt");
    //Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/J2038/2038S_STOKES_RR_IF_0.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}