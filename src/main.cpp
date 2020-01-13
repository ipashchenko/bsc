#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"


using namespace DNest4;


int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/time_machine/bsc/reals/RA/tests/BLLAC_STOKES_0_IF_0_amp120_phase60.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}