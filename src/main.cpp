#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"


using namespace DNest4;


int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/time_machine/bsc/2038X_unsc_IF0_RR.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}