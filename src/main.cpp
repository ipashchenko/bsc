#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"


using namespace DNest4;


int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/time_machine/bsc/data_only.txt");
    //auto times = Data::get_instance().get_times_amp()[2];
    //for (auto t: times) {
    //    std::cout << t << std::endl;
    //}
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}