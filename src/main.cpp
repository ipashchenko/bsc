#include <valarray>
#include <iostream>
#include "Gain.h"
#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"

using namespace DNest4;


void check_random() {
    std::valarray<double> r = make_normal_random(10);
    for (int i=0; i<r.size(); i++) {
        std::cout << r[i] << ", ";
    }
}


void check_dnestmodel() {
    DNestModel model = DNestModel();
    std::cout << "Done" << std::endl;
}


//int main() {
//    //check_random();
//    //check_dnestmodel();
//    return 0;
//}


int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/CLionProjects/ve2/uv_data.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}