#include <valarray>
#include <iostream>
#include "Gain.h"
#include "DNestModel.h"


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


int main() {
    //check_random();
    check_dnestmodel();
    return 0;
}
