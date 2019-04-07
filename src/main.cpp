#include <valarray>
#include <iostream>
#include "Gain.h"


void check_random() {
    std::valarray<double> r = make_normal_random(10);
    for (int i=0; i<r.size(); i++) {
        std::cout << r[i] << ", ";
    }
}

int main() {
    check_random();
    return 0;
}
