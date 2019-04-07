#include "Component.h"
#include <cmath>
#include <Component.h>
#include <iostream>



DFComponent::DFComponent() : dx_(0.0), dy_(0.0), flux_(0.0)
{

}

void DFComponent::ft(std::valarray<double> u, std::valarray<double> v)
{
    std::valarray<double> theta;
    // Phase shift due to not being in a phase center
    theta = 2*PI*mas_to_rad*(u*dx_+v*dy_);

    // Prediction of visibilities
    mu_real = flux_*cos(theta);
    mu_imag = flux_*sin(theta);
}

void DFComponent::set_param_vec(std::valarray<double> param)
{
    dx_ = param[0];
    dy_ = param[1];
    flux_ = param[2];
}

void DFComponent::print(std::ostream &out) const
{
    out << "x=" << dx_ << ", " << "y=" << dy_ << ", " << "flux=" << flux_ << "\n";

}

EGComponent::EGComponent() : dx_(0.0), dy_(0.0), flux_(0.0), bmaj_(0.0), e_(1.0), bpa_(0.0)
{

}

void EGComponent::ft(std::valarray<double> u, std::valarray<double> v)
{
    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Phase shift due to not being in a phase center
    theta = 2*PI*mas_to_rad*(u*dx_+v*dy_);
    // Calculate FT of a Gaussian in a phase center
    c = pow(PI*bmaj_*mas_to_rad, 2)/(4.*log(2.));
    //b = e_*e_*u*u + v*v;
    b = e_*e_*pow((u*cos(bpa_) - v*sin(bpa_)), 2) + pow((u*sin(bpa_)+v*cos(bpa_)), 2);
    ft = flux_*exp(-c*b);

    // Prediction of visibilities
    mu_real = ft*cos(theta);
    mu_imag = ft*sin(theta);
}

void EGComponent::set_param_vec(std::valarray<double> param)
{
    dx_ = param[0];
    dy_ = param[1];
    flux_ = param[2];
    bmaj_ = param[3];
    e_ = param[4];
    bpa_ = param[5];
}

void EGComponent::print(std::ostream &out) const
{
    out << "x=" << dx_ << ", " << "y=" << dy_ << ", " << "flux=" << flux_ << ", "<< "size=" << bmaj_ << ", "<< "e=" << e_ << ", "<< "bpa=" << bpa_ << "\n";

}


CGComponent::CGComponent() : EGComponent() {
    dx_ = 0.0;
    dy_ = 0.0;
    flux_ = 0.0;
    bmaj_ = 0.0;
}

void CGComponent::set_param_vec(std::valarray<double> param) {
    dx_ = param[0];
    dy_ = param[1];
    flux_ = param[2];
    bmaj_ = param[3];
}

void CGComponent::print(std::ostream &out) const
{
    out << "x=" << dx_ << ", " << "y=" << dy_ << ", " << "flux=" << flux_ << ", "<< "size=" << bmaj_ << "\n";
}