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


void CGComponent::from_prior(DNest4::RNG &rng) {
     // Normal diffuse prior for x & y
     dx_ = 5*rng.randn();
     dy_ = 5*rng.randn();
     // Log-normal prior for flux and bmaj
     flux_ = exp(-7. + 2.5*rng.randn());
     bmaj_ = exp(-7. + 2.5*rng.randn());
}


double CGComponent::perturb(DNest4::RNG &rng) {
    double log_H = 0.;
    // Proposals explore the prior
    // For normal priors I usually use the hastings factor to do it
    int which = rng.rand_int(4);
    if(which%4 == 0)
    {
        log_H -= -0.5*pow(dx_/5, 2);
        dx_ += 5*rng.randh();
        log_H += -0.5*pow(dx_/5, 2);

        // Example of the proposal that keeps components sorted
        //double old = params[which];
        //log_H -= -0.5*pow(params[which]/5, 2);
        //
        //do {
        //    double perturbation = 5*rng.randh();
        //    params[which] = old + perturbation;
        //    old = params[which] - perturbation;
        //} while (!are_sorted());
        //log_H += -0.5*pow(params[which]/5, 2);
    }
    else if(which%4 == 1)
    {
        log_H -= -0.5*pow(dy_/5, 2);
        dy_ += 5*rng.randh();
        log_H += -0.5*pow(dy_/5, 2);
    }
    else if(which%4 == 2)
    {
        //// Usual log-uniform prior trick
        //params[which] = log(params[which]);
        //params[which] += 8.*rng.randh();
        //wrap(params[which], -7., 1.);
        //params[which] = exp(params[which]);

        double logflux = log(flux_);
        log_H -= -0.5*pow((logflux+7.0)/2.5, 2);
        logflux += 2.5*rng.randh();
        log_H += -0.5*pow((logflux+7.0)/2.5, 2);
        flux_ = exp(logflux);
    }
    else
    {
        //// Usual log-uniform prior trick
        //params[which] = log(params[which]);
        //params[which] += 12.*rng.randh();
        //wrap(params[which], -10., 2.);
        //params[which] = exp(params[which]);

        double logbmaj = log(bmaj_);
        log_H -= -0.5*pow((logbmaj+7.0)/2.5, 2);
        bmaj_ += 2.5*rng.randh();
        log_H += -0.5*pow((logbmaj+7.0)/2.5, 2);
        bmaj_ = exp(logbmaj);
    }

    //// This naively skip unsorted perturbations
    //if (!are_sorted())
    //    return -1E300;
    //    // Pre-reject
    //else if(rng.rand() >= exp(log_H))
    //    return -1E300;
    //
    //else
    //    log_H = 0.0;

    // Without sorting
    // Pre-reject
     if(rng.rand() >= exp(log_H))
          return -1E300;
     else
          log_H = 0.0;


    //// Calculate model visibilities again if Gaussian parameters changed
    //// Here we can save calculation
    //if(which%4 == 0 || which%4 == 1 || which%4 == 2 || which%4 == 3)
    //    calculate_mu();

    return log_H;
}
