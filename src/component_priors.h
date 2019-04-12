#ifndef COMPONENTS_PRIORS_H
#define COMPONENTS_PRIORS_H

#include "DNest4.h"


DNest4::Gaussian *dx_prior = new DNest4::Gaussian(0, 10);
DNest4::Gaussian *dy_prior = new DNest4::Gaussian(0, 10);
DNest4::Gaussian *logflux_prior = new DNest4::Gaussian(-2.0, 1.0);
DNest4::Gaussian *logbmaj_prior = new DNest4::Gaussian(-5.0, 2.5);

#endif