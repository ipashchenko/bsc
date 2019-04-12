#ifndef GAINS_GP_HP_PRIORS_H
#define GAINS_GP_HP_PRIORS_H

#include "DNest4.h"


DNest4::Gaussian *logamp_prior_of_amp = new DNest4::Gaussian(-3.0, 1.0);
DNest4::Gaussian *logamp_prior_of_phase = new DNest4::Gaussian(-3.0, 1.0);
DNest4::Gaussian *logscale_prior_of_amp = new DNest4::Gaussian(6.0, 1.0);
DNest4::Gaussian *logscale_prior_of_phase = new DNest4::Gaussian(5.0, 1.0);

#endif