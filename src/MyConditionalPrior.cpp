#include "MyConditionalPrior.h"


using namespace DNest4;

MyConditionalPrior::MyConditionalPrior(double x_min, double x_max,
                                       double y_min, double y_max)
    :x_min(x_min)
    ,x_max(x_max)
    ,y_min(y_min)
    ,y_max(y_max)
{

}

void MyConditionalPrior::from_prior(RNG& rng)
{
    // A Cauchy distribution
    const DNest4::Cauchy cauchy(0.0, 1.0);

    // Truncate to (-6, 6)
    do
    {
        typical_flux = cauchy.generate(rng);
    }while(std::abs(typical_flux) >= 6.0);
    // Comment out because we use log of flux
    //typical_flux = exp(typical_flux);

    dev_log_flux = 3.0*rng.rand();

    // Comment out because we use log of size
    //typical_radius = exp(2.0*rng.randn());
    typical_radius = 2.0*rng.randn();
    dev_log_radius = 3.0*rng.rand();
}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(4);

    if(which == 0)
    {
        // A Cauchy distribution
        const DNest4::Cauchy cauchy(0.0, 1.0);

        //typical_flux = log(typical_flux);
        logH += cauchy.perturb(typical_flux, rng);
        if(std::abs(typical_flux) >= 6.0)
        {
            typical_flux = 1.0;
            return -1E300;
        }
        //typical_flux = exp(typical_flux);
    }
    else if(which == 1)
    {
        dev_log_flux += 3.0*rng.randh();
        DNest4::wrap(dev_log_flux, 0.0, 3.0);
    }
    else if(which == 2)
    {
        //typical_radius = log(typical_radius);
        logH -= -0.5*pow(typical_radius/2.0, 2);
        typical_radius += 2.0*rng.randh();
        logH += -0.5*pow(typical_radius/2.0, 2);
        //typical_radius = exp(typical_radius);
    }
    else if(which == 3)
    {
        dev_log_radius += 3.0*rng.randh();
        DNest4::wrap(dev_log_radius, 0.0, 3.0);
    }

    return logH;
}

// x, y, flux, radius
double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    double logp = 0.0;

    // Position
    if(vec[0] < x_min || vec[0] > x_max || vec[1] < y_min || vec[1] > y_max)
        return -1E300;

    // Flux
    //DNest4::Laplace laplace1(log(typical_flux), dev_log_flux);
    DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    //logp += -log(vec[2]) + laplace1.log_pdf(log(vec[2]));
    logp += -vec[2] + laplace1.log_pdf(vec[2]);

    // Radius
    //DNest4::Laplace laplace2(log(typical_radius), dev_log_radius);
    DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    //logp += -log(vec[3]) + laplace2.log_pdf(log(vec[3]));
    logp += -vec[3] + laplace2.log_pdf(vec[3]);

    return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    // Position
    vec[0] = x_min + (x_max - x_min)*vec[0];
    vec[1] = y_min + (y_max - y_min)*vec[1];

    // Flux
    DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    vec[2] = laplace1.cdf_inverse(vec[2]);

    // Radius
    DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    vec[3] = laplace2.cdf_inverse(vec[3]);

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    // Position
    vec[0] = (vec[0] - x_min)/(x_max - x_min);
    vec[1] = (vec[1] - y_min)/(y_max - y_min);

    // Flux
    DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    vec[2] = laplace1.cdf(vec[2]);

    // Radius
    DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    vec[3] = laplace2.cdf(vec[3]);

}

void MyConditionalPrior::print(std::ostream& out) const
{
    out << typical_flux << ' ' << dev_log_flux << ' ';
    out << typical_radius << ' ' << dev_log_radius << ' ';
}

