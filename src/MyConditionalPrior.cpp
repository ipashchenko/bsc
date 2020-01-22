#include "MyConditionalPrior.h"


using namespace DNest4;

MyConditionalPrior::MyConditionalPrior(double std)
    :std(std)
{

}

void MyConditionalPrior::from_prior(RNG& rng)
{
    const DNest4::Gaussian gauss1(0.00, 0.05);
    typical_flux = gauss1.generate(rng);
    const DNest4::Gaussian gauss2(1.00, 0.10);
    dev_log_flux = gauss2.generate(rng);

    const DNest4::Gaussian gauss3(-2.0, 0.05);
    typical_radius = gauss3.generate(rng);
    const DNest4::Gaussian gauss4(2.00, 0.10);
    dev_log_radius = gauss4.generate(rng);

    const DNest4::Gaussian gauss5(-1.0, 0.05);
    typical_alpha = gauss5.generate(rng);
    const DNest4::Gaussian gauss6(0.5, 0.03);
    dev_alpha = gauss6.generate(rng);
}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(6);

    if(which == 0)
    {
        const DNest4::Gaussian gauss1(0.0, 0.05);
        logH += gauss1.perturb(typical_flux, rng);
    }
    else if(which == 1)
    {
        const DNest4::Gaussian gauss2(1.0, 0.10);
        logH += gauss2.perturb(dev_log_flux, rng);
    }
    else if(which == 2)
    {
        const DNest4::Gaussian gauss3(-2.0, 0.05);
        logH += gauss3.perturb(typical_radius, rng);
    }
    else if(which == 3)
    {
        const DNest4::Gaussian gauss4(2.00, 0.10);
        logH += gauss4.perturb(dev_log_radius, rng);
    }
    else if(which == 4)
    {
        const DNest4::Gaussian gauss5(-1.0, 0.05);
        logH += gauss5.perturb(typical_alpha, rng);
    }
    else if(which == 5)
    {
        const DNest4::Gaussian gauss6(0.50, 0.03);
        logH += gauss6.perturb(dev_alpha, rng);
    }

    return logH;
}


double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    double logp = 0.0;

    // Position
    DNest4::Gaussian gaussx(0.0, std);
    DNest4::Gaussian gaussy(0.0, std);
    logp += gaussx.log_pdf(vec[0]);
    logp += gaussy.log_pdf(vec[1]);

    // Flux
    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    DNest4::Gaussian gauss1(typical_flux, dev_log_flux);
    //logp += -log(vec[2]) + laplace1.log_pdf(log(vec[2]));
    logp += gauss1.log_pdf(vec[2]);

    // Radius
    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
    //logp += -log(vec[3]) + laplace2.log_pdf(log(vec[3]));
    logp += gauss2.log_pdf(vec[3]);

    // Spectral index
    DNest4::Gaussian gauss3(typical_alpha, dev_alpha);
    //logp += -log(vec[3]) + laplace2.log_pdf(log(vec[3]));
    logp += gauss3.log_pdf(vec[4]);

    return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    // Position
    DNest4::Gaussian gaussx(0.0, std);
    vec[0] = gaussx.cdf_inverse(vec[0]);
    DNest4::Gaussian gaussy(0.0, std);
    vec[1] = gaussy.cdf_inverse(vec[1]);

    // Flux
    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    DNest4::Gaussian gauss1(typical_flux, dev_log_flux);
    vec[2] = gauss1.cdf_inverse(vec[2]);

    // Radius
    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
    vec[3] = gauss2.cdf_inverse(vec[3]);

    // Spectral index
    DNest4::Gaussian gauss3(typical_alpha, dev_alpha);
    vec[4] = gauss3.cdf_inverse(vec[4]);

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    // Position
    DNest4::Gaussian gaussx(0.0, std);
    vec[0] = gaussx.cdf(vec[0]);
    DNest4::Gaussian gaussy(0.0, std);
    vec[1] = gaussy.cdf(vec[1]);

    // Flux
    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
    DNest4::Gaussian gauss1(typical_flux, dev_log_flux);
    vec[2] = gauss1.cdf(vec[2]);

    // Radius
    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
    vec[3] = gauss2.cdf(vec[3]);

    // Spectral index
    DNest4::Gaussian gauss3(typical_alpha, dev_alpha);
    vec[4] = gauss3.cdf(vec[4]);

}

void MyConditionalPrior::print(std::ostream& out) const
{
    out << typical_flux << ' ' << dev_log_flux << ' ';
    out << typical_radius << ' ' << dev_log_radius << ' ';
    out << typical_alpha << ' ' << dev_alpha << ' ';
}
