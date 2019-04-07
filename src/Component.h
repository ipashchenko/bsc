#ifndef BSC_COMPONENT_H
#define BSC_COMPONENT_H

#include <valarray>
#include <vector>
#include <complex>

const double mas_to_rad = 4.84813681109536e-09;
const double PI = 3.1415926;


class Component {
    public:
        virtual void ft(std::valarray<double> u, std::valarray<double> v) = 0;
        virtual const size_t size() const = 0;
        virtual void set_param_vec(std::valarray<double> param) = 0;
        std::valarray<double> get_mu_real() const {
            return mu_real;
        }

        std::valarray<double> get_mu_imag() const {
            return mu_imag;
        }
        virtual void print(std::ostream& out) const = 0;

    protected:
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
};


class DFComponent : public  Component {
    public:
        DFComponent();
        void ft(std::valarray<double> u, std::valarray<double> v) override;
        void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 3;
        }
        void print(std::ostream& out) const override;

    private:
        // Parameters of a single Delta Function
        double dx_, dy_, flux_;
};


class EGComponent : public Component {

    public:
        EGComponent();
        void ft(std::valarray<double> u, std::valarray<double> v) override;
        void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 6;
        }
        void print(std::ostream& out) const override;

    private:
        // Parameters of a single Gaussian
        double dx_, dy_, flux_, bmaj_, e_, bpa_;
};


class CGComponent : public EGComponent {

    public:
        CGComponent();
        void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 4;
        }
        void print(std::ostream& out) const override;

private:
        // Parameters of a single Gaussian
        double dx_, dy_, flux_, bmaj_;
};


#endif //BSC_COMPONENT_H
