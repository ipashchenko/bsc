#ifndef BSC_SRC_MYEXCEPTIONS_H_
#define BSC_SRC_MYEXCEPTIONS_H_

#include <exception>


class BadReferenceAntenna : public std::exception {
        const char * what () const noexcept override {
            return "No such antenna number in data!";
        }
};

#endif //BSC_SRC_MYEXCEPTIONS_H_