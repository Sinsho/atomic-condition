#ifndef FPINTERFACE_H
#define FPINTERFACE_H

#include "communicator.h"
#include <functional>
#include <gsl/gsl_sf.h>
#include "target.h"


template <typename R, typename ... Types> constexpr size_t getArgumentCount( R(*f)(Types ...))
{
    return sizeof...(Types);
}

template<int... Is>
struct seq { };

template<int N, int... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };

template<int... Is>
struct gen_seq<0, Is...> : seq<Is...> { };


template<class... Args, int... Is>
double callByVector(std::function<double(Args...)> f, std::vector<double>& arguments, seq<Is...>)
{
    return f(arguments[Is]...);
}

template<class... Args>
double callByVector(std::function<double(Args...)> f, std::vector<double>& arguments)
{
    return callByVector(f, arguments, gen_seq<sizeof...(Args)>());
}

template<class... Args, int... Is>
int gslCallByVector(std::function<int(Args...)> f, std::vector<double>& arguments, gsl_sf_result* gslSfResult, seq<Is...>)
{
    return f(arguments[Is]..., gslSfResult);
}

template<class... Args>
int gslCallByVector(std::function<int(Args...)> f, std::vector<double>& arguments, gsl_sf_result* gslSfResult)
{
    return gslCallByVector(f, arguments, gslSfResult, gen_seq<sizeof...(Args) - 1>());
}

struct vecFunc {
    std::function<double(std::vector<double>)> func;
    int argCount;
};

struct gslVecFunc {
    std::function<int(std::vector<double>, gsl_sf_result*)> func;
    int argCount;
};

template<class... Args>
vecFunc vectorizeFunction(std::function<double(Args...)> f, int argCount){
    vecFunc vecStruct;
    vecStruct.func = [f](std::vector<double> args){return callByVector(f, args);};
    vecStruct.argCount = argCount;

    return vecStruct;
};

template<class... Args>
gslVecFunc gslVectorizeFunction(std::function<int(Args...)> f, int argCount){
    gslVecFunc vecStruct;

    vecStruct.func = [f](std::vector<double> args, gsl_sf_result* gslSfResult){return gslCallByVector(f, args, gslSfResult);};
    // Deduct arguments by one because gsl_sf_result is passed
    vecStruct.argCount = argCount - 1;

    return vecStruct;
};

const std::vector<vecFunc> simpleFuncList = {
        vectorizeFunction(std::function<double(double, double, double)>(foo), getArgumentCount(foo)),
};

const std::vector<gslVecFunc> GSLFuncList = {
        gslVectorizeFunction (std::function<int(double, double, double, double, gsl_sf_result*)> (gsl_sf_hyperg_2F1_e), getArgumentCount(gsl_sf_hyperg_2F1_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_bessel_Ynu_e), getArgumentCount(gsl_sf_bessel_Ynu_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_bessel_Jnu_e), getArgumentCount(gsl_sf_bessel_Jnu_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_bessel_Inu_e), getArgumentCount(gsl_sf_bessel_Inu_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_bessel_Inu_scaled_e), getArgumentCount(gsl_sf_bessel_Inu_scaled_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_coulomb_CL_e), getArgumentCount(gsl_sf_coulomb_CL_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_gegenpoly_2_e), getArgumentCount(gsl_sf_gegenpoly_2_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_gegenpoly_3_e), getArgumentCount(gsl_sf_gegenpoly_3_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_laguerre_1_e), getArgumentCount(gsl_sf_laguerre_1_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_laguerre_2_e), getArgumentCount(gsl_sf_laguerre_2_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_laguerre_3_e), getArgumentCount(gsl_sf_laguerre_3_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_poch_e), getArgumentCount(gsl_sf_poch_e)),
        gslVectorizeFunction (std::function<int(double, double, gsl_sf_result*)> (gsl_sf_beta_e), getArgumentCount(gsl_sf_beta_e)),
};

class FloatingPointFunction {
public:
    virtual void call(std::vector<double> x) = 0;
    virtual double getResult() = 0;
    virtual double callAndGetResult(std::vector<double> x) = 0;
    virtual bool isSuccess() = 0;
    virtual int getArgCount() = 0;
    std::vector<InstInfo> getInstInfoList() {
        return comm.getInstInfoList();
    }
    FloatingPointFunction() : comm( Communicator::getInstance() ) {
        comm.initComm();
    }
protected:
    double out;
    std::vector<double> in;
    int status;
    Communicator &comm;
};


class GSLFunction : public FloatingPointFunction {
public:
    GSLFunction(int index) {
        gsl_set_error_handler_off();

        if (index < 0 || index >= GSLFuncList.size()) {
            std::cout << "Invalid index in [GSLFunction]: " << index << '\n';
            GSLFuncRef = GSLFuncList[0];
            return;
        }
        GSLFuncRef = GSLFuncList[index];
        return;
    }

    void call(std::vector<double> x) {
        comm.clear();
        in = x;
        status = GSLFuncRef.func(in, &gslres);
        out = gslres.val;
    }

    double callAndGetResult(std::vector<double> x) {
        call(x);
        return out;
    }

    int getArgCount(){
        return GSLFuncRef.argCount;
    }

    double getResult() { return out; }
    bool isSuccess() { return (status == GSL_SUCCESS && !std::isnan(out) && std::isfinite(out)); }

private:
    gslVecFunc GSLFuncRef;
    gsl_sf_result gslres;
};

class SimpleFunction : public FloatingPointFunction {
public:
    SimpleFunction(int index) {
        if (index < 0 || index >= simpleFuncList.size()) {
            std::cout << "Invalid index in [SimpleFunction]: " << index << '\n';
            funcRef = simpleFuncList[0];
            return;
        }
        funcRef = simpleFuncList[index];
    }
    void call(std::vector<double> x) {
        comm.clear();
        in = x;
        out = funcRef.func(x);
        if (std::isnan(out)) {
            status = -1;
        }
        else {
            status = 0;
        }
    }

    double testFuncResult(std::vector<double> test){
        return funcRef.func(test);
    }

    int getArgCount(){
        return funcRef.argCount;
    }

    double callAndGetResult(std::vector<double> x) {
        call(x);
        return out;
    }

    double getResult() { return out; }
    bool isSuccess() { return (status == 0); }

private:
    vecFunc funcRef;
};

#endif