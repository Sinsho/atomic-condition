#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>
#include <functional>
#include <chrono>
#include <memory>

#include "fpUtil.h"
#include "fpInterface.h"

// libcmaes library
#include "cmaes.h"

class InputFitnessPair {
public:
    double input;
    double fitness;

    InputFitnessPair() {
        input = 0;
        fitness = -std::numeric_limits<double>::max();
    }
    InputFitnessPair(double in, double fit) {
        input = in;
        fitness = fit;
    }
    bool operator < (const InputFitnessPair &rhs) const {
        return fitness < rhs.fitness;
    }
    bool operator > (const InputFitnessPair &rhs) const {
        return fitness > rhs.fitness;
    }
};
class InputFitnessPairGreater {
public:
    bool operator () (const InputFitnessPair& lhs, const InputFitnessPair& rhs) const {
        return lhs.fitness > rhs.fitness;
    }
}pairGreater;
class InputFitnessPairLess {
public:
    bool operator () (const InputFitnessPair& lhs, const InputFitnessPair& rhs) const {
        return lhs.fitness < rhs.fitness;
    }
}pairLess;

struct inputAtomicPair{
    std::vector<double> inputs;
    double atomicCond;
};

class InstructionInfoMultVar{
private:
    uint64_t instID;
    uint64_t opcode;
    inputAtomicPair candidate = {{}, -std::numeric_limits<double>::max()};
    // For recording prioritization results.
    uint32_t topInputCountToEnd;
    double   topInputConditionToEnd;
public:
    InstructionInfoMultVar() {
        instID = 0;
        opcode = 0;
    }
    InstructionInfoMultVar(uint64_t id, uint64_t op) {
        instID = id;
        opcode = op;
    }
    uint64_t getInstID() const { return instID; }
    uint64_t getOpCode() const { return opcode; }
    inputAtomicPair getCandidatePair() {return candidate;}
    void setCandidatePair(inputAtomicPair pair){candidate = pair;}

    void setTopInputCountToEnd(uint32_t c) { topInputCountToEnd = c; }
    void setTopInputConditionToEnd(double d) { topInputConditionToEnd = d; }
    uint32_t getTopInputCountToEnd() const { return topInputCountToEnd; }
    double getTopInputConditionToEnd() const { return topInputConditionToEnd; }
};

class InstructionInfo {
private:
    uint64_t instID;
    uint64_t opcode;
    uint32_t recordSize = 100;
    std::vector<InputFitnessPair> topInputsRandom;
    std::vector<InputFitnessPair> topInputsEvolution;
    // For recording prioritization results.
    uint32_t topInputCountToEnd;
    double   topInputConditionToEnd;

private:
    void pushHeap(double input, double fitness) {
        topInputsRandom.push_back(InputFitnessPair(input, fitness));
        std::push_heap(topInputsRandom.begin(), topInputsRandom.end(), pairGreater);
    }
    void popHeap() {
        std::pop_heap(topInputsRandom.begin(), topInputsRandom.end(), pairGreater);
        topInputsRandom.pop_back();
    }
    InputFitnessPair topHeap() {
        return topInputsRandom[0];
    }
public:
    InstructionInfo() {
        instID = 0;
        opcode = 0;
    }
    InstructionInfo(uint64_t id, uint64_t op) {
        instID = id;
        opcode = op;
    }
    uint64_t getInstID() const { return instID; }
    uint64_t getOpCode() const { return opcode; }
    uint32_t getRecordSize() const { return recordSize; }
    void setRecordSize(uint32_t k) { recordSize = k; }

    std::vector<InputFitnessPair> getInputsRandom() const { return topInputsRandom; }
    std::vector<InputFitnessPair> getInputsEvolution() const { return topInputsEvolution; }
    int getInputsRandomSize() const { return topInputsRandom.size(); }
    int getInputsEvolutionSize() const { return topInputsEvolution.size(); }
    void setInputsEvolution(std::vector<InputFitnessPair> l) { topInputsEvolution = l; }
    InputFitnessPair getInputsEvolutionTop() const { return topInputsEvolution[0]; }

    void setTopInputCountToEnd(uint32_t c) { topInputCountToEnd = c; }
    void setTopInputConditionToEnd(double d) { topInputConditionToEnd = d; }
    uint32_t getTopInputCountToEnd() const { return topInputCountToEnd; }
    double getTopInputConditionToEnd() const { return topInputConditionToEnd; }

    void pushInputFitness(double input, double fitness) {
        if (!std::isfinite(fitness))
            return;
        if (topInputsRandom.size() < recordSize)
        {
            pushHeap(input, fitness);
            return;
        }
        if (fitness > topHeap().fitness)
        {
            popHeap();
            pushHeap(input, fitness);
        }
    }

    void printRandomInfo(std::unique_ptr<FloatingPointFunction> & funcPtr) {
        std::sort(topInputsRandom.begin(), topInputsRandom.end(), pairGreater);

        std::cout << "************************************************\n";
        std::cout << "Instruction ID: " << instID << '\n';
        std::cout << "Instruction OP: " << opcode << '\n';
        std::cout.precision(5);
        uint32_t inputSize = topInputsRandom.size();
        std::cout << "Valid inputs in our records: " << topInputsRandom.size() << '\n';
        std::cout << "Largest fitness after random: " << std::scientific << topInputsRandom[0].fitness << '\n';
        std::cout << "50th    fitness after random: " << std::scientific << topInputsRandom[inputSize/2].fitness << '\n';
        std::cout.precision(16);
        std::cout << "    Largest's input:  " << std::setw(25) << std::scientific << topInputsRandom[0].input << '\n';
        std::vector<double> inputs = {topInputsRandom[0].input};
        std::cout << "    Largest's output: " << std::setw(25) << std::scientific << funcPtr->callAndGetResult(inputs) << '\n';

        std::make_heap(topInputsRandom.begin(), topInputsRandom.end(), pairGreater);
    }
    void printEvolutionInfo(std::unique_ptr<FloatingPointFunction> & funcPtr) {
        std::cout << "------------------\n";
        std::cout << "Instruction ID: " << instID << '\n';
        std::cout << "Instruction OP: " << opcode << '\n';
        std::cout.precision(5);
        uint32_t inputSize = topInputsEvolution.size();
        std::cout << "Valid inputs in our records: " << topInputsEvolution.size() << '\n';
        std::cout << "Largest fitness after evo: " << std::scientific << topInputsEvolution[0].fitness << '\n';
        std::cout << "50th    fitness after evo: " << std::scientific << topInputsEvolution[inputSize/2].fitness << '\n';
        std::cout.precision(16);
        std::cout << "    Largest's input:  " << std::setw(25) << std::scientific << topInputsEvolution[0].input << '\n';
        std::vector<double> inputs = {topInputsEvolution[0].input};
        std::cout << "    Largest's output: " << std::setw(25) << std::scientific << funcPtr->callAndGetResult(inputs) << '\n';
    }

    void printBriefInfo(int funcIndex) {
        std::cout << "2. Searching on Operations.\nFunction Index: " << funcIndex << ", ";
        std::cout << "Analyzing Operation: " << instID << ", OPCode: " << opcode << '\n';
        std::cout << "  Largest's input: " << std::setw(15) << std::scientific << topInputsEvolution[0].input << ", ";
        std::cout << "  Fitness: " << std::scientific << topInputsEvolution[0].fitness << '\n';
    }
};

class EvoSolver {
private:
    std::unique_ptr<FloatingPointFunction> funcUnderTest;
    std::map<uint64_t, InstructionInfo> instMap;
    std::map<uint64_t, InstructionInfoMultVar> instMapMultVar;
    uint32_t unstableInstCount = 0;
    int32_t GSLFuncIndex;
    // For record execution time.
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point finishTime;
    std::chrono::duration<double> elapsedTime;

    // Params:
    // For _1RandomSearch
    double initExpRange = 20;
    double initCenterRate = 0.15;
    uint32_t randomIteration = 100000;
    // For _2EvolutionSearch
    double evoGeometricP = 0.25;
    double evoNormalFactorStart = 1e-2;
    double evoNormalFactorEnd   = 1e-13;
    uint32_t evoIterations = 100;
    uint32_t evoPerIterationSize = 100;

    std::mt19937 mtGenerator;
    std::uniform_real_distribution<double> uni01;

    std::string outPath = "tempOutput.out";

public:
    EvoSolver() :
        mtGenerator(0xdeadbeef),
        uni01(0.0, 1.0)
    {
        // Clean File.
        std::ofstream myfile;
        myfile.open(outPath, std::ofstream::out);
        myfile.close();
    }

    void run(std::unique_ptr<FloatingPointFunction> & funcPtr, int index) {
        startTime = std::chrono::high_resolution_clock::now();
        _init(funcPtr, index);
        if (funcUnderTest->getArgCount() == 1) {
            _1RandomSearch();
            _2EvolutionSearch();
            _3Prioritize();
        } else {
            _1RandomSearchMultVar();
            _2CMAESearch();
        }
        finishTime = std::chrono::high_resolution_clock::now();
        elapsedTime = finishTime - startTime;

        _writeToFile();
        _printInfo();
    }

    void setRandomIteration(uint32_t iter) { randomIteration = iter; }
    void setEvoIteration(uint32_t iter) {evoIterations = iter; }

private:
    void _init(std::unique_ptr<FloatingPointFunction> & funcPtr, int index) {
        GSLFuncIndex = index;
        funcUnderTest = std::move(funcPtr);

        // Clear member values
        unstableInstCount = 0;
        // Fast clear instMap.
        std::map<uint64_t, InstructionInfo> emptyMap;
        std::swap(emptyMap, instMap);
        std::cout << "Current Analyzing Function Index: " << index << std::endl;
    }

    double _initDist() {
        double x = fpUtil::randDouble();
        double p01 = uni01(mtGenerator);
        if (p01 < initCenterRate) {
            p01 = uni01(mtGenerator);
            uint64_t xsign = fpUtil::getDoubleSign(x);
            uint64_t xexpo = 1023 + (uint64_t)(initExpRange*2*p01-initExpRange);
            uint64_t xfrac = fpUtil::getDoubleFrac(x);
            x = fpUtil::buildDouble(xsign, xexpo, xfrac);
        }
        return x;
    }

    void _1RandomSearchMultVar(){
        int dim = funcUnderTest->getArgCount();

        for (int i = 1; i <= randomIteration; i++) {
            std::vector<double> inputs;
            for (int j = 0; j < dim; j++){
                inputs.push_back(_initDist());
            }

            funcUnderTest->call(inputs);

            // print the progress.
            if (i % (randomIteration/100) == 0)
                std::cout << "1. Random Search Processing: " << i / (randomIteration/100) << "%, " << i << "..." << "\r" << std::flush;

            // Only log when the execution is success.
            if (!funcUnderTest->isSuccess())
                continue;
//            funcUnderTest->getResult();

            std::vector<InstInfo> infoList = funcUnderTest->getInstInfoList();


            for (const auto &info : infoList) {
                if (instMapMultVar.count(info.instID) == 0) {
                    instMapMultVar[info.instID] = InstructionInfoMultVar(info.instID, info.opcode);
                }

                double atomicCond = fpUtil::revisedCondition(info.opcode, info.op1, info.op2);
//                std::cout << "Atomic cond: " << atomicCond << std::endl;
                InstructionInfoMultVar &curInst = instMapMultVar[info.instID];
                if (curInst.getCandidatePair().atomicCond <  atomicCond && std::isfinite(atomicCond)){
                    curInst.setCandidatePair({inputs, atomicCond});
                }
            }
        }

        _printInstMapMultVar();
    }

    void _2CMAESearch(){
        std::vector<InstInfo> infoList = funcUnderTest->getInstInfoList();
        for (const auto &info : infoList) {
            libcmaes::FitFunc optFunc = [this, &info] (const double *x, const int N) mutable {
                std::vector<double> inputs;
                for (int i=0; i < N; i++) {
                    inputs.push_back(x[i]);
                }
                uint64_t instID = info.instID;
//                std::cout << infoList == this->funcUnderTest->getInstInfoList() << std::endl;
//                double res = this->funcUnderTest->callAndGetResult(inputs);
//                std::cout << "I hate c++ " << res << std::endl;
//                u->call(inputs);
//                bool blub = infoList == this->funcUnderTest->getInstInfoList();

                this->funcUnderTest->call(inputs);
                std::vector<InstInfo> tmpInfoList = this->funcUnderTest->getInstInfoList();
//                std::cout << infoList ==  << std::endl;
//                double bla = -fpUtil::revisedCondition(info.opcode, info.op1, info.op2);
//                std::cout << bla << std::endl;

                double res = 0;
                for (const auto &tmpInfo : tmpInfoList) {
                    if (instID != tmpInfo.instID)
                        continue;
                    res = fpUtil::revisedCondition(tmpInfo.opcode, tmpInfo.op1, tmpInfo.op2);
                }

//                std::cout << res << std::endl;
                return -res;
            };

            std::vector<double> in = instMapMultVar[info.instID].getCandidatePair().inputs;
            std::vector<double> x0 = in;

//            double bla [funcUnderTest->getArgCount()];
            double bla [funcUnderTest->getArgCount()] = {-0.180807,-209.222,  -37.8287};
            std::vector<double> blub = {-0.0322845, 6.90581e-08, 0.00534196};
            funcUnderTest->call(blub);
            std::vector<InstInfo> infoList2 = funcUnderTest->getInstInfoList();
            for (const auto &info2 : infoList2) {
                double k = fpUtil::revisedCondition(info2.opcode, info2.op1, info2.op2);
                std::cout << "Is: " << k << std::endl;
            }

//            for (int i = 0; i < x0.size(); i++){
//                bla[i] = x0[i];
//            }
//
//            double test = optFunc(bla, funcUnderTest->getArgCount());

            double test = optFunc(bla, funcUnderTest->getArgCount());
            std::cout << "Test : " << test << std::endl;
            double sigma = 0.1;
            libcmaes::CMAParameters<> cmaparams(x0,sigma);
            libcmaes::CMASolutions cmasols = libcmaes::cmaes<>(optFunc,cmaparams);
            cmaparams.set_algo(aCMAES);
            std::cout << "best solution: " << cmasols << std::endl;
            std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
            std::cout << cmasols.run_status() << std::endl;
        }


//        int dim = 3; // problem dimensions.
//        std::vector<double> x0 = {3.0, 2.1, 4.2};
//        double sigma = 0.1;
        //int lambda = 100; // offsprings at each generation.

        //cmaparams.set_algo(BIPOP_CMAES);
    }

    // The results stored in instMap.
    void _1RandomSearch() {
        double x, y;

        for (int i = 1; i <= randomIteration; i++) {
            x = _initDist();
            std::vector<double> inputs = {x};
            funcUnderTest->call(inputs);

            // print the progress.
            if (i % (randomIteration/100) == 0)
                std::cout << "1. Random Search Processing: " << i / (randomIteration/100) << "%, " << i << "..." << "\r" << std::flush;

            // Only log when the execution is success.
            if (!funcUnderTest->isSuccess())
                continue;
            y = funcUnderTest->getResult();

            std::vector<InstInfo> infoList = funcUnderTest->getInstInfoList();
            for (const auto &info : infoList) {
                if (instMap.count(info.instID) == 0) {
                    instMap[info.instID] = InstructionInfo(info.instID, info.opcode);
                }

                double fit = fpUtil::revisedCondition(info.opcode, info.op1, info.op2);

                InstructionInfo &curInst = instMap[info.instID];
                curInst.pushInputFitness(x, fit);
            }
        }
        std::cout << "\n1. Random Search Done.\n";
    }

    // The results stored in instMap.
    void _2EvolutionSearch() {
        std::geometric_distribution<int> geometricDist(evoGeometricP);
        std::normal_distribution<double> normalDist(0, 1);

        double x, y;

        for (auto &kv : instMap) {
            uint64_t instid = kv.first;
            InstructionInfo & curInst = kv.second;

            if (curInst.getInputsRandomSize() == 0)
                continue;

            std::vector<InputFitnessPair> inputList = curInst.getInputsRandom();
            // Sort the inputList, with the largest fitness first.
            std::sort(inputList.begin(), inputList.end(), pairGreater);

            // Filter the actual stable "suspicious" inst.
            double largestFitness = inputList[0].fitness;
            if (largestFitness < 1e1)
                continue;

            unstableInstCount++;

            // Print the info before evolution search.
            // curInst.printRandomInfo(funcUnderTest);

            // Generates new inputs based on current input list.
            for (int i = 0; i < evoIterations; i++) {
                double evoNormalFactor = exp(log(evoNormalFactorStart) +
                    ((double)i / evoIterations) * log(evoNormalFactorEnd / evoNormalFactorStart));
                for (int j = 0; j < evoPerIterationSize; j++) {
                    // 2.1 select.
                    int seletedIndex;
                    // make sure the index is smaller than input list size.
                    while ((seletedIndex = geometricDist(mtGenerator)) >= inputList.size());
                    double curInput = inputList[seletedIndex].input;

                    // 2.2 mutation.
                    double mutation = evoNormalFactor * fabs(curInput) * normalDist(mtGenerator);
                    double newInput = curInput + mutation;
                    if (!std::isfinite(newInput))
                        continue;

                    // 2.3 run with the new input.
                    std::vector<double> inputs = {newInput};
                    funcUnderTest->call(inputs);
                    if (!funcUnderTest->isSuccess())
                        continue;
                    y = funcUnderTest->getResult();

                    std::vector<InstInfo> infoList = funcUnderTest->getInstInfoList();
                    for (const auto &info : infoList) {
                        // We only care about the current instruction.
                        if (info.instID != curInst.getInstID())
                            continue;
                        double newFit = fpUtil::revisedCondition(info.opcode, info.op1, info.op2);
                        if (!std::isfinite(newFit))
                            continue;
                        inputList.push_back(InputFitnessPair(newInput, newFit));
                    }
                }
                // After each iteration, only reserve the top 100 inputs.
                std::sort(inputList.begin(), inputList.end(), pairGreater);
                inputList.resize(curInst.getRecordSize());
            }
            // After all iterations, store the inputList into curInst.
            curInst.setInputsEvolution(inputList);

            // curInst.printEvolutionInfo(funcUnderTest);
            curInst.printBriefInfo(GSLFuncIndex);
        }
    }

    // The results stored in prioritizedInput.
    void _3Prioritize() {

        double x, y;
        int status;
        gsl_sf_result res;

        std::vector<InputFitnessPair> inputsWithCountToEnd;
	std::map<double, double> io;
        // The prior score of an input is calculated from
        // the correspounding instruction to end.
        for (auto &kv : instMap) {
            uint64_t instid = kv.first;
            InstructionInfo & curInst = kv.second;

            if (curInst.getInputsEvolutionSize() == 0)
                continue;

            double curInput = curInst.getInputsEvolutionTop().input;
            double curFit = curInst.getInputsEvolutionTop().fitness;

            std::vector<double> inputs = {curInput};
            funcUnderTest->call(inputs);
            y = funcUnderTest->getResult();

            bool started = false;
            double logConditionToEnd = 0;
            uint32_t countToEnd = 0;
            std::vector<InstInfo> infoList = funcUnderTest->getInstInfoList();
            for (const auto &info : infoList) {
                // we need to find the curInput running to the curInst with the curFit.
                // Considering a inst maybe executed multiple times.
                if (started == false && info.instID == instid) {
                    double tmpFit = fpUtil::revisedCondition(info.opcode, info.op1, info.op2);
                    if (tmpFit == curFit) {
                        started = true;
                    }
                }
                if (started == true) {
                    double rawCond = fpUtil::rawCondition(info.opcode, info.op1, info.op2);
                    // std::cout << "Raw Cond: " << rawCond;
                    // std::cout << ' ' << info.op1 << ' ' << info.op2 << ' ' << info.opcode << "\n";
                    // Condition
                    if ( instMap.count(info.instID) != 0 && std::isfinite(rawCond) ) {
                        logConditionToEnd += log(rawCond);
                    }
                    // Count
                    if ( instMap.count(info.instID) != 0) {
                        countToEnd += 1;
                    }
                }
            }
            curInst.setTopInputCountToEnd(countToEnd);
            curInst.setTopInputConditionToEnd(logConditionToEnd);

            inputsWithCountToEnd.push_back(InputFitnessPair(curInput, countToEnd));
	    io[curInput] = y;
            // std::cout << curInst.getInputsRandomSize() << ' ' << curInst.getInputsEvolutionSize() << '\n';
        }
        // Prioritize based on condition to end.
        sort(inputsWithCountToEnd.begin(), inputsWithCountToEnd.end());

        std::cout << "***********Results after Prioritize***********\n";
        std::cout << "Most suspicious input first: \n";
        for (auto & i : inputsWithCountToEnd) {
            std::cout << std::setprecision(16) << std::scientific << i.input << "\t\t";
	    std::cout << "Output: " << io[i.input] << '\n';
        }
        std::cout << "End of suspicious input list.\n";
    }

    // The evolution search info.
    void _printInfo() {
        std::cout << "************************************************\n";
        std::cout << "Actual    Unstable Instructions Size: " << unstableInstCount << '\n';
        std::cout << "Potential Unstable Instructions Size: " << instMap.size() << '\n';
        std::cout.unsetf(std::ios_base::floatfield);
        std::cout << "Execution Time: " << elapsedTime.count() << " sec.\n";
        std::cout << "************************************************\n";
    }

    void _printInstMapMultVar(){
        for (const auto &entry : instMapMultVar){
            uint64_t index = entry.first;
            InstructionInfoMultVar value = entry.second;

            std::cout << "************************************************\n";
            std::cout << "Index: " << index << " atomicCond: " << value.getCandidatePair().atomicCond << " inputs: ";
            for(double d : value.getCandidatePair().inputs){
                std::cout << d << ", ";
            }
            std::cout << std::endl;
            std::cout << "************************************************\n";
        }
    }

    void _writeToFile() {
        std::ofstream myfile;
        myfile.open(outPath, std::ofstream::out | std::ofstream::app);
        // Write Timestamp.
        std::time_t stamp = time(NULL);
        myfile << "---------------------------------------\n";
        myfile << "Time: " << std::ctime(&stamp) << '\n';
        // Write Function Index.
        myfile << "Function Index: " << GSLFuncIndex << '\n';
        myfile << "Actual    Unstable Instructions Size: " << unstableInstCount << '\n';
        myfile << "Potential Unstable Instructions Size: " << instMap.size() << '\n';
        myfile << "Execution Time: " << elapsedTime.count() << " sec.\n";
        myfile << "Format: InstID, OpCode, MaxAtomCond, Input, Output, CountToEnd, ConditionToEnd\n";
        // Write Computation Results.
        for (const auto &kv : instMap) {
            uint64_t instid = kv.first;
            const InstructionInfo & curInst = kv.second;

            if (curInst.getInputsEvolutionSize() == 0)
                continue;
            InputFitnessPair pair = curInst.getInputsEvolutionTop();
            myfile << "Data: ";
            myfile << curInst.getInstID() << ", ";
            myfile << curInst.getOpCode() << ", ";
            myfile << std::setprecision(5) << std::scientific << pair.fitness << ", ";
            myfile << std::setprecision(16) << std::scientific << pair.input << ", ";
            std::vector<double> inputs = {pair.input};
            myfile << std::setprecision(16) << std::scientific << funcUnderTest->callAndGetResult(inputs) << ", ";
            myfile << curInst.getTopInputCountToEnd() << ", ";
            myfile << curInst.getTopInputConditionToEnd() << "\n";
        }
        myfile << "***************************************\n";
        myfile.close();
    }
};

int main(int argc, char *argv[]) {
    // // Init communicator.
    // Communicator &comm = Communicator::getInstance();

    EvoSolver es;
    std::unique_ptr<FloatingPointFunction> funcPtr;
    //

    if (argc == 1 || (argc > 1 && strcmp(argv[1], "example") == 0)) {
        int index = 0;
        if (argc > 2)
            index = atoi(argv[2]);
        if (index >= simpleFuncList.size()) {
            std::cout << "Invalid index in simpleFuncList\n";
            return 0;
	    }
        funcPtr.reset(new SimpleFunction(index));
        es.run(funcPtr, index);
    }
    else if (argc > 2 && strcmp(argv[1], "gsl") == 0) {
	if (!strcmp(argv[2],"all")) {
	    for (int i = 0; i < GSLFuncList.size(); i++) {
//            funcPtr.reset(new GSLFunction(i));
		    es.run(funcPtr, i);
	    }
	}
	else {
	    int index = atoi(argv[2]);
	    if (index >= GSLFuncList.size()) {
		std::cout << "Invalid index in GSLFuncList.\n";
		return 0;
	    }
//            funcPtr.reset(new GSLFunction(index));
            es.run(funcPtr, index);
	}
    }
    else {
        std::cout << "Invalid argument." << std::endl;
	std::cout << "Valid Example: \n\tbin/gslSolver.out gsl 73\n\tor\n\tbin/gslSolver.out example\n";
    }

    return 0;
}
