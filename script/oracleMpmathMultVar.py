import mpmath
import io


def func_jnu(nu, x):
    return mpmath.besselj(nu, x)


def func_hypergeom2f1(a, b, c, z):
    return mpmath.hyp2f1(a, b, c, z)


def func_ynu(nu, x):
    return mpmath.bessely(nu, x)


def func_inu(nu, x):
    return mpmath.besseli(nu, x)


def func_inu_scaled(nu, x):
    scale = mpmath.exp(-abs(x))
    return mpmath.besseli(nu, x) * scale


def func_coulomb_cl(l, eta):
    return mpmath.coulombc(l, eta)


def func_gegenbauer2(a, z):
    return mpmath.gegenbauer(2, a, z)


def func_gegenbauer3(a, z):
    return mpmath.gegenbauer(3, a, z)


def func_laguerre1(a, z):
    return mpmath.laguerre(1, a, z)


def func_laguerre2(a, z):
    return mpmath.laguerre(2, a, z)


def func_laguerre3(a, z):
    return mpmath.laguerre(3, a, z)


def func_pochhammer(x, n):
    return mpmath.rf(x, n)


def func_beta(x, y):
    return mpmath.beta(x, y)


funcDict = {
    0: func_jnu,
    1: func_hypergeom2f1,
    2: func_ynu,
    3: func_inu,
    4: func_inu_scaled,
    5: func_coulomb_cl,
    6: func_gegenbauer2,
    7: func_gegenbauer3,
    8: func_laguerre1,
    9: func_laguerre2,
    10: func_laguerre3,
    11: func_pochhammer,
    12: func_beta,
}


class GetOracle:
    def __init__(self):
        # Set mpmath precision
        mpmath.mp.prec = 128

    def getOracleValue(self, funcIndex, inputX):
        funcIndex = int(funcIndex)
        # print("Func: ", funcIndex, "Input: ", inputX)
        if funcIndex not in funcDict:
            print("Function Index Error")
            return None
        mpf_list = []

        for d in inputX:
            mpf_list.append(mpmath.mpf(d))

        func = funcDict[funcIndex]
        outY = func(*mpf_list)
        if outY == None:
            return None
        outY = outY.real
        return outY

    def calcRelativeError(self, oracle, output):
        oracle = mpmath.mpf(oracle)
        output = mpmath.mpf(output)
        delta = abs(output - oracle)
        if delta == 0:
            return 0
        try:
            return abs(delta / oracle)
        except:
            return float("inf")


class OutputParser:
    def __init__(self):
        self.outputFile = "tempOutput.out"
        self.writeToFile = "Output.out"
        with io.open(self.outputFile) as f:
            self.rawData = f.readlines()
        self.data = {}
        self.oracleMod = GetOracle()

    def readAndCalculate(self):
        valid = True
        for line in self.rawData:
            # Determine whether valid
            if line.startswith('------------------'):
                valid = True
                continue
            if valid == False:
                continue

            # Read only for valid index.
            if line.startswith("Function Index"):
                curIndex = int(line.split(":")[1])
                self.data[curIndex] = {}
                continue
            if line.startswith("Execution Time"):
                curTime = float(line.split(":")[1].split()[0])
                self.data[curIndex]['time'] = curTime
                continue
            if line.startswith("**************"):
                self.printInfo(curIndex)
                # Wait for input, system('pause')
                # input()
                continue
            if line.startswith("Format:"):
                self.data[curIndex]['input_list'] = []
            # Format: InstID, OpCode, MaxAtomCond, Input, Output, CountToEnd, ConditionToEnd
            if line.startswith("Data:"):
                inputInfo = {}
                terms = line.split(":")[1]
                terms = terms.split(',')
                inputInfo['inst_id'] = int(terms[0])
                inputInfo['op_code'] = int(terms[1])
                inputInfo['atom_cond'] = float(terms[2])
                inputInfo['arg_count'] = int(terms[3])
                inputInfo['inputs'] = []

                for d in range(inputInfo['arg_count']):
                    inputInfo['inputs'].append(float(terms[4 + d]))
                # Input / Output
                inputInfo['output'] = float(terms[4 + inputInfo['arg_count']])
                # print(inputX, outputY)
                # Data for priorization
                countToEnd = int(terms[5 + inputInfo['arg_count']])
                conditionToEnd = float(terms[6 + inputInfo['arg_count']])
                inputInfo['count_to_end'] = countToEnd
                inputInfo['condition_to_end'] = conditionToEnd

                # Calculate Oracle
                oracleVal = self.oracleMod.getOracleValue(curIndex, inputInfo['inputs'])
                relativeErr = self.oracleMod.calcRelativeError(oracleVal, inputInfo['output'])
                inputInfo['oracle'] = oracleVal
                inputInfo['relative_err'] = relativeErr
                self.data[curIndex]['input_list'].append(inputInfo)
                continue

    def printInfo(self, index):
        print("------------------------")
        info = self.data[index]
        print("Function Index:", index)
        inputinfo = info['input_list']
        relList = sorted(inputinfo, key=lambda x: x['relative_err'], reverse=True)
        print("Max Relative Error:")
        for item in relList[:1]:

            print("  Input: ", end="")
            for d in item['inputs']:
                print(str.format("{}, ", format(float(d), '.16e')), end="")
            print("")
            print("  Output:", format(float(item['output']), '.16e'))
            print("  Oracle:", format(float(item['oracle']), '.16e'))
            print("        Relative Error:", format(float(item['relative_err']), '.5e'))
        print("------------------------")


if __name__ == '__main__':
    op = OutputParser()
    op.readAndCalculate()
