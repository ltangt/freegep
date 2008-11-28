import sys
import os
import random
import math
import pdb

class freeGEP() :
    _rand = random.Random()
    _dataset = []
    _head = 0
    _tail = 0
    _chromosomeLen =0
    _funcMaxArg = 0
    _funcList = []
    _funcNameList = []
    _funcArgNumList= []
    _terminalList = []
    _chromosomeList = []
    _bestChromosome =''
    _bestFitness = 0
    _preExpList =[]
    _fitnesses = {}
    _M = 1000.0
    _rateMutation = 0.3
    _rateOneRecombine=0.7
    _logfile = sys.stdout
    
    def load(self, filename) :
        fh = file(filename,'r')
        lines = fh.readlines()
        self._dataset = []
        # load the dataset from the file
        for line in lines:
            if len(line) > 1:
                rowset = line.split()                
                introwset = map(int,rowset)
                self._dataset.append(introwset)
        fh.close()
    
    def addFunc(self, func, funcName, numArg) :
        if (numArg > self._funcMaxArg):
            self._funcMaxArg = numArg
        self._funcList.append(func)
        self._funcNameList.append(funcName)
        self._funcArgNumList.append(numArg)

    def addTerminal(self, term) :
        self._terminalList.append(term)

    def initChromosomes(self, head, num):
       self._head = head
       self._tail = head*(self._funcMaxArg-1)+1
       self._chromosomeLen = self._head + self._tail
       self._chromosomeList = []
       for i in range(num) :
           chromosome = ''
           #pdb.set_trace()
           for h in range(self._head):
               chromosome += self._randGenFuncAndTerm()
           for t in range(self._tail):
               chromosome += self._randGenTerminal()
           self._chromosomeList.append(chromosome)
    
    def printChromosomes(self) :
        print self._chromosomeList

    def decode(self):
        self._preExpList = map(self._decode, self._chromosomeList)

    def evaluate(self) :
        self._fitnesses = {}
        index = 0
        for exp in self._preExpList:
            fit = 0.0
            for row in self._dataset:
                expCopy = exp[:]
                y = self._evaluate(expCopy,row)
                fit += self._M - float(abs(row[-1] - y))
            fit /= float(len(self._dataset))
            chromosome = self._chromosomeList[index]
            if self._fitnesses.has_key(fit) :
                self._fitnesses[fit].append(chromosome)
            else:
                self._fitnesses[fit] = [chromosome]
            index += 1
            
    def printFitnesses(self):
        for fit in self._fitnesses:
            print >>self._logfile, str(fit) + '  ' + str(self._fitnesses[fit])

    def evolution(self, generation):
        for i in range(generation) :
            self.decode()
            print >>self._logfile,'evaluate --------------------------'
            self.evaluate()
            self.printFitnesses()
            print >>self._logfile,'select ----------------------------'
            self._select_replicate()
            print >>self._logfile,'mutation --------------------------'
            self._mutation()
            print >>self._logfile,'one-point recombination -----------'
            self._onePointRecombination()

            self._chromosomeList.append(self._bestChromosome)
            print >>self._logfile,'best '+self._bestChromosome+' '+str(self._bestFitness)
            print 'best '+self._bestChromosome+' '+str(self._bestFitness)
            print 'generating %d' % (i)
            if self._bestFitness == self._M: break
        print >>self._logfile,self._bestChromosome +' '+ str(self._bestFitness)

    def openLogFile(self, filepath):
        self._logfile = file(filepath, 'w')

    def closeLogFile(self):
        self._logfile.close()
        
    # ++abb
    def _evaluate(self, exp, row) :
        op = exp[0]
        del exp[0]
        if op in self._funcNameList:
            funcIndex = self._funcNameList.index(op)
            func = self._funcList[funcIndex]
            funcNumArg = self._funcArgNumList[funcIndex]
            args = [];
            for i in range(funcNumArg):
                args.append(self._evaluate(exp,row))
            return func(args)
        else:
            termIndex = self._terminalList.index(op)
            return row[termIndex]

    def _decode(self,chromosome) :
        chrlist = list(chromosome)
        layers = self._decodeSplitLayers(chrlist)
        preExp = self._decodeGenPreExp(layers,0)
        #print preExp
        return list(preExp)

    def _decodeSplitLayers(self, chrlist) :
        layers =[]
        layer = [chrlist[0]]
        layers.append(layer)
        i = 1
        numArg = self._getNumArg(chrlist[0])
        #print chrlist
        while numArg > 0 :
            layer = []
            nextNumArg = 0
            endIndex = i + numArg
            while i < endIndex :
                nextNumArg += self._getNumArg(chrlist[i])
                layer.append(chrlist[i])
                i += 1
            numArg = nextNumArg
            layers.append(layer)
        return layers

    def _decodeGenPreExp(self, layers, layerindex) :
        if layerindex >= len(layers):
            return ''
        layer = layers[layerindex];
        if len(layer) == 0:
            return ''
        ret = layer[0]
        for i in range(self._getNumArg(layer[0])):
            ret += self._decodeGenPreExp(layers, layerindex+1)
        del layer[0]
        return ret

    def _getNumArg(self, op):
        if op in self._funcNameList:
            funcIndex = self._funcNameList.index(op)
            return self._funcArgNumList[funcIndex]
        else:
            return 0

    def _randGenFuncAndTerm(self):
        totalLen = len(self._funcList)+len(self._terminalList)
        r = self._rand.randint(0, totalLen-1)
        if r < len(self._funcList) :
            return self._funcNameList[r]
        else:
            return self._terminalList[r-len(self._funcList)]

    def _randGenTerminal(self):
        r = self._rand.randint(0, len(self._terminalList)-1 )
        return self._terminalList[r] 
    
    def _select_replicate(self):
        fitkeys = self._fitnesses.keys()
        fitkeys.sort()
        count = len(fitkeys)
        newchromosomeList = []
        self._bestChromosome = self._fitnesses[fitkeys[-1]][0]
        self._bestFitness = fitkeys[-1]
        print >>self._logfile,self._bestChromosome
        for i in range(len(self._chromosomeList)-1) :
            r = self._rand.randint(1, count*(count+1)/2-1)
            idx = int(math.sqrt(r)) -1
            chromosome = self._fitnesses[fitkeys[idx]][0]
            newchromosomeList.append(chromosome)
            print >>self._logfile,chromosome
        self._chromosomeList = newchromosomeList

    def _mutation(self):
        index = 0
        newchromosomeList = []
        while index < len(self._chromosomeList):
            if self._randomProb(self._rateMutation) == False:
                newchromosomeList.append(self._chromosomeList[index])
            else:
                pos = self._rand.randint(0,self._chromosomeLen-1)
                chrlist = list(self._chromosomeList[index])
                if pos < self._head:
                    chrlist[pos] = self._randGenFuncAndTerm()
                else:
                    chrlist[pos] = self._randGenTerminal()
                newchromosomeList.append(''.join(chrlist)) 
            print >>self._logfile, newchromosomeList[-1]
            index += 1
        self._chromosomeList = newchromosomeList

    def _onePointRecombination(self):
        index = 0
        newchromosomeList =[]
        while index < len(self._chromosomeList)-1 :
            if self._randomProb(self._rateOneRecombine) == False:
                newchromosomeList.append(self._chromosomeList[index])
                newchromosomeList.append(self._chromosomeList[index+1])
            else:
                chrlist1 = list(self._chromosomeList[index])
                chrlist2 = list(self._chromosomeList[index+1])
                pos = self._rand.randint(1, self._chromosomeLen-2)
                newchrlist1 = chrlist1[:pos]+chrlist2[pos:]
                newchrlist2 = chrlist2[:pos]+chrlist1[pos:]
                newchromosomeList.append(''.join(newchrlist1))
                newchromosomeList.append(''.join(newchrlist2))
            print >>self._logfile, newchromosomeList[-2]
            print >>self._logfile, newchromosomeList[-1]
            index += 2
        while index < len(self._chromosomeList) :
            newchromosomeList.append(self._chromosomeList[index])
            index += 1
        self._chromosomeList = newchromosomeList

    def _randomProb(self, prob):
        val = self._rand.randint(0,1000)
        if prob <= float(val) / 1000.0:
            return True
        else:
            return False

# Add function
def add(args) :
    # print str(args[0])+' + '+str(args[1]) + ' = '+str(args[0]+args[1])
    return args[0]+args[1]

def sub(args):
    # print str(args[0])+' - '+str(args[1]) + ' = '+str(args[0]-args[1])
    return args[0]-args[1]

def mul(args):
    # print str(args[0])+' * '+str(args[1]) + ' = '+str(args[0]*args[1])
    return args[0]*args[1]

if __name__ == '__main__' :
    myGep = freeGEP()
    myGep.load('dataset.txt')
    myGep.addFunc(add,'+',2)
    myGep.addFunc(sub,'-',2)
    #myGep.addFunc(mul,'*',2)
    myGep.addTerminal('a')
    myGep.addTerminal('b')
    myGep.addTerminal('c')
    myGep.initChromosomes(9,100)
    myGep.printChromosomes()
    myGep.openLogFile('log.txt')
    #myGep._decode('++ab+aabbabaa')
    print '\n'
    myGep.evolution(100)
    print '\n'
    myGep.closeLogFile()
    print myGep._bestChromosome + '  '+str(myGep._bestFitness)

