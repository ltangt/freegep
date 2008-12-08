import sys
import os
import random
import math
import pdb

class freeGEP:
    _rand = random.Random()
    _dataset = []
    _datadimension = 0
    _head = 0
    _tail = 0
    _chromosomeLen =0
    _evalfunc = None
    _evalfunc_reserved = None
    _funcMaxArg = 0
    _funcList = []
    _funcNameList = []
    _funcArgNumList= []
    _constantList = []
    _chromosomeList = []
    _bestChromosome =''
    _bestFitness = 0
    _preExpList =[]
    _fitnesses = {}
    _M = 1000.0
    _rateMutation = 0.3
    _rateInversion = 0.2
    _rateOneRecombine=0.4
    _rateTwoRecombine=0.2
    _rateIS = 0.1
    _rateRIS = 0.1
    _rateConstant = 0.3
    _lenMutation = 2
    _lenInversion = 4
    _lenTransposition = 4
    _logfile = sys.stdout
    _listeners=[]
    _listenerparams=[]
    
    # load the dataset filename
    # @param filename: the file name of the dataset
    def loadfromfile(self, filename) :
        fh = file(filename,'r')
        lines = fh.readlines()
        dataset = []
        for line in lines:
            if len(line) > 1:
                rowset = line.split()
                introwset = map(float,rowset)
                dataset.append(introwset)
        fh.close()
        self.load(dataset)
    
    # load dataset
    # @param  dataset: the dataset to load
    def load(self, dataset):
        self._dataset=dataset
        if len(self._dataset) == 0 :
           return
        self._datadimension = len(self._dataset[0])
                
    # add custom function( or operation) to gene expression    
    # @param func : the function to add, like 'def add(args) : return args[0]+args[1]'
    # @param funcName: the name of the function
    # @param numArg: the number of arguments of the function
    def addFunc(self, func, funcName, numArg) :
        if (numArg > self._funcMaxArg):
            self._funcMaxArg = numArg
        self._funcList.append(func)
        self._funcNameList.append(funcName)
        self._funcArgNumList.append(numArg)
    
    def addConstant(self, constant):
        self._constantList.append(constant)
        
    def addConstants(self, constants):
        map(self.addConstant, constants)
    
    def addListener(self, listener, param):
        self._listeners.append(listener)
        self._listenerparams.append(param)
        
    def setEvalfunc(self, evalfunc, reservedparam):
        self._evalfunc=evalfunc
        self._evalfunc_reserved=reservedparam
    
    # initialize the chromosomes
    # @param head: the length of the head of chromosomes            
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
               chromosome += '.'
           for t in range(self._tail):
               chromosome += self._randGenTerminal()
               if t < self._tail-1 :
                   chromosome += '.'         
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
            evalcount=0
            rowindex=0
            for row in self._dataset:
                expCopy = exp[:]
                y = self._evaluate(expCopy,row)
                if self._evalfunc != None :
                    err = self._evalfunc(index,rowindex,y,row[-1],self._evalfunc_reserved)
                else:
                    err = float(abs(y-row[-1]))
                if err != None :
                    evalcount+=1
                    fit += self._M - err
                rowindex+=1
            fit /= float(evalcount)
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
            self.evaluate()
            #self.printFitnesses()
            self._select_replicate()
            self._mutation()
            self._inversion()
            self._transIS()
            self._transRIS()
            self._onePointRecombination()
            self._twoPointsRecombination()

            self._chromosomeList.append(self._bestChromosome)
            print 'best '+self._bestChromosome+' '+str(self._bestFitness)
            print 'generating %d' % (i)
            if self._bestFitness == self._M: 
                break
            if self._listeners != []:
                listenerindex=0
                for listener in self._listeners:
                    listener(self._chromosomeList, self._bestFitness, 
                             self._listenerparams[listenerindex])
                    listenerindex+=1

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
        elif op[0] == '#' :
            constIndex = int(op[1:])
            return self._constantList[constIndex]
        else :
            termIndex = int(op[1:])
            return row[termIndex]

    def _decode(self,chromosome) :
        chrlist = chromosome.split('.')
        layers = self._decodeSplitLayers(chrlist)
        preExp = self._decodeGenPreExp(layers,0)
        return preExp.split('.')

    def _decodeSplitLayers(self, chrlist) :
        layers =[]
        layer = [chrlist[0]]
        layers.append(layer)
        i = 1
        numArg = self._getNumArg(chrlist[0])
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
            ret += '.'
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
        totalLen = len(self._funcList)+self._datadimension-1;
        r = self._rand.randint(0, totalLen-1)
        if r < len(self._funcList) :
            return self._funcNameList[r]
        else:
            return self._randGenTerminal()

    def _randGenTerminal(self):
        totalLen = self._datadimension-1+len(self._constantList)
        r = self._rand.randint(0, totalLen-1 )
        if r < totalLen*self._rateConstant and len(self._constantList) > 0:
            return '#'+str(r%len(self._constantList))
        else:
            return 'x'+str(r%(self._datadimension-1)) 
    
    def _select_replicate(self):
        fitkeys = self._fitnesses.keys()
        fitkeys.sort()
        count = len(fitkeys)
        newchromosomeList = []
        self._bestChromosome = self._fitnesses[fitkeys[-1]][0]
        self._bestFitness = fitkeys[-1]
        for i in range(len(self._chromosomeList)-1) :
            r = self._rand.randint(1, count*(count+1)/2-1)
            idx = int(math.sqrt(r)) -1
            chromosome = self._fitnesses[fitkeys[idx]][0]
            newchromosomeList.append(chromosome)
        self._chromosomeList = newchromosomeList

    def _mutation(self):
        index = 0
        newchromosomeList = []
        while index < len(self._chromosomeList):
            if self._randomProb(self._rateMutation) == False:
                newchromosomeList.append(self._chromosomeList[index])
            else:
                chrlist = self._chromosomeList[index].split('.')
                for time in range(self._lenMutation) :
                    pos = self._rand.randint(0,self._chromosomeLen-1)
                    if pos < self._head:
                        chrlist[pos] = self._randGenFuncAndTerm()
                    else:
                        chrlist[pos] = self._randGenTerminal()
                newchromosomeList.append('.'.join(chrlist)) 
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
                chrlist1 = self._chromosomeList[index].split('.')
                chrlist2 = self._chromosomeList[index+1].split('.')
                pos = self._rand.randint(1, self._chromosomeLen-2)
                newchrlist1 = chrlist1[:pos]+chrlist2[pos:]
                newchrlist2 = chrlist2[:pos]+chrlist1[pos:]
                newchromosomeList.append('.'.join(newchrlist1))
                newchromosomeList.append('.'.join(newchrlist2))
            index += 2
        while index < len(self._chromosomeList) :
            newchromosomeList.append(self._chromosomeList[index])
            index += 1
        self._chromosomeList = newchromosomeList
        
    def _twoPointsRecombination(self):
        index = 0
        newchromosomeList =[]
        while index < len(self._chromosomeList)-1 :
            if self._randomProb(self._rateTwoRecombine) == False:
                newchromosomeList.append(self._chromosomeList[index])
                newchromosomeList.append(self._chromosomeList[index+1])
            else:
                chrlist1 = self._chromosomeList[index].split('.')
                chrlist2 = self._chromosomeList[index+1].split('.')
                pos1 = self._rand.randint(1, self._chromosomeLen-2)
                pos2 = self._rand.randint(1, self._chromosomeLen-2)
                pos1 = min(pos1,pos2)
                pos2 = max(pos1,pos2)                
                newchrlist1 = chrlist1[:pos1]+chrlist2[pos1:pos2]+chrlist1[pos2:]
                newchrlist2 = chrlist2[:pos1]+chrlist1[pos1:pos2]+chrlist2[pos2:]
                newchromosomeList.append('.'.join(newchrlist1))
                newchromosomeList.append('.'.join(newchrlist2))
            index += 2
        while index < len(self._chromosomeList) :
            newchromosomeList.append(self._chromosomeList[index])
            index += 1
        self._chromosomeList = newchromosomeList
    
    def _inversion(self):
        index = 0
        newchromosomeList = []
        while index < len(self._chromosomeList):
            if self._randomProb(self._rateInversion) == False:
                newchromosomeList.append(self._chromosomeList[index])
            else:
                pos = self._rand.randint(0,self._chromosomeLen-self._lenInversion-1)
                inverselen=self._lenInversion
                chrlist = self._chromosomeList[index].split('.')
                if pos >= self._head-inverselen and pos < self._head :
                    if pos >= self._head-inverselen/2 :
                        pos=self._head
                    else:
                        pos=self._head-inverselen-1
                for i in range(inverselen/2) : # Reverse the fragment of chromosome
                    t = chrlist[pos+i]
                    chrlist[pos+i] = chrlist[pos+inverselen-i]
                    chrlist[pos+inverselen-i]=t
                newchromosomeList.append('.'.join(chrlist)) 
            index += 1
        self._chromosomeList = newchromosomeList
    
    def _transIS(self):
        index = 0
        newchromosomeList = []
        while index < len(self._chromosomeList):
            if self._randomProb(self._rateIS) == False:
                newchromosomeList.append(self._chromosomeList[index])
            else:
                srcpos = self._rand.randint(0,self._chromosomeLen-self._lenTransposition-1)
                destpos = self._rand.randint(1,self._head-1)
                tranlen = self._lenTransposition
                chrlist = self._chromosomeList[index].split('.')
                for i in range(tranlen) :
                    chrlist.insert(destpos+i,chrlist[srcpos+i])
                    if destpos <= srcpos:
                        srcpos-=1
                del chrlist[self._head:self._head+tranlen]
                newchromosomeList.append('.'.join(chrlist))
            index += 1
        self._chromosomeList = newchromosomeList
    
    def _transRIS(self):
        index = 0
        newchromosomeList = []
        while index < len(self._chromosomeList):
            if self._randomProb(self._rateRIS) == False:
                newchromosomeList.append(self._chromosomeList[index])
            else:
                srcpos = self._rand.randint(0,self._head-1)
                tranlen= self._lenTransposition
                chrlist = self._chromosomeList[index].split('.')
                while srcpos < self._head :
                    if chrlist[srcpos] in self._funcNameList:
                        break
                    srcpos += 1
                if srcpos < self._head :
                    for i in range(tranlen) :
                        chrlist.insert(i,chrlist[srcpos+i])
                        srcpos-=1
                    del chrlist[self._head:self._head+tranlen]
                newchromosomeList.append('.'.join(chrlist))
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
    return args[0]+args[1]

def sub(args):
    return args[0]-args[1]

def mul(args):
    return args[0]*args[1]

def div(args):
    if args[1] == 0:
        return 0
    else :
        return args[0]/args[1]

if __name__ == '__main__' :
    myGep = freeGEP()
    myGep.loadfromfile('dataset.txt')
    myGep.addFunc(add,'+',2)
    myGep.addFunc(sub,'-',2)
    myGep.addFunc(mul,'*',2)
    myGep.addFunc(div,'/',2)
    myGep.initChromosomes(9,100)
    myGep.printChromosomes()
    myGep.openLogFile('log.txt')
    #myGep._decode('++ab+aabbabaa')
    print '\n'
    myGep.evolution(200)
    print '\n'
    myGep.closeLogFile()
    print myGep._bestChromosome + '  '+str(myGep._bestFitness)

