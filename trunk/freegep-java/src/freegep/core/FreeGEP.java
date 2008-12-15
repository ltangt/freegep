package freegep.core;
import java.util.*;
import java.io.*;


public class FreeGEP {
	private Random _rand = null;
	private List<float[]> _dataset = null;
	private List<String> _funcNameList = null;
	private Map<String,GEPFunc> _funcMap = null;
	private Map<String,Integer> _funcNumArgs= null;
	private List<String> _chromosomeList = null;
	private List<ArrayList<String>> _preExpList = null;
	private TreeMap<Float, ArrayList<String>> _fitnesses = null;
	private float[] _constants = null;
	private int _datadimension = 0;
	private int _maxNumArgsOfFunc=0;
	private int _head=0;
	private int _tail=0;
	private int _population=0;
	private int _maxGeneration = 0;
	private int _currentGeneration = 0;
	private int _chromosomeLen = 0;
	private String _bestChromosome ="";
	private float _bestFitness= -1f;
	private GEPEvalFunc _evalfunc = null;
	private List<GEPListener> _listeners= null;
	public float _M = 1000.0f;
	public float _rateConstant = 0.3f;
	public float _rateMutation = 0.3f;
	public float _rateOneRecombine=0.4f;
	public float _rateTwoRecombine=0.2f;
	public float _rateInversion = 0.2f;
	public float _rateIS = 0.1f;
	public float _rateRIS = 0.1f;
	public int _lenMutation = 2;
	public int _lenInversion = 4;
	public int _lenTransposition = 4;
	
	
	public interface GEPFunc {
		float cal(float[] args);
	}
	
	public interface GEPEvalFunc {
		float eval(int chrindex, float[] expectvals, float[] actualvals);
	}
	
	public interface GEPListener {
		void listen(List<String> chromosomeList, float bestfitness, String bestChromosome);
	}
	
	public FreeGEP() {
		_rand = new java.util.Random();
		_rand.setSeed(System.currentTimeMillis());
		_dataset = new ArrayList<float[]>();
		_funcNameList = new ArrayList<String>();
		_funcMap = new HashMap<String,GEPFunc>();
		_funcNumArgs = new HashMap<String, Integer>();
		_chromosomeList = new ArrayList<String>();
		_preExpList = new ArrayList<ArrayList<String>>();
		_fitnesses = new TreeMap<Float, ArrayList<String>>();
		_listeners = new ArrayList<GEPListener>();
	}
	
	public boolean loadfromfile(String filename) throws IOException {
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
		String line = null;
		float[] rowvalues = null;
		List<float[]> dataset = new ArrayList<float[]>();
		while ((line=reader.readLine()) != null) {
			String[] rowset=line.split(" ");
			if (rowvalues == null) {
				rowvalues = new float[rowset.length];
			}
			for (int i=0;i<rowvalues.length;i++) {
				rowvalues[i] = Float.parseFloat(rowset[i]);
			}
			dataset.add(rowvalues);
		}
		reader.close();
		fis.close();
		return load(dataset);
	}
	
	public boolean load(List<float[]> dataset) {
		_dataset = dataset;
		if (_dataset.size() > 0) {
			_datadimension = _dataset.get(0).length;
		}
		return true;
	}
	
	public void addFunc(GEPFunc func, String funcName, int numArgs) {
		if (_maxNumArgsOfFunc < numArgs)
			_maxNumArgsOfFunc = numArgs;
		_funcNameList.add(funcName);
		_funcMap.put(funcName, func);
		_funcNumArgs.put(funcName, numArgs);
	}
	
	public void setConstants(float[] values) {
		if (values == null)
			_constants = null;
		else {
			_constants = new float[values.length];
			System.arraycopy(values, 0, _constants, 0, values.length);
		}
	}
	
	public void addListener(GEPListener listener) {
        _listeners.add(listener);
	}
        
    public void setEvalfunc(GEPEvalFunc evalfunc) {
        _evalfunc=evalfunc;
    }
    
    public int getPopulation() {
    	return _population;
    }
    
    public int getCurrentGeneration() {
    	return _currentGeneration;
    }
    
    public int getMaxGeneration() {
    	return _maxGeneration;
    }
    
    public String getChromosome(int chrindex) {
    	return _chromosomeList.get(chrindex);
    }
    
    public String getBestChromosome() {
    	return _bestChromosome;
    }
    
    public float getBestFitness() {
    	return _bestFitness;
    }
    
	public void initChromosomes(int head, int num) {
		_head = head;
		_tail = head * (_maxNumArgsOfFunc - 1) + 1;
		_chromosomeLen = _head + _tail;
		_chromosomeList.clear();
		for (int i = 0; i < num; i++) {
			String chromosome = "";
			for (int h = 0; h < _head; h++) {
				chromosome += randGenFuncAndTerm();
				chromosome += '.';
			}
			for (int t = 0; t < _tail; t++) {
				chromosome += randGenTerminal();
				if (t < _tail - 1)
					chromosome += '.';
			}
			_chromosomeList.add(chromosome);
		}
		_population = num;
	}
	
	public void evolution(int maxGeneration) throws Exception {
		_maxGeneration = maxGeneration;
		checkParameters();
        for (_currentGeneration=0; _currentGeneration<_maxGeneration; _currentGeneration++) {
            decode();
            evaluate();
            select_replicate();
            mutation();
            inversion();
            transIS();
            transRIS();
            onePointRecombination();
            twoPointsRecombination();
            _chromosomeList.add(_bestChromosome);
            System.out.println("best "+_bestChromosome+" "+_bestFitness);
            System.out.println("generating "+_currentGeneration);
            if (_listeners != null && _listeners.size() > 0) {                
            	int index;
            	for (index=0; index<_listeners.size(); index++){
                    GEPListener listener = _listeners.get(index);
            		listener.listen(_chromosomeList, _bestFitness, _bestChromosome);
            	}
            }
            if (_bestFitness == _M)
                break;
        }
	}
	
	private void decode() {
		_preExpList.clear();
		for (int i=0;i<_chromosomeList.size(); i++) {
			_preExpList.add(decode(_chromosomeList.get(i)));
        }
	}
	
	private void checkParameters() throws Exception {
		if (_lenInversion > _head -1) 
			throw new Exception("the length of inversion is larger than the head!");
		if (_lenMutation > _chromosomeLen)
			throw new Exception("the length of mutation is larger than the length of chromosome");
		if (_lenTransposition > _head)
			throw new Exception("the length of transposition is larger than the head!");
	}
	
	private ArrayList<String> decode(String chromosome) {
        String[] chrlist = chromosome.split("[.]");
        List<ArrayList<String>> layers = decodeSplitLayers(chrlist);
        String preExp = decodeGenPreExp(layers,0);
        String[] preExpElems= preExp.split("[.]");
        ArrayList<String> preExpList = new ArrayList<String>(preExpElems.length);
        for (int i=0;i<preExpElems.length;i++)
        	preExpList.add(preExpElems[i]);
        return preExpList;
	}
	
	private List<ArrayList<String>> decodeSplitLayers(String[] chrlist) {
		List<ArrayList<String>> layers= new ArrayList<ArrayList<String>>();
        ArrayList<String> layer = new ArrayList<String>();
        layer.add(chrlist[0]);
        layers.add(layer);
        int i = 1;
        int numArg = getNumArg(chrlist[0]);
        while (numArg > 0) {
            layer = new ArrayList<String>();
            int nextNumArg = 0;
            int endIndex = i + numArg;
            while (i < endIndex) {
                nextNumArg += getNumArg(chrlist[i]);
                layer.add(chrlist[i]);
                i += 1;
            }
            numArg = nextNumArg;
            layers.add(layer);
        }
        return layers;
	}

    public String decodeGenPreExp(List<ArrayList<String>> layers, int layerindex) {
        if (layerindex >= layers.size())
            return "";
        ArrayList<String> layer = layers.get(layerindex);
        if (layer.size() == 0)
            return "";
        String ret = layer.get(0);
        for (int i=0; i<getNumArg(layer.get(0)); i++) {
            ret += '.';
            ret += decodeGenPreExp(layers, layerindex+1);
        }
        layer.remove(0);
        return ret;
    }
	
	public String translateChromosome(String chromosome) {
        String[] chrlist = chromosome.split("[.]");
        List<ArrayList<String>> layers = decodeSplitLayers(chrlist);
        return translateGenMidExp(layers,0);
	}
	
	public float evalute(String chromosome, List<float[]> dataset) {
		float fit = 0.0f;
        float err;
        float[] expects = new float[dataset.size()];
        float[] targets = new float[dataset.size()];
		ArrayList<String> exp = decode(chromosome);
        for (int rowj=0; rowj<dataset.size(); rowj++) {                
        	ArrayList<String> expCopy = new ArrayList<String>(exp);
        	float[] row = dataset.get(rowj);
        	targets[rowj] = row[row.length-1];
            expects[rowj] = evaluate(expCopy, row);
        }
        if (_evalfunc == null) {
        	for (int rowj=0; rowj<dataset.size(); rowj++) {
        		float[] row = dataset.get(rowj);
            	float target = row[row.length-1];
            	float y = expects[rowj];
                err = Math.abs(y-target);
                fit += _M - err;
        	}
            fit /= (float)(dataset.size());
        }else {
        	fit = _evalfunc.eval(0, expects, targets);
        }
        return fit;
	}
	
	private String translateGenMidExp(List<ArrayList<String>> layers, int layerindex) {
        if (layerindex >= layers.size())
            return "";
        ArrayList<String> layer = layers.get(layerindex);
        if (layer.size() == 0)
            return "";
        String midterm = layer.get(0);
        int numarg= getNumArg(layer.get(0));
        StringBuffer ret = new StringBuffer();
        if (numarg == 2) {
            ret.append('(');
            ret.append(translateGenMidExp(layers, layerindex+1));
            ret.append(midterm);
            ret.append(translateGenMidExp(layers, layerindex+1));
            ret.append(')');
        }else if (numarg == 1) {
            ret.append('(');
            ret.append(midterm);
            ret.append(translateGenMidExp(layers, layerindex+1));
            ret.append(')');
        }else {
            if (midterm.charAt(0) == '#')
                ret.append(""+_constants[Integer.parseInt(midterm.substring(1))]);
            else
                ret.append(midterm);
        }
        layer.remove(0);
        return ret.toString();
	}

    private void evaluate() {
        _fitnesses.clear();
        float[] expects = new float[_dataset.size()];
        float[] targets = new float[_dataset.size()];
        for (int expIndex=0; expIndex<_preExpList.size(); expIndex++) {
            float fit = 0.0f;
            float err;
            ArrayList<String> exp=_preExpList.get(expIndex);
            for (int rowj=0; rowj<_dataset.size(); rowj++) {                
            	ArrayList<String> expCopy = new ArrayList<String>(exp);
            	float[] row = _dataset.get(rowj);
            	targets[rowj] = row[row.length-1];
                expects[rowj] = evaluate(expCopy, row);
            }
            if (_evalfunc == null) {
            	for (int rowj=0; rowj<_dataset.size(); rowj++) {
            		float[] row = _dataset.get(rowj);
                	float target = row[row.length-1];
                	float y = expects[rowj];
                    err = Math.abs(y-target);
                    fit += _M - err;
            	}
                fit /= (float)(_dataset.size());
            }else {
            	fit = _evalfunc.eval(expIndex, expects, targets);
            }
            String chromosome = _chromosomeList.get(expIndex);
            if (_fitnesses.containsKey(fit) == false) 
                _fitnesses.put(fit, new ArrayList<String>());
            _fitnesses.get(fit).add(chromosome);
        }
    }
	
    private float evaluate(ArrayList<String> exp, float[] row) {
        String op = exp.get(0);
        exp.remove(0);
        if (_funcMap.containsKey(op)) {
        	GEPFunc func = _funcMap.get(op);
        	int funcNumArg = _funcNumArgs.get(op);
        	float[] args = new float[funcNumArg];
        	for (int i=0;i<funcNumArg; i++)
        		args[i] = evaluate(exp, row);
        	return func.cal(args);
        }else if (op.charAt(0) == '#') {
        	int constIndex = Integer.parseInt(op.substring(1));
        	return _constants[constIndex];
        }else {
        	int termIndex = Integer.parseInt(op.substring(1));
        	return row[termIndex];
        }
    }
    
	private String randGenFuncAndTerm() {
        int totalLen = _funcNameList.size()+_datadimension-1;
        int r = rand(0,totalLen-1);
        if (r < _funcNameList.size())
            return _funcNameList.get(r);
        else
            return randGenTerminal();
	}

    private String randGenTerminal() {
    	int totalLen=0;
    	if (_constants == null)
    		totalLen = _datadimension-1;
    	else
    		totalLen = _datadimension-1+_constants.length;
        int r = rand(0,totalLen);
        if (r < totalLen*_rateConstant && _constants != null && _constants.length > 0)
            return "#"+(r%_constants.length);
        else
            return "x"+r%(_datadimension-1); 
    }
    
	private int getNumArg(String op) {
        if (_funcNumArgs.containsKey(op)) 
        	return _funcNumArgs.get(op);
        else
            return 0;
	}
	
	private void select_replicate() {
		Set<Float> fitkeys=_fitnesses.keySet();
		Float[] keys = new Float[fitkeys.size()];
		fitkeys.toArray(keys);
		float lastKey=_fitnesses.lastKey();
        int count = fitkeys.size();
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        _bestChromosome = _fitnesses.get(lastKey).get(0);
        _bestFitness = lastKey;
        for (int i=0; i<_chromosomeList.size()-1; i++) {
            int r = rand(1,count*(count+1)/2-1);
            int idx = (int)(Math.sqrt(r)) -1;
            int subidx = idx % _fitnesses.get(keys[idx]).size();
            String chromosome = _fitnesses.get(keys[idx]).get(subidx);
            newchromosomeList.add(chromosome);
        }
        _chromosomeList = newchromosomeList;        
	}
	
	private void mutation() {
        int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        while (index < _chromosomeList.size()) {
            if (randomProb(_rateMutation) == false)
                newchromosomeList.add(_chromosomeList.get(index));
            else {
                String[] chrlist = _chromosomeList.get(index).split("[.]");
                for (int time =0; time< _lenMutation; time++) {
                    int pos = rand(0,_chromosomeLen-1);
                    if (pos < _head)
                        chrlist[pos] = randGenFuncAndTerm();
                    else
                        chrlist[pos] = randGenTerminal();
                }
                newchromosomeList.add(joinString('.',chrlist)); 
            }
            index += 1;
        }
        _chromosomeList = newchromosomeList;
	}
	
	private void onePointRecombination() {
        int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        String[] templist = new String[_chromosomeLen];
        while (index < _chromosomeList.size()-1) {
            if (randomProb(_rateOneRecombine) == false) {
                newchromosomeList.add(_chromosomeList.get(index));
                newchromosomeList.add(_chromosomeList.get(index+1));
            }
            else {
                String[] chrlist1 = _chromosomeList.get(index).split("[.]");
                String[] chrlist2 = _chromosomeList.get(index+1).split("[.]");
                int pos = rand(1,_chromosomeLen-2);
                int length = chrlist1.length;
                System.arraycopy(chrlist2, 0, templist, 0, length);
                System.arraycopy(chrlist1, pos, chrlist2, pos, length-pos);
                System.arraycopy(templist, pos, chrlist1, pos, length-pos);
                newchromosomeList.add(joinString('.',chrlist1));
                newchromosomeList.add(joinString('.',chrlist2));
            }
            index += 2;
        }
        while (index < _chromosomeList.size()) {
            newchromosomeList.add(_chromosomeList.get(index));
            index += 1;
        }
        _chromosomeList = newchromosomeList;
	}
    
	private void twoPointsRecombination() {
		int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        String[] templist = new String[_chromosomeLen];
        while (index < _chromosomeList.size()-1) {
            if (randomProb(_rateTwoRecombine) == false) {
                newchromosomeList.add(_chromosomeList.get(index));
                newchromosomeList.add(_chromosomeList.get(index+1));
            }
            else {
            	String[] chrlist1 = _chromosomeList.get(index).split("[.]");
                String[] chrlist2 = _chromosomeList.get(index+1).split("[.]");
                int pos1 = rand(1,_chromosomeLen-2);
                int pos2 = rand(1,_chromosomeLen-2);
                pos1 = Math.min(pos1,pos2);
                pos2 = Math.max(pos1,pos2);
                if (pos2==pos1)
                	pos2=pos1+1;
                int length = chrlist1.length;
                System.arraycopy(chrlist2, 0, templist, 0, length);
                System.arraycopy(chrlist1, pos1, chrlist2, pos1, pos2-pos1);
                System.arraycopy(templist, pos1, chrlist1, pos1, pos2-pos1);
                newchromosomeList.add(joinString('.',chrlist1));
                newchromosomeList.add(joinString('.',chrlist2));
            }
            index += 2;
        }
        while (index < _chromosomeList.size()) {
            newchromosomeList.add(_chromosomeList.get(index));
            index += 1;
        }
        _chromosomeList = newchromosomeList;
	}
	
	private void inversion() {
		int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        while (index < _chromosomeList.size()) {
            if (randomProb(_rateInversion)== false)
                newchromosomeList.add(_chromosomeList.get(index));
            else {
                int pos = rand(0,_chromosomeLen-_lenInversion-1);
                int inverselen=_lenInversion;
                String[] chrlist = _chromosomeList.get(index).split("[.]");
                if (pos >= _head-inverselen && pos < _head ) {
                    if (pos >= _head-inverselen/2)
                        pos=_head;
                    else
                        pos=_head-inverselen-1;
                }
                for (int i=0; i<inverselen/2; i++) {
                    String t = chrlist[pos+i];
                    chrlist[pos+i] = chrlist[pos+inverselen-i];
                    chrlist[pos+inverselen-i]=t;
                }
                newchromosomeList.add(joinString('.',chrlist)); 
            }
            index += 1;
        }
        _chromosomeList = newchromosomeList;
	}
	
	private void transIS() {
		int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        String[] templist = new String[_lenTransposition];
        while (index < _chromosomeList.size()) {
            if (randomProb(_rateIS)== false)
                newchromosomeList.add(_chromosomeList.get(index));
            else {
                int srcpos = rand(0,_chromosomeLen-_lenTransposition-1);
                int destpos = rand(1,_head-1);
                int tranlen = _lenTransposition;
                String[] chrlist = _chromosomeList.get(index).split("[.]");
                System.arraycopy(chrlist, srcpos, templist, 0, tranlen);
                if (destpos+tranlen >= _head) {
                	System.arraycopy(templist, 0, chrlist, destpos, _head-destpos);
                }else {
                	for (int i=0; i<_head-destpos-tranlen; i++) {
                    	chrlist[_head-1-i]=chrlist[_head-tranlen-i];
                    }
                	System.arraycopy(templist,0, chrlist, destpos, tranlen);
                }
                newchromosomeList.add(joinString('.',chrlist)); 
            }
            index += 1;
        }
        _chromosomeList = newchromosomeList;
	}
    
    private void transRIS() {
    	int index = 0;
        ArrayList<String> newchromosomeList = new ArrayList<String>();
        String[] templist = new String[_lenTransposition];
        while (index < _chromosomeList.size()) {
            if (randomProb(_rateIS)== false)
                newchromosomeList.add(_chromosomeList.get(index));
            else {
                int srcpos = rand(0,_head-1);
                int tranlen = _lenTransposition;
                String[] chrlist = _chromosomeList.get(index).split("[.]");
                System.arraycopy(chrlist, srcpos, templist, 0, tranlen);
                while (srcpos < _head) {
                    if (_funcMap.containsKey(chrlist[srcpos]))
                        break;
                    srcpos += 1;
                }
                if (srcpos < _head ) {
                	for (int i=0; i<_head-tranlen; i++) {
                    	chrlist[_head-1-i]=chrlist[_head-tranlen-i];
                    }
                	System.arraycopy(templist,0, chrlist, 0, tranlen);
                }
                newchromosomeList.add(joinString('.',chrlist)); 
            }
            index += 1;
        }
        _chromosomeList = newchromosomeList;
    }
	
    private boolean randomProb(float prob) {
        int val = rand(0,1000);
        if (prob <= (float)(val) / 1000.0f)
            return true;
        else
            return false;
    }
        
    private static String joinString(char link, String[] strs) {
    	StringBuffer buffer = new StringBuffer(strs[0]);
        for (int i=1;i<strs.length;i++) {
        	buffer.append(link);
        	buffer.append(strs[i]);
        }
        return buffer.toString();
    }
    
    private int rand(int min, int max) {
    	if (max==min)
    		return max;
    	int val = _rand.nextInt();
    	val = val < 0? -val:val;
    	return min+(val%(max-min));
    }

	
	
    public static class add implements GEPFunc {
		public float cal(float[] args) {
			return args[0]+args[1];
		}
	}
    
    public static class sub implements GEPFunc {
		public float cal(float[] args) {
			return args[0]-args[1];
		}
	}
    
    public static class mul implements GEPFunc {
		public float cal(float[] args) {
			return args[0]*args[1];
		}
	}
    
    public static class div implements GEPFunc {
		public float cal(float[] args) {
			if (args[1] == 0)
				return Float.MAX_VALUE;
			else
				return args[0]/args[1];
		}
	}
    
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			FreeGEP myGep = new FreeGEP();
		    myGep.loadfromfile("dataset.txt");
		    myGep.addFunc(new FreeGEP.add(),"+",2);
		    myGep.addFunc(new FreeGEP.sub(),"-",2);
		    myGep.addFunc(new FreeGEP.mul(),"*",2);
		    myGep.addFunc(new FreeGEP.div(),"/",2);
		    myGep.setConstants(new float[]{0.2f,2f});
		    myGep.initChromosomes(9,100);
		    myGep.evolution(2000);
		}catch(Exception e) {
			e.printStackTrace();
		}
	}

}

