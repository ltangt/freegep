package freegep.ga;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Random;

public class GA {
	
	private Random random = new Random();
	
	private int generation = 100;
	private int individualLength = 10;	//	the length of each individual
	
	private int individualSize;
	
	private float fraction;
	private float probMutation;
	private float probSingleCrossover;
	private float probDoubleCrossover;
	
	private BitSet[] individuals;
	
	private List<BitSet> trainingSet;
	
	public GA(int size, float select, float mutation, float singleCross, float doubleCross){
		individualSize =  size;
		fraction = select;
		probMutation = mutation;
		probSingleCrossover = singleCross;
		probDoubleCrossover = doubleCross;
		
		individuals = new BitSet[individualSize];
		for(BitSet bs : individuals){
			bs = new BitSet(individualLength);
		}
		
		trainingSet = new ArrayList<BitSet>();
	}
	
	/*
	 * Read training data.
	 */
	public void readTrainingData(String filepath) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filepath));
		String line = null;
		while((line = br.readLine()) != null){
			if(line.length() != individualLength){
				throw new Error("invalid format for training data.");
			}
			else{
				BitSet individual = new BitSet(individualLength);
				for(int i = 0; i < line.length(); ++i){
					if(line.charAt(i) == '1'){
						individual.set(i);
					}
				}
				trainingSet.add(individual);
			}
		}
	}
	
	public void training(){
		if(trainingSet.size() == 0){
			throw new Error("No training data.");
		}
		
		int gen = 0;
		BitSet[] newIndividuaals = new BitSet[individualSize];
		
		while(gen++ < generation){
			//	calculate probability based on fitness
			float[] probs = calculateFitnessProb(individuals);
			
			//	select (1-probSelect)*individualSize
			int portion = (int)(1 - fraction) * individualSize;
			for(int i = 0; i < portion; ++i){
				for(int j = 0; j < individualSize; ++j){
					float rate = random.nextFloat();
					if(rate < probs[j]){
						newIndividuaals[i] = (BitSet)individuals[j].clone();
						break;
					}
				}
			}
			
			//	put mutated individuals into new generations
			for(int i = 0; i < fraction * individualSize / 2; ++i){
				
			}
			
			//	put single crossover individuals into new generation
			
			//	put double crossover individuals into new generation
			
		}
		
	}
	
	private float[] calculateFitnessProb(BitSet[] individuals){
		float[] probs = new float[individualSize];
		Arrays.fill(probs, 0.0f);
		
		float sumFitness = 0.0f;
		for(int i = 0; i < individuals.length; ++i){
			BitSet individual = individuals[i];
			probs[i] = fitness(individual);
			sumFitness += probs[i];
		}
		
		for(int i = 0; i < probs.length; ++i){
			probs[i] /= sumFitness;
		}
		
		return probs;
	}
	
	/*
	 * Mutation operation.
	 */
	private BitSet mutation(BitSet individual){
		int pos = random.nextInt() % individualLength;
		BitSet clone = (BitSet)individual.clone();
		boolean value = clone.get(pos);
		clone.set(pos, !value);
		return clone;
	}
	
	/*
	 * Single crossover operation.
	 */
	private BitSet singleCrossover(BitSet individual1, BitSet individual2){
		int pos = random.nextInt() % individualLength;
		BitSet clone = (BitSet)individual1.clone();
		for(int i = pos; i < individualLength; ++i){
			clone.set(i, individual2.get(i));
		}
		return clone;
	}
	
	/*
	 * Double crossover operation.
	 */
	private BitSet doubleCrossover(BitSet individual1, BitSet individual2){
		int posStart = random.nextInt() % individualLength;
		int posEnd = random.nextInt() % individualLength;
		if(posStart > posEnd){
			int tmp = posStart;
			posStart = posEnd;
			posEnd = tmp;
		}
		BitSet clone = (BitSet)individual1.clone();
		for(int i = posStart; i < posEnd; ++i){
			clone.set(i, individual2.get(i));
		}
		return clone;
	}
	
	/*
	 * Calculate the fitness of a single individual.
	 */
	private float fitness(BitSet individual){
		float fitness = 0;
		
		for(BitSet instance : trainingSet){
			BitSet clone = (BitSet)individual.clone();
			clone.or(instance);
			if(individual.equals(clone)){
				++fitness;
			}
		}
		
		return fitness / trainingSet.size();
	}

	public static void main(String[] args) throws Exception{
		int size = 100;
		
		float probS = 0.11f;
		float probM = 0.02f;
		float probSC = 0.10f;
		float probDC = 0.12f;
		
		String filepath = "./traindata";
		
		GA ga = new GA(size, probS, probM, probSC, probDC);
		ga.readTrainingData(filepath);
	}
}
