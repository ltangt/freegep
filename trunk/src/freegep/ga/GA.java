package freegep.ga;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class GA {
	
	private Random random = new Random();
	
	private int generation = 5;
	private int individualLength = 11;	//	the length of each individual
	
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
		for(int i = 0; i < individualSize; ++i){
			individuals[i] = new BitSet(individualLength);
		}
		initIndividual();
		
		trainingSet = new ArrayList<BitSet>();
	}
	
	public void setGeneration(int gen){
		generation = gen;
	}
	
	private void initIndividual(){
		for(int i = 0; i < individualSize; ++i){
			for(int j = 0; j < individualLength; ++j){
				float ran = random.nextFloat();
				if(ran < 0.5){
					individuals[i].set(j);
				}
			}
		}
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
			System.out.println("Generation " + gen);
			float[] cumuProbs = new float[individualSize];	//	cumulative probabilities
			cumuProbs[0] = probs[0];
			for(int i = 1; i < individualSize; ++i){
				cumuProbs[i] = cumuProbs[i - 1] + probs[i];
			}

			//	directly put the best one into next generation
			float maxFitProb = 0.0f;
			for(int i = 0; i < probs.length; ++i){
				if(probs[i] > maxFitProb){
					newIndividuaals[0] = individuals[i];
					maxFitProb = probs[i];
				}
			}
			
			//	select (1-probSelect)*individualSize
			int portion = (int)((1 - fraction) * individualSize);
			int tournamentSize = 3;
			for(int i = 1; i < portion; ++i){
				int winner = random.nextInt(individualSize);
				float fitness = fitness(individuals[winner]);
				for(int j = 0; j < tournamentSize; ++j){
					int pker = random.nextInt(individualSize);
					float pkerFit = fitness(individuals[pker]);
					if(pkerFit > fitness){
						winner = pker;
						fitness = pkerFit;
					}
				}
				newIndividuaals[i] = (BitSet)individuals[winner].clone();
//				System.out.println("Winner is:" + newIndividuaals[i]);
			}
			
//			for(int i = 0; i < portion; ++i){
//				float ran = random.nextFloat();
//				for(int j = 0; j < individualSize; ++j){
//					if(ran < cumuProbs[j]){
//						newIndividuaals[i] = (BitSet)individuals[j].clone();
//						break;
//					}
//				}
//			}
			
			
			//	put single crossover individuals into new generation
			for(int i = portion; i < portion + fraction * individualSize / 2; ++i){
				BitSet individual1 = (BitSet)individuals[0].clone();
				BitSet individual2 = (BitSet)individuals[0].clone();
				
				float ran1 = random.nextFloat();
				float ran2 = random.nextFloat();
				for(int j = 0; j < individualSize; ++j){
					if(ran1 < cumuProbs[j]){
						individual1 = (BitSet)individuals[j].clone();
//						System.out.println("Crossover select " + j + "th individual for individual1");
						break;
					}
				}
				for(int j = 0; j < individualSize; ++j){
					if(ran2 < cumuProbs[j]){
						individual2 = (BitSet)individuals[j].clone();
//						System.out.println("Crossover select " + j + "th individual for individual2");
						break;
					}
				}
				
				newIndividuaals[i] = singleCrossover(individual1, individual2);
			}
			
			//	put double crossover individuals into new generation
			int doubleCrossoverStart = (int)(portion + fraction * individualSize / 2);
			for(int i = doubleCrossoverStart; i < individualSize; ++i){
				BitSet individual1 = (BitSet)individuals[0].clone();
				BitSet individual2 = (BitSet)individuals[0].clone();
				
				float ran1 = random.nextFloat();
				float ran2 = random.nextFloat();
				for(int j = 0; j < individualSize; ++j){
					if(ran1 < cumuProbs[j]){
						individual1 = (BitSet)individuals[j].clone();
						break;
					}
				}
				
				for(int j = 0; j < individualSize; ++j){
					if(ran2 < cumuProbs[j]){
						individual2 = (BitSet)individuals[j].clone();
						break;
					}
				}
				
				newIndividuaals[i] = doubleCrossover(individual1, individual2);
			}
			
			//	put mutated individuals into new generations
			for(int i = 0; i < individualSize; ++i){
				float ran = random.nextFloat();
				if(ran < probMutation){
					newIndividuaals[i] = mutation(newIndividuaals[i]);
				}
			}
			
			individuals = newIndividuaals;
			
		}
		
	}
	
	private void report(float[] fitness){
		int idx = 0;
		float max = 0.0f;
		for(int i = 0; i < fitness.length; ++i){
			if(fitness[i] > max){
				max = fitness[i];
				idx = i;
			}
		}
		System.out.println("Max fitness is:" + max);
	}
	
	private float[] calculateFitnessProb(BitSet[] individuals){
		float[] probs = new float[individualSize];
		Arrays.fill(probs, 0.0f);
		
		float maxFit = 0.0f;
		int idx = 0;
		
		float sumFitness = 0.0f;
		for(int i = 0; i < individuals.length; ++i){
			BitSet individual = individuals[i];
			probs[i] = fitness(individual);
			if(probs[i] > maxFit){
				maxFit = probs[i];
				idx = i;
			}
			sumFitness += probs[i];
		}
		
		for(int i = 0; i < probs.length; ++i){
			probs[i] /= sumFitness;
		}
		
		System.out.println("Max fitness:" + maxFit + " with bitstring:" + individuals[idx]);
		
		return probs;
	}
	
	/*
	 * Mutation operation.
	 */
	private BitSet mutation(BitSet individual){
		int pos = random.nextInt(individualLength);
		BitSet clone = (BitSet)individual.clone();
		boolean value = clone.get(pos);
		clone.set(pos, !value);
		return clone;
	}
	
	/*
	 * Single crossover operation.
	 */
	private BitSet singleCrossover(BitSet individual1, BitSet individual2){
		int pos = random.nextInt(individualLength);
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
		int posStart = random.nextInt(individualLength);
		int posEnd = random.nextInt(individualLength);
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
		float fitness = 0.0f;
		
		for(BitSet instance : trainingSet){
			boolean classLabel = instance.get(individualLength - 1);
			if(true == classLabel){	//	individual accept positive class
				BitSet clone = (BitSet)individual.clone();
				clone.or(instance);
				if(true == individual.equals(clone)){
					++fitness;
				}
			}
			else{	//	individual reject negative class
				BitSet clone = (BitSet)individual.clone();
				clone.or(instance);
				if(false == individual.equals(clone)){
					++fitness;
				}
			}
		}
		fitness /= trainingSet.size();
//		System.out.println("\t" + individual + "\tFitness:" + fitness);
		
		return fitness;
	}

	public static void main(String[] args) throws Exception{
		int size = 100;
		int gen = 100;
		
		float probS = 0.11f;
		float probM = 0.18f;
		float probSC = 0.10f;
		float probDC = 0.12f;
		
		String filepath = "./traindata";
		
		GA ga = new GA(size, probS, probM, probSC, probDC);
		ga.setGeneration(gen);
		ga.readTrainingData(filepath);
		ga.training();
	}
}
