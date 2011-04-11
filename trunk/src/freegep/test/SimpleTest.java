package freegep.test;
import java.util.List;
import freegep.core.*;

public class SimpleTest {
	
	private static class TestListener implements FreeGEP.GEPListener {
		private FreeGEP _gep = null;

		public TestListener(FreeGEP gep) {
			_gep = gep;
		}
		
		public void listen(List<String> chromosomeList, float bestfitness,
				String bestChromosome) {
			// TODO Auto-generated method stub
			System.out.println("translated: "+_gep.translateChromosome(bestChromosome));
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
		    //myGep.setConstants(new float[]{0.2f,2f});
		    myGep._lenInversion = 2;
		    myGep._lenTransposition = 2;
		    myGep._rateMutation = 0.45f;
		    myGep.initChromosomes(7,200);
		    myGep.addListener(new SimpleTest.TestListener(myGep));
		    myGep.evolution(700);
		}catch(Exception e) {
			e.printStackTrace();
		}

	}

}
