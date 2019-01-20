package expoimpl;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class ExponentialMechanism {
	
	// with BigDecimal
	public static String getRandomCandidateBD(List<String> candidates, List<BigDecimal> probabilities) {
		BigDecimal sum = new BigDecimal(0.0);
		for (int i = 0; i < probabilities.size(); i++) {
			sum = sum.add(probabilities.get(i));
		}
		List<BigDecimal> normalizedProbabilities = new ArrayList<BigDecimal>();
		for (int i = 0; i < probabilities.size(); i++) {
			normalizedProbabilities.add(probabilities.get(i).divide(sum, MathContext.DECIMAL128));
		}
		Random r = new Random();
		double random = r.nextDouble();
		BigDecimal seenSoFar = new BigDecimal(0.0);
		for (int i = 0; i < normalizedProbabilities.size(); i++) {
			seenSoFar = seenSoFar.add(normalizedProbabilities.get(i));
			if (seenSoFar.compareTo(new BigDecimal(random)) >= 0) {
				return candidates.get(i);
			}
		}
		//System.out.println("expo mech returning null");
		//System.out.println(candidates);
		//System.out.println(probabilities);
		//System.out.println(normalizedProbabilities);
		//System.out.println(seenSoFar);
		return null;
	}
	
	
	// with double
	public static String getRandomCandidate(List<String> candidates, List<Double> probabilities) {
		// find sum of probabilities, normalize to [0,1]
		double sum = 0.0;
		for (int i = 0; i < probabilities.size(); i++) {
			sum = sum + probabilities.get(i);
		}
		List<Double> normalizedProbabilities = new ArrayList<Double>();
		for (int i = 0; i < probabilities.size(); i++) {
			normalizedProbabilities.add(((double)probabilities.get(i)/sum));
		}
		Random r = new Random();
		double random = r.nextDouble();
		double seenSoFar = 0.0;
		for (int i = 0; i < normalizedProbabilities.size(); i++) {
			seenSoFar = seenSoFar + normalizedProbabilities.get(i);
			if (seenSoFar >= random) {
				return candidates.get(i);
			}
		}
		//System.out.println("expo mech returning null");
		//System.out.println(candidates);
		//System.out.println(probabilities);
		//System.out.println(normalizedProbabilities);
		//System.out.println(seenSoFar);
		return null;
	}

}
