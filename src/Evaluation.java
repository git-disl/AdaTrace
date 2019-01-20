import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;


public class Evaluation {
	
	public static double calculateTripError(List<Trajectory> original, List<Trajectory> 
			synthetic, int gridCellCount) {
		
		// find data boundaries so that an appropriate grid can be laid
		List<Double> boundaries = Main.getDataBoundaries(original);
		double minX = boundaries.get(0);
		double maxX = boundaries.get(1);
		double minY = boundaries.get(2);
		double maxY = boundaries.get(3);
		
		// lay grid on original and synthetic databases
		Grid ug = new Grid(gridCellCount, minX, maxX, minY, maxY);
		List<GridTrajectory> origDBgrid = Main.convertTrajToGridTraj(original, ug, true);
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(synthetic, ug, true);
		
		// get trip distributions (a high epsilon budget of 10000 yields noiseless results)
		TripDistribution origTD = new TripDistribution(origDBgrid, ug, 10000.0);
		TripDistribution synTD = new TripDistribution(synDBgrid, ug, 10000.0);
		
		// calculate JSD between the two trip distributions
		List<Double> origDist = origTD.getTripProbsAsList();
		List<Double> synDist = synTD.getTripProbsAsList();
		return calcJSD(origDist, synDist);
	}
	
	public static double calculateDiameterError(List<Trajectory> original, List<Trajectory> 
			synthetic, int numberOfBuckets) {
		// get diameters of original database
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		List<Double> originalDiameters = new ArrayList<Double>();
		for (Trajectory t : original) {
			double dia = t.getDiameter();
			originalDiameters.add(dia);
			if (dia < min)
				min = dia;
			if (dia > max)
				max = dia;
		}
		// get diameters of synthetic database
		List<Double> syntheticDiameters = new ArrayList<Double>();
		for (Trajectory t : synthetic) {
			double dia = t.getDiameter();
			syntheticDiameters.add(dia);
		}
		// bucketize the diameters to find diameter counts in each bucket
		double bucketSize = (max-min)/numberOfBuckets;
		List<Integer> origDiaCount = new ArrayList<Integer>(); 
		List<Integer> synDiaCount = new ArrayList<Integer>();
		for (int i = 0; i < numberOfBuckets; i++) {
			double currStart = i * bucketSize;
			double currFinish = (i+1) * bucketSize;
			int cnt = 0;
			//System.out.print("[" + currStart + "," + currFinish + "]  ");
			for (Double d : originalDiameters) {
				if (d >= currStart && d < currFinish)
					cnt++;
			}
			origDiaCount.add(cnt);
			cnt = 0;
			for (Double d : syntheticDiameters) {
				if (d >= currStart && d < currFinish) 
					cnt++;
			}
			synDiaCount.add(cnt);
		}
		// find the diameter distribution of original and synthetic database
		List<Double> origDiaProbs = new ArrayList<Double>();
		List<Double> synDiaProbs = new ArrayList<Double>();
		for (int i = 0; i < numberOfBuckets; i++) {
			origDiaProbs.add(((double)origDiaCount.get(i)/((double)originalDiameters.size())));
			synDiaProbs.add(((double)synDiaCount.get(i)/((double)syntheticDiameters.size())));
		}
		// calculate and return JSD
		return calcJSD(origDiaProbs, synDiaProbs);
	}
	
	public static double calculateTravelDistanceError(List<Trajectory> original, 
			List<Trajectory> synthetic, int numberOfBuckets) {
		// get travel lengths of original database
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		List<Double> originalLengths = new ArrayList<Double>();
		for (Trajectory t : original) {
			double dia = t.getDistanceTravelled();
			originalLengths.add(dia);
			if (dia < min)
				min = dia;
			if (dia > max)
				max = dia;
		}
		// get travel lengths of synthetic database
		List<Double> syntheticLengths = new ArrayList<Double>();
		for (Trajectory t : synthetic) {
			double dia = t.getDistanceTravelled();
			syntheticLengths.add(dia);
		}
		// bucketize the lengths to find counts in each bucket
		double bucketSize = (max-min)/numberOfBuckets;
		List<Integer> origCount = new ArrayList<Integer>(); 
		List<Integer> synCount = new ArrayList<Integer>();
		for (int i = 0; i < numberOfBuckets; i++) {
			double currStart = i * bucketSize;
			double currFinish = (i+1) * bucketSize;
			int cnt = 0;
			//System.out.print("[" + currStart + "," + currFinish + "]  ");
			for (Double d : originalLengths) {
				if (d >= currStart && d < currFinish)
					cnt++;
			}
			origCount.add(cnt);
			cnt = 0;
			for (Double d : syntheticLengths) {
				if (d >= currStart && d < currFinish) 
					cnt++;
			}
			synCount.add(cnt);
		}
		// find the distributions
		List<Double> origProbs = new ArrayList<Double>();
		List<Double> synProbs = new ArrayList<Double>();
		for (int i = 0; i < numberOfBuckets; i++) {
			origProbs.add(((double)origCount.get(i)/((double)originalLengths.size())));
			synProbs.add(((double)synCount.get(i)/((double)syntheticLengths.size())));
		}
		// calculate and return JSD
		return calcJSD(origProbs, synProbs);
	}
	
	public static Query[] generateQueries (List<Trajectory> db, int queryCount, String fn) {
		List<Double> boundaries = Main.getDataBoundaries(db);
		double dataMinX = boundaries.get(0);
		double dataMaxX = boundaries.get(1);
		double dataMinY = boundaries.get(2);
		double dataMaxY = boundaries.get(3);
		
		Query[] tbr = new Query[queryCount];
		for (int i = 0; i < queryCount; i++) {
			tbr[i] = new Query(dataMinX, dataMaxX, dataMinY, dataMaxY, fn);
		}
		return tbr;
	}
	
	public static double calculateQueryError(List<Trajectory> originalDB, 
			List<Trajectory> syntheticDB, Query[] queries, double sanityBound) {
		// evaluate queries on the two databases
		int[] actualAnswers = new int[queries.length];
		int[] syntheticAnswers = new int[queries.length];
		for (int i = 0; i < queries.length; i++) {
			Query q = queries[i];
			actualAnswers[i] = q.evaluateQueryOnDatabase(originalDB);
			syntheticAnswers[i] = q.evaluateQueryOnDatabase(syntheticDB);
		}
		//System.out.println(Arrays.toString(actualAnswers));
		//System.out.println(Arrays.toString(syntheticAnswers));
		// measure relative errors using the sanity bound
		double[] error = new double[queries.length];
		for (int i = 0; i < queries.length; i++) {
			double numerator = Math.abs(syntheticAnswers[i]-actualAnswers[i]);
			double denominator = actualAnswers[i];
			if (originalDB.size() * sanityBound > denominator)
				denominator = originalDB.size() * sanityBound;
			error[i] = ((double)numerator)/((double)denominator);
		}
		// measure average relative error
		double sumErr = 0.0;
		for (int i = 0; i < error.length; i++) {
			sumErr = sumErr + error[i];
		}
		return sumErr/((double)error.length);
	}
	
	public static Map<Pattern, Integer> minePatterns (List<GridTrajectory> db,
			int minSize, int maxSize) {
		Map<Pattern, Integer> tbr = new HashMap<Pattern, Integer>();
		for (int currentSize = minSize; currentSize <= maxSize; currentSize++) {
			//System.out.println("Currently mining patterns of size: " + currentSize);
			for (GridTrajectory t : db) {
				// get sub-lists of length currentSize from t. check if we hit t.size
				for (int i = 0; i <= t.getCells().size()-currentSize; i++) {
					Pattern p = new Pattern(t.getCells().subList(i, i+currentSize));
					if (tbr.containsKey(p)) {
						tbr.put(p,  tbr.get(p)+1);
					} else {
						tbr.put(p,  1);
					}
				}
			}
		}
		return tbr;
	}
	
	// from: http://stackoverflow.com/questions/109383/sort-a-mapkey-value-by-values-java
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue( Map<K, V> map ) {
		List<Map.Entry<K, V>> list = new LinkedList<>( map.entrySet() );
		Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
			@Override
			public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
			{
				return ( o2.getValue() ).compareTo( o1.getValue() );  // reverse order 
			}
		} );
		Map<K, V> result = new LinkedHashMap<>();
		for (Map.Entry<K, V> entry : list) {
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
    }
	
	// same function as above, sorts in reverse order: small to large (ascending)
	public static <K, V extends Comparable<? super V>> Map<K, V> sortReverseOrder
			( Map<K, V> map ) {
		List<Map.Entry<K, V>> list = new LinkedList<>( map.entrySet() );
		Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
			@Override
			public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
			{
				return ( o1.getValue() ).compareTo( o2.getValue() );  // reverse order 
			}
		} );
		Map<K, V> result = new LinkedHashMap<>();
		for (Map.Entry<K, V> entry : list) {
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
	private static boolean checkMapSorting(Map<Pattern, Integer> candidate) {
		Integer previous = null;
		for(Map.Entry<Pattern, Integer> entry : candidate.entrySet()) {
            if (entry.getValue() == null) {
            	System.out.println("null obtained in minePatterns function?!!?");
            	return false;
            }
            if (previous != null) {
                if (entry.getValue() > previous) {
                	System.out.println("sorted map is not sorted !!!");
                	return false;
                }
            }
            previous = entry.getValue();
        }
		return true;
	}
	
	
	public static double patternMiningF1Error (Map<Pattern, Integer> originalPatterns,
			int topK, List<Trajectory> syntheticDatabase, Grid grid, 
			int minFPsize, int maxFPsize) {
		// find synthetic patterns
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(
				syntheticDatabase, grid, true);
		Map<Pattern, Integer> synPatterns = Evaluation.minePatterns(synDBgrid, 
				minFPsize, maxFPsize);
		// sort originalPatterns and synthetic patterns maps
		Map<Pattern, Integer> sortedOrigPatterns = Evaluation.sortByValue(originalPatterns);
		Map<Pattern, Integer> sortedSynPatterns = Evaluation.sortByValue(synPatterns);
		// debug: check if these are actually sorted
		if (!Evaluation.checkMapSorting(sortedOrigPatterns))
			System.out.println("Sorting error @ FP mining.");
		if (!Evaluation.checkMapSorting(sortedSynPatterns))
			System.out.println("Sorting error @ FP mining.");
		// find top-k patterns in original and synthetic
		List<Pattern> origTopK = new ArrayList<Pattern>();
		List<Pattern> synTopK = new ArrayList<Pattern>();
		for (Map.Entry<Pattern, Integer> entry : sortedOrigPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
            if (origTopK.size() < topK) {
            	origTopK.add(entry.getKey());
            	//System.out.println("Pattern: " + entry.getKey().toString() + 
            	//		", supp: " + entry.getValue());
            }
            else {
            	break; // out of the for loop
            }
        }
		for (Map.Entry<Pattern, Integer> entry : sortedSynPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
            if (synTopK.size() < topK) {
            	synTopK.add(entry.getKey());
            }
            else {
            	break; // out of the for loop
            }
        }
		return calcF1(origTopK, synTopK);
	}
	
	public static double patternMiningJaccardSim (Map<Pattern, Integer> originalPatterns,
			int topK, List<Trajectory> syntheticDatabase, Grid grid, 
			int minFPsize, int maxFPsize) {
		// find synthetic patterns
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(
				syntheticDatabase, grid, true);
		Map<Pattern, Integer> synPatterns = Evaluation.minePatterns(synDBgrid, 
				minFPsize, maxFPsize);
		// sort originalPatterns and synthetic patterns maps
		Map<Pattern, Integer> sortedOrigPatterns = Evaluation.sortByValue(originalPatterns);
		Map<Pattern, Integer> sortedSynPatterns = Evaluation.sortByValue(synPatterns);
		// debug: check if these are actually sorted
		if (!Evaluation.checkMapSorting(sortedOrigPatterns))
			System.out.println("Sorting error @ FP mining.");
		if (!Evaluation.checkMapSorting(sortedSynPatterns))
			System.out.println("Sorting error @ FP mining.");
		// find top-k patterns in original and synthetic
		List<Pattern> origTopK = new ArrayList<Pattern>();
		List<Pattern> synTopK = new ArrayList<Pattern>();
		for (Map.Entry<Pattern, Integer> entry : sortedOrigPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
            if (origTopK.size() < topK) {
            	origTopK.add(entry.getKey());
            	//System.out.println("Pattern: " + entry.getKey().toString() + 
            	//		", supp: " + entry.getValue());
            }
            else {
            	break; // out of the for loop
            }
        }
		for (Map.Entry<Pattern, Integer> entry : sortedSynPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
            if (synTopK.size() < topK) {
            	synTopK.add(entry.getKey());
            }
            else {
            	break; // out of the for loop
            }
        }
		return JaccardSimilarity(origTopK, synTopK);
	}
	
	public static double patternMiningSupport (Map<Pattern, Integer> originalPatterns,
			int topK, List<Trajectory> syntheticDatabase, Grid grid, 
			int minFPsize, int maxFPsize) {
		// find synthetic patterns
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(
				syntheticDatabase, grid, true);
		Map<Pattern, Integer> synPatterns = Evaluation.minePatterns(synDBgrid, 
				minFPsize, maxFPsize);
		// sort originalPatterns in descending order of support
		Map<Pattern, Integer> sortedOrigPatterns = Evaluation.sortByValue(originalPatterns);
		// debug: check if these are actually sorted
		if (!Evaluation.checkMapSorting(sortedOrigPatterns))
			System.out.println("Sorting error @ FP mining.");
		// find top-k patterns in original
		List<Pattern> origTopK = new ArrayList<Pattern>();
		for (Map.Entry<Pattern, Integer> entry : sortedOrigPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
            if (origTopK.size() < topK) {
            	origTopK.add(entry.getKey());
            	//System.out.println("Pattern: " + entry.getKey().toString() + 
            	//		", supp: " + entry.getValue());
            }
            else {
            	break; // out of the for loop
            }
        }
		double sumErr = 0.0;
		for (int i = 0; i < origTopK.size(); i++) {
			Pattern freqPat = origTopK.get(i);
			int origSupport = originalPatterns.get(freqPat);
			int synSupport = 0;
			if (synPatterns.containsKey(freqPat))
				synSupport = synPatterns.get(freqPat);
			//System.out.println("origSup: " + origSupport + ", synSup: " + synSupport);
			sumErr = sumErr + Math.abs(origSupport-synSupport)/((double)origSupport);
		}
		return sumErr/((double)origTopK.size());
	}
	
	// uses Kendall's tau coefficient for ranking correlation 
	public static double patternMiningRanking (Map<Pattern, Integer> originalPatterns,
			int rankingN, List<Trajectory> syntheticDatabase, Grid grid, 
			int minFPsize, int maxFPsize) {
		// find synthetic patterns
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(
				syntheticDatabase, grid, true);
		Map<Pattern, Integer> synPatterns = Evaluation.minePatterns(synDBgrid, 
				minFPsize, maxFPsize);
		// sort originalPatterns in descending order of support
		Map<Pattern, Integer> sortedOrigPatterns = Evaluation.sortByValue(originalPatterns);
		// debug: check if these are actually sorted
		if (!Evaluation.checkMapSorting(sortedOrigPatterns))
			System.out.println("Sorting error @ FP mining.");
		// find top-k patterns in original
		List<Pattern> origTopK = new ArrayList<Pattern>();
		for (Map.Entry<Pattern, Integer> entry : sortedOrigPatterns.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained in minePatterns function?!");
			if (origTopK.size() < rankingN) {
            	origTopK.add(entry.getKey());
            	//System.out.println("Pattern: " + entry.getKey().toString() + 
            	//		", supp: " + entry.getValue());
            }
            else {
            	break; // out of the for loop
            }
		}
		// calculate Kendall's tau for topK patterns
		int concordantPairs = 0;
		int reversedPairs = 0;
		for (int i = 0; i < origTopK.size(); i++) {
			for (int j = i+1; j < origTopK.size(); j++) {
				Pattern moreFreq = origTopK.get(i);
				Pattern lessFreq = origTopK.get(j);
				if (!synPatterns.containsKey(moreFreq) && !synPatterns.containsKey(lessFreq)) {
					// tie, do nothing
				} else if (synPatterns.containsKey(moreFreq) && !synPatterns.containsKey(lessFreq)) {
					concordantPairs++;
				} else if (!synPatterns.containsKey(moreFreq) && synPatterns.containsKey(lessFreq)) {
					reversedPairs++;
				} else { // contains both 
					if (synPatterns.get(moreFreq) > synPatterns.get(lessFreq))
						concordantPairs++;
					else if (synPatterns.get(lessFreq) > synPatterns.get(moreFreq))
						reversedPairs++;
				}
			}
		}
		int n = origTopK.size();
		double denom = n*(n-1)/((double)2.0);
		//System.out.println("Concordant: " + concordantPairs + ", reversed: " + reversedPairs);
		return ((double)concordantPairs-reversedPairs)/((double)denom);
	}
	
	// Measure the error in location coverage thru Kendall's tau for rank correlation
	public static double locationCoverageKendallTau(List<Trajectory> original, List<Trajectory>
			synthetic, int gridCell) {
		// find data boundaries so that an appropriate grid can be laid
		List<Double> boundaries = Main.getDataBoundaries(original);
		double minX = boundaries.get(0);
		double maxX = boundaries.get(1);
		double minY = boundaries.get(2);
		double maxY = boundaries.get(3);
			
		// lay grid on original and synthetic databases
		Grid ug = new Grid(gridCell, minX, maxX, minY, maxY);
		List<GridTrajectory> origDBgrid = Main.convertTrajToGridTraj(original, ug, true);
		List<GridTrajectory> synDBgrid = Main.convertTrajToGridTraj(synthetic, ug, true);
	
		// how many trajectories visited each cell?
		List<Cell> cells = ug.getCells();
		int[] actualCounts = new int[cells.size()];
		int[] syntheticCounts = new int[cells.size()];
		for (int i = 0; i < cells.size(); i++) {
			actualCounts[i] = 0;
			syntheticCounts[i] = 0;
			Cell c = cells.get(i);
			for (GridTrajectory t : origDBgrid) {
				if (t.passesThrough(c))
					actualCounts[i] = actualCounts[i]+1;
			}
			for (GridTrajectory t : synDBgrid) {
				if (t.passesThrough(c))
					syntheticCounts[i] = syntheticCounts[i]+1;
			}
		}
		
		// measure kendall tau with actualCounts and syntheticCounts
		int concordantPairs = 0;
		int reversedPairs = 0;
		for (int i = 0; i < actualCounts.length; i++) {
			for (int j = i+1; j < actualCounts.length; j++) {
				if (actualCounts[i] > actualCounts[j]) {
					if (syntheticCounts[i] > syntheticCounts[j])
						concordantPairs++;
					else 
						reversedPairs++;
				}
				if (actualCounts[i] < actualCounts[j]) {
					if (syntheticCounts[i] > syntheticCounts[j])
						reversedPairs++;
					else 
						concordantPairs++;
				}
			}
		}
		int n = cells.size();
		double denom = n*(n-1)/((double)2.0);
		//System.out.println("Concordant: " + concordantPairs + ", reversed: " + reversedPairs);
		return ((double)concordantPairs-reversedPairs)/((double)denom);
		
	}

	public static double calcJSD(List<Double> origProb, List<Double> synProb) {
		List<Double> avgProb = new ArrayList<Double>();
		for (int i = 0; i < origProb.size(); i++) {
			avgProb.add((origProb.get(i)+synProb.get(i))/2.0);
		}
		//System.out.println(origProb);
		//System.out.println(avgProb);
		//System.out.println(synProb);
		return 0.5*calcKL(origProb,avgProb) + 0.5*calcKL(synProb,avgProb);
	}
	
	private static double calcKL(List<Double> prob1, List<Double> prob2) {
		double KL = 0.0;
		for (int i = 0; i < prob1.size(); i++) {
			double Qofi = prob2.get(i);
			double Pofi = prob1.get(i);
			if (Pofi != 0.0)
				KL = KL + Math.log(Pofi/Qofi)*Pofi;
		}
		return KL;
	}
	
	private static double calcF1(List<Pattern> original, List<Pattern> synthetic) {
		double noCorrectPositives = 0.0;
		for (Pattern p : synthetic) {
			for (Pattern p2 : original) {
				if (p.isEqual(p2)) {
					noCorrectPositives++;
					break;
				}
			}
		}
		//System.out.println("# of correct positives: " + noCorrectPositives);
		double precision = noCorrectPositives / synthetic.size();
		double recall = noCorrectPositives / original.size();
		return (2*precision*recall)/(precision+recall);
	}
	
	private static double JaccardSimilarity(List<Pattern> original, List<Pattern> synthetic) {
		List<Pattern> intersection = new ArrayList<Pattern>();
		for (Pattern p : original) {
			for (int i = 0; i < synthetic.size(); i++) {
				if (synthetic.get(i).isEqual(p)) {
					intersection.add(p);
					break;
				}
			}
		}
		List<Pattern> union = new ArrayList<Pattern>();
		for (Pattern p : original) {
			union.add(p);
		}
		for (Pattern p : synthetic) {
			if (!intersection.contains(p))
				union.add(p);
		}
		return ((double)intersection.size())/((double)union.size());
	}
	
	public static double calculateDTW (List<Point> list1, List<Point> list2) {
		double[][] dtwMatrix = new double[list1.size()+1][list2.size()+1];
		
		// fill first row and column
		for (int i = 1; i <= list1.size(); i++) {
			dtwMatrix[i][0] = Double.POSITIVE_INFINITY;
		}
		for (int i = 1; i < list2.size(); i++) {
			dtwMatrix[0][i] = Double.POSITIVE_INFINITY;
		}
		dtwMatrix[0][0] = 0;
		
		// do O(mn) computations
		for (int i = 1; i <= list1.size(); i++) {
			for (int j = 1; j <= list2.size(); j++) {
				Point p1 = list1.get(i-1);
				Point p2 = list2.get(j-1);
				double euclidDist = p1.euclideanDistTo(p2);
				
				// find min
				double min = Double.POSITIVE_INFINITY;
				if (dtwMatrix[i-1][j] < min) 
					min = dtwMatrix[i-1][j];
				if (dtwMatrix[i][j-1] < min)
					min = dtwMatrix[i][j-1];
				if (dtwMatrix[i-1][j-1] < min)
					min = dtwMatrix[i-1][j-1];
				
				// enter into dtw matrix
				dtwMatrix[i][j] = euclidDist + min;
			}
		}
		
		// end
		return dtwMatrix[list1.size()][list2.size()];
	}

}
