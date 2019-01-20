import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.distribution.LaplaceDistribution;
import org.apache.commons.math3.ml.distance.EarthMoversDistance;

import edu.wlu.cs.levy.CG.KDTree;


public class Main {
	
	// normalizes the whole markov transition matrix instead of row by row
	private static List<Double> getMarkovList(List<GridTrajectory> trajs, Grid g) {
		List<Cell> cells = g.getCells();
		// fill actual counts
		double sumAll = 0.0;
		double[][] actualCounts = new double[cells.size()][cells.size()];
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				actualCounts[i][j] = 0.0;
			}
		}
		for (GridTrajectory t : trajs) {
			List<Cell> trajCells = t.getCells();
			for (int i = 0; i < trajCells.size()-1; i++) {
				Cell thisCell = trajCells.get(i);
				Cell nextCell = trajCells.get(i+1);
				int thisCellIndex = cells.indexOf(thisCell);
				int nextCellIndex = cells.indexOf(nextCell);
				actualCounts[thisCellIndex][nextCellIndex] = 
						actualCounts[thisCellIndex][nextCellIndex] + 1.0;
				sumAll = sumAll + 1.0;
			}
		}
		// normalize and return
		List<Double> tbr = new ArrayList<Double>();
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				tbr.add(((double)actualCounts[i][j])/sumAll);
			}
		}
		// return
		return tbr;
	}
	
	public static List<Trajectory> readTrajectories(File f) throws Exception {
		List<Trajectory> tbR = new ArrayList<Trajectory>();
		int trajCount = 0;
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = "";
		String cvsSplitBy = ";";
		while ((line = br.readLine()) != null) {
			if (line.startsWith("#")) {
				trajCount++;
			} else {
				line = line.substring(3); // first 3 chars: ">0:"
				String[] coordPairs = line.split(cvsSplitBy);  // coordPair: x,y
				Trajectory newTraj = new Trajectory();
				for (String s : coordPairs) {
					String[] twocoords = s.split(",");
					double xCoord = Double.parseDouble(twocoords[0]);
					double yCoord = Double.parseDouble(twocoords[1]);
					newTraj.addCoordinates(xCoord, yCoord);					
				}
				tbR.add(newTraj);
				//System.out.println("Trajectory had " + newTraj.getSize() + " locations.");
			}
		}
		br.close();
		System.out.println("Read " + trajCount + " trajectories from file.");
		return tbR;
	}
	
	public static List<Double> getDataBoundaries (List<Trajectory> data) {
		double dataMinX = Double.POSITIVE_INFINITY;
		double dataMinY = Double.POSITIVE_INFINITY;
		double dataMaxX = Double.NEGATIVE_INFINITY;
		double dataMaxY = Double.NEGATIVE_INFINITY;
		for (Trajectory t : data) {
			if (t.getMinXCoord() < dataMinX)
				dataMinX = t.getMinXCoord();
			if (t.getMaxXCoord() > dataMaxX)
				dataMaxX = t.getMaxXCoord();
			if (t.getMinYCoord() < dataMinY)
				dataMinY = t.getMinYCoord();
			if (t.getMaxYCoord() > dataMaxY)
				dataMaxY = t.getMaxYCoord();
		}
		List<Double> tbR = new ArrayList<Double>();
		tbR.add(dataMinX);
		tbR.add(dataMaxX);
		tbR.add(dataMinY);
		tbR.add(dataMaxY);
		return tbR;
	}
	
	public static List<GridTrajectory> convertTrajToGridTraj(List<Trajectory> origDB,
			Grid ug, boolean interp) {
		List<GridTrajectory> tbr = new ArrayList<GridTrajectory>();
		for (Trajectory t : origDB) {
			tbr.add(new GridTrajectory(t, ug, interp));
		}
		return tbr;
	}
	
	public static List<GridTrajectory> convertTrajToGridTraj(List<Trajectory> origDB,
			Grid g, boolean interp, double epsilon, double unusedEpsilon) {
		// perform conversion
		List<GridTrajectory> tbr = new ArrayList<GridTrajectory>();
		for (Trajectory t : origDB) {
			tbr.add(new GridTrajectory(t, g, interp));
		}
		// adaptive grid capability
		// calculate actual densities
		HashMap<Cell, Double> cellDensities = new HashMap<Cell, Double>();
		for (Cell c : g.getCells())
			cellDensities.put(c, 0.0);
		for (GridTrajectory t : tbr) {
			for (Cell c : t.getCells()) {
				cellDensities.put(c, cellDensities.get(c) + 1.0/t.getCells().size());
			}
		}
		// add noise to densities
		LaplaceDistribution ld = new LaplaceDistribution(0, 1.0/epsilon);
		for (Cell c : g.getCells()) {
			double noisydensity = cellDensities.get(c) + ld.sample();
			if (noisydensity < 0.0001)
				noisydensity = 0;
			cellDensities.put(c, noisydensity);
		}
		// divide into further cells
		for (Cell c : g.getCells()) { 
			// efficiency - start
			//List<Trajectory> relevantTrajs = new ArrayList<Trajectory>(); 
			//for (int i = 0; i < origDB.size(); i++) {
			//	if (tbr.get(i).getCells().contains(c)) {
			//		relevantTrajs.add(origDB.get(i));
			//	}
			//}
			// efficiency - end
			c.divideFurther(cellDensities.get(c), epsilon+unusedEpsilon, origDB);
			//System.out.println("div");
		}
		// return
		return tbr;
	}
	
	public static void writeToFile(List<Trajectory> data, String filename) throws Exception {
		PrintWriter writer = new PrintWriter(filename);
		for (int i = 0; i < data.size(); i++) {
			writer.println("#" + i + ":");
			writer.print(">0:");
			Trajectory t = data.get(i);
			for (int j = 0; j < t.getSize(); j++) {
				Point thisp = t.getPoint(j);
				writer.print(thisp.getX() + "," + thisp.getY() + ";");
			}
			writer.print("\r\n");
		}
		writer.close();
	}
	
	public static List<Trajectory> convertGridTrajToTraj(List<GridTrajectory> input, Grid g) {
		List<Trajectory> tbr = new ArrayList<Trajectory>();
		for (GridTrajectory t : input) {
			Trajectory out = new Trajectory();
			List<Cell> cellsForConversion = t.getCells();
			for (Cell c : cellsForConversion) {
				out.addPoint(c.sampleRandomPoint());
			}		
			tbr.add(out);		
		}
		return tbr;
	}
	
	public static double[][] extractMarkovProbs(List<GridTrajectory> origDBgrid, 
			Grid g, double privacyBudget) {
		List<Cell> cells = g.getCells();
		// get actual counts 
		double[][] actualCounts = new double[cells.size()][cells.size()];
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				actualCounts[i][j] = 0.0;
			}
		}
		// fill actualCounts
		for (GridTrajectory t : origDBgrid) {
			List<Cell> trajCells = t.getCells();
			for (int i = 0; i < trajCells.size()-1; i++) {
				Cell thisCell = trajCells.get(i);
				Cell nextCell = trajCells.get(i+1);
				int thisCellIndex = cells.indexOf(thisCell);
				int nextCellIndex = cells.indexOf(nextCell);
				actualCounts[thisCellIndex][nextCellIndex] = 
						actualCounts[thisCellIndex][nextCellIndex] + 
						1.0/((double)(trajCells.size()-1)); 
			}
		}
		// DEBUG - START
		/*Main.printMatrix(actualCounts);
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				if (actualCounts[i][j] > 0) {
					Cell startc = cells.get(i);
					Cell nextc = cells.get(j);
					System.out.println("Nonzero w/ start: " + startc + " next: " + nextc);
				}
			}
		}*/
		// DEBUG-END
		// add noise to obtain noisyCounts
		LaplaceDistribution ld = new LaplaceDistribution(0, 1.0/privacyBudget);
		double[][] noisyCounts = new double[cells.size()][cells.size()];
		for (int i = 0; i < noisyCounts.length; i++) {
			for (int j = 0; j < noisyCounts.length; j++) {
				double newval = actualCounts[i][j] + ld.sample();
				boolean ijAdjacent = false;
				 // you can move to adjacent cell - not current, and no jumps allowed
				Cell c1 = cells.get(i);
				Cell c2 = cells.get(j);
				if (g.areAdjacent(c1, c2))
					ijAdjacent = true;
				if (newval >= 0 && ijAdjacent) { // check for non-negativity and adjacency
					noisyCounts[i][j] = newval;
				}
				else {
					noisyCounts[i][j] = 0.0;
				}
			}
		}
		// DEBUG - START
		/*Main.printMatrix(noisyCounts);
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				if (noisyCounts[i][j] > 0) {
					Cell startc = cells.get(i);
					Cell nextc = cells.get(j);
					//System.out.println("Nonzero w/ start: " + startc + " next: " + nextc);
					//check equality with actualcounts
					//if (Math.abs(noisyCounts[i][j] - actualCounts[i][j]) < 2)
					//	System.out.println("Same w/ start: " + startc + " next: " + nextc);
				}
			}
		}
		*/
		// DEBUG-END
		// obtain and return noisyProbs, i.e., differentially private (noisy) Markov 
		// transition matrix
		double[][] noisyProbs = new double[cells.size()][cells.size()];
		for (int i = 0; i < cells.size(); i++) {
			double sumOfRow = 0.0;
			for (int j = 0; j < cells.size(); j++) {
				sumOfRow = sumOfRow + noisyCounts[i][j];
			}
			if (sumOfRow == 0.0) {
				for (int j = 0; j < cells.size(); j++) {
					noisyProbs[i][j] = 0.0;
				}
			} else {
				for (int j = 0; j < cells.size(); j++) {
					noisyProbs[i][j] = noisyCounts[i][j] / sumOfRow;
				}
			}
		}
		// DEBUG - START
		//Main.printMatrix(noisyProbs);
		for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				if (Double.isNaN(noisyProbs[i][j]))
					System.out.println("NaN observed in 1-step Markov matrix!");
			}
		}
		/*for (int i = 0; i < cells.size(); i++) {
			for (int j = 0; j < cells.size(); j++) {
				if (noisyCounts[i][j] > 0) {
					Cell startc = cells.get(i);
					Cell nextc = cells.get(j);
					System.out.println("Nonzero w/ start: " + startc + " next: " + nextc);
					//check equality with actualcounts
					//if (Math.abs(noisyCounts[i][j] - actualCounts[i][j]) < 2)
					//	System.out.println("Same w/ start: " + startc + " next: " + nextc);
				}
			}
		}*/
		// DEBUG-END
		return noisyProbs;		
	}
	
	public static List<double[][]> precomputeMarkov(double[][] oneStepMatrix, int maxSteps) {
		List<double[][]> tbR = new ArrayList<double[][]>();
		tbR.add(oneStepMatrix); // placeholder dummy, 0-step matrix
		tbR.add(oneStepMatrix); // 1-step matrix
		for (int i = 2; i <= maxSteps; i++) {
			double[][] prevStepMatrix = tbR.get(i-1);
			double[][] currStepMatrix = Main.matrixMultiply(prevStepMatrix, oneStepMatrix);
			tbR.add(currStepMatrix);
		}
		return tbR;
	}
	
	// from the web (stackoverflow??)
	public static double[][] matrixMultiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }
	
	private static Cell findWithMarkov(Cell prevCell, Cell eventualCell, 
			double[][] step1matrix, double[][] stepNmatrix, Grid g) {
		List<Cell> candidateCells = g.getAdjacentCells(prevCell);
		double[] candidateProbs = new double[candidateCells.size()];
		for (int i = 0; i < candidateProbs.length; i++) {
			Cell thisCandidate = candidateCells.get(i);
			// find P[thisCandidate | prevCell] according to step1matrix
			int rowNo = g.getPosInListForm(prevCell);
			int colNo = g.getPosInListForm(thisCandidate);
			double prob1 = step1matrix[rowNo][colNo];
			// find P[eventualCell | thisCandidate] according to stepNmatrix
			rowNo = g.getPosInListForm(thisCandidate);
			colNo = g.getPosInListForm(eventualCell);
			double prob2 = stepNmatrix[rowNo][colNo];
			// find candidate probability
			double candProb = prob1 * prob2;
			candidateProbs[i] = candProb;
		}
		// DEBUG-START (debug: check sum of candidateProbs)
		//System.out.println("Previous cell: " + prevCell);
		//System.out.println("Candidates for next move: " + candidateCells.toString());
		//System.out.println("Candidates probabilities: " + Arrays.toString(candidateProbs));
		// DEBUG-END
		
			// random sampling w/ candidates: candidateCells, weights: candidateProbs
		// find sum of probabilities, normalize to [0,1]
		double sum = 0.0;
		for (int i = 0; i < candidateCells.size(); i++) {
			sum = sum + candidateProbs[i];
		}
		// special (rare) case: sum == 0
		if (sum < 0.00001) {
			// add the first cell in the most direct route 
			List<Cell> directRoute = g.giveInterpolatedRoute(prevCell, eventualCell);
			if (directRoute.size() > 1) {
				return directRoute.get(1);
			} else {
				return eventualCell;
			}
		}
		// handling of special case complete, proceed with regular case
		List<Double> normalizedProbabilities = new ArrayList<Double>();
		for (int i = 0; i < candidateProbs.length; i++) {
			normalizedProbabilities.add(((double)candidateProbs[i]/sum));
		}
		double randomVal = new Random().nextDouble();
		double seenSoFar = 0.0;
		for (int i = 0; i < normalizedProbabilities.size(); i++) {
			seenSoFar = seenSoFar + normalizedProbabilities.get(i);
			if (seenSoFar >= randomVal) {
				return candidateCells.get(i);
			}
		}
		
		// we should never be reaching this line, but if we do:  
		// return the cell in the most direct route between prevCell and eventualCell
		List<Cell> directRoute = g.giveInterpolatedRoute(prevCell, eventualCell);
		if (directRoute.size() > 1) {
			return directRoute.get(1);
		} else {
			return eventualCell;
		}
	}
	
	public static Cell getNthDensestCell(List<GridTrajectory> origDBgrid, Grid grid, int N) {
		List<Cell> cellList = grid.getCells();
		HashMap<Cell, Integer> cellDensities = new HashMap<Cell, Integer>();
		for (Cell c : cellList) 
			cellDensities.put(c, 0);
		for (GridTrajectory t : origDBgrid) {
			for (Cell c : t.getCells()) {
				cellDensities.put(c, cellDensities.get(c)+1);
			}
		}
		Map<Cell, Integer> sortedDensities = Evaluation.sortByValue(cellDensities);
		int currIndex = 1;
		
		for (Map.Entry<Cell, Integer> entry : sortedDensities.entrySet()) {
			if (currIndex != N) {
				currIndex++;
			}
			else {
				return entry.getKey();
			}
		}
		
		// cannot be found? must be error
		return null;		
	}
	
	// from stack overflow
	private static double[] convertListDoubleToArray(List<Double> doubles) {
		 double[] target = new double[doubles.size()];
		 for (int i = 0; i < target.length; i++) {
		    target[i] = doubles.get(i).doubleValue();  // java 1.4 style
		 }
		 return target;
	}
	
	public static boolean evaluateBayesianAttack(List<Double> TripPrior, List<Double> MarkovPrior, 
			List<GridTrajectory> relevantSubset, Grid grid, double VARTHETA) {
		// find posterior belief of relevantSubset (D_Z)
		List<Double> TripPosterior = 
			(new TripDistribution(relevantSubset, grid, 100000.0)).getTripProbsAsList();
		List<Double> MarkovPosterior = Main.getMarkovList(relevantSubset, grid);
		
		// compute difference between prior and posterior	
		double[] TripPriorArray = convertListDoubleToArray(TripPrior);
		double[] MarkovPriorArray = convertListDoubleToArray(MarkovPrior);
		double[] TripPosteriorArray = convertListDoubleToArray(TripPosterior);
		double[] MarkovPosteriorArray = convertListDoubleToArray(MarkovPosterior);
		//EarthMoversDistance emd = new EarthMoversDistance();
		//double tripEMD = emd.compute(TripPrior, TripPosterior);
		//double markovEMD = emd.compute(MarkovPrior, MarkovPosterior);
		double tripEMD = Evaluation.calcJSD(TripPrior, TripPosterior);
		double markovEMD = Evaluation.calcJSD(MarkovPrior, MarkovPosterior);
		
		// check if EMD > vartheta
		if (tripEMD > VARTHETA) {
			return false;
		} 
		if (markovEMD > VARTHETA) {
			return false;
		}
		return true;
	}
	
	public static List<GridTrajectory> generateSyntheticTrajs(Grid g, TripDistribution td, 
			LengthDistribution ld, double[][] markovProbs, int desiredTrajCount) {
		List<GridTrajectory> tbr = new ArrayList<GridTrajectory>();
		// pre-compute Markov transition matrices: 1-step, 2-step, 3-step, ..., 100-step
		List<double[][]> transitionMatrices = new ArrayList<double[][]>(); 
		transitionMatrices = Main.precomputeMarkov(markovProbs, 100);
		//System.out.println("Pre-computation of Markov transition probabilities complete.");
		// synthesize trajectories
		for (int cnt = 0; cnt < desiredTrajCount; cnt++) {
			// determine trip
			Cell[] c = td.sampleStartEndCells();
			Cell startCell = c[0];
			Cell endCell = c[1];
			// determine length
			int minLen = g.findShortestLengthBetween (startCell, endCell);
			int extraLen = ld.getLengthSample(startCell, endCell);
			int totalLen = minLen + extraLen;
			//System.out.println("start: " + startCell + ", end: " + endCell + ", totalLen: " +
			//		totalLen);
			// determine intermediate cells
			Cell[] newTrajCells = new Cell[totalLen];
			newTrajCells[0] = startCell;
			newTrajCells[totalLen-1] = endCell;
			for (int i = 1; i < totalLen-1; i++) { // loops intermediate cells
				Cell prevCell = newTrajCells[i-1];
				Cell eventualCell = newTrajCells[newTrajCells.length-1];
				int stepsToEventual = totalLen-i-1;
				double[][] step1matrix = transitionMatrices.get(1);
				double[][] stepNmatrix = transitionMatrices.get(stepsToEventual);
				newTrajCells[i] = findWithMarkov(prevCell, eventualCell, step1matrix, 
						stepNmatrix, g);
			}
			// initialize and add to synthetic DB
			GridTrajectory synTraj = new GridTrajectory(newTrajCells);
			tbr.add(synTraj);
		}
		return tbr;
	}
	
	public static List<Trajectory> getOutlierTrajs (List<Trajectory> trajs, double[] scores, 
			int number) {
		List<Trajectory> tbr = new ArrayList<Trajectory>();
		
		double[] copy = Arrays.copyOf(scores, scores.length);
		Arrays.sort(copy); // in ascending order
		List<Double> scoresList = new ArrayList<Double>();
		for (int i = 0; i < scores.length; i++) {
			scoresList.add(scores[i]);
		}
		for (int i = 0; i < number; i++) {
			int index = copy.length-1-i; // since copy is in ascending order, begin from last
			double searchForVal = copy[index];
			int foundInTrajs = scoresList.indexOf(searchForVal);
			tbr.add(trajs.get(foundInTrajs));
		}
		
		return tbr;
	}
	
	public static List<Trajectory> generateAttackResilientTrajectories (Grid grid, double[][] markovProbs,
			TripDistribution td, LengthDistribution ld, int desiredNumberOfTrajectories, boolean attacksON,
			List<GridTrajectory> origDBgrid, List<Trajectory> origDBxy) 
			throws Exception {
		
		// Part 1: Differential privacy (Section 5.1)
		List<GridTrajectory> synDB = generateSyntheticTrajs(grid, td, ld, markovProbs, 
				desiredNumberOfTrajectories);
		if (!attacksON) {
			List<Trajectory> synDBxy = Main.convertGridTrajToTraj(synDB, grid);
			System.out.println("Only differential privacy (no attacks)...");
			return synDBxy;
		}
		
		// Part 2: Attacks (Section 5.2)
		
		// ATTACK #1: BAYESIAN INFERENCE
		double VARTHETA = 0.1; // privacy parameter in Defense 1
		// set one cell as the privacy zone - last parameter denotes which cell will be the sensitive zone
		// e.g., if parameter = N, then Nth densest cell will be privacy zone
		Cell SENSITIVEZONE = Main.getNthDensestCell(origDBgrid, grid, 10);
		// set prior belief of adversary		
		List<Double> TripPrior = td.getTripProbsAsList();
		List<Double> MarkovPrior = Main.getMarkovList(origDBgrid, grid);
		// evaluate the defense
		boolean bayesianAttackPass = false;
		while (!bayesianAttackPass) {
			// find subset of synthetic database that passes through privacy zone
			List<GridTrajectory> relevantSubset = new ArrayList<GridTrajectory>();
			for (GridTrajectory t : synDB) {
				if (t.getCells().contains(SENSITIVEZONE)) {
					relevantSubset.add(t);
				}
			}
			//System.out.println("relevant subset size: " + relevantSubset.size()); // debug msg
			if (relevantSubset.size() == 0 || evaluateBayesianAttack(TripPrior, MarkovPrior, 
					relevantSubset, grid, VARTHETA)) {
				bayesianAttackPass = true;
				break;
			} else {
				synDB.removeAll(relevantSubset);
				int newcountrequired = relevantSubset.size();
				List<GridTrajectory> newtrajs = generateSyntheticTrajs(grid, td, ld, markovProbs,
						newcountrequired);
				synDB.addAll(newtrajs);
			}			
		}
		System.out.println("Bayesian defense - OK!");
		// Attack 1 ends
		
		// convert grid trajectories to actual trajectories
		List<Trajectory> synDBxy = Main.convertGridTrajToTraj(synDB, grid);
		
		// Attack 2: PARTIAL SNIFFING
		// choose a sniff region. for simplicity, a sniff region = a cell.
		List<Trajectory> attack2FailingSynthetics = new ArrayList<Trajectory>();
		Cell SNIFFREGION = getNthDensestCell(origDBgrid, grid, 12);
		List<Trajectory> sniffedTrajectories = new ArrayList<Trajectory>();
		for (int i = 0; i < origDBxy.size(); i++) {
			GridTrajectory t = origDBgrid.get(i);
			if (t.getCells().contains(SNIFFREGION)) {
				sniffedTrajectories.add(origDBxy.get(i));
			}
		}
		// for each sniffed trajectory : find most similar synthetic traj
		for (Trajectory sniffedTraj : sniffedTrajectories) {
			List<Point> sniffedActualPts = sniffedTraj.getSniffedPoints(SNIFFREGION);
			double minDistSoFar = Double.POSITIVE_INFINITY;
			int minIndex = -1;
			for (int i = 0; i < synDBxy.size(); i++) {
				Trajectory synTraj = synDBxy.get(i);
				List<Point> synSniffed = synTraj.getSniffedPoints(SNIFFREGION);
				if (synSniffed.size() < 1) {
					continue;
				} else {
					double dist = Evaluation.calculateDTW(sniffedActualPts, synSniffed);
					if (dist < minDistSoFar) {
						minDistSoFar = dist;
						minIndex = i;
					}
				}
			}
			if (minIndex < 0) {
				continue; // no match
			}
			Trajectory mostSimilarSynthetic = synDBxy.get(minIndex);
			// intersection between mostSimSyn and sniffedTraj < threshold
			int intersect = mostSimilarSynthetic.calcIntersectionCount(sniffedTraj);
			int INTERSECTION_THRESHOLD = 5; // this is a parameter
			if (intersect >= INTERSECTION_THRESHOLD) { 
				attack2FailingSynthetics.add(mostSimilarSynthetic);
			}
			// bound sensitive location disclosure for mostSimSyn
			List<Trajectory> toSend = new ArrayList<Trajectory>();
			toSend.add(mostSimilarSynthetic);
			GridTrajectory mostSimGrid = Main.convertTrajToGridTraj(toSend, grid, true).get(0);
			List<Cell> mostSimCells = mostSimGrid.getCells();
			if (mostSimCells.contains(SENSITIVEZONE))
				attack2FailingSynthetics.add(mostSimilarSynthetic);
		}
		// handle attack2FailingSynthetics - remove them from synthetic database
		// generate and add new synthetics that dont pass through sniffZone
		synDBxy.removeAll(attack2FailingSynthetics);
		int newTrajsNeeded = desiredNumberOfTrajectories - synDBxy.size();
		for (int i = 0; i < newTrajsNeeded; i++) {
			GridTrajectory newtraj = generateSyntheticTrajs(grid, td, ld, markovProbs, 1).get(0);
			if (newtraj.getCells().contains(SNIFFREGION) || newtraj.getCells().contains(SENSITIVEZONE)) {
				i--;
			} else {
				List<GridTrajectory> toSend = new ArrayList<GridTrajectory>();
				toSend.add(newtraj);
				Trajectory toAdd = Main.convertGridTrajToTraj(toSend, grid).get(0);
				synDBxy.add(toAdd);
			}
		}
		System.out.println("Sniffing defense - OK!");
		// Attack 2 ends
		
		// ATTACK 3 BEGIN: OUTLIER LEAKAGE
		// PARAMETERS FOR ATTACK 3 - Parameterized so we can run experiments easier
		int TOP_K_NEIGHBORS = 50; // since notion is k-NN distance, determine k
		int N_OUTLIERS = 200; // determine # of outliers (200)
		double CLOSEST_THRESHOLD = 0.1; // in step3a
		int PLAUSIBLE_DENIABILITY_KAPPA = 100;
		double PLAUSIBLE_DENIABILITY_BETA = 0.05; // 0.025
		int outlierCnt = 0;
		
		// (i) Find trip start outliers: trip start deviates significantly from others
		KDTree<Integer> kd = new KDTree<Integer>(2); // 2-dimensional kd-tree
		double [][] elts = new double [synDBxy.size()][2]; // elements in kd-tree
		for (int i = 0; i < synDBxy.size(); i++) {
			elts[i][0] = synDBxy.get(i).getPoint(0).getX();
			elts[i][1] = synDBxy.get(i).getPoint(0).getY();
			kd.insert(elts[i], i);
		} // the kd-tree is now filled
		double[] scores = new double[synDBxy.size()];
		for (int i = 0; i < synDBxy.size(); i++) {
			double[] targ = new double[2];
			targ[0] = synDBxy.get(i).getPoint(0).getX();
			targ[1] = synDBxy.get(i).getPoint(0).getY();
			List<Integer> nbrIndexes = kd.nearest(targ, TOP_K_NEIGHBORS);
			double sum = 0.0;
			for (int j : nbrIndexes) {
				Point thispt = synDBxy.get(i).getPoint(0);
				Point nbrpt = synDBxy.get(j).getPoint(0);
				sum = sum + thispt.euclideanDistTo(nbrpt);
			}
			scores[i] = sum/((double) nbrIndexes.size());
		}
		//System.out.println("ready to get outliers");
		List<Trajectory> outliers = getOutlierTrajs(synDBxy, scores, N_OUTLIERS);
		//System.out.println("Trip start outliers: " + outliers.size());
		//System.out.println(outliers);
				
		// (ii) Plausible deniability for trip start outliers
		for (Trajectory Tout : outliers) {
			// find Tclosest
			Trajectory Tclosest = null;
			double minDist = Double.POSITIVE_INFINITY;
			for (Trajectory Tcand : origDBxy) {
				Point p1 = Tout.getPoint(0);
				Point p2 = Tcand.getPoint(0);
				double currDist = p1.euclideanDistTo(p2);
				if (currDist < minDist) {
					minDist = currDist;
					Tclosest = Tcand;
				}
			}
			// if dist(Tclosest,Tout) > threshold ==> cant be linked, pass privacy test
			List<Double> boundaries = getDataBoundaries(origDBxy);
			double minX = boundaries.get(0);
			double maxX = boundaries.get(1);
			double minY = boundaries.get(2);
			double maxY = boundaries.get(3);
			double maxPossibleDistance = Math.sqrt(
					(maxX-minX)*(maxX-minX) + (maxY-minY)*(maxY-minY));
			double threshold = maxPossibleDistance * CLOSEST_THRESHOLD;
			if (minDist > threshold) { }
			else {
				// otherwise, test for plausible deniability
				double secondthr = maxPossibleDistance * PLAUSIBLE_DENIABILITY_BETA;
				int currk = 0;
				boolean satisfied = false; // plausible deniability is satisfied?
				for (Trajectory Tclosecand : origDBxy) {
					Point p1 = Tout.getPoint(0);
					Point p2 = Tclosecand.getPoint(0);
					double distv = p1.euclideanDistTo(p2);
					if (Math.abs(distv-minDist) <= secondthr)
						currk++;
					if (currk >= PLAUSIBLE_DENIABILITY_KAPPA) {
						satisfied = true;
						break;
					}
				}
				if (!satisfied) {
					synDBxy.remove(Tout);
					outlierCnt++;
				}
			}
		}
		
		//System.out.println("Dataset size: " + synDBxy.size());
		
		// (iii) Find trip end outliers
		kd = new KDTree<Integer>(2); // 2-dimensional kd-tree
		elts = new double [synDBxy.size()][2]; // elements in kd-tree
		for (int i = 0; i < synDBxy.size(); i++) {
			elts[i][0] = synDBxy.get(i).getPoint(synDBxy.get(i).getSize()-1).getX();
			elts[i][1] = synDBxy.get(i).getPoint(synDBxy.get(i).getSize()-1).getY();
			kd.insert(elts[i], i);
		} // the kd-tree is now filled
		scores = new double[synDBxy.size()];
		for (int i = 0; i < synDBxy.size(); i++) {
			double[] targ = new double[2];
			targ[0] = synDBxy.get(i).getPoint(0).getX();
			targ[1] = synDBxy.get(i).getPoint(0).getY();
			List<Integer> nbrIndexes = kd.nearest(targ, TOP_K_NEIGHBORS);
			double sum = 0.0;
			for (int j : nbrIndexes) {
				Point thispt = synDBxy.get(i).getPoint(synDBxy.get(i).getSize()-1);
				Point nbrpt = synDBxy.get(j).getPoint(synDBxy.get(j).getSize()-1);
				sum = sum + thispt.euclideanDistTo(nbrpt);
			}
			scores[i] = sum/((double) nbrIndexes.size());
		}
		//System.out.println("ready to get outliers");
		outliers = getOutlierTrajs(synDBxy, scores, N_OUTLIERS);
		//System.out.println("Outlier count: " + outlierCnt);
		
		// (iv) Fix trip end outliers
		for (Trajectory Tout : outliers) {
			// find Tclosest
			Trajectory Tclosest = null;
			double minDist = Double.POSITIVE_INFINITY;
			for (Trajectory Tcand : origDBxy) {
				Point p1 = Tout.getPoint(Tout.getSize()-1);
				Point p2 = Tcand.getPoint(Tcand.getSize()-1);
				double currDist = p1.euclideanDistTo(p2);
				if (currDist < minDist) {
					minDist = currDist;
					Tclosest = Tcand;
				}
			}
			// if dist(Tclosest,Tout) > threshold ==> cant be linked, pass privacy test
			List<Double> boundaries = getDataBoundaries(origDBxy);
			double minX = boundaries.get(0);
			double maxX = boundaries.get(1);
			double minY = boundaries.get(2);
			double maxY = boundaries.get(3);
			double maxPossibleDistance = Math.sqrt(
					(maxX-minX)*(maxX-minX) + (maxY-minY)*(maxY-minY));
			double threshold = maxPossibleDistance * CLOSEST_THRESHOLD;
			if (minDist > threshold) { }
			else {
				// otherwise, test for plausible deniability
				double secondthr = maxPossibleDistance * PLAUSIBLE_DENIABILITY_BETA;
				int currk = 0;
				boolean satisfied = false; // plausible deniability is satisfied?
				for (Trajectory Tclosecand : origDBxy) {
					Point p1 = Tout.getPoint(Tout.getSize()-1);
					Point p2 = Tclosecand.getPoint(Tclosecand.getSize()-1);
					double distv = p1.euclideanDistTo(p2);
					if (Math.abs(distv-minDist) <= secondthr)
						currk++;
					if (currk >= PLAUSIBLE_DENIABILITY_KAPPA) {
						satisfied = true;
						break;
					}
				}
				if (!satisfied) {
					synDBxy.remove(Tout);
					outlierCnt++;
				}
			}
		}
		//System.out.println("Dataset size: " + synDBxy.size());
		
		// (v) Find length outliers
		kd = new KDTree<Integer>(1); // 1-dimensional kd-tree
		elts = new double [synDBxy.size()][1]; // elements in kd-tree
		for (int i = 0; i < synDBxy.size(); i++) {
			elts[i][0] = synDBxy.get(i).getDistanceTravelled();
			kd.insert(elts[i], i);
		} // the kd-tree is now filled
		scores = new double[synDBxy.size()];
		for (int i = 0; i < synDBxy.size(); i++) {
			double[] targ = new double[1];
			targ[0] = synDBxy.get(i).getDistanceTravelled();
			List<Integer> nbrIndexes = kd.nearest(targ, TOP_K_NEIGHBORS);
			double sum = 0.0;
			for (int j : nbrIndexes) {
				sum = sum + Math.abs(synDBxy.get(i).getDistanceTravelled() -
						synDBxy.get(j).getDistanceTravelled());
			}
			scores[i] = sum/((double) nbrIndexes.size());
		}
		//System.out.println("ready to get outliers");
		outliers = getOutlierTrajs(synDBxy, scores, N_OUTLIERS/2);
		
		// (vi) Fix length outliers
		double mindd = Double.POSITIVE_INFINITY;
		double maxdd = Double.NEGATIVE_INFINITY;
		for (Trajectory t1 : origDBxy) {
			double thistrajdist = t1.getDistanceTravelled();
			if (thistrajdist > maxdd)
				maxdd = thistrajdist;
			if (thistrajdist < mindd)
				mindd = thistrajdist;
		}
		double maxPossibleDistance = maxdd-mindd;
		for (Trajectory Tout : outliers) {
			// find Tclosest
			Trajectory Tclosest = null;
			double minDist = Double.POSITIVE_INFINITY;
			for (Trajectory Tcand : origDBxy) {
				double currDist = Tcand.getDistanceTravelled()-Tout.getDistanceTravelled();
				if (currDist < minDist) {
					minDist = currDist;
					Tclosest = Tcand;
				}
			}
			// if dist(Tclosest,Tout) > threshold ==> cant be linked, pass privacy test
			double threshold = maxPossibleDistance * CLOSEST_THRESHOLD;
			if (minDist > threshold) {}
			else {
				// otherwise, test for plausible deniability
				double secondthr = maxPossibleDistance * PLAUSIBLE_DENIABILITY_BETA;
				int currk = 0;
				boolean satisfied = false; // plausible deniability is satisfied?
				for (Trajectory Tclosecand : origDBxy) {
					double distv = Math.abs( 
							Tclosecand.getDistanceTravelled() - Tout.getDistanceTravelled());
					if (Math.abs(distv-minDist) <= secondthr) 
						currk++;
					if (currk >= PLAUSIBLE_DENIABILITY_KAPPA) {
						satisfied = true;
						break;
					}
				}
				if (!satisfied) {
					synDBxy.remove(Tout);
				}
			}
		}
		
		// (vii) Find location visit outliers
		List<Double> synBoundaries = getDataBoundaries(synDBxy);
		double minX = synBoundaries.get(0);
		double maxX = synBoundaries.get(1);
		double minY = synBoundaries.get(2);
		double maxY = synBoundaries.get(3);
		Grid syngrid = new Grid(25, minX, maxX, minY, maxY); // very detailed
		List<GridTrajectory> synDetailed = convertTrajToGridTraj(synDBxy, syngrid, true);
		// find cells that have fewest visits
		List<Cell> cells = syngrid.getCells();
		HashMap<Cell, Integer> cellVisits = new HashMap<Cell, Integer>();
		for (Cell c : cells) {
			cellVisits.put(c, 0);
		}
		for (GridTrajectory t : synDetailed) {
			for (Cell c : t.getCells()) {
				cellVisits.put(c, cellVisits.get(c)+1);
			}
		}
		Map<Cell, Integer> ascendingCellDensities = Evaluation.sortByValue(cellVisits);
		List<Cell> rareCells = new ArrayList<Cell>();
		for (Map.Entry<Cell, Integer> entry : ascendingCellDensities.entrySet()) {
			if (entry.getValue() == null)
				System.out.println("null obtained?");
			else
				rareCells.add(entry.getKey());
			if (entry.getValue() > N_OUTLIERS)
				break;
        }
		outliers = new ArrayList<Trajectory>();
		for (int i = 0; i < synDetailed.size(); i++) {
			GridTrajectory cand = synDetailed.get(i);
			for (Cell c : cand.getCells()) {
				if (rareCells.contains(c))
					outliers.add(synDBxy.get(i));
			}
			if (outliers.size() >= N_OUTLIERS)
				break;
		}
		//System.out.println("location visiting outlier cnt: " + outliers.size());

		// (viii) Fix location visit outliers
		for (Trajectory Tout : outliers) {
			// find Tclosest
			Trajectory Tclosest = null;
			double minDist = Double.POSITIVE_INFINITY;
			double maxDist = Double.NEGATIVE_INFINITY;
			for (Trajectory Tcand : origDBxy) {
				double currDist = Tout.calculateDTWto(Tcand);
				if (currDist < minDist) {
					minDist = currDist;
					Tclosest = Tcand;
				}
				if (currDist > maxDist)
					maxDist = currDist;
			}
			// if dist(Tclosest,Tout) > threshold ==> cant be linked, pass privacy test
			double threshold = maxDist * CLOSEST_THRESHOLD;
			if (minDist > threshold) {}
			else {
				// otherwise, test for plausible deniability
				double secondthr = maxDist * PLAUSIBLE_DENIABILITY_BETA;
				int currk = 0;
				boolean satisfied = false; // plausible deniability is satisfied?
				for (Trajectory Tclosecand : origDBxy) {
					double distv = Tout.calculateDTWto(Tclosecand);
					if (Math.abs(distv-minDist) <= secondthr)
						currk++;
					if (currk >= PLAUSIBLE_DENIABILITY_KAPPA) {
						satisfied = true;
						break;
					}
				}
				if (!satisfied) {
					synDBxy.remove(Tout);
				}
			}
		}

		// (ix) Add new trajectories in place of the deleted trajectories
		int newcountrequired = desiredNumberOfTrajectories - synDBxy.size();
		for (int i = 0; i < newcountrequired; i++) {
			GridTrajectory newtraj = generateSyntheticTrajs(grid, td, ld, markovProbs, 1).get(0);
			if (newtraj.getCells().contains(SNIFFREGION) || newtraj.getCells().contains(SENSITIVEZONE)) {
				i--;
			} else {
				List<GridTrajectory> toSend = new ArrayList<GridTrajectory>();
				toSend.add(newtraj);
				Trajectory toAdd = Main.convertGridTrajToTraj(toSend, grid).get(0);
				synDBxy.add(toAdd);
			}
		}
		System.out.println("Outlier leakage defense - OK!");
		// Attack 3 ends
		
		if (synDBxy.size() != desiredNumberOfTrajectories) 
			System.out.println("Desired # of trajectories could not be generated!");
		return synDBxy;
	}
	
	/* This is the main procedure for synthetic trajectory generation. 
	 * 
	 * Inputs:
	 * List<Trajectory> database of actual trajectories
	 * Privacy budget epsilon
	 * A boolean to denote whether defenses to attacks should be turned on/off
	 * 
	 * The following parameters are not taken as inputs, but are rather hardcoded within
	 * this function:
	 * Initial cell size for adaptive grid first-level (NxN)
	 * Interpolation (set to true)
	 * Privacy budget distribution (see budgetDistnWeights)
	 */
	public static List<Trajectory> Synthesize_Trajectories (List<Trajectory> originalDatabase,
			double totalEpsilon, boolean attacksON) throws Exception {
		// Hardcoded parameters - BEGIN
		boolean interp = true; // interpolate cells so that every move is to adjacent cell
		int cellCount = 6;
		double[] budgetDistnWeights = {0.05, 0.35, 0.50, 0.10}; // grid, Markov, trip, length
		// Hardcoded parameters - END
		
		// Trajectory synthesis - core components - BEGIN
		List<Double> boundaries = Main.getDataBoundaries(originalDatabase);
		double minX = boundaries.get(0);
		double maxX = boundaries.get(1);
		double minY = boundaries.get(2);
		double maxY = boundaries.get(3);
		Grid grid = new Grid(cellCount, minX, maxX, minY, maxY);
		System.out.print("Adaptive grid....");
		List<GridTrajectory> origDBgrid = convertTrajToGridTraj(originalDatabase, grid, interp, 
				budgetDistnWeights[0]*totalEpsilon, (1.0-budgetDistnWeights[0])*totalEpsilon);
		System.out.print("OK!\r\n");
		System.out.print("Markov mobility model....");
		double[][] markovTransitionProbs = extractMarkovProbs(origDBgrid, grid, 
				budgetDistnWeights[1]*totalEpsilon);
		System.out.print("OK!\r\n");
		System.out.print("Trip distribution....");
		TripDistribution td = new TripDistribution(origDBgrid, grid, budgetDistnWeights[2]*totalEpsilon);
		System.out.print("OK!\r\n");
		System.out.print("Length distribution....");
		LengthDistribution ld = new LengthDistribution(origDBgrid, grid, budgetDistnWeights[3]*totalEpsilon);
		System.out.print("OK!\r\n");
		
		System.out.println("Generating synthetic trajectory database...");
		List<Trajectory> syntheticDB = generateAttackResilientTrajectories(grid, markovTransitionProbs,
				td, ld, origDBgrid.size(), attacksON, origDBgrid, originalDatabase);
		
		return syntheticDB;
	}
	
	
	public static void main(String[] args) throws Exception {
		
		// PART 0 - PARAMETERS
		String inputFilename = "brinkhoff.dat";  // file name/path for actual trajectory database
		double totalEpsilon = 1.0;  // total privacy budget (epsilon)
		boolean attacksON = false;  // want to defend against attacks? (Section 3.3)
		// End of Part 0
		
		// Part 1: Read actual trajectory database
		List<Trajectory> originalDatabase = readTrajectories(new File(inputFilename));
		// End of part 1
		
		// Part 2: Generate synthetic trajectory database - repeat N times
		for (int i = 0; i < 5; i++) {
			List<Trajectory> syntheticDatabase = 
					Synthesize_Trajectories(originalDatabase, totalEpsilon, attacksON);
			String outputFileName = inputFilename + "-eps" + totalEpsilon + "-iteration" + i + ".dat";
			Main.writeToFile(syntheticDatabase, outputFileName);
			System.out.println("Done! Wrote trajectories to file: " + outputFileName);	
		}
		
	}

}
