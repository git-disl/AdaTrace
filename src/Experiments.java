import java.io.File;
import java.util.List;
import java.util.Map;


public class Experiments {
	
	public static void main(String[] args) throws Exception {
		
		// Part 0: Parameters for utility experimentation
		int numberOfQueries = 200;  // query error
		double queryErrorSanityBound = 0.01; // query error
		int gridSizeForTripErr = 6; // trip error
		int bucketCountForDiameterErr = 20; // diameter & length error
		int gridSizeForFPErr = 6; // frequent pattern mining error
		int minFPsize = 2; // frequent pattern mining error - min pattern length
		int maxFPsize = 8; // frequent pattern mining error - max pattern length
		int topK = 100; // frequent pattern mining error 
		int locCoverageGridSize = 15; // location coverage kendall tau
		// end of Part 0
		
		// Part 1: Read original data, prepare it for experimentation
		String inputFilename = "brinkhoff.dat";
		List<Trajectory> originalDatabase = Main.readTrajectories(new File(inputFilename));
		Query[] generatedQueries = Evaluation.generateQueries(originalDatabase, 
				numberOfQueries, inputFilename);
		List<Double> boundaries = Main.getDataBoundaries(originalDatabase);
		double minX = boundaries.get(0);
		double maxX = boundaries.get(1);
		double minY = boundaries.get(2);
		double maxY = boundaries.get(3);
		Grid ug = new Grid(gridSizeForFPErr, minX, maxX, minY, maxY);
		List<GridTrajectory> origDBgrid = Main.convertTrajToGridTraj(originalDatabase, 
				ug, true);
		Map<Pattern, Integer> origPatterns = Evaluation.minePatterns(origDBgrid, minFPsize, 
				maxFPsize);
		System.out.println("Generated queries and mined patterns for original database.");
		
		// Part 2: Read synthetic datasets and assess their quality
		File SYNfolder = new File("SYNTHETIC-DATASETS");
		File[] listOfFiles = SYNfolder.listFiles();
		for (File file : listOfFiles) {
			System.out.println("***********");
			List<Trajectory> syntheticDatabase = Main.readTrajectories(file);
			System.out.println("Filename: " + file.getName());
			System.out.println("Query AvRE:\t\t" +
					Evaluation.calculateQueryError(originalDatabase, syntheticDatabase,
					generatedQueries, queryErrorSanityBound));
			System.out.println("Location coverage kendall-tau:\t\t" +
					Evaluation.locationCoverageKendallTau(originalDatabase, 
					syntheticDatabase, locCoverageGridSize));
			System.out.println("Frequent pattern F1:\t\t" +
					Evaluation.patternMiningF1Error(origPatterns, topK, 
					syntheticDatabase, ug, minFPsize, maxFPsize));
			System.out.println("Frequent pattern support:\t\t" +
					Evaluation.patternMiningSupport(origPatterns, topK, 
					syntheticDatabase, ug, minFPsize, maxFPsize));
			System.out.println("Trip error:\t\t" + 
					Evaluation.calculateTripError(originalDatabase, syntheticDatabase, 
					gridSizeForTripErr));
			System.out.println("Diameter error:\t\t" + 
					Evaluation.calculateDiameterError(originalDatabase, syntheticDatabase, 
					bucketCountForDiameterErr));
			System.out.println("Length error:\t\t" +
					Evaluation.calculateTravelDistanceError(originalDatabase, syntheticDatabase, 
					bucketCountForDiameterErr));
		}
	}

}
