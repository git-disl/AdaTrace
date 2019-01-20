import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class Cell {
	
	private double minX;
	private double maxX;
	private double minY;
	private double maxY;
	private String name = "";
	private List<Cell> level2cells; // for adaptive grid
	private List<Double> level2densities; // for adaptive grid
	
	public Cell (double minx, double miny, double xIncrement, double yIncrement, String nm) {
		this.minX = minx;
		this.minY = miny;
		this.maxX = minx + xIncrement;
		this.maxY = miny + yIncrement;
		this.name = nm;
		level2cells = new ArrayList<Cell>();
		level2densities = new ArrayList<Double>();
	}
	
	public boolean inCell(double xCoord, double yCoord) {
		if (xCoord >= this.minX && xCoord <= this.maxX && yCoord >= this.minY && yCoord <= this.maxY)
			return true;
		else
			return false;
	}
	
	public boolean inCell (Point p) {
		double xcoord = p.getX();
		double ycoord = p.getY();
		return inCell(xcoord, ycoord);
	}
	
	public String toString() {
		return this.name;
	}
	
	public String getName() {
		return this.name;
	}
	
	/*
	// samples and returns a random point in this cell 
	public Point sampleRandomPoint() {
		double xcoord = this.minX + (new Random().nextDouble())*(this.maxX-this.minX);
		double ycoord = this.minY + (new Random().nextDouble())*(this.maxY-this.minY);
		// rounding to 2 digits after decimal place for producing more readable output
		DecimalFormat df = new DecimalFormat("#.##");
		df.setRoundingMode(RoundingMode.HALF_UP);
		xcoord = Double.parseDouble(df.format(xcoord));
		ycoord = Double.parseDouble(df.format(ycoord));
		return new Point(xcoord, ycoord);
	} 
	*/
	
	
	// samples and returns a random point in this cell 
	public Point sampleRandomPoint() {
		if (level2cells.size() == 0 || level2cells.size() == 1) {
			double xcoord = this.minX + (new Random().nextDouble())*(this.maxX-this.minX);
			double ycoord = this.minY + (new Random().nextDouble())*(this.maxY-this.minY);
			// rounding to 2 digits after decimal place for producing more readable output
			DecimalFormat df = new DecimalFormat("#.##");
			df.setRoundingMode(RoundingMode.HALF_UP);
			xcoord = Double.parseDouble(df.format(xcoord));
			ycoord = Double.parseDouble(df.format(ycoord));
			return new Point(xcoord, ycoord);
		} else {
			Random r = new Random();
			double random = r.nextDouble();
			double seenSoFar = 0.0;
			for (int i = 0; i < this.level2densities.size(); i++) {
				seenSoFar = seenSoFar + this.level2densities.get(i);
				if (seenSoFar >= random) {
					// found it !
					Cell myAdaptiveCell = level2cells.get(i);
					double xcoord = myAdaptiveCell.minX + (new Random().nextDouble())*(this.maxX-this.minX);
					double ycoord = myAdaptiveCell.minY + (new Random().nextDouble())*(this.maxY-this.minY);
					DecimalFormat df = new DecimalFormat("#.##");
					df.setRoundingMode(RoundingMode.HALF_UP);
					xcoord = Double.parseDouble(df.format(xcoord));
					ycoord = Double.parseDouble(df.format(ycoord));
					return new Point(xcoord, ycoord);
				}
			}
			// !?!?
			double xcoord = this.minX + (new Random().nextDouble())*(this.maxX-this.minX);
			double ycoord = this.minY + (new Random().nextDouble())*(this.maxY-this.minY);
			// rounding to 2 digits after decimal place for producing more readable output
			DecimalFormat df = new DecimalFormat("#.##");
			df.setRoundingMode(RoundingMode.HALF_UP);
			xcoord = Double.parseDouble(df.format(xcoord));
			ycoord = Double.parseDouble(df.format(ycoord));
			return new Point(xcoord, ycoord);
		}
	}
	
	
	public String getBoundaries() {
		String n = "minX: " + this.minX + ", maxX: " + this.maxX + ", minY: " + this.minY +
				", maxY: " + this.maxY;
		return n;
	}
	
	public void divideFurther (double noisydensity, double EpsLeft, List<Trajectory> db) {
		//int lvl2cell = (int) Math.ceil(Math.sqrt(noisydensity*unusedEpsilon/200.0));
		int lvl2cell = (int) Math.ceil( (5*noisydensity/(db.size()*EpsLeft)) );		
		if (lvl2cell < 0)
			lvl2cell = 1;
		double xIncrement = (this.maxX-this.minX)/((double)lvl2cell);
		double yIncrement = (this.maxY-this.minY)/((double)lvl2cell);
		List<Double> densitiess = new ArrayList<Double>();
		for (int i = 0; i < lvl2cell; i++) {
			for (int j = 0; j < lvl2cell; j++) {
				Cell lvl2cellSelf = new Cell(this.minX+xIncrement*i, minY+yIncrement*j, xIncrement,
						yIncrement, ""+i+","+j);
				double lvl2cellDensity = 0.0;
				for (Trajectory t : db) {
					for (int k = 0; k < t.getSize(); k++) {
						if (lvl2cellSelf.inCell(t.getPoint(k))) {
							lvl2cellDensity = lvl2cellDensity + 1.0;
							break;
						}
					}
				}
				level2cells.add(lvl2cellSelf);
				densitiess.add(lvl2cellDensity);
			}
		}
		// normalize densities
		double totaldensity = 0.0;
		for (Double d : densitiess)
			totaldensity = totaldensity + d;
		for (int i = 0; i < densitiess.size(); i++) {
			level2densities.add(densitiess.get(i)/totaldensity);
		}
	}

}
