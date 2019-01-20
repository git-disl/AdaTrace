import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.Random;


public class Query {
	
	// CIRCULAR QUERY 

	public double centerX;
	public double centerY;
	public double radius;
	
	public Query (double minX, double maxX, double minY, double maxY, String filename) {
		// generate uniformly random query center between (minX-maxX) and (minY-maxY)
		Random r = new Random();
		this.centerX = r.nextDouble()*(maxX-minX)+minX; //<-- for uniform
		this.centerY = r.nextDouble()*(maxY-minY)+minY; //<-- for uniform
		double smallerDiff = -1.0;
		double largerDiff = -1.0;
		if (maxX-minX > maxY-minY) {
			smallerDiff = maxY-minY;
			largerDiff = maxX-minX;
		}
		else {
			smallerDiff = maxX-minX;
			largerDiff = maxY-minY;
		}
		
		// set the radius of query (dataset-dependent, based on observations)
		if (filename.contains("brinkhoff")) {
			this.radius = largerDiff/2.8;
		} else if (filename.contains("taxi")) {
			this.radius = largerDiff/3.0;
		} else if (filename.contains("geolife")) {
			this.radius = largerDiff/2.5;
		} else {
			//System.out.println("File not recognized when generating queries..");
			this.radius = largerDiff/3.0;
		}
		
	}
	
	private boolean evaluateQueryOnTraj(Trajectory t) {
		for (int i = 0; i < t.getSize()-1; i++) {
			double p1x = t.getPoint(i).getX();
			double p2x = t.getPoint(i+1).getX();
			double p1y = t.getPoint(i).getY();
			double p2y = t.getPoint(i+1).getY();
			if (Query.eucDistance(p1x, p1y, centerX, centerY) <= radius)
				return true;
			if (Query.eucDistance(p2x, p2y, centerX, centerY) <= radius) 
				return true;
			if (p1x == p2x && p1y == p2y) { return true; }  
			else if (Query.distanceToSegment(centerX, centerY, p1x, p1y, p2x, p2y) <= radius) { 
				//System.out.println("line in circle");
				return true;
			}
		}
		return false;
	}
	
	public int evaluateQueryOnDatabase(List<Trajectory> dataset) {
		int tbr = 0;
		for (Trajectory t : dataset) {
			if (this.evaluateQueryOnTraj(t) == true)
				tbr++;
		}
		return tbr;
	}
	
	
	// HELPER FUNCTIONS:
	
	public static double eucDistance (double p1x, double p1y, double p2x, double p2y) {
		return Math.sqrt((p2x-p1x)*(p2x-p1x) + (p2y-p1y)*(p2y-p1y));
	}
		
	public static double distanceToSegment(double x3, double y3, double x1, double y1, 
			double x2, double y2) {
		final Point2D p3 = new Point2D.Double(x3, y3);
		final Point2D p1 = new Point2D.Double(x1, y1);
		final Point2D p2 = new Point2D.Double(x2, y2);
		return distanceToSegment(p1, p2, p3);
	}

	
	 public static double distanceToSegment(Point2D p1, Point2D p2, Point2D p3) {
		 final double xDelta = p2.getX() - p1.getX();
		 final double yDelta = p2.getY() - p1.getY();
		 if ((xDelta == 0) && (yDelta == 0)) {
		    throw new IllegalArgumentException("p1 and p2 cannot be the same point");
		 }
		 final double u = ((p3.getX() - p1.getX()) * xDelta + 
				 (p3.getY() - p1.getY()) * yDelta) / (xDelta * xDelta + yDelta * yDelta);
			 final Point2D closestPoint;
		 if (u < 0) {
		    closestPoint = p1;
		 } else if (u > 1) {
		    closestPoint = p2;
		 } else {
		    closestPoint = new Point2D.Double(p1.getX() + u * xDelta, p1.getY() + u * yDelta);
		 }
		 return closestPoint.distance(p3);
	 }
	
	 
	 
	 /*
		// RECTANGULAR QUERY
		private Rectangle2D queryRegion;
		
		// Parameters to the constructor are the boundaries of dataset
		public Query(double minX, double maxX, double minY, double maxY) {
			
			Random r = new Random();
			double upperLeftX = r.nextDouble()*(maxX-minX) + minX;
			double upperLeftY = r.nextDouble()*(maxY-minY) + minY;
			
			// set query boundaries
			double xdiff = maxX-minX;
			double ydiff = maxY-minY;
			double width = xdiff/2.0;
			double height = ydiff/2.0;
			
			this.queryRegion = new Rectangle2D.Double(upperLeftX, upperLeftY, width, height);
			
			// Emre's note: This strategy of choosing the upper left corner yields to very bad
			// results, due to queries with very low counts
			
		}
		
		private boolean evaluateQueryOnTraj (Trajectory t) {
			for (int i = 0; i < t.getSize()-1; i++) {
				double p1x = t.getPoint(i).getX();
				double p2x = t.getPoint(i+1).getX();
				double p1y = t.getPoint(i).getY();
				double p2y = t.getPoint(i+1).getY();
				// if trajectory has a point in the rectangle, return true
				if (queryRegion.contains(p1x, p1y))
					return true;
				if (queryRegion.contains(p2x, p2y))
					return true;
				// if trajectory's line segment between p1-p2 passes thru rectangle, return true
				if (queryRegion.intersectsLine(p1x, p1y, p2x, p2y))
					return true;
			}
			return false;
		}
		
		public int evaluateQueryOnDatabase(List<Trajectory> dataset) {
			int tbr = 0;
			for (Trajectory t : dataset) {
				if (this.evaluateQueryOnTraj(t) == true)
					tbr++;
			}
			return tbr;
		}
	*/	
	 
}
