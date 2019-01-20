import java.util.ArrayList;
import java.util.List;


public class Trajectory {
	
	private List<Point> points;
	
	public Trajectory() {
		this.points = new ArrayList<Point>();
	}
	
	public Trajectory(List<Double> xs, List<Double> ys) {
		if (xs.size() != ys.size()) {
			System.out.println("Couldnt create trajectory since xcoordinates size != y.");
			return;
		}
		this.points = new ArrayList<Point>();
		for (int i = 0; i < xs.size(); i++) {
			Point p = new Point(xs.get(i), ys.get(i));
			this.points.add(p);
		}
	}
	
	public void addCoordinates(double xc, double yc) {
		Point p = new Point(xc, yc);
		this.points.add(p);
	}
	
	public void addPoint (Point inpP) {
		this.points.add(inpP);
	}
	
	public int getSize() {
		return points.size();
	}
	
	public String getPointString(int pos) {
		return this.points.get(pos).toString();
	}
	
	public Point getPoint (int pos) {
		return this.points.get(pos);
	}
	
	public double getMinXCoord() {
		double currMin = this.points.get(0).getX();
		for (Point p : this.points) {
			if (p.getX() < currMin) {
				currMin = p.getX();
			}
		}
		return currMin;
	}
	
	public double getMinYCoord() {
		double currMin = this.points.get(0).getY();
		for (Point p : this.points) {
			if (p.getY() < currMin) {
				currMin = p.getY();
			}
		}
		return currMin;
	}
	
	public double getMaxXCoord() {
		double currMax = this.points.get(0).getX();
		for (Point p : this.points) {
			if (p.getX() > currMax) {
				currMax = p.getX();
			}
		}
		return currMax;
	}
	
	public double getMaxYCoord() {
		double currMax = this.points.get(0).getY();
		for (Point p : this.points) {
			if (p.getY() > currMax) {
				currMax = p.getY();
			}
		}
		return currMax;
	}
	
	public double getDiameter() {
		double maxDist = 0.0;
		for (int i = 0; i < this.getSize(); i++) {
			for (int j = i+1; j < this.getSize(); j++) {
				Point p1 = this.points.get(i);
				Point p2 = this.points.get(j);
				double dist = p1.euclideanDistTo(p2);
				if (dist > maxDist)
					maxDist = dist;
			}
		}
		return maxDist;
	}
	
	public double getDistanceTravelled() {
		double tbr = 0.0;
		for (int i = 0; i < this.points.size()-1; i++) {
			Point thisPoint = this.points.get(i);
			Point nextPoint = this.points.get(i+1);
			tbr = tbr + thisPoint.euclideanDistTo(nextPoint);
		}
		return tbr;
	}
	
	public List<Point> getSniffedPoints (Cell sniffZone) {
		List<Point> tbr = new ArrayList<Point>();
		int pos = 0;
		while (pos < this.getSize()) {
			if (sniffZone.inCell(this.points.get(pos))) {
				// in the sniff zone now
				int sniffpos = pos;
				while (sniffpos < this.getSize() && 
						sniffZone.inCell(this.points.get(sniffpos))) {
					tbr.add(this.points.get(sniffpos));
					sniffpos++;
				}
				return tbr;
			}
			pos++;
		}
		return tbr;
	}
	
	public List<Point> calcIntersectionWith (Trajectory t2) {
		List<Point> intersected = new ArrayList<Point>();
		for (Point p : t2.points) {
			for (int i = 0; i < this.points.size(); i++) {
				Point cand = this.points.get(i);
				if (p.euclideanDistTo(cand) < 0.001) {
					intersected.add(p);
					break;
				}
			}
		}
		return intersected;
	}
	
	public int calcIntersectionCount (Trajectory t2) {
		return this.calcIntersectionWith(t2).size();
	}

	public double calculateDTWto(Trajectory t2) {
		return Evaluation.calculateDTW(this.points, t2.points);
	}
	
}
