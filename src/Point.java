
public class Point {
	
	private double xCoord; // longitude
	private double yCoord; // latitude
	
	public Point(double xc, double yc) {
		this.xCoord = xc;
		this.yCoord = yc;
	}
	
	public double getX() {
		return this.xCoord;
	}
	
	public double getY() {
		return this.yCoord;
	}
	
	public String toString() {
		return "(" + this.xCoord + "," + this.yCoord + ")";
	}
	
	public boolean isEqualTo(Point p2) {
		if (this.getX() != p2.getX())
			return false;
		if (this.getY() != p2.getY())
			return false;
		return true;
	}
	
	public void setX (double xc) {
		this.xCoord = xc;
	}
	
	public void setY (double yc) {
		this.yCoord = yc;
	}
	
	public double euclideanDistTo(Point p2) {
		double sum = (this.getY()-p2.getY())*(this.getY()-p2.getY()) +
				(this.getX()-p2.getX())*(this.getX()-p2.getX());
		return Math.sqrt(sum);
	}



}
