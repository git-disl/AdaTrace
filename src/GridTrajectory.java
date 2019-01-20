import java.util.ArrayList;
import java.util.List;


public class GridTrajectory {
	
private List<Cell> trajCells;
	
	public GridTrajectory (Trajectory t, Grid g, boolean interpWanted) {
		trajCells = new ArrayList<Cell>();
		// convert trajectory t to a gridtrajectory using grid g
		List<Cell> gridCells = g.getCells();
		for (int i = 0; i < t.getSize(); i++) {
			Point p = t.getPoint(i);
			boolean found = false;
			for (Cell c : gridCells) {
				if (c.inCell(p)) {
					trajCells.add(c);
					found = true;
					break;
				}
			}
			if (!found) {
				// this happens due to double precision loss in java. if it happens this 
				// coordinate belongs to a corner cell (on the rightmost column 
				// or uppermost row) of the grid
				double xcoord = p.getX() - 0.001;
				double ycoord = p.getY() - 0.001;
				for (Cell c : gridCells) {
					if (c.inCell(xcoord,  ycoord)) {
						trajCells.add(c);
						found = true;
						break;
					}
				}
			}
			if (!found) {
				// this happens on occasion for the synthetic trajectories published by 
				// [He et al 2015], where the coordinates end up being (0.0, 0.0) or 
				// (1.0, 0.0). solution: add to bottom-left corner of the grid
				trajCells.add(g.getCellMatrix()[0][0]);
				found = true;
			}
		}
		// removal of duplicates phase
		List<Cell> newTrajCells = new ArrayList<Cell>();
		newTrajCells.add(this.trajCells.get(0));
		for (int i = 1; i < this.trajCells.size()-1; i++) {
			if (this.trajCells.get(i).equals(newTrajCells.get(newTrajCells.size()-1))) {
				// do nothing
			} else {
				newTrajCells.add(this.trajCells.get(i));
			}
		}
		try {
			if (this.trajCells.get(this.trajCells.size()-1) != 
					this.trajCells.get(this.trajCells.size()-2)) { // dont duplicate last cell
				newTrajCells.add(this.trajCells.get(this.trajCells.size()-1));
			}
		} catch (Exception e) {
			// do nothing - this try-catch is for DPT, as I got an error during exprmntation
		}
		if (newTrajCells.size() == 1) { // to fix trip distribution
			newTrajCells.add(this.trajCells.get(this.trajCells.size()-1));
		}
		this.trajCells = newTrajCells;
		
		// END - interpolate trajectory cells if wanted&required
		if (interpWanted) {
			List<Cell> finTrajCells = new ArrayList<Cell>();
			for (int i = 0; i < this.trajCells.size()-1; i++) {
				Cell current = this.trajCells.get(i);
				Cell next = this.trajCells.get(i+1);
				if (current == next || g.areAdjacent(current,next)) { 
					// no interp needed
					finTrajCells.add(current);
					//finTrajCells.add(next);
				} 
				else {
					// interp needed
					//System.out.println("Interp needed btwn " + current.toString() + " and " +
					//		next.toString());
					finTrajCells.addAll(g.giveInterpolatedRoute(current,next));
				}
			}
			finTrajCells.add(this.trajCells.get(this.trajCells.size()-1));
			// DEBUG-START
			//if (!finTrajCells.equals(this.trajCells)) {
			//	System.out.println("Before: " + this.trajCells);
			//	System.out.println("After: " + finTrajCells);
			//}
			// DEBUG-END
			this.trajCells = finTrajCells;
		}
	}
	
	// for ShokriConverter.java - no removal of duplicates
	public GridTrajectory (Trajectory t, Grid g, boolean interpWanted, boolean removal) {
		trajCells = new ArrayList<Cell>();
		// convert trajectory t to a gridtrajectory using grid g
		List<Cell> gridCells = g.getCells();
		for (int i = 0; i < t.getSize(); i++) {
			Point p = t.getPoint(i);
			boolean found = false;
			for (Cell c : gridCells) {
				if (c.inCell(p)) {
					trajCells.add(c);
					found = true;
					break;
				}
			}
			if (!found) {
				// this happens due to double precision loss in java. if it happens this 
				// coordinate belongs to a corner cell (on the rightmost column 
				// or uppermost row) of the grid
				double xcoord = p.getX() - 0.001;
				double ycoord = p.getY() - 0.001;
				for (Cell c : gridCells) {
					if (c.inCell(xcoord,  ycoord)) {
						trajCells.add(c);
						found = true;
						break;
					}
				}
			}
			if (!found) {
				// this happens on occasion for the synthetic trajectories published by 
				// [He et al 2015], where the coordinates end up being (0.0, 0.0) or 
				// (1.0, 0.0). solution: add to bottom-left corner of the grid
				trajCells.add(g.getCellMatrix()[0][0]);
				found = true;
			}
		}
				
		// END - interpolate trajectory cells if wanted&required
		if (interpWanted) {
			List<Cell> finTrajCells = new ArrayList<Cell>();
			for (int i = 0; i < this.trajCells.size()-1; i++) {
				Cell current = this.trajCells.get(i);
				Cell next = this.trajCells.get(i+1);
				if (current == next || g.areAdjacent(current,next)) { 
					// no interp needed
					finTrajCells.add(current);
					//finTrajCells.add(next);
				} 
				else {
					// interp needed
					//System.out.println("Interp needed btwn " + current.toString() + " and " +
					//		next.toString());
					finTrajCells.addAll(g.giveInterpolatedRoute(current,next));
				}
			}
			finTrajCells.add(this.trajCells.get(this.trajCells.size()-1));
			// DEBUG-START
			//if (!finTrajCells.equals(this.trajCells)) {
			//	System.out.println("Before: " + this.trajCells);
			//	System.out.println("After: " + finTrajCells);
			//}
			// DEBUG-END
			this.trajCells = finTrajCells;
		}
	}
	
	public GridTrajectory () {
		this.trajCells = new ArrayList<Cell>();
	}
		
	
	public GridTrajectory (Cell startcell, Cell endcell) {
		trajCells = new ArrayList<Cell>();
		trajCells.add(startcell);
		trajCells.add(endcell);
	}
	
	public GridTrajectory (Cell[] cellerino) {
		trajCells = new ArrayList<Cell>();
		for (Cell s : cellerino) 
			trajCells.add(s);
	}
	
	public GridTrajectory (List<Cell> cellerino) {
		this.trajCells = cellerino;
	}
	
	public void addCell (Cell c1) {
		this.trajCells.add(c1);
	}
	
	public String toString() {
		String tbr = "";
		for (Cell c : this.trajCells)
			tbr = tbr + " -> " + c.toString();
		return tbr;
	}

	public List<Cell> getCells() {
		return this.trajCells;
	}
	
	public boolean passesThrough(Cell c) {
		for (int i = 0; i < this.trajCells.size(); i++) {
			if (c == this.trajCells.get(i))
				return true;
		}
		return false;
	}

}
