using System;

namespace betaBarrelProgram.AtomParser
{
	/// <summary>
	/// Summary description for Coordinate.
	/// </summary>
	public class Coordinate
	{
		public double X = 0.0;
		public double Y = 0.0;
		public double Z = 0.0;

		public Coordinate ()
		{
		}

		public Coordinate(double xx, double yy, double zz)
		{
			X = xx;
			Y = yy;
			Z = zz;
		}
	}
}
