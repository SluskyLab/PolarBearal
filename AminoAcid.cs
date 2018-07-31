using System.Collections.Generic;

public class AminoAcid
{
    public List<double> Tilt { get; set; }
    public List<double> TiltfromStrand { get; set; }
    public List<double> TiltfromBarrel { get; set; }
    public List<double> Tilt_even { get; set; }
    public List<double> TiltfromStrand_even { get; set; }
    public List<double> TiltfromBarrel_even { get; set; }
    public List<double> Tilt_odd { get; set; }
    public List<double> TiltfromStrand_odd { get; set; }
    public List<double> TiltfromBarrel_odd { get; set; }
    public List<int> seqIDList { get; set; }
    public List<List<double>> phiPsiList { get; set; }
    
    public List<double> caThetaList { get; set; }
    public List<double> Twist { get; set; }

    public AminoAcid()
	{
        this.Tilt = new List<double>();
        this.TiltfromStrand = new List<double>();
        this.TiltfromBarrel = new List<double>();
        this.Tilt_even = new List<double>();
        this.TiltfromStrand_even = new List<double>();
        this.TiltfromBarrel_even = new List<double>();
        this.Tilt_odd = new List<double>();
        this.TiltfromStrand_odd = new List<double>();
        this.TiltfromBarrel_odd = new List<double>();
        this.seqIDList = new List<int>();
        this.phiPsiList = new List<List<double>>();
        this.caThetaList = new List<double>();
        this.Twist = new List<double>();

	}
    public List<double> GetProperty(string myString)
    {
        if (myString == "Tilt") return Tilt;
        if (myString == "TiltfromStrand") return TiltfromStrand;
        if (myString == "TiltfromBarrel") return TiltfromBarrel;
        if (myString == "caThetaList") return caThetaList;
        if (myString == "Twist") return Twist;
        else
        {
            List<double> myList = new List<double>();
            return myList;
        }

    }
}
