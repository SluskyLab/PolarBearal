/*
**  File: Atom.cs
**  Started: 
**  Contributors: Joanna Slusky, Meghan Franklin, Ryan Feehan
**  Overview: 
**
**  About: 
**
**  Last Edited: 
*/

using System;
using System.Collections.Generic;
using System.Windows.Media.Media3D;


namespace betaBarrelProgram
{
    public class Atom
    {
        public string AtomType { get; set; }
        public string AtomName { get; set; }
        public Vector3D Coords { get; set; }
        public int AtomNum { get; set; }
        public int ResSeqID { get; set; }
        public Vector3D Hydrogen { get; set; }
        public Vector3D e1 { get; set; }
        public Vector3D e2 { get; set; }
        public double partialcharge { get; set; }
        public List<Atom> SCSCNeighAtoms { get; set; }
        public List<Atom> SCBBNeighAtoms { get; set; }
        public List<Atom> BBNeighAtoms { get; set; }

        public Atom(string resName, Vector3D coords, int atomNum, string atomName, string atomType)
        {
            this.AtomNum = atomNum;
            this.Coords = coords;
            this.AtomName = atomName;
            this.AtomType = atomType;

            Tuple<string, string> key = new Tuple<string, string>(resName, atomName);
            if (Global.partialChargesDict.ContainsKey(key) == true) this.partialcharge = Global.partialChargesDict[key];
            else if (atomName == "H") this.partialcharge = 0.31;
            else if (atomName == "N") this.partialcharge = -0.47;
            else if (atomName == "C") this.partialcharge = 0.51;
            else if (atomName == "O") this.partialcharge = -0.51;
            else this.partialcharge = 0;

            this.SCSCNeighAtoms = new List<Atom>();
            this.SCBBNeighAtoms = new List<Atom>();
            this.BBNeighAtoms = new List<Atom>();

        }

        public void translate(Vector3D translationVec)
        {

            this.Coords += translationVec;

        }
    }
}