/*
**  File: Res.cs
**  Started: 
**  Contributors: Joanna Slusky, Meghan Franklin, Ryan Feehan
**  Overview: 
**
**  About: 
**
**  Last Edited: 
*/

using System;
using System.Windows.Media.Media3D;
using System.Collections.Generic;
using System.Collections;

namespace betaBarrelProgram
{

    public class Res : IEnumerable<Atom>
    {
        public int h_bonder { get; set; }
        public int h_bonderID { get; set; }
        public double aDist { get; set; }
        public double bDist { get; set; }
        public int ResNum { get; set; }
        public double Phi { get; set; }
        public double Psi { get; set; }
        public double Omega { get; set; }
        public string SSType { get; set; }
        public string ThreeLetCode { get; set; }
        public string OneLetCode { get; set; }
        public int SeqID { get; set; } //formerly authSeqID
        public double Tilt { get; set; }
        public Res previousStrandRes { get; set; }
        public Res nextStrandRes { get; set; }
        public Res ShearNumNeigh { get; set; }
        public List<int> Neighbors { get; set; }
        public List<Res> SideChainNeighbors { get; set; }
        public List<Res> BackboneNeighbors { get; set; }
        public int ChainNum { get; set; }
        public string ChainName { get; set; }
        public int StrandNum { get; set; }
        public double Dihedral { get; set; }
        public double Twist_next { get; set; }
        public double Twist_prev { get; set; }
        public int ResStrandNum { get; set; }// what number residue in the strand is it
        public double Coil2 { get; set; } // angle of 3 consecutive ca for the ca
        public double Coil1 { get; set; } // for the first ca
        public double Coil3 { get; set; } // for the third ca
        public double Radius { get; set; }  //distance from carbon alpha to axis
        public bool Inward { get; set; } // inward facing or outward facing
        public double Z { get; set; }

        public List<Atom> Atoms { get; set; }
        public Dictionary<string, Vector3D> BackboneCoords { get; set; } // N C CA
        public Vector3D Direction { get; set; } // the coordinates of this carbon alpha - coordinates of the previous carbon alpha
        public Vector3D SideChainDir { get; set; }

        //constructor
        public Res(ref AtomParser.AtomCategory _myAtomCat, int chainNum, int resNum, int FirstAtomNum)
        {
	        public int h_bonder { get; set; }
	        public int h_bonderID { get; set; }
	        public double aDist { get; set; }
	        public double bDist { get; set; }
	        public int ResNum { get; set; }
	        public double Phi { get; set; }
	        public double Psi { get; set; }
	        public double Omega { get; set; }
	        public double BFacCA { get; set; } //Added 6-2-17 for loop purposes
	        public double RelBFac { get; set; }
	        public string SSType { get; set; }
	        public string DSSP { get; set; }
	        public string ThreeLetCode { get; set; }
	        public string OneLetCode { get; set; }
	        public int SeqID { get; set; } //formerly authSeqID
	        public double Tilt { get; set; }
	        public Res previousStrandRes { get; set; }
	        public Res nextStrandRes { get; set; }
	        public Res ShearNumNeigh { get; set; }
	        public List<int> Neighbors { get; set; }
	        public List<Res> SideChainNeighbors { get; set; }
	        public List<Res> BackboneNeighbors { get; set; }
	        public int ChainNum { get; set; }
	        public string ChainName { get; set; }
	        public int StrandNum { get; set; }
	        public double Dihedral { get; set; }
	        public double Twist_next { get; set; }
	        public double Twist_prev { get; set; }
	        public int ResStrandNum { get; set; }// what number residue in the strand is it
	        public double Coil2 { get; set; } // angle of 3 consecutive ca for the ca
	        public double Coil1 { get; set; } // for the first ca
	        public double Coil3 { get; set; } // for the third ca
	        public double Radius { get; set; }  //distance from carbon alpha to axis
	        public bool Inward { get; set; } // inward facing or outward facing
	        public double Z { get; set; }

	        public List<Atom> Atoms { get; set; }
	        public Dictionary<string, Vector3D> BackboneCoords { get; set; } // N C CA
	        public Vector3D Direction { get; set; } // the coordinates of this carbon alpha - coordinates of the previous carbon alpha
	        public Vector3D SideChainDir { get; set; }

	        //constructor
	        public Res(ref AtomParser.AtomCategory _myAtomCat, int chainNum, int resNum, int FirstAtomNum)
	        {
	            this.Neighbors = new List<int>();
	            this.Atoms = new List<Atom>();
	            this.SideChainNeighbors = new List<Res>();
	            this.BackboneNeighbors = new List<Res>();
	            // the coordinates of this carbon alpha -coordinates of the previous carbon alpha
	            this.Direction = new Vector3D();
	            // N C CA
	            this.BackboneCoords = new Dictionary<string, Vector3D>();
            

	            this.SideChainDir = new Vector3D();
	            this.Z = 999;
	            this.ResNum = resNum; //count of residues - residue index
	            this.ChainNum = chainNum;
	            this.ChainName = _myAtomCat.ChainAtomList[chainNum].AuthAsymChain;
	            this.Twist_next = 999;
	            this.Twist_prev = 999;
	            this.Coil2 = 999;
	            this.Coil1 = 999;
	            this.Coil3 = 999;
	            this.Tilt = 999;
	            this.Radius = 999;
	            this.SSType = "X";
	            this.StrandNum = 99;
	            this.ThreeLetCode = _myAtomCat.ChainAtomList[chainNum].CartnAtoms[FirstAtomNum].residue;
	            if (this.ThreeLetCode == "MSE") this.ThreeLetCode = "MET";
	            this.OneLetCode = threeToOne(ThreeLetCode);
	            this.SeqID = Convert.ToInt32(_myAtomCat.ChainAtomList[chainNum].CartnAtoms[FirstAtomNum].authSeqId); //formerly authSeqID; resnum in the PDB

	            //for (int atomctr = FirstAtomNum; atomctr < _myAtomCat.ChainAtomList[chainNum].cartnAtoms.Count; atomctr++) //changed to a while loop on 6-14-16 - works EXPONENTIALLY faster
	            int atomctr = FirstAtomNum;
	            while (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].seqId == _myAtomCat.ChainAtomList[chainNum].CartnAtoms[FirstAtomNum].seqId)
	            {
	                if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].seqId == _myAtomCat.ChainAtomList[chainNum].CartnAtoms[FirstAtomNum].seqId)
	                {
	                    Vector3D coords = new Vector3D(_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].xyz.X, _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].xyz.Y, _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].xyz.Z);
	                    Atom myAtom = new Atom(this.ThreeLetCode, coords, atomctr, _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName, _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomType);
	                    Atoms.Add(myAtom); //collect info from myAtomCat

	                    if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName == "CA" || _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName == "C" || _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName == "N" || _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName == "O")
	                    {
	                        if (BackboneCoords.ContainsKey(_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName) != true)
	                        {
	                            this.BackboneCoords.Add(_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName, myAtom.Coords);
	                        }
	                        else
	                        {
	                            Console.WriteLine("IMPORTANT!!!!! Residue {0} has multiple ocupancy", this.SeqID);
	                        }

	                    }
	                    //Added 6-2-17 for loop purposes
	                    if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].atomName == "CA")
	                    {
	                        this.BFacCA = _myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomctr].bFac;
	                    }

	                }
	                atomctr++;
	                if (atomctr >= _myAtomCat.ChainAtomList[chainNum].cartnAtoms.Count) break;
	            }

            //Added 5-9-16 to determine which direction sidechain is pointing from barrel; also sets the ResSeqID for using with neighbor functions
            foreach (Atom atom1 in this.Atoms)
            {
                atom1.ResSeqID = this.SeqID;
                if (atom1.AtomName == "CB" && this.ThreeLetCode != "Gly")
                {
                    this.SideChainDir = atom1.Coords - this.Atoms[1].Coords;
                    //Console.WriteLine("{0}\t{1}", this.ThreeLetCode, this.SideChainDir);
                }
                else if (this.ThreeLetCode == "Gly")
                {
                    this.SideChainDir = this.Atoms[2].Coords - this.Atoms[0].Coords;
                    //Console.WriteLine("{0}\t{1}", this.ThreeLetCode, this.SideChainDir);
                }
                
            }

        }

        public void translate(Vector3D translationCoords)
        {
            for (int atomCtr = 0; atomCtr < this.Atoms.Count; atomCtr++)
            {
                this.Atoms[atomCtr].translate(translationCoords);
                if (this.Atoms[atomCtr].AtomName == "N" || this.Atoms[atomCtr].AtomName == "CA" || this.Atoms[atomCtr].AtomName == "C" || this.Atoms[atomCtr].AtomName == "O")
                {
                    this.BackboneCoords[this.Atoms[atomCtr].AtomName] += translationCoords;
                }
            }
        }

        public IEnumerator<Atom> GetEnumerator()
        {
            return Atoms.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }

        public string threeToOne(string threeLetters)
        {
            int threeIndex = three2OneTable.IndexOf(threeLetters);
            if (threeIndex == -1) // not found
                return "X";
            else
                return three2OneTable.Substring(threeIndex + 5, 1);

        }

        // member variables
        const string three2OneTable = "ALA -A CYS -C ASP -D GLU -E PHE -F GLY -G " +
    "HIS -H ILE -I LYS -K LEU -L MET -M ASN -N " +
    "PRO -P GLN -Q ARG -R SER -S THR -T VAL -V " +
    "TRP -W TYR -Y ASX -N GLX -Q UNK -X INI -K " +
    "AAR -R ACE -X ACY -G AEI -T AGM -R ASQ -D " +
    "AYA -A BHD -D CAS -C CAY -C CEA -C CGU -E " +
    "CME -C CMT -C CSB -C CSD -C CSE -C CSO -C " +
    "CSP -C CSS -C CSW -C CSX -C CXM -M CYG -C " +
    "CYM -C DOH -D EHP -F FME -M FTR -W GL3 -G " +
    "H2P -H HIC -H HIP -H HTR -W HYP -P KCX -K " +
    "LLP -K LLY -K LYZ -K M3L -K MEN -N MGN -Q " +
    "MHO -M MHS -H MIS -S MLY -K MLZ -K MSE -M " +
    "NEP -H NPH -C OCS -C OCY -C OMT -M OPR -R " +
    "PAQ -Y PCA -Q PHD -D PRS -P PTH -Y PYX -C " +
    "SEP -S SMC -C SME -M SNC -C SNN -D SVA -S " +
    "TPO -T TPQ -Y TRF -W TRN -W TRO -W TYI -Y " +
    "TYN -Y TYQ -Y TYS -Y TYY -Y YOF -Y FOR -X";

    }

    //These are referenced in neighbor-finding functions using energies.
    public class HDonor
    {
        public string resName { get; set; }
        public string atomName { get; set; }
        public string Connector { get; set; }
        public double vDWrad { get; set; }
        public double EpsMin { get; set; }

        public HDonor(List<string> splitLine)
        {
            this.atomName = splitLine[1];
            this.Connector = splitLine[2];
            this.vDWrad = 0.224500; //From CHARMM36 parameter file for proteins
            this.EpsMin = -0.046; //From CHARMM36 parameter file for proteins; all hydrogens added are labeled as type H
        }
    }

    public class HAcceptor
    {
        public string nextAtom { get; set; }
        public string dihedralAtom { get; set; }
        public double eAngle { get; set; }
        public double dihed1 { get; set; }
        public double dihed2 { get; set; }
        public double vDWrad { get; set; }
        public double EpsMin { get; set; }

        public HAcceptor(List<string> splitLine)
        {
            if (splitLine[1].Contains("O"))
            {
                //From CHARMM36 parameter file for proteins; no difference allowed in O vdw radii or Emin values
                this.vDWrad = 1.70;
                this.EpsMin = -0.12;
            }
            if (splitLine[1].Contains("N"))
            {
                //From CHARMM36 parameter file for proteins; no difference observed in N vdw radii
                this.vDWrad = 1.85;
                this.EpsMin = -0.2;
            }
            this.nextAtom = splitLine[2];
            this.dihedralAtom = splitLine[3];
            this.eAngle = Convert.ToDouble(splitLine[5]);
            this.dihed1 = Convert.ToDouble(splitLine[6]);
            this.dihed2 = Convert.ToDouble(splitLine[7]);
        }
    }

  }