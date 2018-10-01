/*
**  File: Poly.cs
**  Started: 
**  Contributors: Meghan Franklin, Ryan Feehan
**  Overview: 
**
**  About: 
**
**  Last Edited: 
*/

using System;
using System.Collections.Generic;
using betaBarrelProgram.BarrelStructures;
using System.Linq;
using System.Windows.Media.Media3D;
using System.Collections;

namespace betaBarrelProgram
{
    namespace Poly
    {
        public class PolyProtein : Protein
        {


            public List<Chain> Chains { get; set; }
            public int ChainCount { get; set; }
            public int totalResNum { get; set; }
            public string PdbName { get; set; }


            public PolyProtein(ref AtomParser.AtomCategory _myAtomCat, string PdbName)
            {
                this.ChainCount = 0;
                this.totalResNum = 0;
                this.PdbName = PdbName.Substring(0, 4).ToUpper();
                this.Chains = new List<BarrelStructures.Chain>();

                //Console.WriteLine("In protein Class");
                for (int chainNum = 0; chainNum < _myAtomCat.chainAtomsList.Count; chainNum++)
                {

                    bool IsItProtein = false;
                    for (int atomNum = 0; atomNum < _myAtomCat.ChainAtomList[chainNum].cartnAtoms.Count; atomNum++)
                    {
                        if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomNum].atomName == "CA") IsItProtein = true; //if there is a CA anywhere in chain, create new Chain
                    }
                    if (IsItProtein == true)
                    {
                        //Console.WriteLine("Creating new chain");
                        Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, false, Global.POLY_DB_DIR);
                        this.Chains.Add(myChain);
                        this.ChainCount++;
                    }
                }
                this.Chains.Sort((x, y) => x.ChainName.CompareTo(y.ChainName)); //This sorts chains by letter - dimers don't necessarily add in order.
                
                for (int i = 0; i < this.ChainCount; i++)
                {
                    this.Chains[i].ChainNum = i;
                    foreach (Res res in this.Chains[i].Residues) res.ChainNum = i;
                }
            }

            public PolyProtein(ref AtomParser.AtomCategory _myAtomCat, int chainNum, string PdbName)
            {
                this.Chains = new List<BarrelStructures.Chain>();
                Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, false, Global.POLY_DB_DIR);
                this.Chains.Add(myChain);
            }

            public void translate(Vector3D translationCoords)
            {
                for (int chainCtr = 0; chainCtr < this.Chains.Count; chainCtr++)
                {
                    this.Chains[chainCtr].translate(translationCoords);
                }
            }
        }

        public class PolyBarrel : Barrel
        {
            public List<Strand> Strands { get; set; }
            public List<Vector3D> NellipseCoords { get; set; }
            public List<Vector3D> CellipseCoords { get; set; }
            public Vector3D Ncentroid { get; set; }
            public Vector3D Ccentroid { get; set; }
            public Vector3D Axis { get; set; }
            public bool Direction { get; set; }
            public double AvgTilt { get; set; }
            public double AvgTilt_even { get; set; }
            public double AvgTilt_odd { get; set; }
            public List<List<List<int>>> protoBarrel { get; set; }
            public double AvgRadius { get; set; }
            public double MinRadius { get; set; }
            public double MaxRadius { get; set; }
            public List<Res> LoopResies { get; set; }
            public string PdbName { get; set; }
            public List<double> StrandLength { get; set; }
            public Vector3D OriginalNcentroid { get; set; }
            public Vector3D OriginalCcentroid { get; set; }
            public Vector3D AxisVector { get; set; }
            public Vector3D NewCaxisPt { get; set; }
            public Vector3D OldCaxisPt { get; set; }
            public int ShearNum { get; set; }
            public List<double> PrevTwists { get; set; }

	        public PolyBarrel(PolyProtein _myProtein, string outputDirectory, string databaseDirectory) // send a protein class instead of chains
	        {
	            this.protoBarrel = new List<List<List<int>>>();
	            this.Strands = new List<Strand>();
	            this.PdbName = _myProtein.PdbName;
	            List<string> outwardStrands = new List<string> { "4MT4", "4K7R", "4MT0", "3D5K", "1EK9", "1YC9", "2LME", "3X2R", "2GR7", "5AZP", "5AZS" }; //This set of pdb ID's have Astrand0 that point outwards of the cell
	            List<string> outside_insertion = new List<string> { "7AHL", "1UUN", "3W9T", "4H56", "4TW1", "3B07", "3O44", "5GAQ", "3J9C" }; //These are the set of exterior inserted barrels for scaling Z-coords appropriately

	            #region Use this to output any info about chains you need before starting to create strands and alter information
	            /*string fileLocation2 = outputDirectory + "PhiPsiAngles\\" + this.PdbName + ".txt";

	                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
	                {
	                    file.WriteLine("the output is \n");
	                    for (int chain = 0; chain < _myProtein.ChainCount; chain++ )
	                    {
	                        file.WriteLine("Chain {0}", _myProtein.Chains[chain].ChainName);
	                        foreach (Res res in _myProtein.Chains[chain])
	                        {
	                            if (res.ResNum != 0)
	                            {
	                                file.WriteLine("Residue {0}({1}), Phi {4}, Psi {5}, SSType {2}, Angle formed at this resiude {3}", res.SeqID, res.ThreeLetCode, res.SSType, (Vector3D.AngleBetween(_myProtein.Chains[chain].Residues[res.ResNum - 1].Direction, res.Direction)), res.Phi, res.Psi);
	                            }
                            
	                        }
	                    }
	                }*/
	            #endregion

	            createAllStrands(ref _myProtein, outputDirectory); //This also changes all SStypes to B for H-bond searching to be more effective; we can reference back to psi/phi angles for orig type if necessary
            
	            if (protoBarrel[0].Count > 4) makeBarrelCircular(ref _myProtein); //remove extra chains if more than 4; known 4x3 proteins have no extra chains, 12-10-15
            
	            //most functions are based on the Strand list, so this region must work properly for any other functions to work properly
	            #region create Strand list 
	            for (int j = 0; j < protoBarrel.Count; j++)
	            {
	                for (int strandCtr = 0; strandCtr < protoBarrel[j].Count; strandCtr++)
	                {
	                    int start = _myProtein.Chains[j].Residues.FindIndex(a => a.SeqID == protoBarrel[j][strandCtr][0]);
	                    int end = _myProtein.Chains[j].Residues.FindIndex(a => a.SeqID == protoBarrel[j][strandCtr].Last());
	                    Strand newStrand = new Strand(_myProtein.Chains[j], start, end, strandCtr);
	                    this.Strands.Add(newStrand);
	                }
	            }
	            #endregion

	            #region Rotate barrel so axis is on z-axis

	            getEllipseCoords(outwardStrands);

	            this.Axis = this.Ccentroid - this.Ncentroid;
                var radii = SharedFunctions.setRadius(this.Strands, this.Axis, this.Ccentroid, this.Ncentroid);
                this.AvgRadius = radii.Item1;
                this.MaxRadius = radii.Item3;
                this.MinRadius = radii.Item2;

	            this.AxisVector = this.Ccentroid - this.Ncentroid; //getNormal(this.NellipseCoords, this.Ncentroid); //This is used for asymmetric barrel ellipses

	            this.OriginalCcentroid = this.Ccentroid;
	            this.OriginalNcentroid = this.Ncentroid;
	            this.NewCaxisPt = ((this.Ncentroid - this.Ccentroid).Length * this.AxisVector) + this.Ncentroid;
	            int axisDirection = 1;
	            if ((this.NewCaxisPt - this.OriginalCcentroid).Length > ((((this.Ncentroid - this.Ccentroid).Length * -1 * this.AxisVector) + this.Ncentroid) - this.OriginalCcentroid).Length)
	            {
	                axisDirection = -1;
	                this.NewCaxisPt = (((this.Ncentroid - this.Ccentroid).Length * -1 * this.AxisVector) + this.Ncentroid);
	            }

	            getEllipseCoords(outwardStrands);
            
	            this.Axis = this.NewCaxisPt - this.Ncentroid;
	            this.OldCaxisPt = this.NewCaxisPt;
	            //this.Axis = this.Ccentroid - this.Ncentroid; july3
	            this.rotateToZ(outwardStrands, _myProtein);
	            centerZ(outside_insertion, outwardStrands, _myProtein);
	            getEllipseCoords(outwardStrands);

	            this.AxisVector = this.Ccentroid - this.Ncentroid;
	            this.Axis = this.Ccentroid - this.Ncentroid;

	            //this.StrandLength = SharedFunctions.getStrandLengths(this.Strands, outputDirectory, this.PdbName);
            
	            #endregion
            
	            /*string fileLocation12 = outputDirectory + "ZCoords\\AllCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation12))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Z-coord";
	                file.WriteLine(newLine);
	                foreach (Strand strand in this.Strands)
	                {
	                    foreach (Res res in strand)
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            if (atom.AtomName == "CA") file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", res.ThreeLetCode, res.SeqID, strand.StrandNum, strand.ChainName, atom.Coords.X, atom.Coords.Y, atom.Coords.Z);
	                        }
	                    }
	                }
	            }*/

	            /*string fileLocation6 = outputDirectory + "\\Renumb\\XYZCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation6))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Chain" + "\t" + "Z-coord";
	                file.WriteLine(newLine);

	                foreach (Chain chain in _myProtein.Chains)
	                {
	                    foreach (Res res in chain)
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            if (atom.AtomName == "CA") file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", res.ThreeLetCode, res.ResNum +1, res.ChainName, atom.Coords.X, atom.Coords.Y, atom.Coords.Z);
	                        }
	                    }
	                }
	            }*/
	            //SharedFunctions.setInOut(this.Strands, outputDirectory, this.PdbName, this.Axis, this.Ccentroid, this.Ncentroid);
	            //SharedFunctions.writePymolScriptForStrands(this.Strands, outputDirectory, Program.MacpolyDBDir, this.PdbName);
	            //writeAminoAcidsTypesToFile(ref _myProtein, outputDirectory);
	            //this.AvgTilt = SharedFunctions.getTiltsByAA(this.Strands, outputDirectory, this.PdbName, this.Axis, ref Program.AADict);
            
	            //Check definitions for loops using strand lists
	            //checkStrandDefnsDSSP(this.Strands, ref _myProtein); //SO MANY PROBLEMS. I can hand-edit these loops.
            
	            for (int chain = 0; chain < _myProtein.ChainCount; chain ++ )
	            {
	                //3b07 and 4tw1 are both heteromers, so need both chains
	                Chain this_chain = _myProtein.Chains[chain];
	                if ((this_chain.ChainName == "E" | this_chain.ChainName == "F") & this.PdbName.ToUpper() == "4TW1")
	                {
	                    Dictionary<string, string> Loops = SharedFunctions.getLoopTurns(this.Strands, ref this_chain, outputDirectory, this.PdbName);
	                    //SharedFunctions.writePymolScriptForLoops(Loops, outputDirectory, Program.MacpolyDBDir, ref this_chain, this.PdbName);
	                    //SharedFunctions.findLoopsHBondingPartnersGeomOnly(Loops, outputDirectory, ref this_chain, this.PdbName, true);
	                    //this.PrevTwists = SharedFunctions.writeTwists(this.Strands, outputDirectory, this.PdbName);
	                }
	                else if ((this_chain.ChainName == "A" | this_chain.ChainName == "B") & this.PdbName.ToUpper() == "3B07")
	                {
	                    Dictionary<string, string> Loops = SharedFunctions.getLoopTurns(this.Strands, ref this_chain, outputDirectory, this.PdbName);
	                    //SharedFunctions.writePymolScriptForLoops(Loops, outputDirectory, Program.MacpolyDBDir, ref this_chain, this.PdbName);
	                    //SharedFunctions.findLoopsHBondingPartnersGeomOnly(Loops, outputDirectory, ref this_chain, this.PdbName, true);
                        //this.PrevTwists = SharedFunctions.writeTwists(this.Strands, outputDirectory, this.PdbName);
	                }
	                else if (this_chain.ChainName == "A")
	                {
	                    Dictionary<string, string> Loops = SharedFunctions.getLoopTurns(this.Strands, ref this_chain, outputDirectory, this.PdbName);
	                    //SharedFunctions.writePymolScriptForLoops(Loops, outputDirectory, Program.MacpolyDBDir, ref this_chain, this.PdbName);
	                    //SharedFunctions.findLoopsHBondingPartnersGeomOnly(Loops, outputDirectory, ref this_chain, this.PdbName, true);
                        //this.PrevTwists = SharedFunctions.writeTwists(this.Strands, outputDirectory, this.PdbName);
	                }

	            }
	            //getTyrVector(outputDirectory, ref _myProtein);

	            sortStrandsAroundEllipse(); //This function allows ease of comparison of strands around the barrel in physical order, rather than in input order, so this function needs to work for any function that finds interactions to properly compare neighboring strands

	            //SharedFunctions.findHBondingPartnersEnergy(this.Strands, outputDirectory, this.PdbName); //SQRWL-like implementation of hydrogen bonding
	            //SharedFunctions.findNearNeighDistOnly(this.Strands, outputDirectory, this.PdbName);

	            /*foreach (Strand strand in this.Strands)
	            {
	                foreach (Res res in strand)
	                {
	                    res.SideChainNeighbors = new List<Res>();
	                    res.BackboneNeighbors = new List<Res>();
	                }
	            }*/
	            //SharedFunctions.findNearestNeighbors(this.Strands, outputDirectory, this.PdbName); //based on CA distances
	            //SharedFunctions.findHBondingPartnersGeomOnly(this.Strands, outputDirectory, this.PdbName, true); //Modelling program-like implementation of hydrogen bonding
                       

	        }//End of class PBarrel

	        public void writeAminoAcidsTypesToFile(ref PolyProtein _myProtein, string outputDirectory)
	        {
	            string fileLocation3 = outputDirectory + "AminoAcids\\AminoAcidTypes" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation3))
	            {
	                file.WriteLine("Sl.No\t*PDB ID*\t*Strands Per Chain*\t*Chains Per Protein*\t*Residue Number*\t*AminoAcid*\t*Chain*\t*Strand Number*\t*Interface*");

	                int sl_no = 1; int res_num1; int res_num2; int range;
	                bool interface_value;
	                for (int j = 0; j < protoBarrel.Count; j++)
	                {
	                    var first_residue = _myProtein.Chains[j].Residues[0].SeqID;

	                    for (int strandCtr = 0; strandCtr < protoBarrel[j].Count; strandCtr++)
	                    {
	                        if (strandCtr == 0 || strandCtr == protoBarrel[j].Count - 1) // to check the residues that are at interface of the different chains
	                        {
	                            interface_value = true;
	                        }
	                        else interface_value = false;

	                        res_num1 = _myProtein.Chains[j].Residues.FindIndex(a => a.SeqID == protoBarrel[j][strandCtr][0]);
	                        res_num2 = _myProtein.Chains[j].Residues.FindIndex(a => a.SeqID == protoBarrel[j][strandCtr].Last());
	                        range = res_num2 - res_num1;
	                        for (int residues = res_num1; residues <= res_num2; residues++)
	                        {
	                            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", sl_no, this.PdbName, protoBarrel[0].Count, protoBarrel.Count, _myProtein.Chains[j].Residues[residues].SeqID, _myProtein.Chains[j].Residues[residues].ThreeLetCode, _myProtein.Chains[j].ChainName, strandCtr, interface_value);
	                            sl_no++;
	                        }
	                    }
	                }
	            }
	        }

	        public void makeBarrelCircular(ref PolyProtein _myProtein) //This is based on the premise that H-bonding between beta sheets in chains only occurs in the barrel
	        {
	            //fill in the barrel
	            int seqnum;
	            int nextchain;
	            int prevchain;
	            double d; double minD = 3.00;
	            int i; //chain count
	            int strandCtr; //strand count

	            for (i = 0; i < protoBarrel.Count; i++)//for each chain
	            {
	                for (strandCtr = 0; strandCtr < protoBarrel[i].Count; )//for each strand strandCtr++
	                {
	                    bool HbondedPrev = false;
	                    bool HbondedNext = false;

	                    if (i == 0)
	                    {
	                        prevchain = _myProtein.Chains.Count - 1;
	                        nextchain = i + 1;
	                    }
	                    else if (i == _myProtein.Chains.Count - 1)
	                    {
	                        prevchain = i - 1;
	                        nextchain = 0;
	                    }
	                    else
	                    {
	                        prevchain = i - 1;
	                        nextchain = i + 1;
	                    }

	                    #region Check for any neighbors in other chains
	                    for (int res1ctr = 0; res1ctr < protoBarrel[i][strandCtr].Count-1; res1ctr++)//for each res in strand
	                    {
	                        seqnum = protoBarrel[i][strandCtr][res1ctr];
	                        Res Residue1 = _myProtein.Chains[i].Residues[_myProtein.Chains[i].Residues.FindIndex(a => a.SeqID == seqnum)];

	                        for (int atomCtr1 = 0; atomCtr1 < Residue1.Atoms.Count; atomCtr1++)
	                        {
	                            if (Residue1.Atoms[atomCtr1].AtomName == "N" && Residue1.ThreeLetCode != "PRO")
	                            {
	                                Atom N = Residue1.Atoms[atomCtr1];

	                                for (int res3ctr = 0; res3ctr < _myProtein.Chains[prevchain].ResidueCount; res3ctr++) //Check previous chain for neighbors
	                                {
	                                    Res Residue2 = _myProtein.Chains[prevchain].Residues[res3ctr];

	                                    for (int atomCtr2 = 0; atomCtr2 < Residue2.Atoms.Count; atomCtr2++)
	                                    {
	                                        if (Residue2.Atoms[atomCtr2].AtomName == "O")
	                                        {
	                                            Atom O = Residue2.Atoms[atomCtr2];
	                                            atomCtr2 = Residue2.Atoms.Count;  // to skip the loop once "O" atom is found, and thus go to next residue number

	                                            d = (O.Coords - N.Hydrogen).Length;
	                                            if (d < minD && (Residue1.SSType == "B" && Residue2.SSType == "B"))
	                                            {
	                                                HbondedPrev = true;
	                                                goto CheckNextChain;
	                                            }
	                                        }
	                                    }
	                                }
	                            CheckNextChain:
	                                for (int res3ctr = 0; res3ctr < _myProtein.Chains[nextchain].ResidueCount; res3ctr++) //Check next chain for neighbors
	                                {
	                                    Res Residue2 = _myProtein.Chains[nextchain].Residues[res3ctr];
	                                    for (int atomCtr2 = 0; atomCtr2 < Residue2.Atoms.Count; atomCtr2++)
	                                    {
	                                        if (Residue2.Atoms[atomCtr2].AtomName == "O")
	                                        {

	                                            Atom O = Residue2.Atoms[atomCtr2];
	                                            atomCtr2 = Residue2.Atoms.Count;  // to skip the loop once "O" atom is found, and thus go to next residue number

	                                            d = (O.Coords - N.Hydrogen).Length;
	                                            if (d < minD && (Residue1.SSType == "B" && Residue2.SSType == "B"))
	                                            {
	                                                HbondedNext = true;
	                                                goto RemoveNonHBondedStrand;
	                                            }
	                                        }
	                                    }
	                                }
	                            }
	                        }

	                    }
	                    #endregion

	                RemoveNonHBondedStrand:
	                    //"3x2r" is an odd structure with strands that stretch into the upper "hair" region
	                    if (HbondedPrev == false && HbondedNext == false && _myProtein.PdbName != "3X2R") 
	                    {
	                        protoBarrel[i].RemoveRange(strandCtr, 1);
	                        //Console.WriteLine("removing strand {0}", strandCtr);
	                    }
	                    else strandCtr++;
	                    HbondedNext = false;
	                    HbondedPrev = false;
	                }
	                if (_myProtein.PdbName == "3X2R" && _myProtein.Chains[i].ChainName != "C")
	                {
	                    protoBarrel[i][3].RemoveRange(9, protoBarrel[i][3].Count - 9);
	                    protoBarrel[i][4].RemoveRange(0, 9);
	                    protoBarrel[i].RemoveRange(0, 2);
                    
	                }
	                if (_myProtein.PdbName == "3X2R" && _myProtein.Chains[i].ChainName == "C")
	                {
	                    protoBarrel[i][2].RemoveRange(9, protoBarrel[i][2].Count - 9);
	                    protoBarrel[i][3].RemoveRange(0, 9);
	                    protoBarrel[i].RemoveRange(0, 1);
                    
	                }
	                if (_myProtein.PdbName == "5GAQ")
	                {
	                    protoBarrel[i].RemoveRange(2, 1);
	                }
	                if (_myProtein.PdbName.Contains("3J9C"))
	                {
	                    protoBarrel[i][0].RemoveRange(0, 8);
	                }
                
	                /*if (_myProtein.PdbName == "7AHL" && (_myProtein.Chains[i].ChainName == "A" || _myProtein.Chains[i].ChainName == "C" || _myProtein.Chains[i].ChainName == "F")) //This stretches into the upper head
	                {
	                    protoBarrel[i][0].RemoveRange(0, 9);
	                }*/
	            } //End of checking each chain
	        }

	        public void createAllStrands(ref PolyProtein _myProtein, string outputDirectory) //This does NOT create the list of Strands
	        {
	            bool strandStart = false;
	            int strandStartRes = 0;
	            int strandEndRes = 0;
	            int CurrentResNum;
	            int strandNum = 0;
	            int betaScore = 0;
	            List<List<int>> StrandsInChain = new List<List<int>>();
	            List<int> myStrand = new List<int>();

	            foreach(Chain chain in _myProtein.Chains)
	            {
	                for (CurrentResNum = 0; CurrentResNum < chain.ResidueCount; CurrentResNum++)
	                {
	                    if (strandStart == false)
	                    {
	                        if (CurrentResNum >= chain.ResidueCount - 5) continue;

	                        if (chain.Residues[CurrentResNum].SSType == "B")
	                        {
	                            strandStartRes = chain.Residues[CurrentResNum].ResNum; //ResNum is index of residue w/in protein; Seq_ID is the PDB residue number
	                            betaScore = 0;
	                            for (int i = 1; i < 5; i++)
	                            {
	                                if (chain.Residues[CurrentResNum + i].SSType == "B") betaScore++;
	                            }
	                            if (betaScore >= 3) strandStart = true;
	                        }
	                    }

	                    else //Strand has been started
	                    {
	                        if ((CurrentResNum == chain.ResidueCount - 1) || (chain.Residues[CurrentResNum + 1].SeqID > (chain.Residues[CurrentResNum].SeqID + 6))) //If we are still forming a chain and reach the end of the residues OR if big gap in seqIDs, end the chain
	                        {
	                            goto EndChain;
	                        }

	                        else
	                        {
	                            if (chain.Residues[CurrentResNum].SSType == "B") continue;
	                            else //if the chain is started and next residue is NOT beta-residue, determine if @ the end of a chain
	                            {
	                                //if there's a big bend in chain, definitely end the chain
	                                if (Vector3D.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction) >= 80)
	                                {
	                                    goto EndChain;
	                                }
	                                //if no bend in chain, check next residue in line 
	                                else
	                                {
	                                    if (chain.Residues[CurrentResNum + 1].SSType == "B" && chain.Residues[CurrentResNum + 2].SSType == "B") continue; // XBB
	                                    else if (chain.Residues[CurrentResNum + 1].SSType == "B" && chain.Residues[CurrentResNum + 2].SSType != "B") //XBX
	                                    {
	                                        double Angle1 = Vector3D.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
	                                        double Angle2 = Vector3D.AngleBetween(chain.Residues[CurrentResNum + 1].Direction, chain.Residues[CurrentResNum + 2].Direction);
	                                        if (Angle1 > 71 || Angle2 > 71 || Angle1 < 45 || Angle2 < 45) goto EndChain; //B-type chain of 7ahl has angle at 71.5
	                                        else continue;
	                                    }
	                                    else //XX
	                                    {
	                                        betaScore = 0;
	                                        for (int i = 2; i < 5; i++)
	                                        {
	                                            try
	                                            {
	                                                if (chain.Residues[CurrentResNum + i].SSType == "B") betaScore++;
	                                            }
	                                            catch (ArgumentOutOfRangeException)
	                                            {
	                                                betaScore = 0;
	                                            }

	                                        }
	                                        if (betaScore >= 3)
	                                        {
	                                            double Angle1 = Vector3D.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction);
	                                            double Angle2 = Vector3D.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
	                                            if (Angle1 + Angle2 > 129 || Angle1 > 80 || Angle2 > 80) //if there's a decent bend at the two non-beta residues between long beta-conf sequences //changed from 135 for 7AHL and 3W9T on 12/3 MWF
	                                            {
	                                                goto EndChain;
	                                            }
	                                        }
	                                        else //this means two res in a row were non-beta AND that less than three res of next 3 were also not beta 
	                                        {
	                                            goto EndChain;
	                                        }

	                                    }
	                                }
	                            } //end of dealing with non-Beta conformation residues
	                        }
	                    } //end of strandStart == true
	                    continue;
	                EndChain:
	                    {
	                        strandEndRes = CurrentResNum - 1;
	                        if (strandEndRes - strandStartRes > 5)
	                        {
	                            for (int j = strandStartRes; j <= strandEndRes; j++)
	                            {
	                                myStrand.Add(chain.Residues[j].SeqID);
	                                chain.Residues[j].SSType = "B";
	                            }
	                            List<int> newList = new List<int>();
	                            newList.AddRange(myStrand);
	                            StrandsInChain.Add(newList);
	                            myStrand.Clear();
	                            strandNum++;
	                        }
	                        strandStart = false;
	                    }

	                }
	                if (_myProtein.PdbName.ToUpper() == "1WP1")
	                {
	                    //this strand has slightly too many non-beta residues in middle to be recognized
	                    strandStartRes = chain.Residues.FindIndex(a => a.SeqID == 88);
	                    strandEndRes = chain.Residues.FindIndex(a => a.SeqID == 99);
	                    for (int j = strandStartRes; j <= strandEndRes; j++)
	                    {
	                        myStrand.Add(chain.Residues[j].SeqID);
	                        chain.Residues[j].SSType = "B";
	                    }
	                    List<int> newList = new List<int>();
	                    newList.AddRange(myStrand);
	                    StrandsInChain.Insert(0, newList);
	                    myStrand.Clear();
	                }
	                List<List<int>> tempList = new List<List<int>>();
	                tempList.AddRange(StrandsInChain);
	                protoBarrel.Add(tempList);
	                StrandsInChain.Clear();
	                strandNum = 0;
	                if (_myProtein.PdbName == "3W9T") //only four strands are recognized in 3W9T, which means it doesn't go through MakeBarrelCircular
	                {
	                    protoBarrel[chain.ChainNum].RemoveRange(2, 1);
	                    //protoBarrel[chain.ChainNum].RemoveRange(0, 1);
	                }
                
	            }
	        }//end of create strands fxn

	        public void extractStrands(ref PolyProtein _myProtein, string outputDirectory) //This doesn't quite work, but it is very close
	        {

	            //what i consider a beta strand based on looking at it...
	            // this code defines a beta strand start as being a residue that is in a Beta cofiguration based on its Phi Psi
	            // and the next 3 out of 4 amino acids are in beta conformation as well.  The strand has a possible end when 
	            // there exists a break of 6 or more residues in the code or 
	            // there exists a non-beta residue and either 4 of the next 5 are not beta or the next 2 in a row are not beta
	            // after that, the end is confirmed by there not being a new start for at least another 11 residues or if the
	            // next start has a direction of greater than 120deg with respect the previous end (or previous end -1)
	            Strands = new List<Strand>();

	            string fileLocation1 = outputDirectory + "color_strand_";
	            string fileLocation3 = ".py";
	            string fileLocation = fileLocation1 + _myProtein.Chains[0].PdbName + fileLocation3;

	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
	            {

	                bool strandStart = false;
	                int strandStartRes = 0;
	                //int strandEndRes = 0;
	                int strandCtr = 0;
	                //bool inLoop = false;
	                for (int chain = 0; chain < _myProtein.Chains.Count; chain++)
	                {
	                    for (int resCtr = 0; resCtr < _myProtein.Chains[chain].Residues.Count; resCtr++)
	                    {
	                        if (strandStart == false)  // if false, find the start
	                        {
	                            if (resCtr + 5 < _myProtein.Chains[chain].Residues.Count)
	                            {
	                                int betaScore = 0;
	                                if (_myProtein.Chains[chain].Residues[resCtr].SSType == "B")
	                                {
	                                    for (int subResCtr = resCtr + 1; subResCtr < resCtr + 5; subResCtr++)
	                                    {
	                                        if (_myProtein.Chains[chain].Residues[subResCtr].SSType == "B")
	                                        {
	                                            betaScore++;
	                                        }
	                                    }
	                                }
	                                if (betaScore > 2) // if at least 3 of the next 4 are in B conformation
	                                {
	                                    strandStart = true;
	                                    strandStartRes = resCtr;
	                                }
	                            }
	                        }
	                        if (strandStart == true) // if true, find the end
	                        {
	                            if (_myProtein.Chains[chain].Residues[resCtr + 1].SeqID - _myProtein.Chains[chain].Residues[resCtr].SeqID > 6)
	                            { // if there is a large break you have found the end of the chain
	                                strandCtr++;
	                                Strand myStrand = new Strand(_myProtein.Chains[chain], strandStartRes, resCtr, strandCtr);
	                                this.Strands.Add(myStrand);
                                
	                                //resCtr++;
	                                strandStart = false;

	                            }
	                            else if (resCtr + 6 < _myProtein.Chains[chain].Residues.Count) // don't want to fall off the end
	                            {
	                                int notBetaScore = 0;
	                                for (int subResCtr = resCtr + 1; subResCtr < resCtr + 6; subResCtr++)
	                                {
	                                    if (_myProtein.Chains[chain].Residues[subResCtr].SSType != "B")
	                                    {
	                                        notBetaScore++;
	                                    }
	                                }
	                                bool notBetaScore2 = false;
	                                if (_myProtein.Chains[chain].Residues[resCtr + 1].SSType != "B" && _myProtein.Chains[chain].Residues[resCtr + 2].SSType != "B") notBetaScore2 = true;

	                                if (_myProtein.Chains[chain].Residues[resCtr].SSType == "B" && (notBetaScore > 3 || notBetaScore2 == true)) //if 2 in a row or if 4 of the next 6 are not beta (tried 4 of 5 and 3 of 4...)
	                                {//Found a possible end of the chain

	                                    for (int subResCtr = resCtr; subResCtr < _myProtein.Chains[chain].Residues.Count - 5; subResCtr++)
	                                    { //keep going until you find a possible beginning

	                                        //find next beginning of strand
	                                        int betaScore = 0;
	                                        if (_myProtein.Chains[chain].Residues[subResCtr].SSType == "B")
	                                        {
	                                            for (int subResCtr2 = subResCtr + 1; subResCtr2 < subResCtr + 5; subResCtr2++)
	                                            {

	                                                if (_myProtein.Chains[chain].Residues[subResCtr2].SSType == "B")
	                                                {
	                                                    betaScore++;
	                                                }
	                                            }
	                                        }

	                                        if (betaScore < 3)
	                                        {
	                                            continue;
	                                        }
	                                        else // if at least 3 of the next 4 are in B conformation
	                                        { //compare strandEndRes to new strand
	                                            //subresiduectr is the possible new beginning
	                                            //resCtr is the possible end

	                                            double angleBetween = Vector3D.AngleBetween(_myProtein.Chains[chain].Residues[subResCtr].Direction, _myProtein.Chains[chain].Residues[resCtr].Direction);
	                                            double angleBetween2 = Vector3D.AngleBetween(_myProtein.Chains[chain].Residues[subResCtr].Direction, _myProtein.Chains[chain].Residues[resCtr - 1].Direction);
	                                            double angleBetween3 = Vector3D.AngleBetween(_myProtein.Chains[chain].Residues[subResCtr].Direction, _myProtein.Chains[chain].Residues[resCtr - 2].Direction);

	                                            file.WriteLine("Vector between res {0} and res {1} is {2} the previous angle is {3} the one before that is {4} ", resCtr, subResCtr, angleBetween, angleBetween2, angleBetween3);

	                                            if (angleBetween < 120 && angleBetween2 < 120 && angleBetween3 < 120 && subResCtr - resCtr < 8)//was 11
	                                            {
	                                                resCtr = subResCtr;
	                                                // inLoop = false;

	                                                break;
	                                            }
	                                            else
	                                            {
	                                                if (resCtr - strandStartRes > 3)
	                                                {
	                                                    strandCtr++;
	                                                    Strand myStrand = new Strand(_myProtein.Chains[chain], strandStartRes, resCtr, strandCtr);
	                                                    this.Strands.Add(myStrand);

	                                                }
	                                                strandStartRes = subResCtr;

	                                                resCtr = subResCtr;
	                                            }

	                                        }
	                                    }

	                                }
	                            }
	                            else
	                            {
	                                for (int subResCtr = _myProtein.Chains[chain].Residues.Count - 1; subResCtr > resCtr; subResCtr--)
	                                {

	                                    if (_myProtein.Chains[chain].Residues[subResCtr].SSType == "B")
	                                    {
	                                        strandCtr++;
	                                        Strand myStrand = new Strand(_myProtein.Chains[chain], strandStartRes, subResCtr, strandCtr);
	                                        this.Strands.Add(myStrand);


	                                        return;
	                                    }
	                                }
	                                strandCtr++;
	                                Strand endStrand = new Strand(_myProtein.Chains[chain], strandStartRes, _myProtein.Chains[chain].Residues.Count - 1, strandCtr);
	                                this.Strands.Add(endStrand);
                                
	                                return;
	                            }

	                        }
	                    } //end of checking each residue
	                }
	            }

	        }

	        public void checkStrandDefnsDSSP(List<Strand> strandlist, ref PolyProtein _myProtein) //Added 6-5-17 for loops
	        {
	            foreach (Strand strand1 in strandlist)
	            {
	                int chainID = strand1.ChainNum;
                 
	                try
	                {
	                    while (_myProtein.Chains[chainID].Residues[strand1.Residues[0].ResNum - 1].DSSP == "E")
	                    {
	                        strand1.Residues.Add(_myProtein.Chains[chainID].Residues[strand1.Residues[0].ResNum - 1]);
	                        strand1.Residues.Sort((x,y) => x.ResNum.CompareTo(y.ResNum));
	                        Console.WriteLine("added res{1} to strand {0}", strand1.StrandNum, _myProtein.Chains[chainID].Residues[strand1.Residues[0].ResNum - 1].ResNum);
	                    }
	                }
	                catch (ArgumentOutOfRangeException)
	                {
	                    continue;
	                }

	                try
	                {
	                    while (_myProtein.Chains[chainID].Residues[strand1.Residues.Last().ResNum + 1].DSSP == "E")
	                    {
	                        strand1.Residues.Add(_myProtein.Chains[chainID].Residues[strand1.Residues.Last().ResNum + 1]);
	                        strand1.Residues.Sort((x, y) => x.ResNum.CompareTo(y.ResNum));
	                        Console.WriteLine("added res{1} to strand {0}", strand1.StrandNum, _myProtein.Chains[chainID].Residues[strand1.Residues.Last().ResNum + 1].ResNum);
	                    }
	                }
	                catch (ArgumentOutOfRangeException)
	                {
	                    continue;
	                }
	            }
	        }

	        public void checkDSSPNeighs(Res Res1, ref Chain _myChain) //Added 6-5-17 for loops
	        {
	            double minD = 2.75;

	            foreach (Res Res2 in _myChain)
	            {
	                if (Math.Abs(Res1.SeqID - Res2.SeqID) > 2)
	                {
	                    double d = (Res1.Atoms[0].Hydrogen - Res2.BackboneCoords["O"]).Length;
	                    if (d < minD && (Res2.SSType == "B" || Res2.DSSP == "E"))
	                    {
	                        if (Res1.Neighbors.Contains(Res2.ResNum) == false) Res1.Neighbors.Add(Res2.ResNum);
	                        if (Res2.Neighbors.Contains(Res1.ResNum) == false) Res2.Neighbors.Add(Res1.ResNum);
	                    }
	                }
	            }
	        }

	        public void getEllipseCoords(List<string> outPDBs)
	        {
	            //this is the top (extracellular) ellipse
	            List<Vector3D> myCEllipse = new List<Vector3D>();
	            //This is the bottom ellipse
	            List<Vector3D> myNEllipse = new List<Vector3D>();
	            List<Vector3D> sortedEllipseCoords = new List<Vector3D>();
	            Vector3D Ncentroid = new Vector3D();
	            Vector3D Ccentroid = new Vector3D();
	            int firstCAIndex; int lastCAIndex;

            
	            foreach (Strand strand in this.Strands)
	            {
	                Vector3D firstCA = new Vector3D();
	                firstCA = strand.Residues[0].BackboneCoords["CA"];
	                firstCAIndex = strand.Residues[0].Atoms.FindIndex(a => a.AtomName == "CA");
	                Vector3D lastCA = new Vector3D();
	                lastCA = strand.Residues.Last().BackboneCoords["CA"];
	                lastCAIndex = strand.Residues.Last().Atoms.FindIndex(a => a.AtomName == "CA");

	                if (outPDBs.Contains(this.PdbName.ToUpper())) //These have numbered strands that run in the opposite direction
	                {
	                    if (strand.StrandNum %2 ==0)
	                    {
	                        myNEllipse.Add(firstCA);
	                        Ncentroid += firstCA;
	                        myCEllipse.Add(lastCA);
	                        Ccentroid += lastCA;
	                        strand.CEllipseCoords = (Point3D) strand.Residues.Last().Atoms[lastCAIndex].Coords;
	                    }
	                    else
	                    {
	                        myCEllipse.Add(firstCA);
	                        myNEllipse.Add(lastCA);
	                        Ncentroid += lastCA;
	                        Ccentroid += firstCA;
	                        strand.CEllipseCoords = (Point3D) strand.Residues[0].Atoms[firstCAIndex].Coords;
	                    }

	                }

	                else
	                { //odd-numbered strands face towards extracellular (C) side; even-numbered strands face downwards towards N side
	                    if (strand.StrandNum %2 ==0)
	                    {
	                        myCEllipse.Add(firstCA);
	                        myNEllipse.Add(lastCA);
	                        Ncentroid += lastCA;
	                        Ccentroid += firstCA;
	                        strand.CEllipseCoords = (Point3D) strand.Residues[0].Atoms[firstCAIndex].Coords;
	                    }
	                    else
	                    {
	                        myNEllipse.Add(firstCA);
	                        myCEllipse.Add(lastCA);
	                        Ccentroid += lastCA;
	                        Ncentroid += firstCA;
	                        strand.CEllipseCoords = (Point3D) strand.Residues.Last().Atoms[lastCAIndex].Coords;
	                    }
	                }
	            }
            
	            this.CellipseCoords = myCEllipse;
	            this.NellipseCoords = myNEllipse;
            
	            Ncentroid = Ncentroid / this.Strands.Count;
	            Ccentroid = Ccentroid / this.Strands.Count;
	            this.Ncentroid = Ncentroid;
	            this.Ccentroid = Ccentroid;
	        }

	        public void rotateToZ(List<String> outPDBs, PolyProtein _myProtein)
	        {
	            double length = this.Axis.Length;

	            double psi = Math.Atan(this.Axis.Y / this.Axis.X);
	            Matrix3D rotationMatrixZT = new Matrix3D(Math.Cos(psi), Math.Sin(psi), 0, 0, -1 * Math.Sin(psi), Math.Cos(psi), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
	            Matrix3D rotationMatrixZ = new Matrix3D(Math.Cos(psi), -1 * Math.Sin(psi), 0, 0, Math.Sin(psi), Math.Cos(psi), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
	            Vector3D axis1 = new Vector3D();
	            axis1 = rotationMatrixZ.Transform(this.Axis);

	            double theta = Math.Atan(axis1.X / axis1.Z);
	            Matrix3D rotationMatrixYT = new Matrix3D(Math.Cos(theta), 0, -1 * Math.Sin(theta), 0, 0, 1, 0, 0, Math.Sin(theta), 0, Math.Cos(theta), 0, 0, 0, 0, 0);
	            Matrix3D rotationMatrixY = new Matrix3D(Math.Cos(theta), 0, Math.Sin(theta), 0, 0, 1, 0, 0, -1 * Math.Sin(theta), 0, Math.Cos(theta), 0, 0, 0, 0, 0);

	            Vector3D axis2 = new Vector3D();
	            axis2 = rotationMatrixY.Transform(axis1);

	                for (int chain = 0; chain < _myProtein.Chains.Count; chain++)
	                {
	                    foreach (Res res in _myProtein.Chains[chain])
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            atom.Coords = rotationMatrixZ.Transform(atom.Coords);
	                            atom.Coords = rotationMatrixY.Transform(atom.Coords);
	                            if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                            {
	                                res.BackboneCoords[atom.AtomName] = atom.Coords;
	                            }
	                        }
	                    }
	                }
                    

	            getEllipseCoords(outPDBs);

	            if (this.Ccentroid.Z < this.Ncentroid.Z) //The getEllipseCoords function handles directionality of strands already, so C ellipse is always the periplasmic edge, N ellipse is the cellular side.
	            {
	                double phi = Math.PI;
	                Matrix3D rotationMatrixY2 = new Matrix3D(Math.Cos(phi), 0, Math.Sin(phi), 0, 0, 1, 0, 0, -1 * Math.Sin(phi), 0, Math.Cos(phi), 0, 0, 0, 0, 0);
	                Vector3D axis3 = new Vector3D();
	                axis3 = rotationMatrixY2.Transform(axis2);

	                for (int chain = 0; chain < _myProtein.Chains.Count; chain++)
	                {
	                    foreach (Res res in _myProtein.Chains[chain])
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            atom.Coords = rotationMatrixY2.Transform(atom.Coords);
	                            if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                            {
	                                res.BackboneCoords[atom.AtomName] = atom.Coords;
	                            }
	                        }
	                    }
	                }

	                /*foreach (Strand strand in this.Strands)
	                {
	                    foreach (Res res in strand)
	                    {
	                        foreach (Atom atom in res)
	                        {
	                            atom.Coords = rotationMatrixY2.Transform(atom.Coords);
	                            if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                            {
	                                res.BackboneCoords[atom.AtomName] = atom.Coords;
	                            }

	                        }
	                    }
	                }*/
	            }

	            Vector3D translate = new Vector3D(this.Ccentroid.X, this.Ccentroid.Y, 0);

	            for (int chain = 0; chain < _myProtein.Chains.Count; chain++)
	            {
	                foreach (Res res in _myProtein.Chains[chain])
	                {
	                    foreach (Atom atom in res)
	                    {
	                        atom.Coords -= translate;
	                        if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                        {
	                            res.BackboneCoords[atom.AtomName] = atom.Coords;
	                        }
	                    }
	                }
	            }

	            /*foreach (Strand strand in this.Strands)
	            {
	                foreach (Res res in strand)
	                {
	                    foreach (Atom atom in res)
	                    {
	                        atom.Coords -= translate;
	                        if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                        {
	                            res.BackboneCoords[atom.AtomName] = atom.Coords;
	                        }
	                    }
	                }
	            }*/

	            getEllipseCoords(outPDBs);
	        }

	        public void centerZ(List<string> out_ins, List<string> out_pdbs, PolyProtein _myProtein)
	        {
	            double factor;
	            double avgZ = (this.Ncentroid.Z + this.Ccentroid.Z) / 2;
	            getEllipseCoords(out_pdbs);
	            if (out_ins.Contains(this.PdbName.ToUpper()) == true) factor = this.Ncentroid.Z + 12.5;
	            else factor = this.Ccentroid.Z - 12.5;

	            for (int chain = 0; chain < _myProtein.Chains.Count; chain++)
	            {
	                foreach (Res res in _myProtein.Chains[chain])
	                {
	                    foreach (Atom atom in res)
	                    {
	                        Vector3D newCoords = new Vector3D(atom.Coords.X, atom.Coords.Y, atom.Coords.Z - factor);
	                        atom.Coords = newCoords;
	                        if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                        {
	                            res.BackboneCoords[atom.AtomName] = atom.Coords;
	                        }
	                    }
	                }
	            }

	            /*foreach (Strand strand in this.Strands)
	            {
	                foreach (Res res in strand)
	                {
	                    foreach (Atom atom in res)
	                    {
	                        //Vector3D newCoords = new Vector3D(atom.Coords.X, atom.Coords.Y, atom.Coords.Z - (avgZ));
	                        Vector3D newCoords = new Vector3D(atom.Coords.X, atom.Coords.Y, atom.Coords.Z - factor);
	                        atom.Coords = newCoords;
	                        if (atom.AtomName == "CA" || atom.AtomName == "C" || atom.AtomName == "N" || atom.AtomName == "O")
	                        {
	                            res.BackboneCoords[atom.AtomName] = atom.Coords;
	                        }
	                    }
	                }
	            }*/
	        }

	        public void sortStrandsAroundEllipse()
	        {
	            List<string> DontSort = new List<string> { "4MT4", "4MT0", "1WP1"}; //1uun doesn't sort properly at chain E
            
	            if (this.PdbName.ToUpper() == "1EK9") //This pdb needs to use the middle of the barrel for sorting
	            {
	                Vector3D axisLine = (this.Ccentroid - this.Ncentroid) / 2;
	                Point3D midPoint = (Point3D)axisLine;
	                Point3D midRes = (Point3D)this.Strands[0].Residues[(this.Strands[0].Residues.Count / 2)].BackboneCoords["CA"];
	                Vector3D upaxis = midRes - midPoint;
	                foreach (Strand strand in this.Strands)
	                {
	                    Point3D midCstrand = (Point3D)strand.Residues[(strand.Residues.Count / 2)].BackboneCoords["CA"];
	                    strand.angle = Vector3D.AngleBetween(midCstrand - midPoint, upaxis);
	                    Vector3D crossP = Vector3D.CrossProduct(upaxis, midCstrand - midPoint);
	                    //Console.WriteLine("{0}{1}\t{2}\t{3}", strand.ChainName, strand.StrandNum, strand.angle, crossP);

	                    if (crossP.Z > 0) strand.angle = 360 - strand.angle; //positive Z means more than halfway around circle from origin
	                }
	            }

	            else
	            {
	                //Console.WriteLine(upaxis);
	                Vector3D upaxis = this.Strands[0].CEllipseCoords - (Point3D)this.Ccentroid;
	                foreach (Strand strand in this.Strands)
	                {
	                    Point3D strandCE = strand.CEllipseCoords;
	                    strand.angle = Vector3D.AngleBetween(strandCE - (Point3D)this.Ccentroid, upaxis);
	                    Vector3D crossP = Vector3D.CrossProduct(upaxis, strandCE - (Point3D)this.Ccentroid);
	                    if (crossP.Z > 0) strand.angle = 360 - strand.angle; //positive Z means more than halfway around circle from origin
	                }
	            }
	            if (DontSort.Contains(this.PdbName.ToUpper()) == false)
	            {
	                this.Strands.Sort((a, b) => a.angle.CompareTo(b.angle));
	            }
	            /*foreach (Strand strand in this.Strands)
	            {
	                Console.WriteLine("{0}{1}", strand.ChainName, strand.StrandNum);
	            }*/
            
	        }

            public IEnumerator<Strand> GetEnumerator()
            {
                return Strands.GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return this.GetEnumerator();
            }
        }//end of PBarrel class
    }
}

  

        

  