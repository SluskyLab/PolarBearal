/*
**  File: Mono.cs
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
using System.Linq;
using System.Windows.Media.Media3D;
using betaBarrelProgram.BarrelStructures;
using System.Collections;


namespace betaBarrelProgram
{
    namespace Mono
    {
        public class MonoProtein : Protein
        {

            public List<Chain> Chains = new List<Chain>();
            public int ChainCount { get; set; }
            public int totalResNum { get; set; }
            public string PdbName { get; set; }


            public MonoProtein(ref AtomParser.AtomCategory _myAtomCat, string PdbName)
            {
                this.Chains = new List<Chain>();
                this.ChainCount = 0;
                this.totalResNum = 0;
                this.PdbName = PdbName;

                for (int chainNum = 0; chainNum < _myAtomCat.chainAtomsList.Count; chainNum++)
                {
                    bool IsItProtein = false;
                    for (int atomNum = 0; atomNum < _myAtomCat.ChainAtomList[chainNum].cartnAtoms.Count; atomNum++)
                    {
                        if (_myAtomCat.ChainAtomList[chainNum].CartnAtoms[atomNum].atomName == "CA") IsItProtein = true;
                    }
                    if (IsItProtein == true)
                    {
                        Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, true, Global.DB_DIR);
                        this.Chains.Add(myChain);
                        this.ChainCount++;
                    }
                }
            }

            public MonoProtein(ref AtomParser.AtomCategory _myAtomCat, int chainNum, string PdbName)
            {
                this.Chains = new List<Chain>();
                Chain myChain = new Chain(ref _myAtomCat, chainNum, PdbName, true, Global.DB_DIR);
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

        public class MonoBarrel : Barrel

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
            public List<List<int>> protoBarrel { get; set; }
            public double AvgRadius { get; set; }
            public List<Res> LoopResies { get; set; }
            public string PdbName { get; set; }
            public string ChainName { get; set; }
            public List<double> StrandLength { get; set; }
            public Vector3D OriginalNcentroid { get; set; }
            public Vector3D OriginalCcentroid { get; set; }
            public Vector3D AxisVector { get; set; }
            public Vector3D NewCaxisPt { get; set; }
            public Vector3D OldCaxisPt { get; set; }
            public int ShearNum { get; set; }
            public List<double> PrevTwists { get; set; }

            public static string path = Global.MONO_OUTPUT_DIR;

            //barrel constructor 
            public MonoBarrel(Chain _myChain, Protein _myProtein)
            {
                this.Strands = new List<Strand>();
                this.ChainName = _myChain.ChainName;
                this.protoBarrel = new List<List<int>>();
                this.PdbName = _myChain.PdbName;
	            
				//Writing info about residues
				/*string fileLocation2 = path + "PhiPsiBBAngle/" + this.PdbName + ".txt";
	                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
	                {
	                        foreach (Res res in _myChain)
	                        {
	                            if (res.ResNum != 0)
	                            {
	                                file.WriteLine("Residue {0}({1}), SSType {2}, phi {4}, psi {5}, Angle formed at this residue {3}", res.SeqID, res.ThreeLetCode, res.SSType, (Vector3D.AngleBetween(_myChain.Residues[res.ResNum - 1].Direction, res.Direction)), res.Phi, res.Psi);
	                            }
	                        }
	                }*/

	            #region makeBarrel
	            //Pattern finding method of barrels
	            /*createStrandsPattern(ref _myChain, ref _myProtein, path);
	            //for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }
	            checkStrandPatt(ref _myChain);
	            makeBarrelCircular(ref _myChain);
	            for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }
	            removeNonBarrelRes(ref _myChain);

	            makeBarrelCircular(ref _myChain);
	            removeNonBarrelRes(ref _myChain);*/
	            //Reset to run diff way
	            /*foreach (Res res1 in _myChain)
	            {
	                res1.Neighbors = new List<int>();
	            }
	            this.protoBarrel = new List<List<int>>();

	            Console.WriteLine("Redefining strands");*/
				//Current method of defining strands
                createStrands(ref _myChain);
				//If DSSP definitions are important
	            /*foreach (Res Res1 in _myChain){
	                checkDSSPNeighs(Res1, ref _myChain);
				}*/
				
                makeBarrelCircular(ref _myChain); // keeps out the beta strands that occlude the barrel
                removeNonBarrelRes(ref _myChain);

                //twice for 2VQI
                makeBarrelCircular(ref _myChain);
                removeNonBarrelRes(ref _myChain);

                //check the ends for inserting loops
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                this.Axis = this.Ccentroid - this.Ncentroid;

                for (int strandCtr = 0; strandCtr < this.protoBarrel.Count; strandCtr++)
                {// gets rid of loops that can bond with each other but go into the middle
                    double angle = (Vector3D.AngleBetween(_myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"] - _myChain.Residues[this.protoBarrel[strandCtr][1]].BackboneCoords["CA"], this.Axis));
                    double direction = Vector3D.AngleBetween(_myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"] - _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"], this.Axis);

                    while ((direction > 90 && angle < 105) || (angle > 75 && direction < 90))
                    {
                        this.protoBarrel[strandCtr].RemoveAt(0);
                        angle = (Vector3D.AngleBetween(_myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"] - _myChain.Residues[this.protoBarrel[strandCtr][1]].BackboneCoords["CA"], this.Axis));
                    }

                    double angle2 = (Vector3D.AngleBetween(_myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"] - _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 2]].BackboneCoords["CA"], this.Axis));

                    while ((angle2 < 105 && direction < 90) || (angle2 > 75 && direction > 90))
                    {
                        this.protoBarrel[strandCtr].RemoveAt(this.protoBarrel[strandCtr].Count - 1);
                        angle2 = (Vector3D.AngleBetween(_myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"] - _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 2]].BackboneCoords["CA"], this.Axis));
                    }
                }

           		//checkStrands is written to use protobarrel list before the strand list is created. 
				//checkStrandDefnsDSSP(ref _myChain);
                
				//redefine axis
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);
                this.Axis = this.Ncentroid - this.Ccentroid;


                for (int strandCtr = 0; strandCtr < protoBarrel.Count; strandCtr++)
                {
                    Strand newStrand = new Strand(_myChain, protoBarrel[strandCtr][0], protoBarrel[strandCtr][protoBarrel[strandCtr].Count - 1], strandCtr);
                    this.Strands.Add(newStrand);
                }

                #endregion

                #region Rotate barrel to z-axis
                this.AvgRadius = SharedFunctions.setRadius(this.Strands, this.Axis, this.Ccentroid, this.Ncentroid);
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                this.AxisVector = SharedFunctions.getNormal(this.NellipseCoords, this.Ncentroid);

                this.OriginalCcentroid = this.Ccentroid;
                this.OriginalNcentroid = this.Ncentroid;
                this.NewCaxisPt = ((this.Ncentroid - this.Ccentroid).Length * this.AxisVector) + this.Ncentroid;
                int axisDirection = 1;
                if ((this.NewCaxisPt - this.OriginalCcentroid).Length > ((((this.Ncentroid - this.Ccentroid).Length * -1 * this.AxisVector) + this.Ncentroid) - this.OriginalCcentroid).Length)
                {
                    axisDirection = -1;
                    this.NewCaxisPt = (((this.Ncentroid - this.Ccentroid).Length * -1 * this.AxisVector) + this.Ncentroid);
                }

                this.Axis = this.NewCaxisPt - this.Ncentroid;
                this.OldCaxisPt = this.NewCaxisPt;
                this.rotateToZ(ref _myChain);

                if (this.Axis.Length > 22.5) setBottom(ref _myChain);
                else centerZ(ref _myChain);

                rotate180(ref _myChain);

                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                this.AxisVector = SharedFunctions.getNormal(this.NellipseCoords, this.Ncentroid);
                this.NewCaxisPt = (axisDirection * (this.Ncentroid - this.Ccentroid).Length * this.AxisVector) + this.Ncentroid;
                this.Axis = this.NewCaxisPt - this.Ncentroid;

                for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                {
                    for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count; resCtr++)
                    {
                        this.Strands[strandCtr].Residues[resCtr].Z = this.Strands[strandCtr].Residues[resCtr].BackboneCoords["CA"].Z;
                    }
                }
                #endregion

                SharedFunctions.setInOut(this.Strands, path, this.PdbName, this.Axis, this.Ccentroid, this.Ncentroid);

                /* ---------output information--------- */
	            this.StrandLength = SharedFunctions.getStrandLengths(this.Strands, path, this.PdbName);
                this.PrevTwists = SharedFunctions.writeTwists(this.Strands, Global.MONO_OUTPUT_DIR, this.PdbName);

	            /*string fileLocation6 = path + "\\ZCoords\\XYZCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation6))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Chain" + "\t" + "XYZcoords";
	                file.WriteLine(newLine);
	                foreach (Res res in _myChain)
	                {
	                    foreach (Atom atom in res)
	                    {
	                        if (atom.AtomName == "CA") file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", res.ThreeLetCode, res.SeqID, res.ChainName, atom.Coords.X, atom.Coords.Y, atom.Coords.Z);
	                    }
	                }
	            }*/

	            string fileLocation15 = path + "\\ZCoords\\AllZCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation15))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Strand" + "\t" + "Z-coord";
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
	            }

	            //DSSP has to be used to define strands before the strands are created; although I could clear and recreate them if necessary...
	            //Dictionary<string, string> Loops = SharedFunctions.getLoopTurns(this.Strands, ref _myChain, path, this.PdbName);
                //SharedFunctions.writePymolScriptForLoops(Loops, Global.MONO_OUTPUT_DIR, Global.MACMONODBDIR, ref _myChain, this.PdbName);
                //SharedFunctions.findLoopsHBondingPartnersGeomOnly(Loops, Global.MONO_OUTPUT_DIR, ref _myChain, this.PdbName, false);

                SharedFunctions.writePymolScriptForStrands(this.Strands, Global.MONO_OUTPUT_DIR, Global.MACMONODBDIR, this.PdbName);
	            //writeAminoAcidsTypesToFile(ref _myChain, path);

	            SharedFunctions.setInOut(this.Strands, path, this.PdbName, this.Axis, this.Ccentroid, this.Ncentroid);

	            /*string fileLocation6 = path + "\\Renumb\\XYZCoords_" + this.PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation6))
	            {
	                string newLine = "Res" + "\t" + "Num" + "\t" + "Chain" + "\t" + "XYZcoords";
	                file.WriteLine(newLine);
	                foreach (Res res in _myChain)
	                {
	                    foreach (Atom atom in res)
	                    {
	                        if (atom.AtomName == "CA") file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", res.ThreeLetCode, res.ResNum + 1, res.ChainName, atom.Coords.X, atom.Coords.Y, atom.Coords.Z);
	                    }
	                }
	            }*/

	            // make list of loop residues
	            //listLoopRes(ref _myChain);

	            //find shear number July 22, 2014
	            //int ShearNum = shearNum(ref _myChain);

	            //SharedFunctions.findNearestNeighbors(this.Strands, path, this.PdbName); //based on CA distances; this is needed for shear number determination now.
	            //getShearNum();
            
	            //SharedFunctions.findHBondingPartnersEnergy(this.Strands, path, this.PdbName); //SQRWL-like implementation of hydrogen bonding
	            //SharedFunctions.findNearNeighDistOnly(this.Strands, path, this.PdbName);
           
	            //Energy adds neighbors to res; geom and dist both add to atom; this section clears both sets of lists to be reformed by a different function.
	            /*foreach (Strand strand in this.Strands)
	            {
	                foreach (Res res in strand)
	                {
	                    res.SideChainNeighbors = new List<Res>();
	                    res.BackboneNeighbors = new List<Res>();
	                    foreach (Atom atom in res)
	                    {
	                        atom.SCSCNeighAtoms = new List<Atom>();
	                        atom.SCBBNeighAtoms = new List<Atom>();
	                        atom.BBNeighAtoms = new List<Atom>();
	                    }
	                }
	            }*/

	            //SharedFunctions.findHBondingPartnersGeomOnly(this.Strands, path, this.PdbName, false); //Modelling program-like implementation of hydrogen bonding
	            //SharedFunctions.findAllNearNeighDistOnly(this.Strands, path, this.PdbName, ref _myChain);
	            //this.Axis = this.Ccentroid - this.Ncentroid; //july3
	            this.AvgTilt = SharedFunctions.getTiltsByAA(this.Strands, path, this.PdbName, this.Axis, ref Global.AADict);

	            //getTyrVector(path, ref _myChain);
            
            }

            public void createStrandsPattern(ref Chain chain, ref Protein _myProtein, string outputDirectory) //This does NOT create the list of Strands
            {
                bool strandStart = false;
                int strandStartRes = 0;
                int strandEndRes = 0;
                int CurrentResNum;
                int UsualStart = 0;
                int strandNum = 0;
                int betaScore = 0;
                List<int> myStrand = new List<int>();

                //These have very large plugs with some short beta strands.
                if (PdbName.ToUpper() == "4EPA" || PdbName.ToUpper() == "4RDR" || PdbName.ToUpper() == "4AIP") UsualStart = 120;
                if (PdbName.ToUpper() == "1XKW" || PdbName.ToUpper() == "3CSL") UsualStart = 125;
                if (PdbName.ToUpper() == "2GSK" || PdbName.ToUpper() == "3V8X" || PdbName.ToUpper() == "4CU4" || PdbName.ToUpper() == "4K3C" || PdbName.ToUpper() == "2HDI") UsualStart = 130;
                if (PdbName.ToUpper() == "3QLB" || PdbName.ToUpper() == "3EFM" || PdbName.ToUpper() == "1FEP" || PdbName.ToUpper() == "1KMO" || PdbName.ToUpper() == "3FHH") UsualStart = 135;
                if (PdbName.ToUpper() == "2IAH" || PdbName.ToUpper() == "4Q35" || PdbName.ToUpper() == "4QKY") UsualStart = 200;
                if (PdbName.ToUpper() == "4C00") UsualStart = 235;
                if (PdbName.ToUpper() == "4K3B") UsualStart = 400;

                for (CurrentResNum = UsualStart; CurrentResNum < chain.ResidueCount; CurrentResNum++)
                {
                    if (strandStart == false)
                    {
                        if (CurrentResNum >= chain.ResidueCount - 3) continue;

                        if (chain.Residues[CurrentResNum].SSType.ToUpper() == "B")
                        {
                            strandStartRes = chain.Residues[CurrentResNum].ResNum; //ResNum is index of residue w/in protein; Seq_ID is the PDB residue number
                            betaScore = 0;
                            for (int i = 1; i < 5; i++)
                            {
                                try
                                {
                                    if (chain.Residues[CurrentResNum + i].SSType.ToUpper() == "B") betaScore++;
                                }

                                catch (ArgumentOutOfRangeException)
                                {
                                    continue;
                                }
                            }
                            if (betaScore >= 2) strandStart = true; //Changed to just 2 to catch short bent upper pieces in esp FpvA
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
                            if (chain.Residues[CurrentResNum].SSType.ToUpper() == "B") continue;
                            else //if the chain is started and next residue is NOT beta-residue, determine if @ the end of a chain
                            {
                                //if there's a big bend in chain, definitely end the chain
                                if (Vector3D.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction) >= 82) //Changed from 80 for 1THQ
                                {
                                    goto EndChain;
                                }
                                //if no bend in chain, check next residue in line 
                                else
                                {
                                    if (chain.Residues[CurrentResNum + 1].SSType.ToUpper() == "B" && chain.Residues[CurrentResNum + 2].SSType.ToUpper() == "B") continue; // XBB
                                    else if (chain.Residues[CurrentResNum + 1].SSType.ToUpper() == "B" && chain.Residues[CurrentResNum + 2].SSType.ToUpper() != "B") //XBX
                                    {
                                        double Angle1 = Vector3D.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
                                        double Angle2 = Vector3D.AngleBetween(chain.Residues[CurrentResNum + 1].Direction, chain.Residues[CurrentResNum + 2].Direction);
                                        if (Angle1 > 71 || Angle2 > 71 || Angle1 < 40 || Angle2 < 40) goto EndChain; //1qj8 needs angle between 40-71 instead of 45-71
                                        else continue;
                                    }
                                    else //XX
                                    {
                                        betaScore = 0;
                                        for (int i = 2; i < 5; i++)
                                        {
                                            try
                                            {
                                                if (chain.Residues[CurrentResNum + i].SSType.ToUpper() == "B") betaScore++;
                                            }
                                            catch (ArgumentOutOfRangeException)
                                            {
                                                betaScore = 0;
                                            }

                                        }
                                        if (betaScore >= 2) //Was beta score of 3/3. Changed to 2/3 for kinky monomerics
                                        {
                                            double Angle1 = Vector3D.AngleBetween(chain.Residues[CurrentResNum - 1].Direction, chain.Residues[CurrentResNum].Direction);
                                            double Angle2 = Vector3D.AngleBetween(chain.Residues[CurrentResNum].Direction, chain.Residues[CurrentResNum + 1].Direction);
                                            if (Angle1 + Angle2 > 136 || Angle1 > 80 || Angle2 > 80) //if there's a decent bend at the two non-beta residues between long beta sequences //changed from 129 for 1T16 and others on 5-11-16
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
                        if (strandEndRes - strandStartRes >= 2) //This results in a strand of 3+ residues; super short strands will be removed later if they don't combine with another strand.
                        {
                            for (int j = strandStartRes; j <= strandEndRes; j++)
                            {
                                myStrand.Add(chain.Residues[j].ResNum);
                                //chain.Residues[j].SSType = "B"; //Removed for monomeric barrels. 
                            }
                            List<int> newList = new List<int>();
                            newList.AddRange(myStrand);
                            protoBarrel.Add(newList);
                            myStrand.Clear();
                            strandNum++;
                        }
                        strandStart = false;
                    }

                }
            }//end of create strands fxn

            public void createStrands(ref Chain _myChain)
            {

                List<int> myStrand1 = new List<int>();
                List<int> myStrand2 = new List<int>();
                List<int> myStrand3 = new List<int>();
                bool strandStart = false;
                int usually0 = 0;
                if (PdbName.ToUpper() == "3KVN") usually0 = 330;
                if (PdbName.ToUpper() == "3RFZ") usually0 = 110;
                if (PdbName.ToUpper() == "4EPA") usually0 = 120;
                if (PdbName.ToUpper() == "4QKY") usually0 = 180;
                if (PdbName.ToUpper() == "4K3B") usually0 = 400;
                if (PdbName.ToUpper() == "4K3C") usually0 = 150;

                int usuallyEnd = _myChain.ResidueCount;
                if (PdbName.ToUpper() == "3RFZ") usuallyEnd = 603;

                for (int res1ctr = usually0; res1ctr < usuallyEnd; res1ctr++)
                {
                    Res Residue1 = _myChain.Residues[res1ctr];
                    for (int atomCtr1 = 0; atomCtr1 < Residue1.Atoms.Count; atomCtr1++)
                    {

                        if (Residue1.Atoms[atomCtr1].AtomName == "N" && Residue1.ThreeLetCode != "PRO")
                        {
                            Atom N = Residue1.Atoms[atomCtr1];
                            for (int res2ctr = usually0; res2ctr < usuallyEnd; res2ctr++)
                            {

                                if (Math.Abs(res2ctr - res1ctr) > 3)
                                {
                                    Res Residue2 = _myChain.Residues[res2ctr];
                                    for (int atomCtr2 = 0; atomCtr2 < Residue2.Atoms.Count; atomCtr2++)
                                    {
                                        if (Residue2.Atoms[atomCtr2].AtomName == "O")
                                        {
                                            Atom O = Residue2.Atoms[atomCtr2];
                                            double d = (O.Coords - N.Hydrogen).Length;
                                            //just using a distance cuttoff for now, not actually calculating hbond angles
                                            //if (d < 2.75)  
                                            //{
                                            //    double w = calculateHydrogenBond(N, O, d);
                                            //}
                                            // establish first two strands

                                            double minD = 2.75; //was 2.75 for 2013 bioinf paper;

                                            if (d < minD && (Residue1.SSType == "B" || Residue2.SSType == "B"))
                                            {
                                                //if (strandStart == false) firstRes = Residue1.SeqID;
                                                if (Residue1.Neighbors.Contains(Residue2.ResNum) == false) Residue1.Neighbors.Add(Residue2.ResNum);
                                                if (Residue2.Neighbors.Contains(Residue1.ResNum) == false) Residue2.Neighbors.Add(Residue1.ResNum);

                                                if (Math.Abs(Residue1.ResNum - Residue2.ResNum) < 75)
                                                {
                                                    if (myStrand1.Count > 0)
                                                    {
                                                        myStrand1.Sort();

                                                        for (int i = myStrand1[0]; i < myStrand1[myStrand1.Count - 1]; i++)
                                                        {
                                                            if (myStrand1.Contains(i) == false) myStrand1.Add(i);
                                                            myStrand1.Sort();
                                                        }
                                                    }
                                                    if (myStrand2.Count > 0)
                                                    {
                                                        myStrand2.Sort();


                                                        for (int i = myStrand2[0]; i < myStrand2[myStrand2.Count - 1]; i++)
                                                        {
                                                            if (myStrand2.Contains(i) == false) myStrand2.Add(i);
                                                            myStrand2.Sort();
                                                        }
                                                    }
                                                    if (myStrand3.Count > 0)
                                                    {
                                                        myStrand3.Sort();

                                                        for (int i = myStrand3[0]; i < myStrand3[myStrand3.Count - 1]; i++)
                                                        {
                                                            if (myStrand3.Contains(i) == false) myStrand3.Add(i);
                                                            myStrand3.Sort();
                                                        }
                                                    }
                                                    int lowerRes;
                                                    int higherRes;
                                                    if (Residue1.ResNum - Residue2.ResNum > 0)
                                                    {
                                                        lowerRes = Residue2.ResNum;
                                                        higherRes = Residue1.ResNum;
                                                    }
                                                    else
                                                    {
                                                        lowerRes = Residue1.ResNum;
                                                        higherRes = Residue2.ResNum;
                                                    }
                                                    if (strandStart == false)
                                                    {
                                                        if ((myStrand2.Contains(lowerRes) == true || myStrand2.Contains(lowerRes + 1) == true) || ((myStrand1.Count > 4 && myStrand2.Count > 4) && (lowerRes > myStrand1[myStrand1.Count - 1] && higherRes > myStrand2[myStrand2.Count - 1])))
                                                        {
                                                            strandStart = true;

                                                        }
                                                    }
                                                    if (strandStart == false)
                                                    {
                                                        //damage control
                                                        if (myStrand1.Count > 0 && (lowerRes > myStrand1[myStrand1.Count - 1] + 6 || lowerRes < myStrand1[0])) // less than 6 misses a bit of strand 1 in 3BS0
                                                        {
                                                            if (myStrand1.Count < 4) myStrand1.Clear();
                                                            else lowerRes = myStrand1[myStrand1.Count - 1];
                                                        }

                                                        if (myStrand2.Count > 0 && (higherRes > myStrand2[myStrand2.Count - 1] + 5 || higherRes < myStrand2[0] - 9)) // 8 is here because of a loop in 1E54
                                                        {
                                                            if (myStrand2.Count < 4) myStrand2.Clear();
                                                            else higherRes = myStrand2[myStrand2.Count - 1];
                                                        }

                                                        //

                                                        if (myStrand1.Contains(lowerRes) == false && myStrand2.Contains(lowerRes) == false) myStrand1.Add(lowerRes);
                                                        if (myStrand2.Contains(higherRes) == false && myStrand1.Contains(higherRes) == false) myStrand2.Add(higherRes);


                                                    }


                                                    if (strandStart == true)
                                                    {
                                                        if (lowerRes > myStrand2[myStrand2.Count - 1] && (myStrand3.Count == 0 || higherRes > myStrand3[myStrand3.Count - 1]))
                                                        {
                                                            if (myStrand3.Count > 0 && lowerRes < myStrand3[0] - 5) break; //some loop getting in the ways
                                                            List<int> newList = new List<int>();

                                                            for (int ctr = 0; ctr < myStrand1.Count; ctr++)
                                                            {
                                                                newList.Add(myStrand1[ctr]);
                                                            }

                                                            protoBarrel.Add(newList);
                                                            myStrand1.Clear();
                                                            for (int ctr = 0; ctr < myStrand2.Count; ctr++)
                                                            {
                                                                myStrand1.Add(myStrand2[ctr]);
                                                            }
                                                            myStrand2.Clear();
                                                            for (int ctr = 0; ctr < myStrand3.Count; ctr++)
                                                            {
                                                                myStrand2.Add(myStrand3[ctr]);
                                                            }
                                                            myStrand3.Clear();

                                                        }

                                                        if (lowerRes > myStrand1[myStrand1.Count - 1] || myStrand3.Count == 0 || higherRes > myStrand3[myStrand3.Count - 1])
                                                        {
                                                            if (myStrand2.Contains(lowerRes) == false && myStrand1.Contains(lowerRes) == false && myStrand3.Contains(lowerRes) == false) myStrand2.Add(lowerRes);
                                                            if (myStrand3.Contains(higherRes) == false && myStrand2.Contains(higherRes) == false && myStrand1.Contains(higherRes) == false) myStrand3.Add(higherRes);
                                                        }
                                                        else
                                                        {

                                                            if (protoBarrel.Count > 0 && lowerRes > protoBarrel[protoBarrel.Count - 1][protoBarrel[protoBarrel.Count - 1].Count - 1])
                                                            {
                                                                if (myStrand1.Contains(lowerRes) == false && myStrand2.Contains(lowerRes) == false && myStrand3.Contains(lowerRes) == false) myStrand1.Add(lowerRes);
                                                                if (myStrand2.Contains(higherRes) == false && myStrand1.Contains(higherRes) == false && myStrand3.Contains(higherRes) == false) myStrand2.Add(higherRes);
                                                            }
                                                        }
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } //end of going through each res and adding to barrel. 

                if (protoBarrel.Count > 4)

                {// go back through from the start of the final strand and make sure we got all of strand 1 and final strand
                    myStrand3.Sort();
                    for (int lastStrndResCtr = myStrand3[0]; lastStrndResCtr < usuallyEnd; lastStrndResCtr++)
                    {
                        if (_myChain.Residues[lastStrndResCtr].Neighbors.Count != 0)
                        {
                            for (int nCtr = 0; nCtr < _myChain.Residues[lastStrndResCtr].Neighbors.Count; nCtr++)
                            {
                                int neighbor = _myChain.Residues[lastStrndResCtr].Neighbors[nCtr];
                                if (neighbor < protoBarrel[1][0])
                                {
                                    if (protoBarrel[0].Contains(neighbor) == false) protoBarrel[0].Add(neighbor);
                                    if (myStrand3.Contains(lastStrndResCtr) == false) myStrand3.Add(lastStrndResCtr);
                                }
                            }
                        }
                    }
                    protoBarrel[0].Sort();


                    // now finish everything       
                    myStrand2.Sort();
                    myStrand1.Sort();
                    myStrand3.Sort();
                    protoBarrel.Add(myStrand1);
                    protoBarrel.Add(myStrand2);
                    protoBarrel.Add(myStrand3);

                }
                //for (int x = 0; x < this.protoBarrel.Count; x++) { Console.WriteLine(x + "\t" + _myChain.Residues[this.protoBarrel[x][0]].SeqID + "\t" + _myChain.Residues[this.protoBarrel[x].Last()].SeqID); }

	            if (PdbName.ToUpper() == "3RBH" || PdbName.ToUpper() == "4AFK") //Added 4afk 6-22-17 as replacement in v4 for 3rbh
	            {
	                List<int> myStrand4 = new List<int>();
	                List<int> myStrand5 = new List<int>();
	                List<int> myStrand6 = new List<int>();
	                myStrand4.AddRange(Enumerable.Range(protoBarrel[4][0]-2, 7));
	                myStrand5.AddRange(Enumerable.Range(protoBarrel[4][20], 8));
	                myStrand6.AddRange(Enumerable.Range(protoBarrel[4][31], 9));
	                protoBarrel.RemoveRange(4, 1);
	                protoBarrel.Insert(4, myStrand6);
	                protoBarrel.Insert(4, myStrand5);
	                protoBarrel.Insert(4, myStrand4);
	            }

	            if (PdbName.ToUpper() == "2YNK")
	            {
	                myStrand1 = new List<int>();
	                myStrand1.AddRange(Enumerable.Range(protoBarrel[1][0], 10));
	                protoBarrel.Insert(1, myStrand1);
	                protoBarrel[2].RemoveRange(0, 10);
	            }
	            if (PdbName.ToUpper() == "4AIP")
	            {
	                myStrand1 = new List<int>();
	                myStrand1.AddRange(Enumerable.Range(protoBarrel[20][0], 12));
	                protoBarrel.Insert(20, myStrand1);
	                protoBarrel[21].RemoveRange(0, 12);
	            }
	            if (PdbName.ToUpper() == "2O4V")
	            {
	                protoBarrel[14].RemoveRange(0, 13);
	            }
	            if (PdbName.ToUpper() == "2X55") //added 6-22-17 for v4
	            {
	                protoBarrel[9].InsertRange(0, Enumerable.Range(protoBarrel[9][0]-4, 4));
	            }
	            if (PdbName.ToUpper() == "4FT6") //added 6-22-17 for v4
	            {
	                protoBarrel[2].RemoveRange(15, 5);
	            }
	            if (PdbName.ToUpper() == "4ZGV") //added 6-22-17 for v4
	            {
	                protoBarrel.RemoveRange(0, 2);
	                protoBarrel[8].RemoveRange(25, 23);
	                protoBarrel[9].RemoveRange(0, 5);

	                myStrand1 = new List<int>();
	                myStrand1.AddRange(Enumerable.Range(protoBarrel[15][0], 17));
	                protoBarrel.Insert(15, myStrand1);
	                protoBarrel[16].RemoveRange(0, 20);
	            }
	            if (PdbName.ToUpper() == "5FP1") //added 6-22-17 for v4
	            {
	                protoBarrel[7].RemoveRange(14, 9);
	                protoBarrel[8].RemoveRange(0, 22);
	            }
	            if (PdbName.ToUpper() == "5FQ8") //added 6-22-17 for v4
	            {
	                protoBarrel[3].RemoveRange(27, 13);
	                protoBarrel[4].RemoveRange(0, 41);
	                protoBarrel[19].RemoveRange(13, 15);
	                protoBarrel.RemoveRange(20, 2);
	            }
	            if (PdbName.ToUpper() == "5FR8") //added 6-22-17 for v4
	            {
	                protoBarrel[19].RemoveRange(17, 32);
	            }
	            if (PdbName.ToUpper() == "3SY7") //added 6-22-17 for v4
	            {
	                protoBarrel[6].RemoveRange(9, 4);
	            }
	            if (PdbName.ToUpper() == "4FSP") //added 6-22-17 for v4
	            {
	                protoBarrel[6].RemoveRange(10, 4);
	            }
	            if (PdbName.ToUpper() == "5LDV") //added 6-22-17 for v4
	            {
	                protoBarrel[6].RemoveRange(11, 16);
	            }
	            if (PdbName.ToUpper() == "5T3R") //added 6-22-17 for v4
	            {
	                protoBarrel.RemoveRange(19, 1);
	            }
	            if (PdbName.ToUpper() == "2MAF") //Added 2-17-17 for creating design scaffolds
	            {
	                protoBarrel[4].RemoveRange(12, 19);
	                protoBarrel[5].RemoveRange(0, 19);
	            }
	            if (PdbName.ToUpper() == "3DZM") //Added 2-17-17 for creating design scaffolds
	            {
	                protoBarrel[7].RemoveRange(0, 23);
               
	            }
	            if (PdbName.ToUpper() == "1FEP") //Added 6-5-17 for loop definitions
	            {
	                protoBarrel[16].RemoveRange(22, 10);
	            }
	            if (PdbName.ToUpper() == "3FIP") //Added 5-25-18 - this should have been fixed YEARS ago.
	            {
	                protoBarrel[5].RemoveRange(13, protoBarrel[5].Count-13);
	                protoBarrel[7].RemoveRange(0, 28);
	                protoBarrel.RemoveRange(6, 1);
	            }
	            if (PdbName.ToUpper() == "5M9B") //Added 5-25-18 for v5DB
	            {
	                protoBarrel[16].RemoveRange(22, 16);
	            }
	            if (PdbName.ToUpper() == "5MDO") //Added 5-25-18 for v5DB
	            {
	                protoBarrel[12].InsertRange(0, Enumerable.Range(protoBarrel[12][0]-3, 3));
	                protoBarrel[13].AddRange(Enumerable.Range(protoBarrel[13].Last(), 4));
	            }
			}//End of CreateStrands

	        public void checkStrandDefnsDSSP(ref Chain _myChain) //Added 6-5-17 for loops
	        {
	            for (int strandNum = 0; strandNum < this.protoBarrel.Count; strandNum ++)
	            {
	                //Remove turn residues
	                while (_myChain.Residues[this.protoBarrel[strandNum][0]].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum][0]].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(0);
	                    Console.WriteLine("Removed res{1} from beg of strand {0}", strandNum, removed_res);
	                }

	                while (_myChain.Residues[this.protoBarrel[strandNum].Last()].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum].Last()].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(this.protoBarrel[strandNum].Count - 1);
	                    Console.WriteLine("Removed res{1} from end of strand {0}", strandNum, removed_res);
	                }


	                //Add to the beginning of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum][0] - 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum][0] - 1 > this.protoBarrel[strandNum - 1].Last() + 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum][0] - 1], ref _myChain);
	                                Console.WriteLine("added res{1} to beg of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break;  }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into prev strand"); break; }
                        
	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }

	                //For adding to the ends of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum].Last() + 1 < this.protoBarrel[strandNum + 1][0] - 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1], ref _myChain);
	                                Console.WriteLine("added res{1} to end of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break;  }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into next strand"); break; }
	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }
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

            // this will calculate the shear number of the beta barrel.  see murzin lesk and chothia 1994 http://www.mrc-lmb.cam.ac.uk/tcb/pdf/chc/97_jmb_236_1369_94.pdf
            public int shearNum(ref Chain _myChain)
            {
                string filelocation = path + "ShearNum_" + PdbName + ".txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(filelocation))
                {
                    // first determine nearest hbonding neighbor with N coming from N-term  (arbitrary)           
                    for (int strandCtr = 0; strandCtr < this.Strands.Count(); strandCtr++)
                    {
                        int strand2 = strandCtr + 1;
                        if (strandCtr == this.Strands.Count() - 1) strand2 = 0;

                        for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count(); resCtr++)
                        {

                            Res firstRes = this.Strands[strandCtr].Residues[resCtr];
                            double minD = 3.4;
                            if (firstRes.Atoms[0].AtomName == "N" && firstRes.ThreeLetCode != "PRO")
                            {
                                Atom firstAtom = firstRes.Atoms[0];

                                for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
                                {
                                    Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[3];
                                    if (secondAtom.AtomName == "O")
                                    {
                                        double d = (secondAtom.Coords - firstAtom.Coords).Length;
                                        if (d < minD)
                                        {
                                            minD = d;
                                            firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
                                            firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
                                            firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length;
                                        }

                                    }

                                }


                            }
                            // looking at o=n as well just to be safe... 
                            if (firstRes.h_bonder == 0 && firstRes.Atoms[3].AtomName == "O")
                            {
                                Atom firstAtom = firstRes.Atoms[3];

                                for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
                                {
                                    Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[0];
                                    if (secondAtom.AtomName == "N")
                                    {
                                        double d = (secondAtom.Coords - firstAtom.Coords).Length;

                                        if (d < minD)
                                        {
                                            minD = d;
                                            firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
                                            firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
                                            firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length;
                                        }
                                    }
                                }
                            }
                            if (firstRes.h_bonder == 0) file.WriteLine("strand {0} res {1} has no partner!", strandCtr, firstRes.SeqID);
                            else file.WriteLine("strand {0} res {1}s partner is {2}; k = {3}", strandCtr, firstRes.SeqID, firstRes.h_bonderID, firstRes.ResStrandNum);
                        }
                    }
                    return 1;
                }
            }

            public void getShearNum() //7-11-16 MWF; produces an integer for shear number.
            {
                int shearNum = 0;
                int k = 0; int l = 0; int j = 0; int m = 0;
                int StrandNum1; int StrandNum2; int StrandNum0 = 0;
                string filelocation = path + "ShearNum/ShearNum_" + PdbName + ".txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(filelocation))
                {
                    int res_num = 0;
                    Res first_res = this.Strands[0].Residues[res_num];
                    Res next_res = this.Strands[0].Residues[res_num];
                    Res prev_res = this.Strands[0].Residues[res_num];
                    StrandNum1 = first_res.StrandNum;
                    if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) k = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
                    else k = first_res.ResStrandNum;
                    while (StrandNum1 < this.Strands.Count - 1)
                    {
                        if (first_res.ShearNumNeigh == null)
                        {
                            file.WriteLine("End of seq, restarting");
                            prev_res = first_res;
                            m = l;

                            //If a residue with no neighbors is reached - strand sticks up and need a new residue to restart; go to opposite end of strand and work back towards this residue
                            bool dir = false;
                            if (first_res.ResStrandNum > this.Strands[StrandNum1].Count() / 2)
                            {
                                res_num = 0;
                                dir = true;
                            }
                            else res_num = this.Strands[StrandNum1].Residues.Count() - 1;

                            while (first_res.ShearNumNeigh == null)
                            {
                                if (dir == true)
                                {
                                    res_num += 1;
                                    first_res = this.Strands[StrandNum1].Residues[res_num];
                                }
                                else
                                {
                                    res_num -= 1;
                                    first_res = this.Strands[StrandNum1].Residues[res_num];
                                }
                            }

                            StrandNum1 = first_res.StrandNum;
                            if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
                            else j = first_res.ResStrandNum;

                            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", prev_res.SeqID, prev_res.StrandNum, m, first_res.SeqID, first_res.StrandNum, j, m - j);
                            shearNum += j - m;
                        }

                        if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
                        else j = first_res.ResStrandNum;

                        next_res = first_res.ShearNumNeigh;
                        StrandNum2 = next_res.StrandNum;
                        if (this.Strands[StrandNum2].Residues.Last().Z < this.Strands[StrandNum2].Residues[0].Z) l = this.Strands[StrandNum2].Count() - next_res.ResStrandNum;
                        else l = next_res.ResStrandNum;

                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", first_res.SeqID, first_res.StrandNum, j, next_res.SeqID, next_res.StrandNum, l);
                        prev_res = first_res;
                        StrandNum0 = prev_res.StrandNum;
                        m = j;
                        first_res = next_res;
                        StrandNum1 = first_res.StrandNum;

                    }
                    if (StrandNum1 == this.Strands.Count - 1)
                    {
                        if (first_res.ShearNumNeigh == null)
                        {
                            file.WriteLine("End of seq, restarting");
                            prev_res = first_res;
                            m = l;

                            //If a residue with no neighbors is reached - strand sticks up and need a new residue to restart; go to opposite end of strand and work back towards this residue
                            bool dir = false;
                            if (first_res.ResStrandNum > this.Strands[StrandNum1].Count() / 2)
                            {
                                res_num = 0;
                                dir = true;
                            }
                            else res_num = this.Strands[StrandNum1].Residues.Count() - 1;

                            while (first_res.ShearNumNeigh == null)
                            {
                                if (dir == true)
                                {
                                    res_num += 1;
                                    first_res = this.Strands[StrandNum1].Residues[res_num];
                                }
                                else
                                {
                                    res_num -= 1;
                                    first_res = this.Strands[StrandNum1].Residues[res_num];
                                }
                            }

                            StrandNum1 = first_res.StrandNum;
                            if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
                            else j = first_res.ResStrandNum;

                            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", prev_res.SeqID, prev_res.StrandNum, m, first_res.SeqID, first_res.StrandNum, j, m - j);
                            shearNum += j - m;
                        }

                        if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
                        else j = first_res.ResStrandNum;

                        next_res = first_res.ShearNumNeigh;
                        StrandNum2 = next_res.StrandNum;
                        if (this.Strands[StrandNum2].Residues.Last().Z < this.Strands[StrandNum2].Residues[0].Z) l = this.Strands[StrandNum2].Count() - next_res.ResStrandNum;
                        else l = next_res.ResStrandNum;

                        file.WriteLine("End of seq");
                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", first_res.SeqID, first_res.StrandNum, j, next_res.SeqID, next_res.StrandNum, l, Math.Abs(l - k));
                        shearNum += Math.Abs(l - k);
                        this.ShearNum = shearNum;
                    }


                }
            }

            public void rotate180(ref Chain _myChain)
            {
                double phi = Math.PI;
                Matrix3D rotationMatrixY2 = new Matrix3D(Math.Cos(phi), 0, Math.Sin(phi), 0, 0, 1, 0, 0, -1 * Math.Sin(phi), 0, Math.Cos(phi), 0, 0, 0, 0, 0);

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = rotationMatrixY2.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords);
                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
            }

            public void setBottom(ref Chain _myChain)
            {
                double adjuster = 12;
                //double adjuster = .0962 * this.Axis.Length + 8.3358;
                //double adjuster = 11.75;
                // if (this.Axis.Length < 30) adjuster = 11;
                // if (this.Axis.Length > 40) adjuster = 12.75;
                double adjustment = adjuster - this.Ncentroid.Z;

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {

                        Vector3D newCoords = new Vector3D(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.X, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Y, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Z + (adjustment));
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = newCoords;


                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
            }

            public void centerZ(ref Chain _myChain)
            {
                /* double avgZ = 0;
                 int totalRes = 0;
                 for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                 {
                     for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count; resCtr++)
                     {
                         avgZ = avgZ + this.Strands[strandCtr].Residues[resCtr].BackboneCoords["CA"].Z;
                         totalRes++;
                     }
                 }
                 avgZ = avgZ / totalRes;
                 */
                double avgZ = (this.Ncentroid.Z + this.Ccentroid.Z) / 2;

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {

                        Vector3D newCoords = new Vector3D(_myChain.Residues[resCtr].Atoms[atomCtr].Coords.X, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Y, _myChain.Residues[resCtr].Atoms[atomCtr].Coords.Z - (avgZ));
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = newCoords;


                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }


            }

            public void listLoopRes(ref Chain _myChain)
            {

                string fileLocation = Global.MONO_OUTPUT_DIR + "\\Loop_out_" + this.PdbName + ".py";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                {

                    file.WriteLine("from pymol import cmd, stored");
                    this.LoopResies = new List<Res>();
                    if (this.Strands.Count > 0 && this.Strands[0].Residues[0].ResNum < 18)
                    {
                        file.WriteLine("cmd.select(\"loopFirst\", \"resi {0}-{1} & chain{2}\")", _myChain.Residues[0].SeqID, this.Strands[0].Residues[0].SeqID, this.ChainName);
                        file.WriteLine("cmd.color ( \"salmon\", \"loopFirst\")");
                        for (int resCtr = _myChain.Residues[0].ResNum; resCtr < this.Strands[0].Residues[0].ResNum; resCtr++)
                        {
                            _myChain.Residues[resCtr].Z = _myChain.Residues[resCtr].BackboneCoords["CA"].Z;
                            _myChain.Residues[resCtr].Inward = false;
                            this.LoopResies.Add(_myChain.Residues[resCtr]);
                        }
                    }

                    if (_myChain.Residues[_myChain.Residues.Count - 1].ResNum - this.Strands[this.Strands.Count - 1].Residues[this.Strands[this.Strands.Count - 1].Residues.Count - 1].ResNum < 18)
                    {
                        file.WriteLine("cmd.select(\"loopLast\", \"resi {0}-{1} & chain{2}\")", this.Strands[this.Strands.Count - 1].Residues[this.Strands[this.Strands.Count - 1].Residues.Count - 1].SeqID, _myChain.Residues[_myChain.Residues.Count - 1].SeqID, this.ChainName);
                        file.WriteLine("cmd.color ( \"salmon\", \"loopLast\")");
                        for (int resCtr = _myChain.Residues[0].ResNum; resCtr < this.Strands[0].Residues[0].ResNum; resCtr++)
                        {
                            _myChain.Residues[resCtr].Z = _myChain.Residues[resCtr].BackboneCoords["CA"].Z;
                            _myChain.Residues[resCtr].Inward = false;
                            this.LoopResies.Add(_myChain.Residues[resCtr]);
                        }
                    }
                    for (int strandCtr = 0; strandCtr < this.Strands.Count - 1; strandCtr++)
                    {

                        int loopLength = this.Strands[strandCtr + 1].Residues[0].ResNum - this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum;
                        //Console.WriteLine("LoopLength = {0}", loopLength);
                        if (loopLength < 18)
                        {

                            file.WriteLine("cmd.select(\"loop{0}\", \"resi {1}-{2} & chain{3}\")", strandCtr, this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].SeqID, this.Strands[strandCtr + 1].Residues[0].SeqID, this.ChainName);
                            file.WriteLine("cmd.color ( \"salmon\", \"loop{0}\")", strandCtr);
                            for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                            {
                                _myChain.Residues[resCtr].Z = _myChain.Residues[resCtr].BackboneCoords["CA"].Z;
                                //double direction = Vector3D.AngleBetween(_myChain.Residues[resCtr].Direction, this.Axis);
                                //Res myRes = _myChain.Residues[resCtr];
                                //double angle = Vector3D.AngleBetween(myRes.BackboneCoords["CA"] - ((myRes.BackboneCoords["N"] + myRes.BackboneCoords["C"]) / 2), this.Axis);
                                //if (direction > 150)
                                //{
                                //    if (angle > 90)
                                //    {
                                //        myRes.Inward = true;
                                //        Console.WriteLine("resi {0} is inward", myRes.SeqID);
                                //    }
                                //    else myRes.Inward = false;
                                //}
                                //else if (direction < 30)
                                //{
                                //    if (angle < 90)
                                //    {
                                //        myRes.Inward = true;
                                //        Console.WriteLine("resi {0} is inward", myRes.SeqID);
                                //    }
                                //    else myRes.Inward = false;
                                //}
                                //else myRes.Inward = false;
                                _myChain.Residues[resCtr].Inward = false;
                                this.LoopResies.Add(_myChain.Residues[resCtr]);

                            }
                        }
                        else if (loopLength < 50)
                        {
                            bool enteringBarrel = false;

                            Vector3D direction = new Vector3D();
                            direction = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[0].BackboneCoords["CA"];
                            double angle = Vector3D.AngleBetween(this.Axis, direction);


                            if (angle < 90 && this.Ccentroid.Z > 0)
                            {
                                for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                                {
                                    double lowestZ = this.Ccentroid.Z;
                                    if (_myChain.Residues[resCtr].BackboneCoords["CA"].Z < lowestZ - 6)
                                    {
                                        enteringBarrel = true;
                                        break;
                                    }
                                }
                            }
                            else if (angle < 90 && this.Ccentroid.Z < 0)
                            {
                                for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                                {
                                    double highestZ = this.Ccentroid.Z;
                                    if (_myChain.Residues[resCtr].BackboneCoords["CA"].Z > highestZ + 6)
                                    {
                                        enteringBarrel = true;
                                        break;
                                    }
                                }
                            }
                            else if (angle > 90 && this.Ncentroid.Z < 0)
                            {
                                for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                                {
                                    double highestZ = this.Ncentroid.Z;
                                    if (_myChain.Residues[resCtr].BackboneCoords["CA"].Z > highestZ + 3)
                                    {
                                        enteringBarrel = true;
                                        break;
                                    }
                                }
                            }
                            else if (angle > 90 && this.Ncentroid.Z > 0)
                            {
                                for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                                {
                                    double lowestZ = this.Ncentroid.Z;
                                    if (_myChain.Residues[resCtr].BackboneCoords["CA"].Z < lowestZ - 6)
                                    {
                                        enteringBarrel = true;
                                        break;
                                    }
                                }
                            }

                            if (enteringBarrel == false)
                            {
                                file.WriteLine("cmd.select(\"loop{0}\", \"resi {1}-{2} & chain{3}\")", strandCtr, this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].SeqID, this.Strands[strandCtr + 1].Residues[0].SeqID, this.ChainName);
                                file.WriteLine("cmd.color ( \"salmon\", \"loop{0}\")", strandCtr);
                                for (int resCtr = this.Strands[strandCtr].Residues[this.Strands[strandCtr].Residues.Count - 1].ResNum + 1; resCtr < this.Strands[strandCtr + 1].Residues[0].ResNum; resCtr++)
                                {
                                    _myChain.Residues[resCtr].Z = _myChain.Residues[resCtr].BackboneCoords["CA"].Z;

                                    _myChain.Residues[resCtr].Inward = false;
                                    this.LoopResies.Add(_myChain.Residues[resCtr]);
                                }
                            }
                        }

                    }
                }
            }

            public void rotateToZ(ref Chain _myChain)
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


                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = rotationMatrixZ.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords);
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = rotationMatrixY.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords);



                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }

                if (this.Strands[1].Residues[0].BackboneCoords["CA"].Z > this.Strands[0].Residues[0].BackboneCoords["CA"].Z)
                {
                    double phi = Math.PI;
                    Matrix3D rotationMatrixY2 = new Matrix3D(Math.Cos(phi), 0, Math.Sin(phi), 0, 0, 1, 0, 0, -1 * Math.Sin(phi), 0, Math.Cos(phi), 0, 0, 0, 0, 0);
                    Vector3D axis3 = new Vector3D();
                    axis3 = rotationMatrixY2.Transform(axis2);


                    for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                    {
                        for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                        {
                            _myChain.Residues[resCtr].Atoms[atomCtr].Coords = rotationMatrixY2.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords);
                            if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                            {

                                _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                            }
                        }
                    }
                }
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);

                Vector3D translate = new Vector3D(this.Ccentroid.X, this.Ccentroid.Y, 0);

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {
                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords -= translate;

                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
                setCEllipseCoords(ref _myChain);
                setNEllipseCoords(ref _myChain);
            }

            public void rotateToZmem(ref Chain _myChain)
            {
                double halfLength = this.Axis.Length * .5;

                double theta = Math.Acos(14.5 / halfLength);
                Matrix3D rotationMatrixY = new Matrix3D(Math.Cos(theta), 0, Math.Sin(theta), 0, 0, 1, 0, 0, -1 * Math.Sin(theta), 0, Math.Cos(theta), 0, 0, 0, 0, 0);
                //Console.WriteLine("theta ={0}", theta * 180 / Math.PI);

                for (int resCtr = 0; resCtr < _myChain.Residues.Count; resCtr++)
                {
                    for (int atomCtr = 0; atomCtr < _myChain.Residues[resCtr].Atoms.Count; atomCtr++)
                    {

                        _myChain.Residues[resCtr].Atoms[atomCtr].Coords = rotationMatrixY.Transform(_myChain.Residues[resCtr].Atoms[atomCtr].Coords);



                        if (_myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "CA" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "C" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "N" || _myChain.Residues[resCtr].Atoms[atomCtr].AtomName == "O")
                        {

                            _myChain.Residues[resCtr].BackboneCoords[_myChain.Residues[resCtr].Atoms[atomCtr].AtomName] = _myChain.Residues[resCtr].Atoms[atomCtr].Coords;
                        }
                    }
                }
            }

            public void makeBarrelCircular(ref Chain _myChain)
            {
                //fill in the barrel
                /*if (PdbName == "4K3C")
                {
                    for (int strndCtr2 = 0; strndCtr2 < this.protoBarrel.Count; strndCtr2++)
                    {
                        if (this.protoBarrel[strndCtr2].Count == 0) this.protoBarrel.RemoveAt(strndCtr2);
                    }
                }
                if (PdbName == "4K3C")
                {
                    for (int strndCtr2 = 0; strndCtr2 < 1; strndCtr2++)
                    {
                        this.protoBarrel.RemoveAt(0);
                    }
                    this.protoBarrel.RemoveAt(0);
                }*/

                for (int strndCtr2 = 0; strndCtr2 < this.protoBarrel.Count; strndCtr2++)
                {
                    if (this.protoBarrel[strndCtr2].Count == 0) this.protoBarrel.RemoveAt(strndCtr2);
                    for (int i = this.protoBarrel[strndCtr2][0]; i < this.protoBarrel[strndCtr2][this.protoBarrel[strndCtr2].Count - 1]; i++)
                    {
                        if (this.protoBarrel[strndCtr2].Contains(i) == false)
                        {
                            this.protoBarrel[strndCtr2].Add(i);
                            //Console.WriteLine("added res to strand {0}", strndCtr2);
                            this.protoBarrel[strndCtr2].Sort();
                        }


                    }
                }


                for (int strndCtr = 0; strndCtr < this.protoBarrel.Count;)
                {
                    bool HbondedPrev = false;
                    bool HbondedNext = false;


                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count; resCtr++)
                    { // for each residue
                        for (int nCtr = 0; nCtr < _myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors.Count; nCtr++)
                        {// for each neighbor
                            if (strndCtr == 0)
                            {
                                if (this.protoBarrel[this.protoBarrel.Count - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedPrev = true;
                                if (this.protoBarrel[1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedNext = true;
                            }
                            else if (strndCtr == this.protoBarrel.Count - 1)
                            {
                                if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedPrev = true;
                                if (this.protoBarrel[0].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedNext = true;
                            }
                            else
                            {
                                if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedPrev = true;
                                if (this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr])) HbondedNext = true;
                            }


                        }


                    }
                    if (HbondedPrev == false || HbondedNext == false)
                    {
                        this.protoBarrel.Remove(this.protoBarrel[strndCtr]);
                        //Console.WriteLine("removing strand {0}", strndCtr);
                    }
                    else strndCtr++;
                    /*
                    //Ryan edit
                    if (PdbName == "4K3C")
                    {
                        if (HbondedPrev == false && HbondedNext == false)
                        {
                            this.protoBarrel.Remove(this.protoBarrel[strndCtr]);
                            Console.WriteLine("removing strand {0}", strndCtr);
                        }
                        else if (this.protoBarrel[strndCtr].Count == 0) this.protoBarrel.Remove(this.protoBarrel[strndCtr]);
                        else strndCtr++;
                    }
                    else
                    {
                        if (HbondedPrev == false || HbondedNext == false)
                        {
                            this.protoBarrel.Remove(this.protoBarrel[strndCtr]);
                            Console.WriteLine("removing strand {0}", strndCtr);
                        }
                        else strndCtr++;
                    }
                    */
                }

            }

            public void makeBarrelCircularPatt(ref Chain chain) //This is based on the premise that H-bonding between beta sheets in chains only occurs in the barrel
            {
                //fill in the barrel
                int nextchain;
                int prevchain;
                double d; double minD = 3.00;
                int i; //strand count
                Vector3D strand1vec;
                Vector3D strand2vec;

                #region Combine sequential strands, i.e. upper and lower
                for (i = 1; i < protoBarrel.Count;)//for each strand
                {
                    strand1vec = chain.Residues[protoBarrel[i - 1].Last()].BackboneCoords["CA"] - chain.Residues[protoBarrel[i - 1][0]].BackboneCoords["CA"];
                    strand2vec = chain.Residues[protoBarrel[i].Last()].BackboneCoords["CA"] - chain.Residues[protoBarrel[i][0]].BackboneCoords["CA"];
                    //Console.WriteLine("{0} {1} {2} {3} {4}", strand1vec, strand1vec.Length, strand2vec, strand2vec.Length, (strand1vec + strand2vec).Length);
                    double longest = strand1vec.Length;
                    if (strand2vec.Length > longest) longest = strand2vec.Length;

                    if ((strand1vec + strand2vec).Length > longest)
                    {
                        protoBarrel[i - 1].AddRange(protoBarrel[i]);
                        protoBarrel.RemoveRange(i, 1);
                        //Console.WriteLine("Combined two strands");
                    }
                    else if (protoBarrel[i - 1].Count <= 3) protoBarrel.RemoveRange(i - 1, 1);
                    else i++;
                }
                #endregion

                for (i = 0; i < protoBarrel.Count;)//for each strand
                {
                    int HbondedPrev = 0;
                    int HbondedNext = 0;

                    if (i == 0)
                    {
                        prevchain = protoBarrel.Count - 1;
                        nextchain = i + 1;
                    }
                    else if (i == protoBarrel.Count - 1)
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
                    for (int res1ctr = 0; res1ctr < protoBarrel[i].Count - 1; res1ctr++)//for each res in strand
                    {
                        Res Residue1 = chain.Residues[protoBarrel[i][res1ctr]];

                        for (int atomCtr1 = 0; atomCtr1 < Residue1.Atoms.Count; atomCtr1++)
                        {
                            if (Residue1.Atoms[atomCtr1].AtomName == "N" && Residue1.ThreeLetCode != "PRO")
                            {
                                Atom N = Residue1.Atoms[atomCtr1];

                                for (int res3ctr = 0; res3ctr < protoBarrel[prevchain].Count; res3ctr++) //Check previous strand for neighbors
                                {
                                    Res Residue2 = chain.Residues[protoBarrel[prevchain][res3ctr]];

                                    for (int atomCtr2 = 0; atomCtr2 < Residue2.Atoms.Count; atomCtr2++)
                                    {
                                        if (Residue2.Atoms[atomCtr2].AtomName == "O")
                                        {
                                            Atom O = Residue2.Atoms[atomCtr2];
                                            atomCtr2 = Residue2.Atoms.Count;  // to skip the loop once "O" atom is found, and thus go to next residue number

                                            d = (O.Coords - N.Hydrogen).Length;
                                            if (d < minD && (Residue1.SSType == "B" && Residue2.SSType == "B"))
                                            {
                                                HbondedPrev += 1;
                                                if (Residue1.Neighbors.Contains(Residue2.ResNum) == false) Residue1.Neighbors.Add(Residue2.ResNum);
                                                if (Residue2.Neighbors.Contains(Residue1.ResNum) == false) Residue2.Neighbors.Add(Residue1.ResNum);

                                            }
                                        }
                                    }
                                }
                            CheckNextChain:
                                for (int res3ctr = 0; res3ctr < protoBarrel[nextchain].Count; res3ctr++) //Check next strand for neighbors
                                {
                                    Res Residue2 = chain.Residues[protoBarrel[nextchain][res3ctr]];
                                    for (int atomCtr2 = 0; atomCtr2 < Residue2.Atoms.Count; atomCtr2++)
                                    {
                                        if (Residue2.Atoms[atomCtr2].AtomName == "O")
                                        {

                                            Atom O = Residue2.Atoms[atomCtr2];
                                            atomCtr2 = Residue2.Atoms.Count;  // to skip the loop once "O" atom is found, and thus go to next residue number

                                            d = (O.Coords - N.Hydrogen).Length;
                                            if (d < minD && (Residue1.SSType == "B" && Residue2.SSType == "B"))
                                            {
                                                HbondedNext += 1;
                                                if (Residue1.Neighbors.Contains(Residue2.ResNum) == false) Residue1.Neighbors.Add(Residue2.ResNum);
                                                if (Residue2.Neighbors.Contains(Residue1.ResNum) == false) Residue2.Neighbors.Add(Residue1.ResNum);
                                            }
                                        }
                                    }
                                }
                            }
                        }


                    }
                    #endregion

                    if ((HbondedPrev == 0 && HbondedNext == 0) || protoBarrel[i].Count == 0)
                    {
                        protoBarrel.RemoveRange(i, 1);
                        //Console.WriteLine("removing strand {0}", i);
                    }
                    else if (protoBarrel[i].Count <= 4 && (HbondedPrev == 0 || HbondedNext == 0))
                    {
                        protoBarrel.RemoveRange(i, 1);
                        //Console.WriteLine("removing strand {0}", i);
                    }
                    else i++;

                    HbondedNext = 0;
                    HbondedPrev = 0;

                } //End of checking each chain
            }

            public void removeNonBarrelRes(ref Chain _myChain)
            {// removes residues that are not h-bonded to the next or the previous strand.
             //bool deletedSomething = false;
                for (int strndCtr = 0; strndCtr < this.protoBarrel.Count; strndCtr++)
                {
                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count;)
                    {
                        bool markForDeletion = true;

                        if (_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors.Count != 0)

                        {
                            for (int nCtr = 0; nCtr < _myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors.Count; nCtr++)
                            {
                                if (strndCtr == 0)
                                {
                                    if (this.protoBarrel[this.protoBarrel.Count - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else if (strndCtr == this.protoBarrel.Count - 1)
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[0].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[this.protoBarrel[strndCtr][resCtr]].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                            }
                        }

                        if (markForDeletion == true)
                        {
                            this.protoBarrel[strndCtr].RemoveAt(resCtr);
                            //Console.WriteLine("removed a residue from Strand {0})", strndCtr);
                        }
                        else resCtr++;
                        //resCtr--;

                    }

                    this.protoBarrel[strndCtr].Sort();

                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count;)
                    {// if there are no residues before or after 5 that are part of the sheet delete that residue too.
                        bool markForDeletion = true;
                        for (int ctr2 = this.protoBarrel[strndCtr][resCtr] - 5; ctr2 < this.protoBarrel[strndCtr][resCtr] + 5; ctr2++)
                        {
                            if (ctr2 > -1 && ctr2 != this.protoBarrel[strndCtr][resCtr])
                            {
                                if (this.protoBarrel[strndCtr].Contains(ctr2) == true)
                                {
                                    markForDeletion = false;
                                }

                            }
                        }
                        if (markForDeletion == true)
                        {
                            this.protoBarrel[strndCtr].RemoveAt(resCtr);
                            //Console.WriteLine("removed a residue from Strand {0} because no nearby residues", strndCtr);
                            //deletedSomething = true;
                        }
                        else resCtr++;
                    }

                    if (this.protoBarrel[strndCtr].Count == 0) protoBarrel.RemoveAt(strndCtr);
                }
            }

            public void removeNonBarrelResPatt(ref Chain _myChain)
            {// removes residues that are not h-bonded to the next or the previous strand.
             //bool deletedSomething = false;
                for (int strndCtr = 0; strndCtr < this.protoBarrel.Count; strndCtr++)
                {
                    int halfway;
                    if (this.protoBarrel[strndCtr].Count > 4) halfway = this.protoBarrel[strndCtr].Count / 2;
                    else halfway = this.protoBarrel[strndCtr].Count - 2;
                    for (int resCtr = this.protoBarrel[strndCtr].Count - 1; resCtr > halfway;)
                    {
                        bool markForDeletion = true;
                        int index2 = protoBarrel[strndCtr][resCtr];
                        if (_myChain.Residues[index2].Neighbors.Count != 0)
                        {
                            for (int nCtr = 0; nCtr < _myChain.Residues[index2].Neighbors.Count; nCtr++)
                            {
                                if (strndCtr == 0)
                                {
                                    if (this.protoBarrel[this.protoBarrel.Count - 1].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else if (strndCtr == this.protoBarrel.Count - 1)
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true || this.protoBarrel[0].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                                else
                                {
                                    if (this.protoBarrel[strndCtr - 1].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true || this.protoBarrel[strndCtr + 1].Contains(_myChain.Residues[index2].Neighbors[nCtr]) == true) markForDeletion = false;

                                }
                            }
                        }

                        if (markForDeletion == true)
                        {
                            this.protoBarrel[strndCtr].RemoveAt(resCtr);
                            //Console.WriteLine("removed a residue from Strand {0})", strndCtr);
                        }
                        //else resCtr++;
                        resCtr--;

                    }

                    this.protoBarrel[strndCtr].Sort();

                    for (int resCtr = 0; resCtr < this.protoBarrel[strndCtr].Count;)
                    {// if there are no residues before or after 5 that are part of the sheet delete that residue too.
                        bool markForDeletion = true;
                        for (int ctr2 = this.protoBarrel[strndCtr][resCtr] - 5; ctr2 < this.protoBarrel[strndCtr][resCtr] + 5; ctr2++)
                        {
                            if (ctr2 > -1 && ctr2 != this.protoBarrel[strndCtr][resCtr])
                            {
                                if (this.protoBarrel[strndCtr].Contains(ctr2) == true)
                                {
                                    markForDeletion = false;
                                }

                            }
                        }
                        if (markForDeletion == true)
                        {
                            //this.protoBarrel[strndCtr].RemoveAt(resCtr);
                            //Console.WriteLine("removed a residue from Strand {0})", strndCtr);
                            //deletedSomething = true;
                        }
                        else resCtr++;
                    }
                }
            }
			
	        public bool check_neighbors(int strand_num, Res Res1)
	        {
	            bool markForDeletion = true;
	            for (int nCtr = 0; nCtr < Res1.Neighbors.Count; nCtr++)
	            {
	                if (strand_num == 0)
	                {
	                    if (this.protoBarrel[this.protoBarrel.Count - 1].Contains(Res1.Neighbors[nCtr]) == true || this.protoBarrel[strand_num + 1].Contains(Res1.Neighbors[nCtr]) == true) markForDeletion = false;

	                }
	                else if (strand_num == this.protoBarrel.Count - 1)
	                {
	                    if (this.protoBarrel[strand_num - 1].Contains(Res1.Neighbors[nCtr]) == true || this.protoBarrel[0].Contains(Res1.Neighbors[nCtr]) == true) markForDeletion = false;

	                }
	                else
	                {
	                    if (this.protoBarrel[strand_num - 1].Contains(Res1.Neighbors[nCtr]) == true || this.protoBarrel[strand_num + 1].Contains(Res1.Neighbors[nCtr]) == true) markForDeletion = false;

	                }
	            }
	            return markForDeletion;
	        }
			

            public void getTilts()
            {
                for (int strandCtr = 0; strandCtr < Strands.Count; strandCtr++)
                {
                    this.Strands[strandCtr].getTilts(this.Axis, strandCtr);
                }
            }

            public void setCEllipseCoords(ref Chain _myChain)
            {
                //this is the top (extracellular) ellipse
                List<Vector3D> myEllipse = new List<Vector3D>();
                Vector3D centroid = new Vector3D();

                for (int strandCtr = 0; strandCtr < this.protoBarrel.Count; strandCtr++)
                {
                    Vector3D firstCA = new Vector3D();
                    firstCA = _myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"];
                    Vector3D lastCA = new Vector3D();
                    lastCA = _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"];



                    if (strandCtr % 2 == 0)
                    {
                        myEllipse.Add(lastCA);
                        centroid += lastCA;

                    }
                    else
                    {
                        myEllipse.Add(firstCA);
                        centroid += firstCA;

                    }

                }
                this.CellipseCoords = myEllipse;
                centroid = centroid / this.protoBarrel.Count;
                this.Ccentroid = centroid;


            }

            public void setNEllipseCoords(ref Chain _myChain)
            {
                //this is the bottom (periplasmic) ellipse
                List<Vector3D> myEllipse = new List<Vector3D>();
                Vector3D centroid = new Vector3D();

                for (int strandCtr = 0; strandCtr < this.protoBarrel.Count; strandCtr++)
                {
                    Vector3D firstCA = new Vector3D();
                    firstCA = _myChain.Residues[this.protoBarrel[strandCtr][0]].BackboneCoords["CA"];
                    Vector3D lastCA = new Vector3D();
                    lastCA = _myChain.Residues[this.protoBarrel[strandCtr][this.protoBarrel[strandCtr].Count - 1]].BackboneCoords["CA"];


                    if (strandCtr % 2 == 0)
                    {
                        myEllipse.Add(firstCA);
                        centroid += firstCA;

                    }
                    else
                    {
                        myEllipse.Add(lastCA);
                        centroid += lastCA;
                    }

                }
                this.NellipseCoords = myEllipse;
                centroid = centroid / this.protoBarrel.Count;
                this.Ncentroid = centroid;


            }

            public void checkStrands()
            {
                int direction = 0;
                for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                {
                    if (strandCtr % 2 == 0)
                    {
                        if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[0].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[2].BackboneCoords["CA"]), this.Axis) < Vector3D.AngleBetween(this.Strands[strandCtr].Residues[2].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[0].BackboneCoords["CA"], this.Axis)) direction--;
                        else direction++;
                    }
                    else
                    {
                        if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[0].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[2].BackboneCoords["CA"]), this.Axis) < Vector3D.AngleBetween(this.Strands[strandCtr].Residues[2].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[0].BackboneCoords["CA"], this.Axis)) direction++;
                        else direction--;
                    }
                }
                if (direction < 0)
                {
                    this.Direction = false;
                }
                else this.Direction = true;

                for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                {
                    if (this.Direction == false)
                    {
                        int markForDeletion = 0;

                        if (strandCtr % 2 == 0)
                        {
                            for (int residueCtr = 1; residueCtr < this.Strands[strandCtr].Residues.Count - 1; residueCtr++)
                            {
                                if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"]), this.Axis) > Vector3D.AngleBetween(this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"], this.Axis)) markForDeletion++;
                                else markForDeletion--;
                            }

                        }
                        else
                        {
                            for (int residueCtr = 1; residueCtr < this.Strands[strandCtr].Residues.Count - 1; residueCtr++)
                            {
                                if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"]), this.Axis) < Vector3D.AngleBetween(this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"], this.Axis)) markForDeletion++;
                                else markForDeletion--;
                            }
                        }
                        if (markForDeletion > 0)
                        {
                            //Console.WriteLine("Removed Strand {0} !!!!!!", strandCtr);
                            this.Strands.Remove(this.Strands[strandCtr]);


                            if (strandCtr == 0) this.Direction = true;
                            else strandCtr--;

                            //break;
                        }
                    }


                    if (this.Direction == true)
                    {

                        int markForDeletion = 0;

                        if (strandCtr % 2 == 0)
                        {
                            for (int residueCtr = 1; residueCtr < this.Strands[strandCtr].Residues.Count - 1; residueCtr++)
                            {
                                double forwardAngle = Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"]), this.Axis);
                                double backwardAngle = Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"]), this.Axis);
                                if (forwardAngle < backwardAngle) markForDeletion++;
                                //if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"]), this.Axis) < Vector3D.AngleBetween(this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"], this.Axis)) markForDeletion++;
                                else markForDeletion--;
                            }

                        }
                        else
                        {
                            for (int residueCtr = 1; residueCtr < this.Strands[strandCtr].Residues.Count - 1; residueCtr++)
                            {
                                if (Vector3D.AngleBetween((this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"]), this.Axis) > Vector3D.AngleBetween(this.Strands[strandCtr].Residues[residueCtr + 1].BackboneCoords["CA"] - this.Strands[strandCtr].Residues[residueCtr - 1].BackboneCoords["CA"], this.Axis)) markForDeletion++;
                                else markForDeletion--;
                            }
                        }
                        if (markForDeletion > 0)
                        {

                            //Console.WriteLine("Removed Strand {0} !!!!!!", strandCtr);
                            this.Strands.Remove(this.Strands[strandCtr]);

                            strandCtr--;
                            if (strandCtr == 0) this.Direction = false;
                            //break;
                        }

                    }
                }
                //make sure strandNum is correct. 
                for (int strandCtr = 0; strandCtr < this.Strands.Count; strandCtr++)
                {
                    this.Strands[strandCtr].StrandNum = strandCtr;
                    for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count; resCtr++)
                    {
                        this.Strands[strandCtr].Residues[resCtr].StrandNum = strandCtr;
                    }
                }

            }

            public void writeAminoAcidsTypesToFile(ref Chain _myChain, string outputDirectory)
            {
                string fileLocation3 = outputDirectory + "AminoAcids\\AminoAcidTypes" + this.PdbName + ".txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation3))
                {
                    file.WriteLine("Sl.No\t*PDB ID*\t*Strands Per Chain*\t*Chains Per Protein*\t*Residue Number*\t*AminoAcid*\t*Chain*\t*Strand Number*\t*Interface*");

                    int sl_no = 1;
                    bool interface_value;
                    for (int strandCtr = 0; strandCtr < protoBarrel.Count; strandCtr++)
                    {
                        if (strandCtr == 0 || strandCtr == protoBarrel.Count - 1) // to check the residues that are at interface of the different chains
                        {
                            interface_value = true;
                        }
                        else interface_value = false;

                        for (int residues = protoBarrel[strandCtr][0]; residues <= protoBarrel[strandCtr].Last(); residues++)
                        {
                            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", sl_no, this.PdbName, protoBarrel.Count, protoBarrel.Count, _myChain.Residues[residues].SeqID, _myChain.Residues[residues].ThreeLetCode, _myChain.ChainName, strandCtr, interface_value);
                            sl_no++;
                        }
                    }
                }
            }

            public void getTyrVector(string outputDirectory, ref Chain myChain)
            {
                string newLine;
                Vector3D strandVect;
                Vector3D ringVect;
                Vector3D Cg;
                Vector3D OH;

                string fileLocation = outputDirectory + "FoldingPatterns/TyrVect" + this.PdbName + ".txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                {
                    newLine = "Res" + "\t" + "ResNum" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Angle Btwn Ring/Strand" + "\t" + "Dihedral";
                    file.WriteLine(newLine);

                    foreach (Strand strand in this.Strands)
                    {
                        for (int resIndex = 0; resIndex < strand.Residues.Count; resIndex++)
                        {
                            Res res1 = strand.Residues[resIndex];
                            if (res1.ThreeLetCode == "TYR")
                            {
                                //Use CA coords of neighboring residues to determine strand vector
                                if (resIndex == strand.Residues.Count - 1) strandVect = strand.Residues[resIndex].BackboneCoords["CA"] - strand.Residues[resIndex - 2].BackboneCoords["CA"];
                                else if (resIndex == 0) strandVect = strand.Residues[resIndex + 2].BackboneCoords["CA"] - strand.Residues[resIndex].BackboneCoords["CA"];
                                else strandVect = strand.Residues[resIndex + 1].BackboneCoords["CA"] - strand.Residues[resIndex - 1].BackboneCoords["CA"];

                                Cg = res1.Atoms[0].Coords;
                                OH = res1.Atoms[0].Coords;
                                foreach (Atom atom1 in res1)
                                {
                                    if (atom1.AtomName == "CG") Cg = atom1.Coords;
                                    if (atom1.AtomName == "OH") OH = atom1.Coords;
                                }
                                ringVect = OH - Cg;
                                double dihedral = SharedFunctions.CalculateTorsion(myChain.Residues[res1.ResNum - 2].BackboneCoords["CA"], res1.BackboneCoords["CA"], Cg, OH);
                                file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", res1.ThreeLetCode, res1.SeqID, strand.StrandNum, res1.ChainName, Vector3D.AngleBetween(ringVect, strandVect), dihedral);
                            }

                        }
                    }
                }
            }

            public IEnumerator<Strand> GetEnumerator()
            {
                return Strands.GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return this.GetEnumerator();
            }
        }
		
	    public class SBarrel : Barrel
	    {
	        public List<Strand> Strands { get; set; }
	        public List<double> PrevTwists { get; set; }
	        public double AvgTilt { get; set; }
	        public double AvgTilt_even { get; set; }
	        public double AvgTilt_odd { get; set; }
	        public List<List<int>> protoBarrel { get; set; }
	        public string PdbName { get; set; }
	        public string ChainName { get; set; }
	        public List<double> StrandLength { get; set; }
	        public int ShearNum { get; set; }

	        public static string path = Global.SOL_OUTPUT_DIR;

	        //barrel constructor 
	        public SBarrel(Chain _myChain, Protein _myProtein)
	        {
	            this.Strands = new List<Strand>();
	            this.ChainName = _myChain.ChainName;
	            this.protoBarrel = new List<List<int>>();
	            this.PdbName = _myChain.PdbName;

	            //Writing info about residues
	            /*string fileLocation2 = path + "PhiPsiBBAngle/" + this.PdbName + ".txt";
	                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
	                {
	                        foreach (Res res in _myChain)
	                        {
	                            if (res.ResNum != 0)
	                            {
	                                file.WriteLine("Residue {0}({1}), SSType {2}, phi {4}, psi {5}, Angle formed at this residue {3}", res.SeqID, res.ThreeLetCode, res.SSType, (Vector3D.AngleBetween(_myChain.Residues[res.ResNum - 1].Direction, res.Direction)), res.Phi, res.Psi);
	                            }
	                        }
	                }*/

	            createStrandsDSSPonly(ref _myChain);
	            for (int strandCtr = 0; strandCtr < protoBarrel.Count; strandCtr++)
	            {
	                Strand newStrand = new Strand(_myChain, protoBarrel[strandCtr][0], protoBarrel[strandCtr][protoBarrel[strandCtr].Count - 1], strandCtr);
	                this.Strands.Add(newStrand);
	            }

	            //SharedFunctions.writePymolScriptForStrands(this.Strands, Program.soloutDirectory, Program.MacsolDBDir, this.PdbName);
	            this.PrevTwists = SharedFunctions.writeTwists(this.Strands, path, this.PdbName);
                Dictionary<string, string> Loops = SharedFunctions.getLoopTurns(this.Strands, ref _myChain, Global.MONO_OUTPUT_DIR, this.PdbName);

	            //find shear number July 22, 2014
	            //int ShearNum = shearNum(ref _myChain);
	            //getShearNum();

	        }

	        public void createStrandsDSSPonly(ref Chain _myChain)
	        {
	            List<int> myStrand1 = new List<int>();
	            int usually0 = 0;
	            int usuallyEnd = _myChain.ResidueCount;

	            #region DefineStrands
	            for (int res1ctr = usually0; res1ctr < usuallyEnd; res1ctr++)
	            {
	                Res Residue1 = _myChain.Residues[res1ctr];

	                if (Residue1.DSSP != "E")
	                {
	                    continue;
	                }
	                else
	                {
	                    //add Residue1.ResNum
	                    if (myStrand1.Count > 0)
	                    {
	                        if (Residue1.ResNum == myStrand1.Last() + 1)
	                        {
	                            myStrand1.Add(Residue1.ResNum);
	                        }
	                        else if (Residue1.ResNum <= myStrand1.Last() + 2)
	                        { //contiguous strands, only one residue missing at most
	                            myStrand1.Add(Residue1.ResNum);
	                        }
	                        else
	                        {
	                            //have reached a break in strands, so add old strand and start a new one
	                            List<int> newList = new List<int>();
	                            for (int i = myStrand1[0]; i <= myStrand1.Last(); i++)
	                            {
	                                newList.Add(i);
	                            }
	                            protoBarrel.Add(newList);
	                            myStrand1.Clear();
	                            //Start new strand from here 
	                            myStrand1.Add(Residue1.ResNum);
	                        }
	                    }
	                    else
	                    {
	                        myStrand1.Add(Residue1.ResNum);
	                    }
	                }

	            } //end of going through each res and adding to barrel. 
	            if (myStrand1.Count > 0) { protoBarrel.Add(myStrand1); }

	            /*if (PdbName.ToUpper() == "1EST") //added 6-22-17 for v4
	            {
	                //protoBarrel.RemoveRange(0, 2);
	                //protoBarrel[1].RemoveRange(25, 23);
	                //protoBarrel[9].RemoveRange(0, 5);

	                myStrand1 = new List<int>();
	                myStrand1.AddRange(Enumerable.Range(protoBarrel[1][9], 10));
	                protoBarrel.Insert(2, myStrand1);
	                protoBarrel[1].RemoveRange(7, 12);
	                protoBarrel[8].AddRange(Enumerable.Range(protoBarrel[9][0], protoBarrel[9].Count));
	                protoBarrel.RemoveAt(9);
	            }
	            if (PdbName.ToUpper() == "2SGA") //added 6-22-17 for v4
	            {
	                myStrand1 = new List<int>();
	                myStrand1.AddRange(Enumerable.Range(protoBarrel[2][4], 5));
	                protoBarrel.Insert(3, myStrand1);
	                protoBarrel[2].RemoveRange(3, 6);
	                //protoBarrel[8].AddRange(Enumerable.Range(protoBarrel[9][0], protoBarrel[9].Count));
	                //protoBarrel.RemoveAt(9);
	            }*/
	            #endregion
	        }

	        public void checkStrandDefnsDSSP(ref Chain _myChain) //Added 6-5-17 for loops
	        {
	            for (int strandNum = 0; strandNum < this.protoBarrel.Count; strandNum++)
	            {
	                //Remove turn residues
	                while (_myChain.Residues[this.protoBarrel[strandNum][0]].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum][0]].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(0);
	                    Console.WriteLine("Removed res{1} from beg of strand {0}", strandNum, removed_res);
	                }

	                while (_myChain.Residues[this.protoBarrel[strandNum].Last()].DSSP == "T")
	                {
	                    int removed_res = _myChain.Residues[this.protoBarrel[strandNum].Last()].ResNum;
	                    this.protoBarrel[strandNum].RemoveAt(this.protoBarrel[strandNum].Count - 1);
	                    Console.WriteLine("Removed res{1} from end of strand {0}", strandNum, removed_res);
	                }


	                //Add to the beginning of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum][0] - 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum][0] - 1 > this.protoBarrel[strandNum - 1].Last() + 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum][0] - 1], ref _myChain);
	                                Console.WriteLine("added res{1} to beg of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum][0] - 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break; }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into prev strand"); break; }

	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }

	                //For adding to the ends of strands
	                try
	                {
	                    while (_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].DSSP == "E")
	                    {
	                        try
	                        {
	                            if (this.protoBarrel[strandNum].Last() + 1 < this.protoBarrel[strandNum + 1][0] - 1)
	                            {
	                                this.protoBarrel[strandNum].Add(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                //checkDSSPNeighs(_myChain.Residues[this.protoBarrel[strandNum].Last() + 1], ref _myChain);
	                                Console.WriteLine("added res{1} to end of strand {0}", strandNum, _myChain.Residues[this.protoBarrel[strandNum].Last() + 1].ResNum);
	                                this.protoBarrel[strandNum].Sort();
	                            }
	                            else { break; }
	                        }
	                        catch (ArgumentOutOfRangeException) { Console.WriteLine("Ran into next strand"); break; }
	                    }
	                }
	                catch (ArgumentOutOfRangeException) { continue; }
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

	        // this will calculate the shear number of the beta barrel.  see murzin lesk and chothia 1994 http://www.mrc-lmb.cam.ac.uk/tcb/pdf/chc/97_jmb_236_1369_94.pdf
	        public int shearNum(ref Chain _myChain)
	        {
	            string filelocation = path + "ShearNum_" + PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(filelocation))
	            {
	                // first determine nearest hbonding neighbor with N coming from N-term  (arbitrary)           
	                for (int strandCtr = 0; strandCtr < this.Strands.Count(); strandCtr++)
	                {
	                    int strand2 = strandCtr + 1;
	                    if (strandCtr == this.Strands.Count() - 1) strand2 = 0;

	                    for (int resCtr = 0; resCtr < this.Strands[strandCtr].Residues.Count(); resCtr++)
	                    {

	                        Res firstRes = this.Strands[strandCtr].Residues[resCtr];
	                        double minD = 3.4;
	                        if (firstRes.Atoms[0].AtomName == "N" && firstRes.ThreeLetCode != "PRO")
	                        {
	                            Atom firstAtom = firstRes.Atoms[0];

	                            for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
	                            {
	                                Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[3];
	                                if (secondAtom.AtomName == "O")
	                                {
	                                    double d = (secondAtom.Coords - firstAtom.Coords).Length;
	                                    if (d < minD)
	                                    {
	                                        minD = d;
	                                        firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
	                                        firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
	                                        firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length;
	                                    }

	                                }

	                            }


	                        }
	                        // looking at o=n as well just to be safe... 
	                        if (firstRes.h_bonder == 0 && firstRes.Atoms[3].AtomName == "O")
	                        {
	                            Atom firstAtom = firstRes.Atoms[3];

	                            for (int resCtr2 = 0; resCtr2 < this.Strands[strand2].Residues.Count(); resCtr2++)
	                            {
	                                Atom secondAtom = this.Strands[strand2].Residues[resCtr2].Atoms[0];
	                                if (secondAtom.AtomName == "N")
	                                {
	                                    double d = (secondAtom.Coords - firstAtom.Coords).Length;

	                                    if (d < minD)
	                                    {
	                                        minD = d;
	                                        firstRes.h_bonder = this.Strands[strand2].Residues[resCtr2].ResNum;
	                                        firstRes.h_bonderID = this.Strands[strand2].Residues[resCtr2].SeqID;
	                                        firstRes.bDist = (this.Strands[strand2].Residues[resCtr2].Atoms[1].Coords - firstRes.Atoms[1].Coords).Length;
	                                    }
	                                }
	                            }
	                        }
	                        if (firstRes.h_bonder == 0) file.WriteLine("strand {0} res {1} has no partner!", strandCtr, firstRes.SeqID);
	                        else file.WriteLine("strand {0} res {1}s partner is {2}; k = {3}", strandCtr, firstRes.SeqID, firstRes.h_bonderID, firstRes.ResStrandNum);
	                    }
	                }
	                return 1;
	            }
	        }

	        public void getShearNum() //7-11-16 MWF; produces an integer for shear number.
	        {
	            int shearNum = 0;
	            int k = 0; int l = 0; int j = 0; int m = 0;
	            int StrandNum1; int StrandNum2; int StrandNum0 = 0;
	            string filelocation = path + "ShearNum/ShearNum_" + PdbName + ".txt";
	            using (System.IO.StreamWriter file = new System.IO.StreamWriter(filelocation))
	            {
	                int res_num = 0;
	                Res first_res = this.Strands[0].Residues[res_num];
	                Res next_res = this.Strands[0].Residues[res_num];
	                Res prev_res = this.Strands[0].Residues[res_num];
	                StrandNum1 = first_res.StrandNum;
	                if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) k = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
	                else k = first_res.ResStrandNum;
	                while (StrandNum1 < this.Strands.Count - 1)
	                {
	                    if (first_res.ShearNumNeigh == null)
	                    {
	                        file.WriteLine("End of seq, restarting");
	                        prev_res = first_res;
	                        m = l;

	                        //If a residue with no neighbors is reached - strand sticks up and need a new residue to restart; go to opposite end of strand and work back towards this residue
	                        bool dir = false;
	                        if (first_res.ResStrandNum > this.Strands[StrandNum1].Count() / 2)
	                        {
	                            res_num = 0;
	                            dir = true;
	                        }
	                        else res_num = this.Strands[StrandNum1].Residues.Count() - 1;

	                        while (first_res.ShearNumNeigh == null)
	                        {
	                            if (dir == true)
	                            {
	                                res_num += 1;
	                                first_res = this.Strands[StrandNum1].Residues[res_num];
	                            }
	                            else
	                            {
	                                res_num -= 1;
	                                first_res = this.Strands[StrandNum1].Residues[res_num];
	                            }
	                        }

	                        StrandNum1 = first_res.StrandNum;
	                        if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
	                        else j = first_res.ResStrandNum;

	                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", prev_res.SeqID, prev_res.StrandNum, m, first_res.SeqID, first_res.StrandNum, j, m - j);
	                        shearNum += j - m;
	                    }

	                    if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
	                    else j = first_res.ResStrandNum;

	                    next_res = first_res.ShearNumNeigh;
	                    StrandNum2 = next_res.StrandNum;
	                    if (this.Strands[StrandNum2].Residues.Last().Z < this.Strands[StrandNum2].Residues[0].Z) l = this.Strands[StrandNum2].Count() - next_res.ResStrandNum;
	                    else l = next_res.ResStrandNum;

	                    file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", first_res.SeqID, first_res.StrandNum, j, next_res.SeqID, next_res.StrandNum, l);
	                    prev_res = first_res;
	                    StrandNum0 = prev_res.StrandNum;
	                    m = j;
	                    first_res = next_res;
	                    StrandNum1 = first_res.StrandNum;

	                }
	                if (StrandNum1 == this.Strands.Count - 1)
	                {
	                    if (first_res.ShearNumNeigh == null)
	                    {
	                        file.WriteLine("End of seq, restarting");
	                        prev_res = first_res;
	                        m = l;

	                        //If a residue with no neighbors is reached - strand sticks up and need a new residue to restart; go to opposite end of strand and work back towards this residue
	                        bool dir = false;
	                        if (first_res.ResStrandNum > this.Strands[StrandNum1].Count() / 2)
	                        {
	                            res_num = 0;
	                            dir = true;
	                        }
	                        else res_num = this.Strands[StrandNum1].Residues.Count() - 1;

	                        while (first_res.ShearNumNeigh == null)
	                        {
	                            if (dir == true)
	                            {
	                                res_num += 1;
	                                first_res = this.Strands[StrandNum1].Residues[res_num];
	                            }
	                            else
	                            {
	                                res_num -= 1;
	                                first_res = this.Strands[StrandNum1].Residues[res_num];
	                            }
	                        }

	                        StrandNum1 = first_res.StrandNum;
	                        if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
	                        else j = first_res.ResStrandNum;

	                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", prev_res.SeqID, prev_res.StrandNum, m, first_res.SeqID, first_res.StrandNum, j, m - j);
	                        shearNum += j - m;
	                    }

	                    if (this.Strands[StrandNum1].Residues.Last().Z < this.Strands[StrandNum1].Residues[0].Z) j = this.Strands[StrandNum1].Count() - first_res.ResStrandNum;
	                    else j = first_res.ResStrandNum;

	                    next_res = first_res.ShearNumNeigh;
	                    StrandNum2 = next_res.StrandNum;
	                    if (this.Strands[StrandNum2].Residues.Last().Z < this.Strands[StrandNum2].Residues[0].Z) l = this.Strands[StrandNum2].Count() - next_res.ResStrandNum;
	                    else l = next_res.ResStrandNum;

	                    file.WriteLine("End of seq");
	                    file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", first_res.SeqID, first_res.StrandNum, j, next_res.SeqID, next_res.StrandNum, l, Math.Abs(l - k));
	                    shearNum += Math.Abs(l - k);
	                    this.ShearNum = shearNum;
	                }


	            }
	        }

	    }

    }

 }
