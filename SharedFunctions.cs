/*
**  File: SharedFunctions.cs
**  Started: 
**  Contributors: Meghan Franklin, Ryan Feehan
**  Overview: 
**
**  About: 
**
**  Last Edited: 
*/


using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using betaBarrelProgram.BarrelStructures;
using betaBarrelProgram.Mono;
using betaBarrelProgram.Poly;
using System.IO;
using System.Windows.Media.Media3D;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace betaBarrelProgram
{
    public class SharedFunctions
    {
        
        public static Vector3D sete1(Atom atomB, Vector3D atomC, Vector3D atomD, double angle, double dihedral)
        {
            Vector3D uCB = -(atomB.Coords - atomC);
            uCB.Normalize();
            Vector3D dDC = -(atomC - atomD);
            double angle2 = 180 - angle;
            double dihe2 = 180 + dihedral;
            double rsin = Math.Sin(angle2);
            double rcos = Math.Cos(angle2);
            double rsinsin = rsin * Math.Sin(dihe2);
            double rsincos = rsin * Math.Cos(dihe2);

            Vector3D dp2 = dDC - (uCB * Vector3D.DotProduct(dDC, uCB));
            dp2.Normalize();
            Vector3D cp1 = Vector3D.CrossProduct(uCB, dDC);
            cp1.Normalize();
            return ((uCB * rcos) + (dp2 * rsincos) + (cp1 * rsinsin) + atomB.Coords);
        }

        public static double calculateEvdw(HAcceptor acceptor, HDonor hAtom, double d)
        {
            double Evdw;
            double DsigmaIJ = d / (acceptor.vDWrad + hAtom.vDWrad);
            double Eij = Math.Sqrt(acceptor.EpsMin * hAtom.EpsMin);

            if (DsigmaIJ <= 0.8254) Evdw = 10;
            else if (DsigmaIJ <= 1) Evdw = 57.273 * (1 - DsigmaIJ);
            else if (DsigmaIJ < 10 / 9) Evdw = Eij * Math.Pow((10 - 9 * DsigmaIJ), 57.273 / (9 * Eij)) - Eij;
            else if (DsigmaIJ < 4 / 3) Evdw = (Eij / 4) * Math.Pow((9 * DsigmaIJ - 10), 2) - Eij;
            else Evdw = 0;

            return Evdw;
        }
        
        public static double CalculateTorsion(Vector3D _atom1Coords, Vector3D _atom2Coords, Vector3D _atom3Coords, Vector3D _atom4Coords)
        {
            Vector3D vec1 = new Vector3D();
            vec1 = _atom2Coords - _atom1Coords;
            Vector3D vec2 = new Vector3D();
            vec2 = _atom3Coords - _atom2Coords;
            Vector3D vec3 = new Vector3D();
            vec3 = _atom4Coords - _atom3Coords;

            Vector3D projVec = new Vector3D(Math.Sqrt(Vector3D.DotProduct(vec2, vec2)) * vec1.X, Math.Sqrt(Vector3D.DotProduct(vec2, vec2)) * vec1.Y,
            Math.Sqrt(Vector3D.DotProduct(vec2, vec2)) * vec1.Z);

            return (Math.Atan2((Vector3D.DotProduct(projVec, Vector3D.CrossProduct(vec2, vec3))),
                Vector3D.DotProduct(Vector3D.CrossProduct(vec1, vec2), Vector3D.CrossProduct(vec2, vec3)))) * 180 / Math.PI;

        }
        
        public static double calculateHydrogenBond(Atom N, Atom O, double d)
        {
            double aMax = 37;
            double BMax = 49;
            double sigmaD = .67;
            double sigmaDSQ = .449;
            double cosAmax = .7986;
            double cosBmax = .561;
            double d0 = 2.08;
            double a = (Math.PI * Vector3D.AngleBetween(O.Coords - N.Hydrogen, N.Hydrogen - N.Coords)) / 180;

            double B = (Math.PI * Vector3D.AngleBetween(O.e1, N.Hydrogen - O.Coords)) / 180;
            if (O.e2!=null)
            {
                double B2 = (Math.PI * Vector3D.AngleBetween(O.e2, N.Hydrogen - O.Coords)) / 180;
                if (B2 < B) B = B2;

            }
            if (B > BMax || a > aMax)
            {
                return 999;
            }
            double denominator = .176;
            double w = Math.Sqrt((sigmaDSQ - Math.Pow(d - d0, 2)) * (Math.Cos(a) - cosAmax) * (Math.Cos(B) - cosBmax)) / denominator;
            return w;
        }

        public static void writePymolScriptForStrands(List<Strand> strandlist, string outputDirectory, string DBdirectory, string pdbName)
        {
            List<string> chain_names = new List<string>();
            string fileLocation = outputDirectory + "strands_" + pdbName + ".py";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                //If you do NOT use MacPymol, uncomment several lines below
                //file.WriteLine("import pymol"); //For PC
                //file.WriteLine("import cmd"); //For PC
                //file.WriteLine("from pymol import stored"); //For PC
                file.WriteLine("from pymol import cmd, stored"); //For MacPyMOL
                string[] colors = { "white", "red", "orange", "purple", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "red", "orange", "yellow", "green", "cyan", "blue", "purple" };
                string pdb_file = DBdirectory + pdbName + ".pdb";
                file.WriteLine("cmd.hide(\"everything\", \"all\")");
                file.WriteLine("cmd.color(\"wheat\",\"all\")");

                foreach (Strand strand in strandlist)
                {
                    file.WriteLine("cmd.select(\"{0}strand{1}\", \"resi {2}-{3} & chain {0} \")", strand.ChainName, strand.StrandNum, strand.Residues[0].SeqID, strand.Residues.Last().SeqID);
                    file.WriteLine("cmd.color (\"{0}\", \"{1}strand{2}\")", colors[strand.StrandNum], strand.ChainName, strand.StrandNum);
                    if (chain_names.Contains(strand.ChainName) == false) chain_names.Add(strand.ChainName);
                    file.WriteLine("\n");
                }

                file.Write("cmd.select(\"barrel\", \"");
                for (int i = 0; i < chain_names.Count; i++)
                {
                    if (i < chain_names.Count - 1) file.Write("{0}strand* or ", chain_names[i]);
                    else file.WriteLine("{0}strand*\")", chain_names[i]);

                }
                file.WriteLine("cmd.show(\"cartoon\", \"barrel\")");
                file.WriteLine("cmd.zoom(\"barrel\")");
            }
        }
		
        public static void writePymolScriptForLoops(Dictionary<string,string> looplist, string outputDirectory, string DBdirectory, ref Chain myChain, string pdbName)
        {
            List<string> chain_names = new List<string>();
            string fileLocation = outputDirectory + "Pymol\\Loops_" + pdbName + ".py";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))
            {
                //If you do NOT use MacPymol, uncomment several lines below
                //file.WriteLine("import pymol"); //For PC
                //file.WriteLine("import cmd"); //For PC
                //file.WriteLine("from pymol import stored"); //For PC
                file.WriteLine("from pymol import cmd, stored"); //For MacPyMOL
                string[] colors = { "white", "red", "orange", "yellow", "forest", "cyan", "blue", "purple", "red", "orange", "yellow", "forest", "cyan", "blue", "purple", "red", "orange", "yellow", "forest", "cyan", "blue", "purple", "red", "orange", "yellow", "forest", "cyan", "blue", "purple", "red", "orange", "yellow", "forest", "cyan", "blue", "purple"};
                string pdb_file = DBdirectory + pdbName + ".pdb";
                file.WriteLine("cmd.load(\"{0}\")", pdb_file);
                file.WriteLine("cmd.hide(\"everything\", \"all\")");
                file.WriteLine("cmd.color(\"grey50\",\"all\")");
                file.WriteLine("cmd.show(\"cartoon\", \"chain {0}\")", myChain.ChainName);

                foreach (KeyValuePair<string,string> loop in looplist)
                {
                    file.WriteLine("cmd.select(\"{0}{1}\", \"resi {2} & chain {0} \")", myChain.ChainName, loop.Key, loop.Value);
                    int loopNum = Convert.ToInt16(loop.Key.Substring(1, loop.Key.Length - 1));
                    file.WriteLine("cmd.color (\"{0}\", \"{2}{1}\")", colors[loopNum], loop.Key, myChain.ChainName);
                    //file.WriteLine("\n");
                }
                file.WriteLine("cmd.center (\"chain {0}\")", myChain.ChainName);
            }
        }

        public static Vector3D getNormal(List<Vector3D> myEllipse, Vector3D myCentroid)
        {
            List<Vector3D> rawPoints = myEllipse;

            //average points to find centroid
            Vector3D centroid = new Vector3D();
            for (int cur = 0; cur < rawPoints.Count; cur++)
            {
                centroid += myEllipse[cur];
            }
            centroid /= rawPoints.Count;

            //move points to around orgin
            for (int cur = 0; cur < rawPoints.Count; cur++) myEllipse[cur] -= centroid;

            Matrix<double> A = DenseMatrix.Build.Random(3, 3);

            double X = myEllipse[0].X;
            double Y = myEllipse[0].Y;
            double Z = myEllipse[0].Z;

            A[0, 0] = (X * X);
            A[0, 1] = (X * Y);
            A[0, 2] = (X * Z);

            A[1, 0] = (Y * X);
            A[1, 1] = (Y * Y);
            A[1, 2] = (Y * Z);

            A[2, 0] = (Z * X);
            A[2, 1] = (Z * Y);
            A[2, 2] = (Z * Z);

            for (int cur = 1; cur < rawPoints.Count; cur++)
            {
                X = myEllipse[cur].X;
                Y = myEllipse[cur].Y;
                Z = myEllipse[cur].Z;

                A[0, 0] = A[0, 0] + (X * X);
                A[0, 1] = A[0, 1] + (X * Y);
                A[0, 2] = A[0, 2] + (X * Z);

                A[1, 0] = A[1, 0] + (Y * X);
                A[1, 1] = A[1, 1] + (Y * Y);
                A[1, 2] = A[1, 2] + (Y * Z);

                A[2, 0] = A[2, 0] + (Z * X);
                A[2, 1] = A[2, 1] + (Z * Y);
                A[2, 2] = A[2, 2] + (Z * Z);
            }


            SortedList<double, Vector3D> myEigenSolution = new SortedList<double, Vector3D>();

            Evd<double> eigen = A.Evd();


            for (int vecCtr = 0; vecCtr < A.RowCount; vecCtr++)
            {
                //double lambda = eigen.EigenValues.At(vecCtr).Real;
                var tmplambda = eigen.EigenValues.At(vecCtr);
                double lambda = tmplambda.Real;

                Vector3D vec3D = new Vector3D(eigen.EigenVectors.At(vecCtr, 0), eigen.EigenVectors.At(vecCtr, 1), eigen.EigenVectors.At(vecCtr, 2));

                Vector3D vec3d = new Vector3D();

                vec3d.X = eigen.EigenVectors.At(0, vecCtr);
                vec3d.Y = eigen.EigenVectors.At(1, vecCtr);
                vec3d.Z = eigen.EigenVectors.At(2, vecCtr);

                myEigenSolution.Add(lambda, vec3d);

            }
            //Vector3D EigenvectorA = new Vector3D();
            //Vector3D EigenvectorB = new Vector3D();
            Vector3D EigenvectorN = new Vector3D();

            double eigenValN = myEigenSolution.Keys.Min();
            EigenvectorN = myEigenSolution[eigenValN];

            return EigenvectorN;

            //Matrix3D covariance = new Matrix3D(sumX2,sumxy,sumxz,0,sumxy,sumY2,sumyz,0, sumxz,sumyz,sumZ2,0,0,0,0,0);

        }

        public static double setRadius(List<Strand> strandlist, Vector3D axis, Vector3D CCentroid, Vector3D NCentroid)
        {
            double AvgRad = 0;
            double totalRes = 0;
            foreach (Strand strand in strandlist)
            {
                foreach (Res myRes in strand)
                {
                    myRes.Radius = (Vector3D.CrossProduct(myRes.BackboneCoords["CA"] - NCentroid, myRes.BackboneCoords["CA"] - CCentroid)).Length / axis.Length;
                    AvgRad += myRes.Radius;
                    totalRes++;
                }
            }
            return AvgRad / totalRes;
        }

        public static void setInOut(List<Strand> strandlist, string outputDirectory, string pdbName, Vector3D axis, Vector3D CCentroid, Vector3D NCentroid)
        {
            double direction; double angleUncertainty; double angle;
                foreach (Strand strand in strandlist)
                {
                    direction = Vector3D.AngleBetween(strand.Residues[0].BackboneCoords["CA"] - strand.Residues[strand.Residues.Count - 1].BackboneCoords["CA"], axis);

                    foreach (Res myRes in strand)
                    {
                        angleUncertainty = 10;
                        angle = Vector3D.AngleBetween(myRes.BackboneCoords["CA"] - ((myRes.BackboneCoords["N"] + myRes.BackboneCoords["C"]) / 2), axis);

                        if ((angle < 90 + angleUncertainty && angle > 90 - angleUncertainty || myRes.ResNum == 0 || myRes.ResNum == strand.Residues.Count - 1) && myRes.ThreeLetCode != "GLY")
                        {
                            //This combines the double-check in/out function
                            bool inward = false;
                            Vector3D myCentroid = new Vector3D();
                            if ((myRes.BackboneCoords["CA"] - NCentroid).Length < (myRes.BackboneCoords["CA"] - CCentroid).Length) myCentroid = NCentroid;
                            else myCentroid = CCentroid;

                            if ((myRes.BackboneCoords["CA"] - myCentroid).Length > (myRes.Atoms[4].Coords - myCentroid).Length) inward = true;

                            if (inward == true)
                            {
                                myRes.Inward = true;
                            }
                            else myRes.Inward = false;
                        }
                        else
                        {
                            if (direction > 90)
                            {
                                if (angle > 90) myRes.Inward = true;
                                else myRes.Inward = false;
                            }
                            else
                            {
                                if (angle < 90) myRes.Inward = true;
                                else myRes.Inward = false;

                            }
                        }
                }
            }
        }

        public static void getInOut(List<Strand> strandlist, string outputDirectory, string pdbName, Vector3D axis, Vector3D CCentroid, Vector3D NCentroid)
        {
            string fileLocation3 = outputDirectory + "InOut\\InOut_" + pdbName + ".txt";
            double direction; double angleUncertainty; double angle;

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation3))
            {
                file.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n", "Res", "SeqID", "Strand", "Chain", "Inward", "ZCoord");

                foreach (Strand strand in strandlist)
                {
                    direction = Vector3D.AngleBetween(strand.Residues[0].BackboneCoords["CA"] - strand.Residues[strand.Residues.Count - 1].BackboneCoords["CA"], axis);

                    foreach (Res myRes in strand)
                    {
                        angleUncertainty = 10;
                        angle = Vector3D.AngleBetween(myRes.BackboneCoords["CA"] - ((myRes.BackboneCoords["N"] + myRes.BackboneCoords["C"]) / 2), axis);

                        if ((angle < 90 + angleUncertainty && angle > 90 - angleUncertainty || myRes.ResNum == 0 || myRes.ResNum == strand.Residues.Count - 1) && myRes.ThreeLetCode != "GLY")
                        {
                            //This combines the double-check in/out function
                            bool inward = false;
                            Vector3D myCentroid = new Vector3D();
                            if ((myRes.BackboneCoords["CA"] - NCentroid).Length < (myRes.BackboneCoords["CA"] - CCentroid).Length) myCentroid = NCentroid;
                            else myCentroid = CCentroid;

                            if ((myRes.BackboneCoords["CA"] - myCentroid).Length > (myRes.Atoms[4].Coords - myCentroid).Length) inward = true;

                            if (inward == true)
                            {
                                myRes.Inward = true;
                            }
                            else myRes.Inward = false;
                        }
                        else
                        {
                            if (direction > 90)
                            {
                                if (angle > 90) myRes.Inward = true;
                                else myRes.Inward = false;
                            }
                            else
                            {
                                if (angle < 90) myRes.Inward = true;
                                else myRes.Inward = false;

                            }
                        }

                        file.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n", myRes.ThreeLetCode, myRes.SeqID, myRes.StrandNum, myRes.ChainName, myRes.Inward, myRes.BackboneCoords["CA"].Z);
                        //file.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n", myRes.ThreeLetCode, myRes.ResNum + 1, myRes.StrandNum, myRes.ChainName, myRes.Inward, myRes.BackboneCoords["CA"].Z);
                    }
                }
            }
        }

        public static List<double> getStrandLengths(List<Strand> strandlist, string outputDirectory, string pdbName)
        {
            string fileLocation3 = outputDirectory + "Tilts\\StrandLengths_" + pdbName + ".txt";
            double height; List<double> all_lengths = new List<double>();
            string fileOfPDBs = Global.MACMONODBDIR;
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation3))
            {
                file.Write("{0}\t{1}\t{2}\n", "Chain", "Strand #", "Length");

                foreach (Strand strand in strandlist)
                {
                    Res res1 = strand.Residues[0];
                    Res res2 = strand.Residues.Last();
                    height = Math.Abs(res1.BackboneCoords["CA"].Z - res2.BackboneCoords["CA"].Z);
                    all_lengths.Add(height);
                    file.Write("{0}\t{1}\t{2}\n", strand.ChainName, strand.StrandNum, height);
                }
            }
            
            return all_lengths;
        }
		
        public static Dictionary<string, string> getLoopTurns(List<Strand> strandlist, ref Chain myChain, string outputDirectory, string pdbName)
        {
            //string fileLocation = outputDirectory + "RosettaLoops\\Loops\\" + pdbName + "_Loops_Test.txt";
            //string fileLocation2 = outputDirectory + "RosettaLoops\\Turns\\" + pdbName + "_Turns_Test.txt";
            string fileLocation = outputDirectory + "LoopData\\v4Loops\\" + pdbName + "_Loops_Hairpin.txt";
            string fileLocation2 = outputDirectory + "TurnData\\v4Turns\\" + pdbName + "_Turns_Hairpin.txt";
            //string fileLocation = outputDirectory + "Loops\\" + pdbName + "_Loops.txt";
            //string fileLocation2 = outputDirectory + "Turns\\" + pdbName + "_Turns.txt";
            int turn_count = 0; int loop_count = 0; bool is_turn; string loop_seq = ""; List<Tuple<double, double>> phipsi = new List<Tuple<double, double>>();
            string newline;
            Dictionary<string, string> all_loops = new Dictionary<string, string>();

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))
            {
                file.Write("PDB\tID\tLength\tSequence\tRes\tDist\tStrand1\tStrand2\n");
                 using (System.IO.StreamWriter file2 = new System.IO.StreamWriter(fileLocation2, true))
                 {
                     file2.Write("PDB\tID\tLength\tSequence\tRes\tDist\tStrand1\tStrand2\n");
                     
                     foreach (Strand strand in strandlist)
                     {
                         if (strand.ChainName == myChain.ChainName && strand.StrandNum != strandlist.Count - 1 && strandlist.IndexOf(strand) + 1 < strandlist.Count && strandlist[strandlist.IndexOf(strand) + 1].ChainName == myChain.ChainName)
                         {
                             int strand_id = strandlist.IndexOf(strand);
                             Res first_res = strand.Residues.Last();
                             Res last_res = strandlist[strand_id + 1].Residues[0];
                             //Determine if loop/turn; turns are on inside so Z-coords will be negative
                             if (first_res.BackboneCoords["CA"].Z - strand.Residues[0].BackboneCoords["CA"].Z < 0) { turn_count += 1; is_turn = true; }
                             else { loop_count += 1; is_turn = false; }
                             loop_seq = "";
                             phipsi.Clear();
                             List<double> avgBFac = new List<double>();
                             double distance = (first_res.Atoms[1].Coords - last_res.Atoms[1].Coords).Length;
                             //Check to see if loop is contiguous between strands
                             if (last_res.ResNum - first_res.ResNum == last_res.SeqID - first_res.SeqID)
                             {
                                 for (int loop_res = first_res.ResNum; loop_res < last_res.ResNum + 1; loop_res++)
                                 {
                                     loop_seq += myChain.Residues[loop_res].OneLetCode;
                                     phipsi.Add(Tuple.Create(myChain.Residues[loop_res].Phi, myChain.Residues[loop_res].Psi));
                                     avgBFac.Add(myChain.Residues[loop_res].RelBFac);
                                 }
                                 //Console.WriteLine(loop_seq);
                                 //Console.WriteLine(phipsi);
                                 newline = loop_seq.Length.ToString() + "\t" + loop_seq + "\t" + first_res.SeqID.ToString() + "-" + last_res.SeqID.ToString() + "\t" + distance.ToString() + "\t" + strand_id.ToString() + "\t" + (strand_id + 1).ToString() + "\t";
                                //newline = loop_seq.Length.ToString() + "\t" + loop_seq + "\t" + (first_res.ResNum + 1).ToString() + "-" + (last_res.ResNum + 1).ToString() + "\t";

                                  //if (loop_seq.Length > 0) { newline += avgBFac.Average().ToString("F4") + "\t"; }
                                 
                                 //Only print the phi/psi angles if less than 16 residues
                                 /*if (loop_seq.Length < 16 && loop_seq.Length > 0)
                                 {
                                     foreach (Tuple<double, double> x in phipsi) { newline += "(" + x.Item1.ToString("F3") + "," + x.Item2.ToString("F3") +"), " ; }
                                 }*/
                                 
                             }
                             else
                             {
                                 if (last_res.ResNum < first_res.ResNum)
                                 {
                                     Console.WriteLine("INCORRECT LOOP DEFINITION");
                                     if (is_turn == true) { turn_count -= 1; is_turn = false; loop_count += 1; }
                                     else { loop_count -= 1; is_turn = true; turn_count += 1; }
                                 }
                                 //incomplete loop
                                 //newline = "INCOMPLETE\t" + first_res.SeqID.ToString() + "-" + last_res.SeqID.ToString();
                                 newline = "INCOMPLETE\t" + (first_res.ResNum + 1).ToString() + "-" + (last_res.ResNum + 1).ToString() + "\t" + distance.ToString() + "\t" + strand_id.ToString() + "\t" + (strand_id + 1).ToString();
                             }

                             //Write data to file and to dictionary
                             if (is_turn == true) 
                             { 
                                 file2.WriteLine("{2}\tT{0}\t{1}", turn_count, newline, pdbName);
                                 string key = "T" + turn_count.ToString();
                                 string value = first_res.SeqID.ToString() + "-" + last_res.SeqID.ToString();
                                 all_loops.Add(key, value); 
                             }
                             else 
                             {
                                 file.WriteLine("{2}\tL{0}\t{1}", loop_count, newline, pdbName);
                                 string key = "L" + loop_count.ToString();
                                 string value = first_res.SeqID.ToString() + "-" + last_res.SeqID.ToString();
                                 all_loops.Add(key, value);
                             }
                         }
                     }

                 }
            }
            return all_loops;
        }

        public static void findLoopsHBondingPartnersGeomOnly(Dictionary<string, string> looplist, string outputDirectory, ref Chain myChain, string pdbName, bool polybarrel)
        {
            string newLine;
            bool inLoop;
            string fileLocation = outputDirectory + "HBonding\\GeomSC" + pdbName + ".txt";
            string fileLocation2 = outputDirectory + "HBonding\\GeomBBone" + pdbName + ".txt";
            string fileLocation3 = outputDirectory + "HBonding\\GeomBBSC" + pdbName + ".txt";
            Char delimiter = '-';

            using (System.IO.StreamWriter test_file = new System.IO.StreamWriter(outputDirectory + "HBonding\\TestFile"+pdbName + ".txt"))
            {
                using (System.IO.StreamWriter SCfile = new System.IO.StreamWriter(fileLocation))
                {
                    newLine = "Atom\tNum\tRes\tChain\tAtom2\tNum\tRes2\tChain\tDist\tLoopID\tPartInLoop\tLoopPos";
                    SCfile.WriteLine(newLine);
                    using (System.IO.StreamWriter backbone = new System.IO.StreamWriter(fileLocation2))
                    {
                        backbone.WriteLine(newLine);
                        using (System.IO.StreamWriter bboneSC = new System.IO.StreamWriter(fileLocation3))
                        {
                            bboneSC.WriteLine(newLine);

                            foreach (KeyValuePair<string, string> loop in looplist)
                            {
                                var loop_ID = loop.Key;
                                int loop_start = Convert.ToInt16(loop.Value.Split(delimiter)[0]);
                                int loop_end = Convert.ToInt16(loop.Value.Split(delimiter)[1]);
                                //Console.WriteLine("{0},{1},{2}", loop_ID, loop_start, loop_end);

                                loop_start = myChain.Residues.FindIndex(a => a.SeqID == loop_start);
                                loop_end = myChain.Residues.FindIndex(a => a.SeqID == loop_end);
                                //Console.WriteLine("{0},{1},{2}", loop_ID, loop_start, loop_end);

                                int loop_pos = -1;
                                for (int i = loop_start; i < loop_end + 1; i++)
                                {
                                    Res Residue1 = myChain.Residues[i];
                                    loop_pos++;

                                    for (int j = loop_start - 3; j < loop_end + 4; j++)
                                    {
                                        Res Residue2 = myChain.Residues[j];
                                        if (j >= loop_start && j <= loop_end) inLoop = true;
                                        else inLoop = false;

                                        foreach (Atom atom1 in Residue1)
                                        {
                                            if ((atom1.AtomName == "N" && Residue1.ThreeLetCode != "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode != "GLY") || atom1.AtomName == "C" || atom1.AtomType == "H") continue;
                                            else if (atom1.AtomType == "C" && atom1.AtomName != "CA") continue; //skip any carbon atoms that aren't CA of gly
                                            else //Compare to chains on either side
                                            {

                                                checkLoopsGeom(Residue1, atom1, Residue2, loop_ID.ToString(), inLoop, loop_pos, backbone, bboneSC, SCfile, test_file);
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

        public static void checkLoopsGeom(Res Residue1, Atom atom1, Res Residue2, String LoopID, bool InLoop, int loop_position, StreamWriter backbone, StreamWriter bboneSC, StreamWriter SCfile, StreamWriter test_file)
        {
            double dist; double angle1; double angle2;
            Atom prevatom1 = atom1; Point3D prevatom1_coords;
            Atom prevatom2; Point3D prevatom2_coords;
            string newLine;
            List<string> bbatoms = new List<string> { "O", "C", "N", "CA", "OXT" };

            foreach (Atom atom2 in Residue2)
            {
                if (atom2.AtomType == "H" || atom2.AtomName == "C") continue; //Skip all hydrogen; C can't h-bond
                else if (atom2.AtomName == atom1.AtomName && Residue1.ResNum == Residue2.ResNum || atom2.AtomName == "CA" && Residue1.ResNum == Residue2.ResNum) continue; //skip self-atom in loop and self-CA bonds as impossible
                else if (bbatoms.Contains(atom2.AtomName) && atom1.AtomName == "O") //O to bbatom
                {
                    if (Residue1.ResNum == Residue2.ResNum || atom2.AtomName == "N" && Residue2.ResNum == Residue1.ResNum + 1 ) continue; //can't self-bb-bond or O-next N
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;

                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "CA")]; //N-CA, CA-N, O-C
                        if (atom2.AtomName == "CA") prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "N")];
                        if (atom2.AtomName == "O") prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = 180- Vector3D.AngleBetween(interatom, intraatom);
                        test_file.WriteLine("{0}_{1} {2} {3}_{4} {5} {6} {7}", atom1.AtomName, atom1.ResSeqID, prevatom1.AtomName, atom2.AtomName, atom2.ResSeqID, prevatom2.AtomName, angle1, angle2);

                        if (angle1 > 90 && angle2 > 90)
                        {
                            if (atom1.BBNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.BBNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, Residue2.ChainName, dist, LoopID, InLoop, loop_position);
                                backbone.WriteLine(newLine);
                            }
                            if (atom2.BBNeighAtoms.Contains(atom1) == false) atom2.BBNeighAtoms.Add(atom2);
                        }
                    }
                }
                else if ((atom1.AtomName == "N" && Residue1.ThreeLetCode == "PRO" && atom2.AtomName == "N" && Residue2.ThreeLetCode == "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode == "GLY" && atom2.AtomName == "CA" && Residue2.ThreeLetCode == "GLY")) //check Gly-Gly or Pro-Pro
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
                        for (var i = Residue2.Atoms.FindIndex(a => a.AtomName == atom2.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue2.Atoms[i].AtomType == "C")
                            {
                                prevatom2 = Residue2.Atoms[i];
                                prevatom2_coords = (Point3D)Residue2.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = 180 - Vector3D.AngleBetween(interatom, intraatom);
                        if (angle1 > 90 && angle2 > 90)
                        {
                            if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCSCNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, Residue2.ChainName, dist, LoopID, InLoop, loop_position);
                                SCfile.WriteLine(newLine);
                            }
                            if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                        }
                    }
                }

                else if (bbatoms.Contains(atom2.AtomName) && bbatoms.Contains(atom1.AtomName) == false) //SC atom to bb atom
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
						//These lines set up in case can't find a prev atom, ie CA
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
						//Work backwards from atom1/atom2 to find closest C as this is the atom to which it is bonded previously in the residue
						//ALso need to ID the atom before this
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
						//Bbatoms have specific geometries
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "CA")]; //N-CA, CA-N, O-C
                        if (atom2.AtomName == "CA") prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "N")];
                        if (atom2.AtomName == "O") prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = 180 - Vector3D.AngleBetween(interatom, intraatom); //180 is because intraatom points towards towards atom2 instead of away from it
                        //test_file.WriteLine("{0}_{1} {2} {3}_{4} {5} {6} {7}", atom1.AtomName, atom1.ResSeqID, prevatom1.AtomName, atom2.AtomName, atom2.ResSeqID, prevatom2.AtomName, angle1, angle2);
                        if (angle1 > 90 && angle2 > 90) //changed to && from || 2/22/18
                        {
                            if (atom1.SCBBNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCBBNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, Residue2.ChainName, dist, LoopID, InLoop, loop_position);
                                bboneSC.WriteLine(newLine);
                            }
                            if (atom2.SCBBNeighAtoms.Contains(atom1) == false) atom2.SCBBNeighAtoms.Add(atom1);
                        }
                    }
                }
                else if (bbatoms.Contains(atom1.AtomName) || bbatoms.Contains(atom1.AtomName)) continue;
                else
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
                        for (var i = Residue2.Atoms.FindIndex(a => a.AtomName == atom2.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue2.Atoms[i].AtomType == "C")
                            {
                                prevatom2 = Residue2.Atoms[i];
                                prevatom2_coords = (Point3D)Residue2.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = 180 - Vector3D.AngleBetween(interatom, intraatom);

                        if (angle1 > 90 && angle2 > 90)
                        {
                            if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCSCNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, Residue2.ChainName, dist, LoopID, InLoop, loop_position);
                                SCfile.WriteLine(newLine);
                            }
                            if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                        }

                    }
                }
            }
        }
		

        public static void findNearestNeighbors(List<Strand> strandlist, string outputDirectory, string pdbName)
        {
            int next_strand; Res neighbor; double dist; double neighDist;
            string fileLocation = outputDirectory + "\\NearestNeighMono\\CANeighbors\\NearestNeigh" + pdbName + ".txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                file.WriteLine("Res\tResNum\tStrand\tNearNeigh\tResNum\tStrand\tDistance Between\tZCoord");
                foreach (Strand strand in strandlist)
                {
                    int strandIndex = strandlist.IndexOf(strand);
                    if (strandIndex == strandlist.Count - 1) next_strand = 0;
                    else next_strand = strandIndex + 1;
                    foreach (Res Residue1 in strand)
                    {
                        neighbor = strandlist[next_strand].Residues[0]; //obviously, this is the neighbor furthest away b/c strands are antiparallel
                        foreach (Res Residue2 in strandlist[next_strand]) //potential neighbors are only sought out on next strand because barrels are circular and this prevents duplicates in reversed order in the output file
                        {
                            dist = (Residue1.BackboneCoords["CA"] - Residue2.BackboneCoords["CA"]).Length;
                            neighDist = (Residue1.BackboneCoords["CA"] - neighbor.BackboneCoords["CA"]).Length;
                            if (dist < neighDist)
                            {
                                neighbor = Residue2;
                            }
                        }
                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t\t{7}", Residue1.ThreeLetCode, Residue1.SeqID, strand.StrandNum, neighbor.ThreeLetCode, neighbor.SeqID, strandlist[next_strand].StrandNum, (Residue1.BackboneCoords["CA"] - neighbor.BackboneCoords["CA"]).Length, Residue1.BackboneCoords["CA"].Z);
                        if ((Residue1.BackboneCoords["CA"] - neighbor.BackboneCoords["CA"]).Length < 6) Residue1.ShearNumNeigh = neighbor;
                    }
                }
            }
        }
    
        //This is the current implementation for partners
        public static void findHBondingPartnersGeomOnly(List<Strand> strandlist, string outputDirectory, string pdbName)
        {
            string newLine; bool IFstatus;
            string fileLocation = outputDirectory + "HBonding\\GeomSC" + pdbName + ".txt";
            string fileLocation2 = outputDirectory + "HBonding\\GeomBBone" + pdbName + ".txt";
            string fileLocation3 = outputDirectory + "HBonding\\GeomBBSC" + pdbName + ".txt";

            using (System.IO.StreamWriter SCfile = new System.IO.StreamWriter(fileLocation))
            {
                newLine = "Atom" + "\t" + "Num" + "\t" + "Res" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Atom2" + "\t" + "Num" + "\t" + "Res2" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Dist";
                SCfile.WriteLine(newLine);
                using (System.IO.StreamWriter backbone = new System.IO.StreamWriter(fileLocation2))
                {
                    backbone.WriteLine(newLine);

                    using (System.IO.StreamWriter bboneSC = new System.IO.StreamWriter(fileLocation3))
                    {
                        bboneSC.WriteLine(newLine);

                        foreach (Strand strand in strandlist)
                        {
                            int strandIndex = strandlist.IndexOf(strand);

                            foreach (Res Residue1 in strand)
                            {
                                foreach (Atom atom1 in Residue1)
                                {
                                    if ((atom1.AtomName == "N" && Residue1.ThreeLetCode != "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode != "GLY") || atom1.AtomName == "C" || atom1.AtomType == "H") continue;
                                    else //Compare to chains on either side
                                    {
                                        int strand_check = 0;
                                        while (strand_check < strandlist.Count())
                                        {
                                            ////Console.WriteLine(strand_check);
                                            if (strand_check != strandIndex)
                                            {
                                                foreach (Res Residue2 in strandlist[strand_check])
                                                {
                                                    if ((strandIndex == 0 && (strand_check == strandlist.Count() - 1 || strand_check == strandlist.Count() - 2)) || (strandIndex == 1 && strand_check == strandlist.Count() - 1)) IFstatus = true;
                                                    else if ((strandIndex == strandlist.Count() - 1 && (strand_check == 0 || strand_check == 1)) || (strandIndex == strandlist.Count() - 2 && strand_check == 0)) IFstatus = true;
                                                    else IFstatus = false;

                                                    checkStrandsGeom(Residue1, atom1, Residue2, strand, strandlist, strand_check, IFstatus, backbone, bboneSC, SCfile);
                                                }
                                            }
                                            strand_check += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        public static void checkStrandsGeom(Res Residue1, Atom atom1, Res Residue2, Strand strand, List<Strand> strandlist, int strand_num_compare, bool IFstatus, StreamWriter backbone, StreamWriter bboneSC, StreamWriter SCfile)
        {
            double dist; double angle1; double angle2;
            Atom prevatom1 = atom1; Point3D prevatom1_coords;
            Atom prevatom2; Point3D prevatom2_coords;
            string newLine;
            List<string> bbatoms = new List<string> { "O", "C", "N", "CA", "OXT" };

            foreach (Atom atom2 in Residue2)
            {
                if (atom2.AtomType == "H") continue; //Skip all hydrogen
                else if (bbatoms.Contains(atom2.AtomName) && atom1.AtomName == "O")
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        if (atom1.BBNeighAtoms.Contains(atom2) == false)
                        {
                            atom1.BBNeighAtoms.Add(atom2);
                            newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                            backbone.WriteLine(newLine);
                        }
                        if (atom2.BBNeighAtoms.Contains(atom1) == false) atom2.BBNeighAtoms.Add(atom2);
                    }
                }
                else if ((atom1.AtomName == "N" && Residue1.ThreeLetCode == "PRO" && atom2.AtomName == "N" && Residue2.ThreeLetCode == "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode == "GLY" && atom2.AtomName == "CA" && Residue2.ThreeLetCode == "GLY")) //bbatoms.Contains(atom2.AtomName) = false
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
                        for (var i = Residue2.Atoms.FindIndex(a => a.AtomName == atom2.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue2.Atoms[i].AtomType == "C")
                            {
                                prevatom2 = Residue2.Atoms[i];
                                prevatom2_coords = (Point3D)Residue2.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = Vector3D.AngleBetween(interatom, intraatom);
                        if (angle1 > 90 || angle2 > 90)
                        {
                            if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCSCNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                                SCfile.WriteLine(newLine);
                            }
                            if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                        }
                    }
                }

                else if (bbatoms.Contains(atom2.AtomName) && bbatoms.Contains(atom1.AtomName) == false) //(atom1.AtomName == "N" && Residue1.ThreeLetCode == "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode == "GLY") ||
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
                        for (var i = Residue2.Atoms.FindIndex(a => a.AtomName == atom2.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue2.Atoms[i].AtomType == "C")
                            {
                                prevatom2 = Residue2.Atoms[i];
                                prevatom2_coords = (Point3D)Residue2.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = Vector3D.AngleBetween(interatom, intraatom);
                        if (angle1 > 90 || angle2 > 90)
                        {
                            if (atom1.SCBBNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCBBNeighAtoms.Add(atom2);
                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                                bboneSC.WriteLine(newLine);
                            }
                            if (atom2.SCBBNeighAtoms.Contains(atom1) == false) atom2.SCBBNeighAtoms.Add(atom1);
                        }
                    }
                }
                else if (bbatoms.Contains(atom1.AtomName) || bbatoms.Contains(atom1.AtomName)) continue;
                else
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        prevatom1 = Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom1_coords = (Point3D)prevatom1.Coords;
                        prevatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == "C")];
                        prevatom2_coords = (Point3D)prevatom2.Coords;
                        for (var i = Residue1.Atoms.FindIndex(a => a.AtomName == atom1.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue1.Atoms[i].AtomType == "C")
                            {
                                prevatom1 = Residue1.Atoms[i];
                                prevatom1_coords = (Point3D)Residue1.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }
                        for (var i = Residue2.Atoms.FindIndex(a => a.AtomName == atom2.AtomName) - 1; i > 0; i--)
                        {
                            if (Residue2.Atoms[i].AtomType == "C")
                            {
                                prevatom2 = Residue2.Atoms[i];
                                prevatom2_coords = (Point3D)Residue2.Atoms[i].Coords; i = 0;
                            }
                            else continue;
                        }

                        Vector3D interatom = ((Point3D)atom1.Coords) - prevatom1_coords;
                        Vector3D intraatom = atom1.Coords - atom2.Coords;
                        angle1 = Vector3D.AngleBetween(interatom, intraatom);
                        interatom = ((Point3D)atom2.Coords) - prevatom2_coords;
                        angle2 = Vector3D.AngleBetween(interatom, intraatom);

                        if (angle1 > 90 || angle2 > 90)
                        {
                            if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                            {
                                atom1.SCSCNeighAtoms.Add(atom2);

                                newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                                SCfile.WriteLine(newLine);
                            }
                            if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                        }

                    }
                }
            }
        }

        public static void findNearNeighDistOnly(List<Strand> strandlist, string outputDirectory, string pdbName)
        {
            //Neighbor residues are added to Strand residues 

            string newLine;
            int next_strand; int prev_strand; int next_next_strand; int prev_prev_strand;
            bool IFstatus = false;
            List<string> bbatoms = new List<string> { "O", "C", "N", "CA" };

            string fileLocation = outputDirectory + "NearestNeigh/NNSideChains4strand" + pdbName + ".txt";
            string fileLocation2 = outputDirectory + "NearestNeigh/NNBackBone4strand" + pdbName + ".txt";
            string fileLocation3 = outputDirectory + "NearestNeigh/NNBBSC4strand" + pdbName + ".txt";
            using (System.IO.StreamWriter SCfile = new System.IO.StreamWriter(fileLocation))
            {
                newLine = "Atom" + "\t" + "Num" + "\t" + "Res" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Atom2" + "\t" + "Num" + "\t" + "Res2" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Dist" + "\t" + "IF" + "\t" + "Res1 In/Out" + "\t" + "Res2 In/Out";
                SCfile.WriteLine(newLine);
                using (System.IO.StreamWriter backbone = new System.IO.StreamWriter(fileLocation2))
                {
                    backbone.WriteLine(newLine);
                    using (System.IO.StreamWriter bboneSC = new System.IO.StreamWriter(fileLocation3))
                    {
                        bboneSC.WriteLine(newLine);

                        foreach (Strand strand in strandlist)
                        {
                            int strandIndex = strandlist.IndexOf(strand);

                            if (strandIndex == strandlist.Count - 1)
                            {
                                next_strand = 0;
                                next_next_strand = 1;
                                prev_strand = strandIndex - 1;
                                prev_prev_strand = strandIndex - 2;
                            }
                            else if (strandIndex == strandlist.Count - 2)
                            {
                                next_strand = strandlist.Count - 1;
                                next_next_strand = 0;
                                prev_strand = strandIndex - 1;
                                prev_prev_strand = strandIndex - 2;
                            }
                            else if (strandIndex == 0)
                            {
                                next_strand = 1;
                                next_next_strand = 2;
                                prev_strand = strandlist.Count - 1;
                                prev_prev_strand = strandlist.Count - 2;
                            }
                            else if (strandIndex == 1)
                            {
                                next_strand = 2;
                                next_next_strand = 3;
                                prev_strand = 0;
                                prev_prev_strand = strandlist.Count - 1;
                            }

                            else
                            {
                                next_strand = strandIndex + 1;
                                next_next_strand = strandIndex + 2;
                                prev_strand = strandIndex - 1;
                                prev_prev_strand = strandIndex - 2;
                            }

                            foreach (Res Residue1 in strand)
                            {
                                foreach (Atom atom1 in Residue1)
                                {
                                    if ((atom1.AtomName == "N" && Residue1.ThreeLetCode != "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode != "GLY") || atom1.AtomName == "C" || atom1.AtomType == "H") continue;
                                    else //Compare to chains on either side
                                    {
                                        foreach (Res Residue2 in strandlist[prev_strand])
                                        {
                                            if (strandIndex == 0) IFstatus = true;
                                            else IFstatus = false;

                                            checkNeighborStrand(Residue1, atom1, Residue2, strand, strandlist, prev_strand, IFstatus, backbone, bboneSC, SCfile);

                                        }

                                        foreach (Res Residue3 in strandlist[next_strand])
                                        {
                                            if (strandIndex == strandlist.Count - 1) IFstatus = true;
                                            else IFstatus = false;

                                            checkNeighborStrand(Residue1, atom1, Residue3, strand, strandlist, next_strand, IFstatus, backbone, bboneSC, SCfile);

                                        }

                                        foreach (Res Residue4 in strandlist[prev_prev_strand])
                                        {
                                            if (strandIndex == 0 || strandIndex == 1) IFstatus = true;
                                            else IFstatus = false;

                                            checkNeighborStrand(Residue1, atom1, Residue4, strand, strandlist, prev_prev_strand, IFstatus, backbone, bboneSC, SCfile);
                                        }

                                        foreach (Res Residue5 in strandlist[next_next_strand])
                                        {
                                            if (strandIndex == strandlist.Count - 1 || strandIndex == strandlist.Count - 2) IFstatus = true;
                                            else IFstatus = false;

                                            checkNeighborStrand(Residue1, atom1, Residue5, strand, strandlist, next_next_strand, IFstatus, backbone, bboneSC, SCfile);

                                        }
                                    }
                                }
                            }
                        }
                        //SCfile.WriteLine();
                    }
                }
            }
        }

        public static void checkNeighborStrand(Res Residue1, Atom atom1, Res Residue2, Strand strand, List<Strand> strandlist, int strand_num_compare, bool IFstatus, StreamWriter backbone, StreamWriter bboneSC, StreamWriter SCfile)
        {
            double dist;
            string newLine;
            List<string> bbatoms = new List<string> { "O", "C", "N", "CA", "OXT" };

            foreach (Atom atom2 in Residue2)
            {
                if (atom2.AtomType == "H") continue; //Skip all hydrogen
                else if (bbatoms.Contains(atom2.AtomName) && atom1.AtomName == "O")
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 4)
                    {
                        if (atom1.BBNeighAtoms.Contains(atom2) == false)
                        {
                            atom1.BBNeighAtoms.Add(atom2);
                            newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                            backbone.WriteLine(newLine);
                        }
                        if (atom2.BBNeighAtoms.Contains(atom1) == false) atom2.BBNeighAtoms.Add(atom2);
                    }
                }
                else if (bbatoms.Contains(atom2.AtomName) == false && ((atom1.AtomName == "N" && Residue1.ThreeLetCode == "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode == "GLY")))
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                        {
                            atom1.SCSCNeighAtoms.Add(atom2);
                            newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                            SCfile.WriteLine(newLine);
                        }
                        if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                    }
                }

                else if (bbatoms.Contains(atom2.AtomName) && ((atom1.AtomName == "N" && Residue1.ThreeLetCode == "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode == "GLY") || bbatoms.Contains(atom1.AtomName) == false))
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        if (atom1.SCBBNeighAtoms.Contains(atom2) == false)
                        {
                            atom1.SCBBNeighAtoms.Add(atom2);
                            newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                            bboneSC.WriteLine(newLine);
                        }
                        if (atom2.SCBBNeighAtoms.Contains(atom1) == false) atom2.SCBBNeighAtoms.Add(atom1);
                    }
                }
                else if (bbatoms.Contains(atom1.AtomName) || bbatoms.Contains(atom1.AtomName)) continue;
                else
                {
                    dist = (atom1.Coords - atom2.Coords).Length;
                    if (dist < 6)
                    {
                        if (atom1.SCSCNeighAtoms.Contains(atom2) == false)
                        {
                            atom1.SCSCNeighAtoms.Add(atom2);

                            newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[strand_num_compare].StrandNum, Residue2.ChainName, dist, IFstatus, Residue1.Inward, Residue2.Inward);
                            SCfile.WriteLine(newLine);
                        }
                        if (atom2.SCSCNeighAtoms.Contains(atom1) == false) atom2.SCSCNeighAtoms.Add(atom1);
                    }
                }
            }
        }

        public static void findHBondingPartnersEnergy(List<Strand> strandlist, string outputDirectory, string pdbName)
        {
            //Neighbor residues are added to Strand residues 
            Dictionary<Tuple<string, string>, HAcceptor> HAcceptorList = new Dictionary<Tuple<string, string>, HAcceptor>();
            Dictionary<Tuple<string, string>, HDonor> HDonorList = new Dictionary<Tuple<string, string>, HDonor>();

            if (File.Exists(Global.parameterFile))
            {
                using (StreamReader sr = new StreamReader(Global.parameterFile))
                {
                    String line;
                    string resName = "ALA";
                    string atomName;
                    // Read and display lines from the file until the end of the file (null) is reached.
                    while ((line = sr.ReadLine()) != null)
                    {
                        ////Console.WriteLine(line);
                        string[] spLine = Array.FindAll<string>(((string)line).Split(
                        new char[] { ' ', '\t', ',', ';' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                        List<string> splitLine = new List<string>(spLine);
                        if (splitLine.Count == 0 || splitLine[0].Contains("!")) continue;
                        else
                        {
                            if (splitLine[0].Contains("RESI"))
                            {
                                resName = splitLine[1];
                                if (resName == "HSE") resName = "HIS";
                            }
                            else if (splitLine[1] == "O" || splitLine[1] == "HN") continue;
                            else if (splitLine[0].Contains("DONOR"))
                            {
                                atomName = splitLine[1];
                                HDonor newDonor = new HDonor(splitLine);
                                HDonorList.Add(new Tuple<string, string>(resName, atomName), newDonor);
                            }
                            else if (splitLine[0].Contains("ACCEPTOR"))
                            {
                                atomName = splitLine[1];
                                HAcceptor acceptor = new HAcceptor(splitLine);
                                HAcceptorList.Add(new Tuple<string, string>(resName, atomName), acceptor);
                            }
                        }
                    }
                }
            }

            double dist;
            string newLine;
            int next_strand; int prev_strand;
            List<string> saltBridgeAA = new List<string> { "ARG", "LYS", "ASP", "GLU" };
            List<string> saltBridgeAtoms = new List<string> { "O1-", "N1+" };

            //This set of variables is needed for the energy calculation
            double B = 35; double d0 = 2.08; double sigmaD = 0.67; double sigmaD2 = sigmaD * sigmaD; double alphaMax = 37; double betaMax = 49;
            double d; double hBondE;
            double cosAM = Math.Cos(alphaMax);
            double cosBM = Math.Cos(betaMax);
            double wdenom = sigmaD * Math.Sqrt((1 - cosAM) * (1 - cosBM));


            string fileLocation = outputDirectory + "HBonding/SCIntEnergy" + pdbName + ".txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                newLine = "Atom" + "\t" + "Num" + "\t" + "Res" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Atom2" + "\t" + "Num" + "\t" + "Res2" + "\t" + "Strand" + "\t" + "Chain" + "\t" + "Energy or SB";
                file.WriteLine(newLine);

                foreach (Strand strand in strandlist)
                {
                    int strandIndex = strandlist.IndexOf(strand);
                    if (strandIndex == strandlist.Count - 1)
                    {
                        next_strand = 0;
                        prev_strand = strandIndex - 1;
                    }
                    else if (strandIndex == 0)
                    {
                        next_strand = 1;
                        prev_strand = strandlist.Count - 1;
                    }
                    else
                    {
                        next_strand = strandIndex + 1;
                        prev_strand = strandIndex - 1;
                    }

                    foreach (Res Residue1 in strand)
                    {
                        foreach (Atom atom1 in Residue1)
                        {
                            if ((atom1.AtomName == "N" && Residue1.ThreeLetCode != "PRO") || (atom1.AtomName == "CA" && Residue1.ThreeLetCode != "GLY") || atom1.AtomName == "C" || atom1.AtomName == "H") continue;
                            else //Compare to chains on either side
                            {
                                foreach (Res Residue2 in strandlist[prev_strand])
                                {
                                    foreach (Atom atom2 in Residue2)
                                    {
                                        if ((atom2.AtomType == "C" && atom2.AtomName != "CA") || atom2.AtomName == "H") continue; //Skip all carbons
                                        else if (atom2.AtomName == "CA" && Residue2.ThreeLetCode != "GLY") continue;
                                        else if (atom2.AtomName == "N" && atom1.AtomName == "O") //N-O is a backbone interaction
                                        {
                                            dist = (atom1.Coords - atom2.Hydrogen).Length;
                                            if (dist < 2.761 && dist >= 1.20)
                                            {
                                                if (Residue1.BackboneNeighbors.Contains(Residue2) == false) Residue1.BackboneNeighbors.Add(Residue2);
                                                if (Residue2.BackboneNeighbors.Contains(Residue1) == false) Residue2.BackboneNeighbors.Add(Residue1);
                                            }
                                        }
                                        else if ((atom1.AtomName == "O" || atom1.AtomName == "CA") && (atom2.AtomName == "CA" || atom2.AtomName == "C" || atom2.AtomName == "O")) continue; //ignore all other backbone interactions, eg O-O
                                        else
                                        {
                                            if ((atom1.Coords - atom2.Coords).Length < 4 && ((atom1.AtomType == "N1+" && atom2.AtomType == "O1-") || (atom2.AtomType == "N1+" && atom1.AtomType == "O1-")))
                                            {
                                                if (Residue1.SideChainNeighbors.Contains(Residue2) == false)
                                                {
                                                    Residue1.SideChainNeighbors.Add(Residue2);
                                                    newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\tSB", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[prev_strand].StrandNum, Residue2.ChainName);
                                                    file.WriteLine(newLine);
                                                }
                                                if (Residue2.SideChainNeighbors.Contains(Residue1) == false) Residue2.SideChainNeighbors.Add(Residue1);
                                            }
                                            else if (atom1.AtomType != "H") continue; //Having considered backbone and salt bridges, only look at H atoms in first residue and non-C atoms in second residue

                                            else if ((atom1.Coords - atom2.Coords).Length < 4.5)//only H atoms in res1 and O or N atoms in res2 are left now
                                            {
                                                Point3D Hdonor; Vector3D HAatom1; Vector3D HAatom2;
                                                Tuple<string, string> keyD = new Tuple<string, string>(Residue1.ThreeLetCode, atom1.AtomName);
                                                Tuple<string, string> keyA = new Tuple<string, string>(Residue2.ThreeLetCode, atom2.AtomName);
                                                if (HAcceptorList.ContainsKey(keyA) == true)
                                                {
                                                    Hdonor = (Point3D)Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == HDonorList[keyD].Connector)].Coords;

                                                    Vector3D e0 = (((Point3D)atom1.Coords) - Hdonor) / (((Point3D)atom1.Coords) - Hdonor).Length; // unit vect D->H
                                                    Vector3D nB = atom1.Coords - atom2.Coords; //A -> H
                                                    Vector3D nA = atom2.Coords - atom1.Coords; //H ->A

                                                    HAatom1 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == HAcceptorList[keyA].nextAtom)].Coords;
                                                    HAatom2 = Residue2.Atoms[Residue2.Atoms.FindIndex(a => a.AtomName == HAcceptorList[keyA].dihedralAtom)].Coords;
                                                    atom2.e1 = SharedFunctions.sete1(atom2, HAatom1, HAatom2, HAcceptorList[keyA].eAngle, HAcceptorList[keyA].dihed1);
                                                    atom2.e2 = SharedFunctions.sete1(atom2, HAatom1, HAatom2, HAcceptorList[keyA].eAngle, HAcceptorList[keyA].dihed2);

                                                    double alpha = Vector3D.AngleBetween(e0, nA);
                                                    double beta = Vector3D.AngleBetween(atom2.e1, nB);
                                                    double beta2 = Vector3D.AngleBetween(atom2.e2, nB);
                                                    d = nA.Length;

                                                    double multiTop = (sigmaD2 - Math.Pow((d - d0), 2)) * (Math.Cos(alpha) - cosAM) * (Math.Cos(beta) - cosBM);
                                                    if (multiTop < 0) multiTop = 0;
                                                    double w = Math.Sqrt(multiTop) / (wdenom);
                                                    hBondE = (1 - w) * SharedFunctions.calculateEvdw(HAcceptorList[keyA], HDonorList[keyD], d) + w * B * Global.partialChargesDict[keyA] * Global.partialChargesDict[keyD];

                                                    if (d < 2.75)
                                                    {
                                                        newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue2.SeqID, Residue2.ThreeLetCode, strandlist[prev_strand].StrandNum, Residue2.ChainName, hBondE, d);
                                                        file.WriteLine(newLine);
                                                        if (Residue1.SideChainNeighbors.Contains(Residue2) == false) Residue1.SideChainNeighbors.Add(Residue2);
                                                        if (Residue2.SideChainNeighbors.Contains(Residue1) == false) Residue2.SideChainNeighbors.Add(Residue1);
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }

                                foreach (Res Residue3 in strandlist[next_strand])
                                {
                                    foreach (Atom atom2 in Residue3)
                                    {
                                        if ((atom2.AtomType == "C" && atom2.AtomName != "CA") || atom2.AtomName == "H") continue; //Skip all carbons
                                        else if (atom2.AtomName == "CA" && Residue3.ThreeLetCode != "GLY") continue;
                                        else if (atom2.AtomName == "N" && atom1.AtomName == "O") //N-O is a backbone interaction
                                        {
                                            dist = (atom1.Coords - atom2.Hydrogen).Length;
                                            if (dist < 2.761 && dist >= 1.20)
                                            {
                                                if (Residue1.BackboneNeighbors.Contains(Residue3) == false) Residue1.BackboneNeighbors.Add(Residue3);
                                                if (Residue3.BackboneNeighbors.Contains(Residue1) == false) Residue3.BackboneNeighbors.Add(Residue1);
                                            }
                                        }
                                        else if ((atom1.AtomName == "O" || atom1.AtomName == "CA") && (atom2.AtomName == "CA" || atom2.AtomName == "C" || atom2.AtomName == "O")) continue; //ignore all other backbone interactions, eg O-O
                                        else
                                        {
                                            if ((atom1.Coords - atom2.Coords).Length < 4 && ((atom1.AtomType == "N1+" && atom2.AtomType == "O1-") || (atom2.AtomType == "N1+" && atom1.AtomType == "O1-")))
                                            {
                                                if (Residue1.SideChainNeighbors.Contains(Residue3) == false)
                                                {
                                                    Residue1.SideChainNeighbors.Add(Residue3);
                                                    newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\tSB", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue3.SeqID, Residue3.ThreeLetCode, strandlist[next_strand].StrandNum, Residue3.ChainName);
                                                    file.WriteLine(newLine);
                                                }
                                                if (Residue3.SideChainNeighbors.Contains(Residue1) == false) Residue3.SideChainNeighbors.Add(Residue1);
                                            }
                                            else if (atom1.AtomType != "H") continue; //Having considered backbone and salt bridges, only look at H atoms in first residue and non-C atoms in second residue

                                            else if ((atom1.Coords - atom2.Coords).Length < 4.5)//only H atoms in res1 and O or N atoms in res2 are left now
                                            {
                                                Point3D Hdonor; Vector3D HAatom1; Vector3D HAatom2;
                                                Tuple<string, string> keyD = new Tuple<string, string>(Residue1.ThreeLetCode, atom1.AtomName);
                                                Tuple<string, string> keyA = new Tuple<string, string>(Residue3.ThreeLetCode, atom2.AtomName);
                                                if (HAcceptorList.ContainsKey(keyA) == true)
                                                {
                                                    Hdonor = (Point3D)Residue1.Atoms[Residue1.Atoms.FindIndex(a => a.AtomName == HDonorList[keyD].Connector)].Coords;

                                                    Vector3D e0 = (((Point3D)atom1.Coords) - Hdonor) / (((Point3D)atom1.Coords) - Hdonor).Length; // unit vect D->H
                                                    Vector3D nB = atom1.Coords - atom2.Coords; //A -> H
                                                    Vector3D nA = atom2.Coords - atom1.Coords; //H ->A

                                                    HAatom1 = Residue3.Atoms[Residue3.Atoms.FindIndex(a => a.AtomName == HAcceptorList[keyA].nextAtom)].Coords;
                                                    HAatom2 = Residue3.Atoms[Residue3.Atoms.FindIndex(a => a.AtomName == HAcceptorList[keyA].dihedralAtom)].Coords;
                                                    atom2.e1 = SharedFunctions.sete1(atom2, HAatom1, HAatom2, HAcceptorList[keyA].eAngle, HAcceptorList[keyA].dihed1);
                                                    atom2.e2 = SharedFunctions.sete1(atom2, HAatom1, HAatom2, HAcceptorList[keyA].eAngle, HAcceptorList[keyA].dihed2);

                                                    double alpha = Vector3D.AngleBetween(e0, nA);
                                                    double beta = Vector3D.AngleBetween(atom2.e1, nB);
                                                    double beta2 = Vector3D.AngleBetween(atom2.e2, nB);
                                                    d = nA.Length;

                                                    double multiTop = (sigmaD2 - Math.Pow((d - d0), 2)) * (Math.Cos(alpha) - cosAM) * (Math.Cos(beta) - cosBM);
                                                    if (multiTop < 0) multiTop = 0;
                                                    double w = Math.Sqrt(multiTop) / (wdenom);
                                                    hBondE = (1 - w) * SharedFunctions.calculateEvdw(HAcceptorList[keyA], HDonorList[keyD], d) + w * B * Global.partialChargesDict[keyA] * Global.partialChargesDict[keyD];

                                                    if (d < 2.75)
                                                    {
                                                        newLine = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", atom1.AtomName, Residue1.SeqID, Residue1.ThreeLetCode, strand.StrandNum, Residue1.ChainName, atom2.AtomName, Residue3.SeqID, Residue3.ThreeLetCode, strandlist[next_strand].StrandNum, Residue3.ChainName, hBondE, d);
                                                        file.WriteLine(newLine);
                                                        if (Residue1.SideChainNeighbors.Contains(Residue3) == false) Residue1.SideChainNeighbors.Add(Residue3);
                                                        if (Residue3.SideChainNeighbors.Contains(Residue1) == false) Residue3.SideChainNeighbors.Add(Residue1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                /*foreach (Res neigh in Residue1.SideChainNeighbors)
                                {
                                    file.WriteLine("Residue {0}({1}) in strand {2} of chain {3} may interact with with residue {4}({5}) in strand {6} of chain {7}", Residue1.SeqID, Residue1.ThreeLetCode, Residue1.StrandNum, Residue1.ChainName, neigh.SeqID, neigh.ThreeLetCode, neigh.StrandNum, neigh.ChainName);
                                }*/
                            }
                        }
                    }
                }
            }
            /*string fileLocation5 = outputDirectory + "HBondingBackbone" + this.PdbName + ".txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation5))
            {
                foreach (Strand strand in strandlist)
                {
                    foreach (Res Residue1 in strand)
                    {
                        foreach (Res neigh in Residue1.BackboneNeighbors)
                        {
                            file.WriteLine("Residue {0}({1}) in strand {2} of chain {3} has a backbone interaction with residue {4}({5}) in strand {6} of chain {7}", Residue1.ResNum, Residue1.ThreeLetCode, Residue1.StrandNum, Residue1.ChainName, neigh.ResNum, neigh.ThreeLetCode, neigh.StrandNum, neigh.ChainName);

                        }
                    }
                    file.WriteLine();
                }
            }*/

        }

        public static void getTwist(List<Strand> strandlist, ref Dictionary<string, AminoAcid> _aaDict)
        {
            // first determine the relevant position on the next strand
            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {
                for (int resCtr = 0; resCtr < strandlist[strandCtr].Residues.Count; resCtr++)
                {
                    double closestDistance = 5.75;
                    Vector3D myResCoords = new Vector3D();
                    myResCoords = strandlist[strandCtr].Residues[resCtr].BackboneCoords["CA"];
                    //set closest previous strand residue
                    if (strandCtr != 0)
                    {
                        for (int prevStrandResCtr = 0; prevStrandResCtr < strandlist[strandCtr - 1].Residues.Count; prevStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandCtr - 1].Residues[prevStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].previousStrandRes = strandlist[strandCtr - 1].Residues[prevStrandResCtr];
                            }
                        }
                    }
                    else
                    {
                        for (int prevStrandResCtr = 0; prevStrandResCtr < strandlist[strandlist.Count - 1].Residues.Count; prevStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandlist.Count - 1].Residues[prevStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].previousStrandRes = strandlist[strandlist.Count - 1].Residues[prevStrandResCtr];
                            }
                        }
                    }

                    //set closest next strand residue
                    closestDistance = 5.75;
                    if (strandCtr != strandlist.Count - 1)
                    {
                        for (int nextStrandResCtr = 0; nextStrandResCtr < strandlist[strandCtr + 1].Residues.Count; nextStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandCtr + 1].Residues[nextStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].nextStrandRes = strandlist[strandCtr + 1].Residues[nextStrandResCtr];
                            }
                        }
                    }
                    else
                    {
                        for (int nextStrandResCtr = 0; nextStrandResCtr < strandlist[0].Residues.Count; nextStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[0].Residues[nextStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].nextStrandRes = strandlist[0].Residues[nextStrandResCtr];
                            }
                        }
                    }
                }
            }
            //then determine the dihedral angle
            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {
                for (int resCtr = 0; resCtr < strandlist[strandCtr].Residues.Count; resCtr++)
                {

                    Res myRes = strandlist[strandCtr].Residues[resCtr];
                    if (resCtr + 1 < strandlist[strandCtr].Residues.Count)
                    {

                        Vector3D c1 = new Vector3D();
                        Vector3D c2 = new Vector3D();
                        Vector3D c3 = new Vector3D();
                        Vector3D c4 = new Vector3D();
                        Vector3D c5 = new Vector3D();
                        Vector3D c6 = new Vector3D();

                        c3 = strandlist[strandCtr].Residues[resCtr].BackboneCoords["CA"];
                        c4 = strandlist[strandCtr].Residues[resCtr + 1].BackboneCoords["CA"];

                        if (strandlist[strandCtr].Residues[resCtr].previousStrandRes != null && strandlist[strandCtr].Residues[resCtr].previousStrandRes != null && myRes.previousStrandRes.ResStrandNum != 0 && resCtr + 1 < strandlist[strandCtr].Residues.Count)
                        {
                            Strand previousStrand = strandlist[strandlist[strandCtr].Residues[resCtr].previousStrandRes.StrandNum];
                            c1 = previousStrand.Residues[myRes.previousStrandRes.ResStrandNum - 1].BackboneCoords["CA"];
                            c2 = previousStrand.Residues[myRes.previousStrandRes.ResStrandNum].BackboneCoords["CA"];
                            strandlist[strandCtr].Residues[resCtr].Twist_prev = CalculateTorsion(c1, c2, c3, c4);
                        }
                        if (strandlist[strandCtr].Residues[resCtr].nextStrandRes != null && strandlist[strandCtr].Residues[resCtr].nextStrandRes != null && myRes.nextStrandRes.ResStrandNum != 0 && resCtr + 1 < strandlist[strandCtr].Residues.Count)
                        {
                            Strand nextStrand = strandlist[strandlist[strandCtr].Residues[resCtr].nextStrandRes.StrandNum];
                            c5 = nextStrand.Residues[myRes.nextStrandRes.ResStrandNum - 1].BackboneCoords["CA"];
                            c6 = nextStrand.Residues[myRes.nextStrandRes.ResStrandNum].BackboneCoords["CA"];
                            strandlist[strandCtr].Residues[resCtr].Twist_next = CalculateTorsion(c4, c3, c6, c5);

                            //_aaDict[strandlist[strandCtr].Residues[resCtr].ThreeLetCode].Twist.Add(strandlist[strandCtr].Residues[resCtr].Dihedral);
                        }
                    }
                    //_aaDict[strandlist[strandCtr].Residues[resCtr].ThreeLetCode].Twist.Add(-170.0); //if I can't get a dihedral let me know

                }
            }
        }

        public static List<double> writeTwists(List<Strand> strandlist, string outputDirectory, string pdbName)//, ref Dictionary<string, AminoAcid> _aaDict)
        {
            //returns theta as definied by Murzin 1994
            List<double> all_prev_twist = new List<double>();
            List<double> all_next_twist = new List<double>();
            // first determine the relevant position on the next strand
            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {
                for (int resCtr = 0; resCtr < strandlist[strandCtr].Residues.Count; resCtr++)
                {
                    double closestDistance = 5.75;
                    Vector3D myResCoords = new Vector3D();
                    myResCoords = strandlist[strandCtr].Residues[resCtr].BackboneCoords["CA"];
                    //set closest previous strand residue
                    if (strandCtr != 0)
                    {
                        for (int prevStrandResCtr = 0; prevStrandResCtr < strandlist[strandCtr - 1].Residues.Count; prevStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandCtr - 1].Residues[prevStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].previousStrandRes = strandlist[strandCtr - 1].Residues[prevStrandResCtr];
                            }
                        }
                    }
                    else
                    {
                        for (int prevStrandResCtr = 0; prevStrandResCtr < strandlist[strandlist.Count - 1].Residues.Count; prevStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandlist.Count - 1].Residues[prevStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].previousStrandRes = strandlist[strandlist.Count - 1].Residues[prevStrandResCtr];
                            }
                        }
                    }

                    //set closest next strand residue
                    closestDistance = 5.75;
                    if (strandCtr != strandlist.Count - 1)
                    {
                        for (int nextStrandResCtr = 0; nextStrandResCtr < strandlist[strandCtr + 1].Residues.Count; nextStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[strandCtr + 1].Residues[nextStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].nextStrandRes = strandlist[strandCtr + 1].Residues[nextStrandResCtr];
                            }
                        }
                    }
                    else
                    {
                        for (int nextStrandResCtr = 0; nextStrandResCtr < strandlist[0].Residues.Count; nextStrandResCtr++)
                        {
                            double currentDistance = (myResCoords - strandlist[0].Residues[nextStrandResCtr].BackboneCoords["CA"]).Length;
                            if (currentDistance < closestDistance)
                            {
                                closestDistance = currentDistance;
                                strandlist[strandCtr].Residues[resCtr].nextStrandRes = strandlist[0].Residues[nextStrandResCtr];
                            }
                        }
                    }
                }
            }
            //then determine the dihedral angle
            string fileLocation3 = outputDirectory + "Twists\\PrevTwist_" + pdbName + ".txt";
            string fileLocation4 = outputDirectory + "Twists\\NextTwist_" + pdbName + ".txt";

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation3, true))
            {
                using (System.IO.StreamWriter file2 = new System.IO.StreamWriter(fileLocation4, true))
                {
                    for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
                    {
                        for (int resCtr = 0; resCtr < strandlist[strandCtr].Residues.Count; resCtr++)
                        {

                            Res myRes = strandlist[strandCtr].Residues[resCtr];
                            if (resCtr + 1 < strandlist[strandCtr].Residues.Count)
                            {

                                Vector3D c1 = new Vector3D();
                                Vector3D c2 = new Vector3D();
                                Vector3D c3 = new Vector3D();
                                Vector3D c4 = new Vector3D();
                                Vector3D c5 = new Vector3D();
                                Vector3D c6 = new Vector3D();

                                c3 = strandlist[strandCtr].Residues[resCtr].BackboneCoords["CA"];
                                c4 = strandlist[strandCtr].Residues[resCtr + 1].BackboneCoords["CA"];

                                if (strandlist[strandCtr].Residues[resCtr].previousStrandRes != null && strandlist[strandCtr].Residues[resCtr].previousStrandRes != null && myRes.previousStrandRes.ResStrandNum != 0 && resCtr + 1 < strandlist[strandCtr].Residues.Count)
                                {
                                    Strand previousStrand = strandlist[strandlist[strandCtr].Residues[resCtr].previousStrandRes.StrandNum];
                                    c1 = previousStrand.Residues[myRes.previousStrandRes.ResStrandNum - 1].BackboneCoords["CA"];
                                    c2 = previousStrand.Residues[myRes.previousStrandRes.ResStrandNum].BackboneCoords["CA"];
                                    strandlist[strandCtr].Residues[resCtr].Twist_prev = CalculateTorsion(c1, c2, c3, c4);
                                    double distance = (c2 - c3).Length;
                                    file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", previousStrand.Residues[myRes.previousStrandRes.ResStrandNum - 1].SeqID, previousStrand.Residues[myRes.previousStrandRes.ResStrandNum].SeqID, strandlist[strandCtr].Residues[resCtr].SeqID, strandlist[strandCtr].Residues[resCtr + 1].SeqID, strandlist[strandCtr].Residues[resCtr].Twist_prev, distance, strandCtr, previousStrand.StrandNum);
                                    all_prev_twist.Add(strandlist[strandCtr].Residues[resCtr].Twist_prev);
                                }
                                if (strandlist[strandCtr].Residues[resCtr].nextStrandRes != null && strandlist[strandCtr].Residues[resCtr].nextStrandRes != null && myRes.nextStrandRes.ResStrandNum != 0 && resCtr + 1 < strandlist[strandCtr].Residues.Count)
                                {
                                    Strand nextStrand = strandlist[strandlist[strandCtr].Residues[resCtr].nextStrandRes.StrandNum];
                                    c5 = nextStrand.Residues[myRes.nextStrandRes.ResStrandNum - 1].BackboneCoords["CA"];
                                    c6 = nextStrand.Residues[myRes.nextStrandRes.ResStrandNum].BackboneCoords["CA"];
                                    double distance = (c6 - c3).Length;
                                    strandlist[strandCtr].Residues[resCtr].Twist_next = CalculateTorsion(c4, c3, c6, c5);
                                    file2.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", strandlist[strandCtr].Residues[resCtr + 1].SeqID, strandlist[strandCtr].Residues[resCtr].SeqID, nextStrand.Residues[myRes.nextStrandRes.ResStrandNum].SeqID, nextStrand.Residues[myRes.nextStrandRes.ResStrandNum - 1].SeqID, strandlist[strandCtr].Residues[resCtr].Twist_next, distance, strandCtr, nextStrand.StrandNum);

                                    all_next_twist.Add(strandlist[strandCtr].Residues[resCtr].Twist_next);
                                    //_aaDict[strandlist[strandCtr].Residues[resCtr].ThreeLetCode].Twist.Add(strandlist[strandCtr].Residues[resCtr].Dihedral);
                                }
                            }
                            //_aaDict[strandlist[strandCtr].Residues[resCtr].ThreeLetCode].Twist.Add(-170.0); //if I can't get a dihedral let me know

                        }
                    }
                }
            }
            return (all_prev_twist);
        }

        public static double getTiltsByAA(List<Strand> strandlist, string outputDirectory, string PdbName, Vector3D axis, ref Dictionary<string, AminoAcid> _AADict)
        {
            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
            int totalResidues = 0;

            List<int> thePts = new List<int>();
            for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
            {
                thePts.Add(_AADict[AAarray[aaCtr]].Tilt.Count);
            }
            List<int[]> myList = new List<int[]>();

            double AvgTilt =0;
            string fileLocation = outputDirectory + "\\Tilts\\TiltAngles" + PdbName + ".txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))
            {
                for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
                {

                    List<int> numOfAA = new List<int>();
                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {
                        numOfAA.Add(_AADict[AAarray[aaCtr]].Tilt.Count);
                    }
                    strandlist[strandCtr].getTiltsByAA(axis, strandCtr, ref _AADict);
                    AvgTilt += strandlist[strandCtr].AvgTilt; //*(strandlist[strandCtr].Residues.Count - 2);
                    totalResidues += (strandlist[strandCtr].Residues.Count - 2);

                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {

                        for (int valueCtr = numOfAA[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt.Count; valueCtr++)
                        {
                            double tiltFromStrand = _AADict[AAarray[aaCtr]].Tilt[valueCtr] - strandlist[strandCtr].AvgTilt;
                            _AADict[AAarray[aaCtr]].TiltfromStrand.Add(tiltFromStrand);

                        }
                    }
                    int[] aminoAcidTally = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

                    //using just the aminoAcids i'm measuring for tilt ie minus one for each end
                    for (int residueCtr = 1; residueCtr < strandlist[strandCtr].Residues.Count - 1; residueCtr++)
                    {
                        for (int idCtr = 0; idCtr < 20; idCtr++)
                        {
                            if (AAarray[idCtr] == strandlist[strandCtr].Residues[residueCtr].ThreeLetCode) aminoAcidTally[idCtr]++;
                        }

                    }
                    aminoAcidTally[20] = strandCtr;
                    myList.Add(aminoAcidTally);

                }
                AvgTilt = AvgTilt / strandlist.Count;
                //AvgTilt = AvgTilt / totalResidues;
                for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                {// for all amino acids
                    for (int valueCtr = thePts[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt.Count; valueCtr++)
                    {
                        double tiltFromBarrel = _AADict[AAarray[aaCtr]].Tilt[valueCtr] - AvgTilt;
                        ////using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Ryan\Dropbox\output\AA_Tilt_special.txt", true))
                        ////{

                        //    if (aaCtr == 8 && tiltFromBarrel < -7.5 && tiltFromBarrel > -12.5)
                        //    {// his in -10 bin
                        //        int seqID = _AADict[AAarray[aaCtr]].seqIDList[valueCtr];
                        //       // file.WriteLine("the HIS {0} in this protein has a tilt of {1}", seqID, tiltFromBarrel);
                        //    }
                        //    if (aaCtr == 13 && tiltFromBarrel < 22.5 && tiltFromBarrel > 17.5)
                        //    {// phe in 20 bin
                        //        int seqID = _AADict[AAarray[aaCtr]].seqIDList[valueCtr];
                        //        //file.WriteLine("the  PHE {0} in this protein has a tilt of {1}", seqID, tiltFromBarrel);
                        //    }
                        //}

                        _AADict[AAarray[aaCtr]].TiltfromBarrel.Add(tiltFromBarrel);

                    }

                }
                for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
                {
                    for (int i = 0; i < 21; i++)
                    {
                        file.Write("{0}\t", myList[strandCtr][i]);
                    }

                    file.Write("{0}\t{1}\n", strandlist[strandCtr].AvgTilt, AvgTilt);
                }
            }
            return AvgTilt;

        }

        public static void getTiltsByAA_divided(List<Strand> strandlist, string outputDirectory, string PdbName, ref Dictionary<string, AminoAcid> _AADict, Vector3D axis, out double AvgTilt_even, out double AvgTilt_odd)
        {
            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
            int totalResidues_odd = 0;
            int totalResidues_even = 0;

            AvgTilt_even = 0; AvgTilt_odd = 0;
            //get the number of data points for each amino acid
            List<int> evenPts = new List<int>();
            List<int> oddPts = new List<int>();
            for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
            {
                evenPts.Add(_AADict[AAarray[aaCtr]].Tilt_even.Count);
                oddPts.Add(_AADict[AAarray[aaCtr]].Tilt_odd.Count);

            }

            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {

                if (strandCtr % 2 == 0)
                {
                    List<int> numOfAA = new List<int>();
                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {
                        numOfAA.Add(_AADict[AAarray[aaCtr]].Tilt_even.Count);
                    }

                    strandlist[strandCtr].getTiltsbyAA_divided(axis, strandCtr, ref _AADict);


                    AvgTilt_even += strandlist[strandCtr].AvgTilt_even * (strandlist[strandCtr].Residues.Count - 2);
                    totalResidues_even += (strandlist[strandCtr].Residues.Count - 2);

                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {

                        for (int valueCtr = numOfAA[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt_even.Count; valueCtr++)
                        {
                            double tiltFromStrand_even = _AADict[AAarray[aaCtr]].Tilt_even[valueCtr] - strandlist[strandCtr].AvgTilt_even;
                            _AADict[AAarray[aaCtr]].TiltfromStrand_even.Add(tiltFromStrand_even);

                        }
                    }
                }
                else
                {
                    List<int> numOfAA = new List<int>();
                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {
                        numOfAA.Add(_AADict[AAarray[aaCtr]].Tilt_odd.Count);
                    }
                    strandlist[strandCtr].getTiltsbyAA_divided(axis, strandCtr, ref _AADict);


                    AvgTilt_odd += strandlist[strandCtr].AvgTilt_odd * (strandlist[strandCtr].Residues.Count - 2);
                    totalResidues_odd += (strandlist[strandCtr].Residues.Count - 2);

                    for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
                    {

                        for (int valueCtr = numOfAA[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt_odd.Count; valueCtr++)
                        {
                            double tiltFromStrand_odd = _AADict[AAarray[aaCtr]].Tilt_odd[valueCtr] - strandlist[strandCtr].AvgTilt_odd;
                            _AADict[AAarray[aaCtr]].TiltfromStrand_odd.Add(tiltFromStrand_odd);

                        }
                    }
                }
            }
            AvgTilt_even = AvgTilt_even / totalResidues_even;
            AvgTilt_odd = AvgTilt_odd / totalResidues_odd;

            for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
            {
                for (int valueCtr = evenPts[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt_even.Count; valueCtr++)
                {
                    double tiltFromBarrel_even = _AADict[AAarray[aaCtr]].Tilt_even[valueCtr] - AvgTilt_even;
                    _AADict[AAarray[aaCtr]].TiltfromBarrel_even.Add(tiltFromBarrel_even);
                }
            }
            for (int aaCtr = 0; aaCtr < _AADict.Count; aaCtr++)
            {
                for (int valueCtr = oddPts[aaCtr]; valueCtr < _AADict[AAarray[aaCtr]].Tilt_odd.Count; valueCtr++)
                {
                    double tiltFromBarrel_odd = _AADict[AAarray[aaCtr]].Tilt_odd[valueCtr] - AvgTilt_odd;
                    _AADict[AAarray[aaCtr]].TiltfromBarrel_odd.Add(tiltFromBarrel_odd);
                }
            }
        }

        public static void runBetaBarrel_RYAN(string pdb, ref Protein m_Protein, ref Barrel m_Barrel)
        {
            try
            {
                Directory.SetCurrentDirectory(Global.MONO_DB_DIR);
                string PDB = pdb.Substring(0, 4).ToUpper();
                string pdbFileName = PDB + ".pdb";

                if (File.Exists(pdbFileName))
                {
                    AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                    ////Console.WriteLine("opened {0}", pdbFileName);
                    myAtomCat = Program.readPdbFile(pdbFileName, ref Global.partialChargesDict);
                    int chainNum = 0;

                    int stop = myAtomCat.ChainAtomList.Count();
                    ////Console.Write("Protein Class {0}", chainNum);

                    m_Protein = new MonoProtein(ref myAtomCat, chainNum, PDB);

                    ////Console.Write("creating barrel class..");

                    m_Barrel = new MonoBarrel(m_Protein.Chains[0], m_Protein);
                }
                else
                {
                    Directory.SetCurrentDirectory(Global.POLY_DB_DIR);

                    pdbFileName = pdb.Substring(0, 4).ToUpper() + ".pdb";
                    if (File.Exists(pdbFileName))
                    {
                        AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                        //Console.WriteLine("opened {0}", pdbFileName);
                        myAtomCat = Program.readPdbFile(pdbFileName, ref Global.partialChargesDict);

                        PolyProtein newProt = new PolyProtein(ref myAtomCat, pdb); //For considering all chains

                        //Console.WriteLine("creating barrel class..");
                        m_Protein = newProt;
                        m_Barrel = new PolyBarrel(newProt, Global.POLY_OUTPUT_DIR, Global.POLY_DB_DIR);
                    }
                    else
                    {
                        pdbFileName = pdb.Substring(0, 4).ToLower() + ".pdb";
                        if (File.Exists(pdbFileName))
                        {
                            AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                            //Console.WriteLine("opened {0}", pdbFileName);
                            myAtomCat = Program.readPdbFile(pdbFileName, ref Global.partialChargesDict);

                            //Console.WriteLine("opened {0}", pdbFileName);
                            myAtomCat = Program.readPdbFile(pdbFileName, ref Global.partialChargesDict);

                            PolyProtein newProt = new PolyProtein(ref myAtomCat, pdb); //For considering all chains

                            //Console.WriteLine("creating barrel class..");
                            m_Protein = newProt;
                            m_Barrel = new PolyBarrel(newProt, Global.POLY_OUTPUT_DIR, Global.POLY_DB_DIR);
                        }
                        else
                        {

                            //Console.WriteLine("could not find {0}", pdbFileName);
                        }
                    }
                }
            }
            catch
            {
                Console.WriteLine("could not run: " + pdb);
            }
        }


        #region JoannasFunctions
        //-----------------------------------Joanna's functions----------------------------

        //outputs amino acids' PhiPsi to a file
        public static void outputAA(ref Dictionary<string, AminoAcid> _aaDict)
        {
            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

            outputThisAngleVariable("Tilt", ref _aaDict, true);
            outputThisAngleVariable("caThetaList", ref _aaDict, true);
            outputThisAngleVariable("TiltfromStrand", ref _aaDict, true);
            outputThisAngleVariable("TiltfromBarrel", ref _aaDict, true);
            outputThisAngleVariable("Twist", ref _aaDict, true);

            string fileLocation = Global.POLY_OUTPUT_DIR + "AA_phiPsi.txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    for (int valueCtr = 0; valueCtr < _aaDict[AAarray[aaCtr]].phiPsiList.Count; valueCtr++)
                    {
                        file.Write("PhiPsi\t{0}\t", AAarray[aaCtr]);
                        file.Write("{0}\t", _aaDict[AAarray[aaCtr]].phiPsiList[valueCtr][0], _aaDict[AAarray[aaCtr]].phiPsiList[valueCtr][1]);
                        file.WriteLine();
                    }
                }
            }
        }

        //outputs amino acid angle to a file
        public static void outputThisAngleVariable(string varName, ref Dictionary<string, AminoAcid> _aaDict, bool isNeg)
        {
            string fileLocation = Global.POLY_OUTPUT_DIR + "AA_" + varName + ".txt";

            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

            List<List<int>> binnedVariable = new List<List<int>>();
            make180List(ref binnedVariable);
            int adder = 0;
            if (isNeg == true) adder = 180;


            for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
            {
                for (int valueCtr = 0; valueCtr < _aaDict[AAarray[aaCtr]].GetProperty(varName).Count; valueCtr++)
                {

                    int variable = Convert.ToInt32(Math.Round((_aaDict[AAarray[aaCtr]].GetProperty(varName)[valueCtr] + adder) / 5));

                    binnedVariable[aaCtr][variable]++;
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {

                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("{0}\t{1}\t", varName, AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedVariable[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }
            }
        }

        //helper function to make bins for outputThisAngleVariable(...) function
        public static void make180List(ref List<List<int>> _myList)
        {
            for (int aaCtr = 0; aaCtr < 20; aaCtr++)
            {
                List<int> OneEightyList = new List<int>();

                for (int binCtr = 0; binCtr < 72; binCtr++)
                {
                    OneEightyList.Add(0);
                }

                _myList.Add(OneEightyList);
            }
        }

        //outputs amino acid tilts
        public static void outputDividedAA(ref Dictionary<string, AminoAcid> _aaDict)
        {
            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

            string fileLocation = Global.POLY_OUTPUT_DIR + "AA_Tilt.txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                for (int AActr = 0; AActr < AAarray.Count(); AActr++)
                {
                    if (_aaDict[AAarray[AActr]].Tilt_even.Count > 0)
                    {
                        file.WriteLine("{0} EVEN instances:{1} avgTilt: {2} from strand {3} from barrel {4}", AAarray[AActr], _aaDict[AAarray[AActr]].Tilt_even.Count, _aaDict[AAarray[AActr]].Tilt_even.Average(), _aaDict[AAarray[AActr]].TiltfromStrand_even.Average(), _aaDict[AAarray[AActr]].TiltfromBarrel_even.Average());
                    }
                    else file.WriteLine("{0} EVEN instances:{1} avgTilt: {2} from strand {3} from barrel {4}", AAarray[AActr], 0, 0, 0, 0);

                    if (_aaDict[AAarray[AActr]].Tilt_odd.Count > 0)
                    {
                        file.WriteLine("{0} ODD instances:{1} avgTilt: {2} from strand {3} from barrel {4}", AAarray[AActr], _aaDict[AAarray[AActr]].Tilt_odd.Count, _aaDict[AAarray[AActr]].Tilt_odd.Average(), _aaDict[AAarray[AActr]].TiltfromStrand_odd.Average(), _aaDict[AAarray[AActr]].TiltfromBarrel_odd.Average());
                    }
                    else file.WriteLine("{0} ODD instances:{1} avgTilt: {2} from strand {3} from barrel {4}", AAarray[AActr], 0, 0, 0, 0);
                    file.WriteLine();
                }

            }
            fileLocation = Global.POLY_OUTPUT_DIR + "AA_Tilt_raw2.txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                for (int AActr = 0; AActr < AAarray.Count(); AActr++)
                {

                    file.Write("EVEN \t {0}\tTilt\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].Tilt_even.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].Tilt_even[instanceCtr]);
                    }
                    file.WriteLine();
                    file.Write("EVEN \t {0}\tTilt from Strand\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].TiltfromStrand_even.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].TiltfromStrand_even[instanceCtr]);
                    }
                    file.WriteLine();
                    file.Write("EVEN \t {0}\tTilt from Barrel\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].TiltfromBarrel_even.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].TiltfromBarrel_even[instanceCtr]);
                    }
                    file.WriteLine();
                }
                for (int AActr = 0; AActr < AAarray.Count(); AActr++)
                {

                    file.Write("ODD \t {0}\tTilt\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].Tilt_odd.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].Tilt_odd[instanceCtr]);
                    }
                    file.WriteLine();
                    file.Write("ODD \t {0}\tTilt from Strand\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].TiltfromStrand_odd.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].TiltfromStrand_odd[instanceCtr]);
                    }
                    file.WriteLine();
                    file.Write("ODD \t {0}\tTilt from Barrel\t", AAarray[AActr]);
                    for (int instanceCtr = 0; instanceCtr < _aaDict[AAarray[AActr]].TiltfromBarrel_odd.Count; instanceCtr++)
                    {
                        file.Write("{0}\t", _aaDict[AAarray[AActr]].TiltfromBarrel_odd[instanceCtr]);
                    }
                    file.WriteLine();
                }
            }


            List<List<int>> binnedTilt_even = new List<List<int>>();
            List<List<int>> binnedTilt_odd = new List<List<int>>();
            List<List<int>> binnedTiltStrand_even = new List<List<int>>();
            List<List<int>> binnedTiltStrand_odd = new List<List<int>>();
            List<List<int>> binnedTiltBarrel_even = new List<List<int>>();
            List<List<int>> binnedTiltBarrel_odd = new List<List<int>>();

            for (int aaCtr = 0; aaCtr < 20; aaCtr++)
            {
                List<int> OneEightyList = new List<int>();

                for (int binCtr = 0; binCtr < 72; binCtr++)
                {
                    OneEightyList.Add(0);
                }
                binnedTilt_even.Add(OneEightyList);
                binnedTilt_odd.Add(OneEightyList);
                binnedTiltStrand_even.Add(OneEightyList);
                binnedTiltStrand_odd.Add(OneEightyList);
                binnedTiltBarrel_even.Add(OneEightyList);
                binnedTiltBarrel_odd.Add(OneEightyList);

            }
            for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
            {
                for (int valueCtr = 0; valueCtr < _aaDict[AAarray[aaCtr]].Tilt_even.Count; valueCtr++)
                {

                    int tilt = Convert.ToInt32(Math.Round(_aaDict[AAarray[aaCtr]].Tilt_even[valueCtr] / 5));
                    int strand = Convert.ToInt32(Math.Round((_aaDict[AAarray[aaCtr]].TiltfromStrand_even[valueCtr] + 180.0) / 5));
                    int barrel = Convert.ToInt32(Math.Round((_aaDict[AAarray[aaCtr]].TiltfromBarrel_even[valueCtr] + 180.0) / 5));


                    binnedTilt_even[aaCtr][tilt]++;
                    binnedTiltStrand_even[aaCtr][strand]++;
                    binnedTiltBarrel_even[aaCtr][barrel]++;
                }
            }
            for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
            {
                for (int valueCtr = 0; valueCtr < _aaDict[AAarray[aaCtr]].Tilt_odd.Count; valueCtr++)
                {
                    int tilt = Convert.ToInt32(Math.Round(_aaDict[AAarray[aaCtr]].Tilt_odd[valueCtr] / 5));
                    int strand = Convert.ToInt32(Math.Round((_aaDict[AAarray[aaCtr]].TiltfromStrand_odd[valueCtr] + 180.0) / 5));
                    int barrel = Convert.ToInt32(Math.Round((_aaDict[AAarray[aaCtr]].TiltfromBarrel_odd[valueCtr] + 180.0) / 5));

                    binnedTilt_odd[aaCtr][tilt]++;
                    binnedTiltStrand_odd[aaCtr][strand]++;
                    binnedTiltBarrel_odd[aaCtr][barrel]++;
                }
            }
            fileLocation = Global.POLY_OUTPUT_DIR + "AA_Tilt_raw5.txt";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {

                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("EVEN\tTilt\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTilt_even[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }

                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("EVEN\tTilt_strand\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTiltStrand_even[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }
                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("EVEN\tTilt_barrel\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTiltBarrel_even[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }


                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("ODD\tTilt\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTilt_odd[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }

                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("ODD\tTilt_strand\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTiltStrand_odd[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }
                for (int aaCtr = 0; aaCtr < _aaDict.Count; aaCtr++)
                {
                    file.Write("ODD\tTilt_barrel\t{0}\t", AAarray[aaCtr]);
                    for (int valueCtr = 0; valueCtr < 72; valueCtr++)
                    {
                        file.Write("{0}\t", binnedTiltBarrel_odd[aaCtr][valueCtr]);

                    }
                    file.WriteLine();
                }


            }
        }

        //returns bottom strand residues
        public static List<Vector3D> getTopEllipseCoords(Barrel myBarrel)
        {
            List<Vector3D> myEllipse = new List<Vector3D>();
            for (int strandCtr = 0; strandCtr < myBarrel.Strands.Count; strandCtr++)
            {
                Vector3D firstCA = new Vector3D();
                firstCA = myBarrel.Strands[strandCtr].Residues[0].BackboneCoords["CA"];
                Vector3D lastCA = new Vector3D();
                lastCA = myBarrel.Strands[strandCtr].Residues[myBarrel.Strands[strandCtr].Residues.Count - 1].BackboneCoords["CA"];

                if (strandCtr % 2 == 0)
                {
                    myEllipse.Add(lastCA);
                }
                else
                {
                    myEllipse.Add(firstCA);
                }
            }
            return myEllipse;

        }

        //returns top strand residues
        public static List<Vector3D> getBottomEllipseCoords(Barrel myBarrel)
        {
            List<Vector3D> myEllipse = new List<Vector3D>();


            for (int strandCtr = 0; strandCtr < myBarrel.Strands.Count; strandCtr++)
            {
                Vector3D firstCA = new Vector3D();
                firstCA = myBarrel.Strands[strandCtr].Residues[0].BackboneCoords["CA"];
                Vector3D lastCA = new Vector3D();
                lastCA = myBarrel.Strands[strandCtr].Residues[myBarrel.Strands[strandCtr].Residues.Count - 1].BackboneCoords["CA"];

                if (strandCtr % 2 == 0)
                {
                    myEllipse.Add(firstCA);
                }
                else
                {
                    myEllipse.Add(lastCA);
                }
            }
            return myEllipse;
        }

        //returns top strand residues
        public static List<Vector3D> getTopEllipseCoords(List<Strand> strandlist)
        {
            List<Vector3D> myEllipse = new List<Vector3D>();
            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {
                Vector3D firstCA = new Vector3D();
                firstCA = strandlist[strandCtr].Residues[0].BackboneCoords["CA"];
                Vector3D lastCA = new Vector3D();
                lastCA = strandlist[strandCtr].Residues[strandlist[strandCtr].Residues.Count - 1].BackboneCoords["CA"];

                if (strandCtr % 2 == 0)
                {
                    myEllipse.Add(lastCA);
                }
                else
                {
                    myEllipse.Add(firstCA);
                }
            }
            return myEllipse;

        }

        //returns bottom strand residues
        public static List<Vector3D> getBottomEllipseCoords(List<Strand> strandlist)
        {
            List<Vector3D> myEllipse = new List<Vector3D>();


            for (int strandCtr = 0; strandCtr < strandlist.Count; strandCtr++)
            {
                Vector3D firstCA = new Vector3D();
                firstCA = strandlist[strandCtr].Residues[0].BackboneCoords["CA"];
                Vector3D lastCA = new Vector3D();
                lastCA = strandlist[strandCtr].Residues[strandlist[strandCtr].Residues.Count - 1].BackboneCoords["CA"];

                if (strandCtr % 2 == 0)
                {
                    myEllipse.Add(firstCA);
                }
                else
                {
                    myEllipse.Add(lastCA);
                }
            }
            return myEllipse;
        }

        //makes a dictionary? of amino acids
        public static Dictionary<string, AminoAcid> makeAADict()
        {
            string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

            Dictionary<string, AminoAcid> AADict = new Dictionary<string, AminoAcid>();

            for (int i = 0; i < AAarray.Count(); i++)
            {
                AminoAcid ALA = new AminoAcid();
                AADict.Add(AAarray[i], ALA);
            }
            return AADict;

        }

        //returns average position of a list of 3D vectors
        public static Vector3D averagePosition(List<Vector3D> AllCaxyzCoord)
        {
            Vector3D averagePosition = new Vector3D();

            foreach (Vector3D position in AllCaxyzCoord)
            {
                averagePosition += position / AllCaxyzCoord.Count;

            }
            return averagePosition;
        }


    }
}


        #endregion