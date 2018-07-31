/*
**  File: BarrelEllipse.cs
**  Started: 7/15/2015
**  Contributors: Ryan Feehan
**  Overview: Outputs data and images for best-fit ellipse at the top and bottom of betaBarrels 
**
**  About:  The member class ellipse can find the normal, semi-major and semi-minor axises 
**          of a best fit a ellipse for a group of points.  The barrelEllipse class makes
**          elipsis for the strand ends of the betaBarrel.  When the whole database is run
**          it writes relevent ellipse data results for each barrel to mono and poly results
**          file.  During this process, each PDB has a file created that conatains
**          a PDB and various CGOs, including top and bottom ellipse CGOs.  Note the 
**          atoms coordinates are adjusted while creating the barrel, so CGOs are only
**          accurate for the outputed PDB.
**
**  Last Edited: 7/15/2017 by Ryan Feehan
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media.Media3D;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;
using betaBarrelProgram.BarrelStructures;
using System.Diagnostics;

namespace betaBarrelProgram
{
    public class BarrelEllipse
    {

        /*    
        ** Will run a user defined PDB through barrelEllipse, createing
        ** PDBs and CGOs helpful for determining if a successful betaBarrel
        ** was created by the program
        */
        static public void testEllipseSinglePDB()
        {
            Stopwatch stopWatch = new Stopwatch();
            string elapsedTime;
            TimeSpan ts;


            Console.WriteLine("Enter pdb:");
            string pdb = Console.ReadLine();

            string IN_FILE = Global.MONO_DB_DIR + pdb + ".pdb";
            OUTPUT_PATH = Global.MONO_OUTPUT_DIR + "ellipse\\";

            if (!File.Exists(IN_FILE))
            { 
                IN_FILE = Global.POLY_DB_DIR + pdb + ".pdb";
                OUTPUT_PATH = Global.POLY_OUTPUT_DIR + "ellipse\\";
            }

            stopWatch.Start();
            
            if (File.Exists(IN_FILE))
            {
                        ts = stopWatch.Elapsed;
                        elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
                        Console.WriteLine("RunTime " + elapsedTime);
                        BarrelEllipse newBarrel = new BarrelEllipse(pdb);

            }
            else
            {
                Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
                Console.WriteLine("could not open {0}", IN_FILE);
                Console.ReadLine();

            }

            ts = stopWatch.Elapsed;
            elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);
        }

        /*
        **  
        */
        static public void runEllipseData()
        {
            Stopwatch stopWatch = new Stopwatch();
            string elapsedTime;
            TimeSpan ts;
            string file;

            Console.WriteLine("Mono or Poly?(m/p)");
            string input = Console.ReadLine();

            if (input == "m")
            {
                DB_FILE = Global.MACMONODBDIR;
                OUTPUT_PATH = Global.MONO_OUTPUT_DIR + "ellipse\\";
            }
            else
            {
                DB_FILE = Global.MACPOLYDBDIR;
                OUTPUT_PATH = Global.POLY_OUTPUT_DIR + "ellipse\\";
            }

            stopWatch.Start();
            file = OUTPUT_PATH + "ELLIPSE_RESULTS.txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(file))
            {
                output.Write("\tProtein:\tAngle Between Ellipsis\tTop Eccentricity:\tBottom Eccentricity:\tTop Deviation:\tBottom Deviation:");
            }
            

            if (File.Exists(DB_FILE))
            {
                using (StreamReader sr = new StreamReader(DB_FILE))
                {
                    String line;

                    while ((line = sr.ReadLine()) != null)
                    {
                        string[] splitLine = line.Split(new char[] { ' ', '\t', ',' });
                        string pdb = splitLine[0];
                        ts = stopWatch.Elapsed;
                        elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
                        Console.WriteLine("RunTime " + elapsedTime);
                        BarrelEllipse newBarrel = new BarrelEllipse(pdb);
                        newBarrel.EllipseResults(pdb);

                    }

                }
            }
            else
            {
                Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
                Console.WriteLine("could not open {0}", DB_FILE);
                Console.ReadLine();

            }

            ts = stopWatch.Elapsed;
            elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        public BarrelEllipse(string PDBid)
        {

            pdbDir = OUTPUT_PATH + PDBid + "\\";
            if (!System.IO.Directory.Exists(pdbDir))
            {
                System.IO.Directory.CreateDirectory(pdbDir);
            }
            _protein = null;
            _barrel = null;

            SharedFunctions.runBetaBarrel(PDBid, ref _protein, ref _barrel);
            strandlist = _barrel.Strands;

            BottomEllipse = new Ellipse(SharedFunctions.getBottomEllipseCoords(strandlist));
            TopEllipse = new Ellipse(SharedFunctions.getTopEllipseCoords(strandlist));

            //draw CGOs for top and bottom of strands 
            DrawEllipse("BottomStrands", SharedFunctions.getBottomEllipseCoords(strandlist));
            DrawEllipse("TopStrands", SharedFunctions.getTopEllipseCoords(strandlist));

            //draw CGOs for top and bottom ellipse 
            DrawSemiMajorMinor("BMM", BottomEllipse.m_Centroid, BottomEllipse.m_MajorAxis, BottomEllipse.m_minorAxis, BottomEllipse.m_Normal);
            DrawSemiMajorMinor("TMM", TopEllipse.m_Centroid, TopEllipse.m_MajorAxis, TopEllipse.m_minorAxis, TopEllipse.m_Normal);

            DrawEllipse("Bot", BottomEllipse.newEllipseABCtoPoints());
            DrawEllipse("Top", TopEllipse.newEllipseABCtoPoints());

            //furthestDistance("Dis", PDBid, ref TopEllipse, ref BottomEllipse);

            WritePDB(PDBid);
            SharedFunctions.writePymolScriptForStrands(strandlist, pdbDir, pdbDir, PDBid);
            ELLIPSE_BARREL_PYMOL_SCRIPT(PDBid);
        }


        /*
        ** Pre:     
        ** Post:      
        ** About: creates pymol CGO of an ellipse
        */
        private void DrawEllipse(string fileName, List<Vector3D> Ellipse)
        {

            string fileLocation = pdbDir + "\\" + fileName + ".py";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                file.Write("\nfrom pymol.cgo import *    # get constants");
                file.Write("\nfrom pymol import cmd");
                file.Write("\n");
                file.Write("\nobj = [");
                file.Write("\n\tBEGIN, LINES,");
                file.Write("\n\tCOLOR, 119/256.0, 119/256.0, 119/256.0,");
                file.Write("\n");

                int count = 0;
                for (int i = 0; i < Ellipse.Count; i++)
                {
                    if (i == 0)
                    {
                        file.Write("\n\tVERTEX," + Ellipse[Ellipse.Count - 1].X + "," + Ellipse[Ellipse.Count - 1].Y + ", " + Ellipse[Ellipse.Count - 1].Z + ",");
                        file.Write("\n\tVERTEX," + Ellipse[i].X + "," + Ellipse[i].Y + ", " + Ellipse[i].Z + ",");
                        file.Write("\n");
                    }
                    else
                    {
                        if (Ellipse[i] != Ellipse[i - 1])
                        {

                            file.Write("\n\tVERTEX," + Ellipse[i].X + "," + Ellipse[i].Y + ", " + Ellipse[i].Z + ",");
                            file.Write("\n\tVERTEX," + Ellipse[i - 1].X + "," + Ellipse[i - 1].Y + ", " + Ellipse[i - 1].Z + ",");
                            file.Write("\n");
                        }

                    }
                    count = i;
                }




                file.Write("\n");
                file.Write("\n\tEND");
                file.Write("\n]");
                file.Write("\n");
                file.Write("\ncmd.load_cgo(obj, '" + fileName + "')");

            }
        }

        /*
        ** Pre:     
        ** Post:      
        ** About: creates pymol CGO for the major and minor axises of an ellipse
        */
        private void DrawSemiMajorMinor(string fileName, Vector3D centroid, Vector3D A, Vector3D B, Vector3D C)
        {
            string fileLocation = pdbDir + "\\" + fileName + ".py";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                file.Write("\nfrom pymol.cgo import *    # get constants");
                file.Write("\nfrom pymol import cmd");
                file.Write("\n");
                file.Write("\nobj = [");
                file.Write("\n\tBEGIN, LINES,");
                file.Write("\n\tCOLOR, 20.0, 0.0, 0.0,");
                file.Write("\n");

                file.Write("\n\tVERTEX," + centroid.X + "," + centroid.Y + ", " + centroid.Z + ",");
                file.Write("\n\tVERTEX," + (centroid.X + A.X) + "," + (centroid.Y + A.Y) + ", " + (centroid.Z + A.Z) + ",");
                file.Write("\n");


                file.Write("\n\tCOLOR, 0.0, 20.0, 0.0,");
                file.Write("\n\tVERTEX," + centroid.X + "," + centroid.Y + ", " + centroid.Z + ",");
                file.Write("\n\tVERTEX," + (centroid.X + B.X) + "," + (centroid.Y + B.Y) + ", " + (centroid.Z + B.Z) + ",");
                file.Write("\n");



                file.Write("\n\tCOLOR, 0.0, 0.0, 20.0,");
                file.Write("\n\tVERTEX," + centroid.X + "," + centroid.Y + ", " + centroid.Z + ",");
                file.Write("\n\tVERTEX," + (centroid.X + C.X) + "," + (centroid.Y + C.Y) + ", " + (centroid.Z + C.Z) + ",");
                file.Write("\n");

                file.Write("\n");
                file.Write("\n\tEND");
                file.Write("\n]");
                file.Write("\n");
                file.Write("\ncmd.load_cgo(obj, 'SemiMajorMinor" + fileName + "')");

            }
        }


        /*
        ** Pre:     
        ** Post:      
        ** About: creates pymol CGO of an ellipse
        */
        private void furthestDistance(string fileName, string PDBID, ref Ellipse topEllipse, ref Ellipse bottomEllipse)
        {
            int count = topEllipse.ellipsePoints360.Count;
            int largestDisNumTop = 0, largestDisNumBot = 0;
            int smallestDisNumTop = 0, smallestDisNumBot = 0;

                        //remove centroid from ellipse points for good calculations
                        int top0 = 0, bot0 = 0;
                        for (int i = 0; i < topEllipse.ellipsePoints360.Count; i++)
                        {
                            topEllipse.ellipsePoints360[i] -= topEllipse.m_Centroid;
                            bottomEllipse.ellipsePoints360[i] -= bottomEllipse.m_Centroid;
                        }

                                    for (int i = 0; i < topEllipse.ellipsePoints360.Count; i++)
                                    {
                                        if (Math.Floor(topEllipse.ellipsePoints360[i].Y) == 0) if (topEllipse.ellipsePoints360[i].X > 0) top0 = i;
                                        if (Math.Floor(bottomEllipse.ellipsePoints360[i].Y) == 0) if (bottomEllipse.ellipsePoints360[i].X > 0) bot0 = i;
                                    }

                                    int topCtr = top0, botCtr = bot0;

                                    double largestDis = Math.Abs(topEllipse.ellipsePoints360[topCtr].Z + topEllipse.m_Centroid.Z) + Math.Abs(bottomEllipse.ellipsePoints360[botCtr].Z + bottomEllipse.m_Centroid.Z);
                                    double smallestDis = Math.Abs(topEllipse.ellipsePoints360[topCtr].Z + topEllipse.m_Centroid.Z) + Math.Abs(bottomEllipse.ellipsePoints360[botCtr].Z + bottomEllipse.m_Centroid.Z);

                                    for (int i = 0; i < count; i++)
                                    {
                                        if (topCtr >= 360) topCtr -= 360;
                                        if (botCtr >= 360) botCtr -= 360;

                                        double dis = Math.Abs(topEllipse.ellipsePoints360[topCtr].Z + topEllipse.m_Centroid.Z) + Math.Abs(bottomEllipse.ellipsePoints360[botCtr].Z + bottomEllipse.m_Centroid.Z);

                                        if (largestDis < dis)
                                        {
                                            largestDisNumTop = topCtr;
                                            largestDisNumBot = botCtr;
                                            largestDis = dis;
                                        }
                                        if (smallestDis > dis)
                                        {
                                            smallestDisNumTop = topCtr;
                                            smallestDisNumBot = botCtr;
                                            smallestDis = dis;
                                        }
                                        topCtr++;
                                        botCtr++;
                                    }
            //add centroid back to ellipse points
            for (int i = 0; i < topEllipse.ellipsePoints360.Count; i++)
            {
                topEllipse.ellipsePoints360[i] += topEllipse.m_Centroid;
                bottomEllipse.ellipsePoints360[i] += bottomEllipse.m_Centroid;
            }


            string fileLocation = pdbDir + "\\" + fileName + ".py";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(fileLocation))
            {
                output.Write("\nfrom pymol.cgo import *    # get constants");
                output.Write("\nfrom pymol import cmd");
                output.Write("\n");
                output.Write("\nobj = [");
                output.Write("\n\tBEGIN, LINES,");
                output.Write("\n\tCOLOR, 1.0, 1.0, 1.0,");
                output.Write("\n");


                output.Write("\n\tVERTEX," + topEllipse.ellipsePoints360[largestDisNumTop].X + "," + topEllipse.ellipsePoints360[largestDisNumTop].Y + ", " + topEllipse.ellipsePoints360[largestDisNumTop].Z + ",");
                output.Write("\n\tVERTEX," + bottomEllipse.ellipsePoints360[largestDisNumBot].X + "," + bottomEllipse.ellipsePoints360[largestDisNumBot].Y + ", " + bottomEllipse.ellipsePoints360[largestDisNumBot].Z + ",");
                output.Write("\n");

                output.Write("\n\tVERTEX," + topEllipse.ellipsePoints360[smallestDisNumTop].X + "," + topEllipse.ellipsePoints360[smallestDisNumTop].Y + ", " + topEllipse.ellipsePoints360[smallestDisNumTop].Z + ",");
                output.Write("\n\tVERTEX," + bottomEllipse.ellipsePoints360[smallestDisNumBot].X + "," + bottomEllipse.ellipsePoints360[smallestDisNumBot].Y + ", " + bottomEllipse.ellipsePoints360[smallestDisNumBot].Z + ",");
                output.Write("\n");
                output.Write("\n\tCOLOR, 0.0, 100.0, 0.0,");
                output.Write("\n\tVERTEX," + topEllipse.ellipsePoints360[top0].X + "," + topEllipse.ellipsePoints360[top0].Y + ", " + topEllipse.ellipsePoints360[top0].Z + ",");
                output.Write("\n\tVERTEX," + bottomEllipse.ellipsePoints360[bot0].X + "," + bottomEllipse.ellipsePoints360[bot0].Y + ", " + bottomEllipse.ellipsePoints360[bot0].Z + ",");
                output.Write("\n");

                output.Write("\n");
                output.Write("\n\tEND");
                output.Write("\n]");
                output.Write("\n");
                output.Write("\ncmd.load_cgo(obj, '" + fileName + "')");

            }
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        public void EllipseResults(string PDBID)
        {
        string file = OUTPUT_PATH + "ELLIPSE_RESULTS.txt";
            using (System.IO.StreamWriter output = File.AppendText(file))
            {
                //this seemed wrong 7/19/2017
                //if (TopEllipse.m_Normal.Z< 0) TopEllipse.m_Normal *= -1;
                //if (BottomEllipse.m_Normal.Z< 0) BottomEllipse.m_Normal *= -1;

                TopEllipse.findDeviation();
                BottomEllipse.findDeviation();

                output.Write("\n\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}", PDBID, Vector3D.AngleBetween(TopEllipse.m_Normal, BottomEllipse.m_Normal), Eccentricity(TopEllipse), Eccentricity(BottomEllipse), TopEllipse.deviation, BottomEllipse.deviation);

            }
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        private double Eccentricity(Ellipse ellip)
        {
            double ecc = Math.Sqrt(1 - ((ellip.m_minorAxis.Length * ellip.m_minorAxis.Length) / (ellip.m_MajorAxis.Length * ellip.m_MajorAxis.Length)));
            return (ecc);
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        public string WritePDB(string _pdbId)
        {
            string fileName = _pdbId;
            fileName += ".pdb";
           
            FileStream fileStream = new FileStream(Path.Combine(pdbDir, fileName), FileMode.Create, FileAccess.Write);
            StreamWriter fileWriter = new StreamWriter(fileStream);
            string header = "HEADER    " + _pdbId + "                    " + DateTime.Now;
            fileWriter.WriteLine(header);

            try
            {
                string line = "";
                int atomCount = 1;

                foreach (Chain chain in _protein.Chains)
                {
                    foreach (Res res in chain)
                    {
                        foreach (Atom atom in res.Atoms)
                        {
                            line = "ATOM  ";

                            string atomIdStr = atomCount.ToString();
                            line += atomIdStr.PadLeft(5, ' ');
                            line += " ";

                            string atomName = atom.AtomName;
                            if (atomName != "" && atom.AtomType != "H" && atomName.Length < 4)
                            {
                                atomName = " " + atomName;
                            }
                            line += atomName.PadRight(4, ' ');

                            line += " ";
                            line += res.ThreeLetCode;

                            line += " ";
                            line += res.ChainName;

                            line += res.SeqID.ToString().PadLeft(4, ' ');
                            line += "    ";
                            line += FormatDoubleString(atom.Coords.X, 4, 3);
                            line += FormatDoubleString(atom.Coords.Y, 4, 3);
                            line += FormatDoubleString(atom.Coords.Z, 4, 3);

                            line += "  1.00"; //(for dummy occupancy)

                            line += "  0.00"; //(for dummy bfactor)

                            line += "    ";
                            line += atom.AtomType;
                            fileWriter.WriteLine(line);
                            atomCount++;
                        }
                    }
                    fileWriter.WriteLine("END");
                }
            }
            catch (Exception ex)
            {
                string errorMsg = ex.Message;
                throw ex;
            }
            finally
            {
                fileWriter.Close();
            }
            return Path.Combine(pdbDir, fileName);
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        private string FormatDoubleString(double val, int numPre, int numPost)
        {
            string valStr = val.ToString();
            int dotIndex = valStr.IndexOf(".");
            if (dotIndex == -1)
            {
                // return the int part, plus ".0  "
                valStr = valStr.PadLeft(numPre, ' ');
                valStr += ".";
                int i = 0;
                while (i < numPost)
                {
                    valStr += "0";
                    i++;
                }
                return valStr;
            }
            string intPartStr = valStr.Substring(0, dotIndex).PadLeft(numPre, ' ');
            int subStrLen = valStr.Length - dotIndex - 1;
            if (subStrLen > numPost)
            {
                subStrLen = numPost;
            }
            string fractStr = valStr.Substring(dotIndex + 1, subStrLen).PadRight(3, '0');
            return intPartStr + "." + fractStr;
        }

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        public void ELLIPSE_BARREL_PYMOL_SCRIPT(string _pdbId)
        {
            string fileLocation = pdbDir + "loadEllipse.py";
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
            {
                file.WriteLine("cmd.reinitialize('everything')");
                file.WriteLine( "cmd.load('" + _pdbId + ".pdb')" );
                file.WriteLine("execfile('BottomStrands.py')");
                file.WriteLine("execfile('TopStrands.py')");
                file.WriteLine("execfile('strands_" + _pdbId + ".py')");
                file.WriteLine("execfile('Top.py')");
                file.WriteLine("execfile('Bot.py')");
                file.WriteLine("cmd.bg_color('white')");
            }
        }

        /*     
        ** About:    
        */
        class Ellipse
        {
            public Ellipse(List<Vector3D> rawPoints)
            {

                //some algorithm steps from http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/EigenPlane.pdf
                // does not normalize correctly
                _RawPoints = rawPoints;
                List<Vector3D> myEllipse = rawPoints;

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
                A = A * (1.0 / rawPoints.Count);

                SortedList<double, Vector3D> myEigenSolution = new SortedList<double, Vector3D>();
                Evd<double> eigen = A.Evd();

                for (int vecCtr = 0; vecCtr < A.RowCount; vecCtr++)
                {
                    double lambda = eigen.EigenValues.At(vecCtr).Real;
                    //Vector3D vec3D = new Vector3D(eigen.EigenVectors.At(vecCtr, 0), eigen.EigenVectors.At(vecCtr, 1), eigen.EigenVectors.At(vecCtr, 2));
                    Vector3D vec3d = new Vector3D();

                    vec3d.X = eigen.EigenVectors.At(0, vecCtr);
                    vec3d.Y = eigen.EigenVectors.At(1, vecCtr);
                    vec3d.Z = eigen.EigenVectors.At(2, vecCtr);

                    myEigenSolution.Add(lambda, vec3d);
                }

                double eigenValN = myEigenSolution.Keys.Min();
                double eigenValA = myEigenSolution.Keys.Max();
                double eigenValB = myEigenSolution.Keys[1];

                Vector3D EigenvectorA = myEigenSolution[eigenValA];
                Vector3D EigenvectorB = myEigenSolution[eigenValB];
                Vector3D EigenvectorN = myEigenSolution[eigenValN];

                //alter eigenvalues associated with major and minor axis into unit vectors (the normal should already be a unit vector)
                Vector3D unitA = EigenvectorA / EigenvectorA.Length;
                Vector3D unitB = EigenvectorB / EigenvectorB.Length;

                double MajorMagnitude = 4*Math.Sqrt(eigenValA / myEllipse.Count);
                double MinorMagnitude = 4*Math.Sqrt(eigenValB / myEllipse.Count);
                //double MajorMagnitude = (4* Math.Sqrt(eigenValA / myEllipse.Count))/2; //from noline pafe for semi?
                //double MinorMagnitude = (4* Math.Sqrt(eigenValB / myEllipse.Count))/2;


                this.m_MajorAxis = (unitA * MajorMagnitude);
                this.m_minorAxis = (unitB * MinorMagnitude);
                this.m_Normal = EigenvectorN;

                if (double.IsNaN(m_Normal.X)) m_Normal = Vector3D.CrossProduct(m_MajorAxis, m_minorAxis);

                this.m_Centroid = centroid;

                ellipsePoints360 = newEllipseABCtoPoints();
            }

            /*
            public Ellipse(List<Vector3D> rawPoints)
            {

                //some algorithm steps from http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/EigenPlane.pdf
                // does not normalize correctly
                _RawPoints = rawPoints;
                List<Vector3D> myEllipse = rawPoints;

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
                    double lambda = eigen.EigenValues.At(vecCtr).Real;

                    Vector3D vec3D = new Vector3D(eigen.EigenVectors.At(vecCtr, 0), eigen.EigenVectors.At(vecCtr, 1), eigen.EigenVectors.At(vecCtr, 2));

                    Vector3D vec3d = new Vector3D();

                    vec3d.X = eigen.EigenVectors.At(0, vecCtr);
                    vec3d.Y = eigen.EigenVectors.At(1, vecCtr);
                    vec3d.Z = eigen.EigenVectors.At(2, vecCtr);

                    myEigenSolution.Add(lambda, vec3d);

                }

                double eigenValN = myEigenSolution.Keys.Min();
                double eigenValA = myEigenSolution.Keys.Max();
                double eigenValB = myEigenSolution.Keys[1];

                Vector3D EigenvectorA = myEigenSolution[eigenValA];
                Vector3D EigenvectorB = myEigenSolution[eigenValB];
                Vector3D EigenvectorN = myEigenSolution[eigenValN];
                //alter eigenvalues associated with major and minor axis into unit vectors (the normal should already be a unit vector)
                Vector3D unitA = EigenvectorA / EigenvectorA.Length;
                Vector3D unitB = EigenvectorB / EigenvectorB.Length;


                myEigenSolution.Keys.Max();

                double MajorMagnitude = Math.Sqrt(2 * eigenValA / myEllipse.Count);
                double MinorMagnitude = Math.Sqrt(2 * eigenValB / myEllipse.Count);
                //double MajorMagnitude = (4* Math.Sqrt(eigenValA / myEllipse.Count))/2; //from noline pafe for semi?
                //double MinorMagnitude = (4* Math.Sqrt(eigenValB / myEllipse.Count))/2;
                double NormalMagnitude = Math.Sqrt(2 * eigenValN / myEllipse.Count);


                this.m_MajorAxis = (unitA * MajorMagnitude);

                this.m_minorAxis = (unitB * MinorMagnitude);

                this.m_Normal = EigenvectorN * NormalMagnitude;
                if (double.IsNaN(m_Normal.X)) m_Normal = Vector3D.CrossProduct(m_MajorAxis, m_minorAxis);

                this.m_Centroid = centroid;

                ellipsePoints360 = newEllipseABCtoPoints();
            }
            */
            //new version in progress 7-10-17
            //Works for all ellipsis
            public List<Vector3D> newEllipseABCtoPoints()
            {
                List<Vector3D> EllipsePts = new List<Vector3D>();
                //unit vector in direction of major axis
                Vector3D U = new Vector3D( (m_MajorAxis.X/ magnitude(m_MajorAxis)), (m_MajorAxis.Y / magnitude(m_MajorAxis)), (m_MajorAxis.Z / magnitude(m_MajorAxis)));
                //unit vector in direction of minor axis
                Vector3D V = new Vector3D( (m_minorAxis.X/ magnitude(m_minorAxis)), (m_minorAxis.Y / magnitude(m_minorAxis)), (m_minorAxis.Z / magnitude(m_minorAxis)) );
                //center of ellipse
                Vector3D C =  m_Centroid;
                //
                double a = magnitude(2* m_MajorAxis);
                //
                double b = magnitude(2 * m_minorAxis);
                //
                double theta;
                Vector3D newVector;
                for (int i = -180; i < 180; i++)
                {
                    //define t
                    theta = 2 * Math.PI * i / 180.0;
                    //simple and works, thank you NASA
                    newVector = C + (Math.Cos(theta) * m_MajorAxis) + (Math.Sin(theta) * m_minorAxis);
                    //add point to list
                    EllipsePts.Add(newVector);
                }
                return (EllipsePts);
            }

            private double magnitude(Vector3D a)
            {
                return ( Math.Sqrt((a.X * a.X) + (a.Y * a.Y) + (a.Z * a.Z)) );
            }

            public void findDeviation()
            {
                double closestDistance;
                double[] closestDistances = new double[_RawPoints.Count];

                //for each strand end
                for (int i = 0; i < _RawPoints.Count; i++)
                {
                    closestDistance = (_RawPoints[i] - ellipsePoints360[0]).Length;

                    //find the ellipse point closest to that strand end
                    foreach (Vector3D vec in ellipsePoints360)
                    {
                        double temp = (_RawPoints[i] - vec).Length;
                        if (temp < closestDistance) closestDistance = temp;
                    }
                    //add to collection
                    closestDistances[i] = closestDistance;
                }

                //find average of closest ellipse distances to each strand end
                double sum = 0.0;
                foreach (double dub in closestDistances) sum += dub;
                double average = sum / _RawPoints.Count;

                //find deviation of closest ellipse distances to each strand end
                for (int i = 0; i < _RawPoints.Count; i++) closestDistances[i] -= average;
                for (int i = 0; i < _RawPoints.Count; i++) closestDistances[i] = closestDistances[i] * closestDistances[i];
                foreach (double dub in closestDistances) deviation += dub;

                deviation = Math.Sqrt(deviation);
            }


            //ellipse member variables
            public Vector3D m_Centroid;
            public Vector3D m_MajorAxis;
            public Vector3D m_minorAxis;
            public Vector3D m_Normal;
            public List<Vector3D> ellipsePoints360;
            public List<Vector3D> _RawPoints;
            public double deviation;

        }

        static string OUTPUT_PATH;
        static string DB_FILE;
        public string pdbDir;
        public Protein _protein;
        public Barrel _barrel;
        public List<Strand> strandlist;
        private Ellipse TopEllipse;
        private Ellipse BottomEllipse;
        public List<Vector3D> adjustBarrelTop;
        public List<Vector3D> adjustBarrelBottom;
    }

}