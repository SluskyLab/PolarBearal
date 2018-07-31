using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using betaBarrel.AtomParser;
using System.Windows.Media.Media3D;

using System.Collections;
using System.IO;
using System.Xml;
using System.Data;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;

using Meta.Numerics;
using Meta.Numerics.Matrices;




namespace betaBarrelProgram
{

    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }
        public void betaBarrel(object sender, MouseButtonEventArgs e)
        {
            

                //report when starting application 
                Console.WriteLine("Starting at {0}", DateTime.Now);


                Dictionary<string, int> pdbBeta = new Dictionary<string, int>();



                Dictionary<string, AminoAcid> AADict = makeAADict();
                
            
                
                //if making an using program to make output for pymol program start here
                //string fileLocation = @"C:\Users\Joanna\output\OmpDraw.py";
                //using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                //{
                //    file.WriteLine("import pymol");
                //    file.WriteLine("import cmd");
                //    file.WriteLine("from pymol import stored");
                //    file.WriteLine("cmd.alter(\"all\", \"ss=L\")");
                //}


              //  string fileLocation = @"C:\Users\Joanna\output\Omp_may24_2014.txt"; // output file
                string fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\Omp_feb17_2015.txt"; // output file
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                {
                    file.WriteLine("res\ttilt\tstrandTiltt\tbarrelTilt\tcoil1\tcoil2\tcoil3\tphi\tpsi\ttwist_prev\ttwist_next\tinward\trad\tavgRad\tZ\tLoop\tStrand\tlastStrand\tAxisL\tpdb\tseqID\tX\tY\tZ\taDist");
                    //  file.WriteLine("res\tatom\tZ\tpdb\tseqID\tzRes");

                }

               // string fileLocation2 = @"C:\Users\Joanna\output\myShearNumbers.txt"; // shear number output file
                string fileLocation2 = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\Omp_feb17_2015.txt"; // shear number output file
               
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
                {
                    file.WriteLine("PDB\tstrandNum\tshear\tavgAlpha\tavgA\tavgB");


                }
     

                string fileOfPDBs = "07_25_2014.txt"; //input file with list of xml files
               // string fileOfPDBs = "02_09_2015.txt"; //make the LIST of xml polymeric

                if (File.Exists(fileOfPDBs))
                //if (File.Exists(workingDirectory + pdbXMLfile)) currently in bin/debug
                {
                    using (StreamReader sr = new StreamReader(fileOfPDBs))
                    {
                        String line;
                        // Read and display lines from the file until the end of
                        // the file is reached.
                        while ((line = sr.ReadLine()) != null)
                        {
                            string[] splitLine = Array.FindAll<string>(((string)line).Split(
                            new char[] { ' ', '\t', ',' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                             string pdb = splitLine[0];
                            if (pdb != "IDs")
                            {
                                string fileName = pdb + ".xml";
                                //string fileName = pdb + ".pdb";
                                runBetaBarrel(fileName, ref AADict);
                            }
                        }
                    }
                }
                else
                {
                    Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
                    Console.WriteLine("could not open {0}", fileOfPDBs);
                    Console.ReadLine();
                }


                Console.WriteLine("Finished at {0}", DateTime.Now);
                Console.ReadLine();

                return;
            }
                // the next line defines the start of a "region": this allows easy formatting with Visual Studio
               

            #region staticFunctions
            //-----------------------------------Joanna's functions----------------------------


            public static void outputAA(ref Dictionary<string, AminoAcid> _aaDict)
            {

                string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

                outputThisAngleVariable("Tilt", ref _aaDict, true);
                outputThisAngleVariable("caThetaList", ref _aaDict, true);
                outputThisAngleVariable("TiltfromStrand", ref _aaDict, true);
                outputThisAngleVariable("TiltfromBarrel", ref _aaDict, true);
                outputThisAngleVariable("Twist", ref _aaDict, true);

               // string fileLocation = @"C:\Users\Joanna\output\AA_phiPsi.txt";
                string fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_phiPsi.txt";
                //string fileLocation = @"C:\Documents and Settings\joanna\My Documents\Output\AA_PhiPsi.txt";

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
            public static void outputThisAngleVariable(string varName, ref Dictionary<string, AminoAcid> _aaDict, bool isNeg)
            {
                //string fileLocation = @"C:\Users\Joanna\output\AA_" + varName + ".txt";
                string fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_" + varName + ".txt";

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
            public static void outputDividedAA(ref Dictionary<string, AminoAcid> _aaDict)
            {
                string[] AAarray = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

                //string fileLocation = @"C:\Users\Joanna\output\AA_Tilt.txt";
                string fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_Tilt.txt";
                //string fileLocation = @"C:\Documents and Settings\joanna\My Documents\Output\AA_Tilt.txt";

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
              //  fileLocation = @"C:\Users\Joanna\output\AA_Tilt_raw2.txt";

                fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_Tilt_raw2.txt";
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
           //     fileLocation = @"C:\Users\Joanna\output\AA_Tilt_raw5.txt";
                fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_Tilt_raw5.txt";
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

            //public static void runBetaBarrel(string fileName, Dictionary<string, int> pdbBeta, ref Dictionary<string, AminoAcid> _aaDict)
            public static void runBetaBarrel(string fileName, ref Dictionary<string, AminoAcid> _aaDict)
            {

                PeptideAnalysis myModule = new PeptideAnalysis();
                string pdbXMLfile = fileName;//"1I78.xml";//"1KMO.xml";//"1QD6.xml";// "1I78.xml";
                string pdbName = pdbXMLfile.Remove(4);

                AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
                if (Convert.ToString(pdbXMLfile[5]) == Convert.ToString("p")) myAtomCat = readPdbFile(pdbXMLfile);
                else
                {
                    AtomParser.XmlAtomParser myXMLparser = new AtomParser.XmlAtomParser();
                    if (File.Exists(pdbXMLfile))
                    {
                        Console.WriteLine("opened {0}", pdbXMLfile);
                       // using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Joanna\output\AA_Tilt_special.txt", true))

                        using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\AA_Tilt_special.txt", true))
                        {
                            file.WriteLine(pdbXMLfile);
                        }
                        string emptyString = ""; // tells it to pass all atoms
                        myXMLparser.ParseXmlFileAndPassStorage(pdbXMLfile, emptyString, ref myAtomCat);

                    }
                }
                //Protein newProt = new Protein(ref myAtomCat, pdbBeta[pdbName]);
                int chainNum = 0;
                if (pdbName == "3RFZ" || pdbName == "3KVN") chainNum = 1;
                //if (pdbName == " ")
                int stop = myAtomCat.ChainAtomList.Count();
                //chainNum = stop;
                Console.Write("Protein Class {0}",chainNum);

                //**uncomment the below line while using monomeric beta barrels
               Protein newProt = new Protein(ref myAtomCat, chainNum, pdbName);

                //** uncomment the below line while using Polymeric B Barrels
               //PolyProtein newProt = new PolyProtein(ref myAtomCat, 0, pdbName); // For considering only One chain at a time in polymeric beta barrel
               //PolyProtein newProt = new PolyProtein(ref myAtomCat,pdbName); //For cinsidering all chains
               
                                    
                    Console.WriteLine("ChainNum {0}",chainNum);

               // while (newProt.Chains[0].Residues.Count < 136)// for QJ8 and 1QJP
               /* while (chainNum > stop-1)
                {

                    chainNum++;
                    if (myAtomCat.ChainAtomList.Count() <= chainNum) break;
                    else
                    {
                        newProt = new Protein(ref myAtomCat, chainNum, pdbName);

                        Console.WriteLine("{0} Residues", newProt.Chains[0].Residues.Count);
                    }

                }
              //  Console.WriteLine("chainNum", chainNum);
                */

                if (myAtomCat.ChainAtomList.Count() <= chainNum) return;

                //using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Joanna\output\chain_list.txt", true))

                //## WRITE A NEW CLASS polyBarrel to extract the beta Barrel only in the membrane. For now exclude the part after that.....
                Console.Write("creating barrel class..");

                //** uncomment the below line while using monomeric beta barrel
                Barrel myBarrel = new Barrel(newProt.Chains[0]);

                //** uncomment the below line while using polymeric beta barrels, created on 17-03-2015
               // PBarrel myBarrel = new PBarrel(newProt); 

                if (myBarrel.Strands.Count < 8) return;
               // myBarrel.getTiltsByAA(ref _aaDict);
                //myBarrel.getTwist(ref _aaDict);


                ////myBarrel.getTiltsByAA_divided(ref _aaDict);
                //string[] AAarray = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY", "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};

                ////calculateEllipses(getTopEllipseCoords(myBarrel));
                ////calculateEllipses(getBottomEllipseCoords(myBarrel));

                //string fileLocation1 = @"C:\Users\Joanna\output\SSType_";// @"C:\Documents and Settings\joanna\My Documents\Output\colorStrand_";
                ////string fileLocation1 = @"C:\Documents and Settings\joanna\My Documents\Output\colorStrand_";
                //string fileLocation3 = ".txt";
                //string fileLocation = fileLocation1 + pdbName + fileLocation3;

                //using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                //{

                //    for (int residueCtr = 0; residueCtr < newProt.Chains[0].Residues.Count; residueCtr++)
                //    {

                //        file.WriteLine(" residue {0} seqID {1} SSType {2}", newProt.Chains[0].Residues[residueCtr].ResNum, newProt.Chains[0].Residues[residueCtr].SeqID, newProt.Chains[0].Residues[residueCtr].SSType);
                //    }
                //    file.WriteLine("pdbName");
                //}

                //tell me about the strands

                //string fileLocation = @"C:\Users\Joanna\output\OmpDraw" + myBarrel.PdbName + ".py";
                string fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\OmpDraw" + myBarrel.PdbName + ".py";
                
                bool lastStrand = false;
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation))
                {
                    file.WriteLine("import pymol");
                    file.WriteLine("import cmd");
                    file.WriteLine("from pymol import stored");
                    file.WriteLine("cmd.alter(\"all\", \"ss='L'\")");
                    for (int strandCtr = 0; strandCtr < myBarrel.Strands.Count; strandCtr++)
                    {

                        file.WriteLine("cmd.alter(\"chain {2} & resi {0}-{1}\", \"ss ='S'\")", myBarrel.Strands[strandCtr].Residues[0].SeqID, myBarrel.Strands[strandCtr].Residues[myBarrel.Strands[strandCtr].Residues.Count - 1].SeqID, myBarrel.ChainName);
                    }
                    file.WriteLine("cmd.pseudoatom(\"Centroid\", name =\"c\", pos=[{0}, {1}, {2}])", myBarrel.OriginalCcentroid.X, myBarrel.OriginalCcentroid.Y, myBarrel.OriginalCcentroid.Z);
                    file.WriteLine("cmd.pseudoatom(\"Centroid\", name =\"n\", pos=[{0}, {1}, {2}])", myBarrel.OriginalNcentroid.X, myBarrel.OriginalNcentroid.Y, myBarrel.OriginalNcentroid.Z);
                    file.WriteLine("cmd.pseudoatom(\"Centroid\", name =\"c2\", pos=[{0}, {1}, {2}])", myBarrel.OldCaxisPt.X, myBarrel.OldCaxisPt.Y, myBarrel.OldCaxisPt.Z);

                    file.WriteLine("cmd.bond(\"Centroid////c\",\"Centroid////n\")");
                    file.WriteLine("cmd.bond(\"Centroid////c2\",\"Centroid////n\")");

                }


                double sumA = 0;
                double totalA = 0;
                double sumB = 0;
                double totalB = 0;


                //fileLocation = @"C:\Users\Joanna\output\Omp_may24_2014.txt";
                fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\Omp_may24_2014.txt";
                
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))


                    for (int strandCtr = 0; strandCtr < myBarrel.Strands.Count; strandCtr++)
                    {
                        if (strandCtr == myBarrel.Strands.Count - 1) lastStrand = true;
                        for (int resCtr = 0; resCtr < myBarrel.Strands[strandCtr].Residues.Count; resCtr++)
                        {
                            //Res myRes = myBarrel.Strands[strandCtr].Residues[resCtr];
                            //if (myRes.Inward == false && myRes.ThreeLetCode != "TRP")
                            //{
                            //    for (int atomCtr = 0; atomCtr < myBarrel.Strands[strandCtr].Residues[resCtr].Atoms.Count; atomCtr++)
                            //    {
                            //        string myAtName = myBarrel.Strands[strandCtr].Residues[resCtr].Atoms[atomCtr].AtomName;
                            //        if (chargedAtomNames.Contains(myAtName))
                            //        {


                            //            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", myRes.ThreeLetCode, myAtName, -1 * myBarrel.Strands[strandCtr].Residues[resCtr].Atoms[atomCtr].Coords.Z, myBarrel.PdbName, myRes.SeqID, myRes.Z);
                            //        }

                            //    }
                            //}
                            Res myRes = myBarrel.Strands[strandCtr].Residues[resCtr];
                            file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\tFalse\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}", myRes.ThreeLetCode, myRes.Tilt, myBarrel.Strands[strandCtr].AvgTilt, myBarrel.AvgTilt, myRes.Coil1, myRes.Coil2, myRes.Coil3, myRes.Phi, myRes.Psi, myRes.Twist_prev, myRes.Twist_next, myRes.Inward, myRes.Radius, myBarrel.AvgRadius, -1 * (myRes.Z), strandCtr, lastStrand, myBarrel.Axis.Length, myBarrel.PdbName, myRes.SeqID, myRes.BackboneCoords["CA"].X, myRes.BackboneCoords["CA"].Y, myRes.BackboneCoords["CA"].Z, myRes.aDist, myRes.bDist);
                            sumA += myRes.aDist;
                            sumB += myRes.bDist;
                            if (myRes.aDist > 0) totalA++;
                            if (myRes.bDist > 0) totalB++;


                        }

                    }
                // tell me about the loops
                int XtraLoops = 0;
                int periloops = 0;
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))
                    for (int resCtr = 0; resCtr < myBarrel.LoopResies.Count; resCtr++)
                    {
                        Res myRes = myBarrel.LoopResies[resCtr];
                        file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\tTrue\t999\tFalse\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}", myRes.ThreeLetCode, 999, 999, myBarrel.AvgTilt, 999, 999, 999, myRes.Phi, myRes.Psi, 999, 999, myRes.Inward, 999, myBarrel.AvgRadius, -1 * (myRes.Z), myBarrel.Axis.Length, myBarrel.PdbName, myRes.SeqID, myRes.BackboneCoords["CA"].X, myRes.BackboneCoords["CA"].Y, myRes.BackboneCoords["CA"].Z);
                        if (myRes.Z > 0) periloops++;
                        if (myRes.Z < 0) XtraLoops++;
                    }
                if (periloops > XtraLoops)
                {
                    Console.WriteLine("periloops are greater!");
                }

               // fileLocation = @"C:\Users\Joanna\output\myShearNumbers.txt";
                fileLocation = @"C:\Users\Pushpa\Documents\Visual Studio 2013\Projects\betaBarrel_release\output\myShearNumbers.txt";
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation, true))
                {
                    double shear = Math.Tan(Math.PI * myBarrel.AvgTilt / 180) * myBarrel.Strands.Count * (sumB / totalB) / (sumA / totalA);
                    file.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", myBarrel.PdbName, myBarrel.Strands.Count, shear, myBarrel.AvgTilt, sumA / totalA, sumB / totalB);
                }


            }
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
                        myEllipse.Add(lastCA);


                    }
                    else
                    {
                        myEllipse.Add(firstCA);

                    }

                }
                return myEllipse;

            }
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
                        myEllipse.Add(firstCA);


                    }
                    else
                    {
                        myEllipse.Add(lastCA);

                    }

                }
                return myEllipse;

            }

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
            public static void calculateEllipses(List<Vector3D> myEllipse)
            {
                //algorithm from http://www.ahinson.com/algorithms/Sections/InterpolationRegression/EigenPlane.pdf
                //get points for two ellipsies

                Vector3D centroid = new Vector3D();


                for (int strandCtr = 0; strandCtr < myEllipse.Count; strandCtr++)
                {

                    centroid += myEllipse[strandCtr] / myEllipse.Count;
                }


                // center the points
                for (int strandCtr = 0; strandCtr < myEllipse.Count; strandCtr++)
                {
                    myEllipse[strandCtr] = myEllipse[strandCtr] - centroid;

                }

                //set up covariance matrix
                SquareMatrix covarMatrix = new SquareMatrix(3);

                for (int strandCtr = 0; strandCtr < myEllipse.Count; strandCtr++)
                {
                    covarMatrix[0, 0] = covarMatrix[0, 0] + (myEllipse[strandCtr].X * myEllipse[strandCtr].X);
                    covarMatrix[0, 1] = covarMatrix[0, 1] + (myEllipse[strandCtr].X * myEllipse[strandCtr].Y);
                    covarMatrix[0, 2] = covarMatrix[0, 2] + (myEllipse[strandCtr].X * myEllipse[strandCtr].Z);

                    covarMatrix[1, 0] = covarMatrix[1, 0] + (myEllipse[strandCtr].Y * myEllipse[strandCtr].X);
                    covarMatrix[1, 1] = covarMatrix[1, 1] + (myEllipse[strandCtr].Y * myEllipse[strandCtr].Y);
                    covarMatrix[1, 2] = covarMatrix[1, 2] + (myEllipse[strandCtr].Y * myEllipse[strandCtr].Z);

                    covarMatrix[2, 0] = covarMatrix[2, 0] + (myEllipse[strandCtr].Z * myEllipse[strandCtr].X);
                    covarMatrix[2, 1] = covarMatrix[2, 1] + (myEllipse[strandCtr].Z * myEllipse[strandCtr].Y);
                    covarMatrix[2, 2] = covarMatrix[2, 2] + (myEllipse[strandCtr].Z * myEllipse[strandCtr].Z);
                }


                ComplexEigensystem myEigensystem = covarMatrix.Eigensystem();

                Vector3D EigenvectorA = new Vector3D();
                Vector3D EigenvectorB = new Vector3D();
                Vector3D EigenvectorN = new Vector3D();
                SortedList<double, Vector3D> myEigenSolution = new SortedList<double, Vector3D>();

                for (int vecCtr = 0; vecCtr < myEigensystem.Dimension; vecCtr++)
                {
                    Complex[] vec = myEigensystem.Eigenvector(vecCtr);
                    double val = myEigensystem.Eigenvalue(vecCtr).Re;

                    Vector3D vec3d = new Vector3D();
                    if (vec.Count<Complex>() == 3)
                    {
                        vec3d.X = vec[0].Re;
                        vec3d.Y = vec[1].Re;
                        vec3d.Z = vec[2].Re;

                    }
                    myEigenSolution.Add(val, vec3d);
                }

                double eigenValN = myEigenSolution.Keys.Min();
                EigenvectorN = myEigenSolution[eigenValN];
                double eigenValA = myEigenSolution.Keys.Max();
                EigenvectorA = myEigenSolution[eigenValA];
                double eigenValB = myEigenSolution.Keys[1];
                EigenvectorB = myEigenSolution[eigenValB];

                Vector3D normalizedA = new Vector3D();
                Vector3D normalizedB = new Vector3D();

                normalizedA = EigenvectorA / EigenvectorA.Length;
                normalizedB = EigenvectorB / EigenvectorB.Length;

                double scaleA = 4 * Math.Sqrt(eigenValA / myEllipse.Count);
                double scaleB = 4 * Math.Sqrt(eigenValB / myEllipse.Count);

                Vector3D axisA = new Vector3D();
                axisA = normalizedA * scaleA;

                Vector3D axisB = new Vector3D();
                axisB = normalizedB * scaleB;

                double lengthA = axisA.Length;
                double lengthB = axisB.Length;

                Console.WriteLine("centroid is {0},{1},{2}, and the major axis is {3},{4},{5}, and the minor axis is {6},{7},{8}", centroid.X, centroid.Y, centroid.Z, axisA.X, axisA.Y, axisA.Z, axisB.X, axisB.Y, axisB.Z);



            }

            public static double DotProduct(Vector3D vec1, Vector3D vec2)
            {
                return System.Windows.Media.Media3D.Vector3D.DotProduct(vec1, vec2);

            }
            public static Vector3D CrossProduct(Vector3D vec1, Vector3D vec2)
            {
                return System.Windows.Media.Media3D.Vector3D.CrossProduct(vec1, vec2);
            }
            public static double AngleBetween(Vector3D vec1, Vector3D vec2)
            {
                return System.Windows.Media.Media3D.Vector3D.AngleBetween(vec1, vec2);
            }
            public static Vector3D Multiply(SquareMatrix Mat1, Vector3D vec1)
            {
                Vector3D product = new Vector3D();
                product.X = (Mat1[0, 0] * vec1.X) + (Mat1[0, 1] * vec1.Y) + (Mat1[0, 2] * vec1.Z);
                product.Y = (Mat1[1, 0] * vec1.X) + (Mat1[1, 1] * vec1.Y) + (Mat1[1, 2] * vec1.Z);
                product.Z = (Mat1[2, 0] * vec1.X) + (Mat1[2, 1] * vec1.Y) + (Mat1[2, 2] * vec1.Z);

                return product;
            }




            public static Vector3D averagePosition(List<Vector3D> AllCaxyzCoord)
            {
                Vector3D averagePosition = new Vector3D();

                foreach (Vector3D position in AllCaxyzCoord)
                {
                    averagePosition += position / AllCaxyzCoord.Count;

                }
                return averagePosition;
            }

            //public static void calculateCaAngles(ref AtomParser.AtomInfo[] _chainYouWant)
            //{
            //    int calphaCtr = 0;
            //    List<int> caList = new List<int>();

            //    for (int atomC = 0; atomC < _chainYouWant.Length; atomC++)
            //    {
            //        if (_chainYouWant[atomC].atomName == "CA")
            //        {
            //            caList.Add(atomC);
            //            if (calphaCtr > 1)
            //            {
            //                int atomB = caList[calphaCtr - 1];
            //                int atomA = caList[calphaCtr - 2];


            //                var xyzCoordA = new DenseVector(coordArray(_chainYouWant, atomA));
            //                var xyzCoordB = new DenseVector(coordArray(_chainYouWant, atomB));
            //                var xyzCoordC = new DenseVector(coordArray(_chainYouWant, atomC));
            //                var AtoBVec = new DenseVector(xyzCoordA - xyzCoordB);
            //                var BtoCVec = new DenseVector(xyzCoordB - xyzCoordC);



            //                double dotABC = AtoBVec.DotProduct(BtoCVec);

            //                //     Console.WriteLine("angle around residue {0} is {1}.", calphaCtr - 1, Math.Acos(AtoBVec.DotProduct(BtoCVec) / Math.Sqrt(AtoBVec.DotProduct(AtoBVec) * BtoCVec.DotProduct(BtoCVec))) * 180 / Math.PI);


            //            }
            //            calphaCtr++;
            //        }

            //    }
            //}
            public static Vector3D coordVector(AtomParser.AtomInfo[] chainYouWant, int atomNum)
            {
                Vector3D coords = new Vector3D();
                coords.X = (chainYouWant[atomNum].xyz.X);
                coords.Y = (chainYouWant[atomNum].xyz.Y);
                coords.Z = (chainYouWant[atomNum].xyz.Z);
                return coords;


            }

            public static AtomCategory readPdbFile(string pdbfilename)
            {
                AtomCategory myAtomCategory = new AtomCategory();
                ChainAtoms myChainAtoms = new ChainAtoms();
                ArrayList mychainAtoms = new ArrayList();

                string currentChain = "A";
                if (File.Exists(pdbfilename))
                {
                    using (StreamReader sr = new StreamReader(pdbfilename))
                    {
                        String line;
                        // Read and display lines from the file until the end of
                        // the file is reached.
                        Int32 seqID = -1;
                        Int32 previousSeqID = -1;
                        while ((line = sr.ReadLine()) != null)
                        {

                            string[] splitLine = Array.FindAll<string>(((string)line).Split(
                            new char[] { ' ', '\t', ',', ';' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                            if (splitLine[0] == "ATOM")
                            {

                                if (splitLine[4] != currentChain) Console.WriteLine("pdb file has more than 1 chain");
                                currentChain = splitLine[4];
                                AtomInfo myAtom = new AtomInfo();
                                myAtom.altConfID = "";
                                myAtom.atomId = Convert.ToInt32(splitLine[1]);
                                myAtom.atomName = splitLine[2];
                                myAtom.residue = splitLine[3];
                                myAtom.authResidue = splitLine[3];
                                myAtom.authSeqId = splitLine[5];
                                if (Convert.ToInt32(myAtom.authSeqId) != previousSeqID)
                                {
                                    // seqID++;
                                    seqID = seqID + 1;
                                }

                                myAtom.seqId = seqID.ToString();
                                previousSeqID = Convert.ToInt32(myAtom.authSeqId);



                                AtomParser.Coordinate coord = new Coordinate(Convert.ToDouble(splitLine[6]), Convert.ToDouble(splitLine[7]), Convert.ToDouble(splitLine[8]));
                                myAtom.xyz = coord;
                                myAtom.atomType = splitLine[10];

                                myAtom.bFac = Convert.ToDouble(splitLine[9]);
                                myChainAtoms.AddAtom(myAtom);
                            }
                        }

                    }
                }

                myChainAtoms.asymChain = currentChain;
                myChainAtoms.authAsymChain = currentChain;
                myAtomCategory.AddChainAtoms(myChainAtoms);
                //myAtomCategory.chainAtomsList = mychainAtoms;

                myAtomCategory.Resolution = 2.5;


                return myAtomCategory;
            }

            public static string ConvertConfToLargeAlphabetGeoSeq(ArrayList _dihArray)
            {
                string bigPValue = "P";
                string littlePValue = "p";
                string geoSeq = "";
                for (int resPosition = 0; resPosition < _dihArray.Count; resPosition = resPosition + 3)
                {
                    double phiValue = (double)_dihArray[resPosition];
                    double psiValue = (double)_dihArray[resPosition + 1];
                    double omgValue = (double)_dihArray[resPosition + 2];
                    // possible values of geo hash:
                    // trans values: lowercase
                    // phi < 0:
                    //  alpha-R (A): -100 <= psi <= 50
                    //   if not: beta (B) if phi < -100, polyproII (P) if phi >= -100
                    // phi > 0:
                    //  alpha-L (L): -50 <= psi <= 100 : L
                    //  beta(B): phi > 150 and not alpha-L : G (gamma) otherwise
                    bool transPeptide = new bool();
                    transPeptide = true;
                    if (Math.Abs(omgValue) < 90)
                    { transPeptide = false; }
                    if (phiValue < 0)
                    {
                        if (psiValue >= -100 && psiValue <= 50)
                        { // alphaR or D
                            if (phiValue < -100)
                            {// D
                                if (transPeptide) { geoSeq += "D"; }
                                else { geoSeq += "d"; }
                            }
                            else
                            {
                                if (transPeptide) { geoSeq += "A"; }
                                else { geoSeq += "a"; }
                            }
                        }
                        else // beta or polyproII
                        {
                            if (phiValue < -100)
                            {
                                if (transPeptide) { geoSeq += "B"; }
                                else { geoSeq += "b"; }
                            }
                            else // polyproII
                            {
                                if (transPeptide) { geoSeq += bigPValue; } // switch to B to merge
                                else { geoSeq += littlePValue; } // switch to b to merge
                            }
                        }
                    }
                    else // phi >= 0
                    {
                        if (psiValue >= -50 && psiValue <= 100)
                        { // alpha-L
                            if (transPeptide) { geoSeq += "L"; }
                            else { geoSeq += "l"; }
                        }
                        else
                        {
                            if (phiValue > 150)
                            {  //extension of beta
                                if (transPeptide) { geoSeq += "B"; }
                                else { geoSeq += "b"; }
                            }
                            else // gamma conformation
                            {
                                if (transPeptide) { geoSeq += "G"; }
                                else { geoSeq += "g"; }
                            }
                        }
                    }
                }
                return geoSeq;
            }

            #endregion
            #region externalClassDefs
            public class TurnListSorter : IComparer<List<string>>
            {
                public int Compare(List<string> _l1, List<string> _l2)
                {
                    if (_l1.Count != _l2.Count)
                    { return _l2.Count.CompareTo(_l1.Count); } // returns in descending order
                    return _l1[0].CompareTo(_l2[0]); // returns in alphabetical order if equal counts
                }
            }
            public class LKandCLIDSorter : IComparer<string> // format: loopKey_clusterID
            {
                public int Compare(string _s1, string _s2)
                {
                    string lk1 = _s1.Substring(0, _s1.IndexOf("_"));
                    int clID1 = Convert.ToInt32(_s1.Substring(_s1.IndexOf("_") + 1));
                    string lk2 = _s2.Substring(0, _s2.IndexOf("_"));
                    int clID2 = Convert.ToInt32(_s2.Substring(_s2.IndexOf("_") + 1));
                    if (lk1 != lk2)
                    { return lk1.CompareTo(lk2); }
                    else
                    { return clID1.CompareTo(clID2); }
                }
            }

            public class SupplDataSorter : IComparer<string>
            {
                public int Compare(string _s1, string _s2)
                {
                    string[] splitIntoWords1 = Array.FindAll<string>(_s1.Split(
                            new char[] { ' ', '\t', ',' }), delegate(string s)
                            {
                                return !String.IsNullOrEmpty(s);
                            });
                    string[] splitIntoWords2 = Array.FindAll<string>(_s2.Split(
                            new char[] { ' ', '\t', ',' }), delegate(string s)
                            {
                                return !String.IsNullOrEmpty(s);
                            });
                    // [0] is loopID with length, [1] is clusterID, [2] is pdb and chainID, [6] is seq
                    if (splitIntoWords1[0] != splitIntoWords2[0]) // loopIDs not equal
                    {
                        if (splitIntoWords1[0].Substring(0, 1) != splitIntoWords2[0].Substring(0, 1))
                        { return splitIntoWords2[0].CompareTo(splitIntoWords1[0]); } // returns L before H
                        if (splitIntoWords1[0].Substring(1, 1) != splitIntoWords2[0].Substring(1, 1))
                        { return splitIntoWords1[0].CompareTo(splitIntoWords2[0]); } // returns L1 before L2
                        // if here, length is only difference
                        return splitIntoWords1[6].Length.CompareTo(splitIntoWords2[6].Length); //returns L1-10 before L1-11
                    }
                    if (splitIntoWords1[1] != splitIntoWords2[1]) // clusterIDs not equal
                    { return splitIntoWords1[1].CompareTo(splitIntoWords2[1]); }
                    return splitIntoWords1[2].CompareTo(splitIntoWords2[2]); // sorting pdbIDs in alphabetical order
                }
            }
            public class myUniqueSequenceContainer
            {
                public myUniqueSequenceContainer()
                {
                    uniqueSeqObject = null;
                    uniqueSeqCount = -1;
                }
                public myUniqueSequenceContainer(string _seq)
                {
                    uniqueSeqObject = (object)_seq;
                    uniqueSeqCount = 1;
                }
                public myUniqueSequenceContainer(string _seq, int _count)
                {
                    uniqueSeqObject = (object)_seq;
                    uniqueSeqCount = _count;
                }
                // memberVariables
                public object uniqueSeqObject = new object();
                public int uniqueSeqCount = new int();
            }
            public class myUniqueSeqContainerSorter : IComparer<myUniqueSequenceContainer>
            {
                public int Compare(myUniqueSequenceContainer _s1, myUniqueSequenceContainer _s2)
                {
                    if (_s1.uniqueSeqCount != _s2.uniqueSeqCount)
                    { return _s2.uniqueSeqCount.CompareTo(_s1.uniqueSeqCount); }
                    else
                    {
                        string seq1 = (string)_s1.uniqueSeqObject;
                        string seq2 = (string)_s2.uniqueSeqObject;
                        return seq1.CompareTo(seq2);
                    }
                }
            }
            public class myGeoSeqPair
            {
                public myGeoSeqPair()
                {
                    geoSeqObject = null;
                    geoSeqCount = -1;
                }
                public myGeoSeqPair(string _geoSeq, int _count)
                {
                    geoSeqObject = (object)_geoSeq;
                    geoSeqCount = _count;
                }
                // member variables
                public object geoSeqObject = new object();
                public int geoSeqCount = new int();
            }
            public class myGeoSeqPairSorter : IComparer<myGeoSeqPair>
            {
                public int Compare(myGeoSeqPair _g1, myGeoSeqPair _g2)
                {
                    if (_g1.geoSeqCount != _g2.geoSeqCount)
                    { return _g2.geoSeqCount.CompareTo(_g1.geoSeqCount); } // descending count order
                    else
                    {
                        string g1str = (string)_g1.geoSeqObject;
                        string g2str = (string)_g2.geoSeqObject;
                        return g1str.CompareTo(g2str); // if equal count, ascending alphabetical order
                    }
                }
            }
            public class myQualityObj
            {
                public myQualityObj()
                {
                    myPdbIDobj = null;
                    myResolution = 999;
                    myBfac = 999;
                    myConfE = 999;
                }
                public myQualityObj(string _pdbID, double _res, double _bfac, double _confE)
                {
                    myPdbIDobj = (object)_pdbID;
                    myResolution = _res;
                    myBfac = _bfac;
                    myConfE = _confE;
                }

                // memberVariables
                public object myPdbIDobj = new object();
                public double myResolution = new double();
                public double myBfac = new double();
                public double myConfE = new double();
            }

            public class myQualitySorter : IComparer<myQualityObj>
            {
                public int Compare(myQualityObj _q1, myQualityObj _q2)
                {
                    // compare resolution, then bfactor, then confE
                    double theTolerance = (double)0.05;
                    if (Math.Abs(_q1.myResolution - _q2.myResolution) > theTolerance)
                    { return _q1.myResolution.CompareTo(_q2.myResolution); } // sort in ascending order
                    if (Math.Abs(_q1.myBfac - _q2.myBfac) > theTolerance)
                    { return _q1.myBfac.CompareTo(_q2.myBfac); } // ditto
                    if (Math.Abs(_q1.myConfE - _q2.myConfE) > theTolerance)
                    { return _q1.myConfE.CompareTo(_q2.myConfE); } // ditto

                    return 0;
                }
            }

            public class myDihPDBpair
            {
                public myDihPDBpair()
                {
                    myPDBID = "XXXX";
                    myDihedralValue = 999;
                }

                public myDihPDBpair(string _pdbID, double _dihValue)
                {
                    myPDBID = _pdbID;
                    myDihedralValue = _dihValue;
                }
                // member variables
                public string myPDBID;
                public double myDihedralValue = new double();
            }

            public class myDihSorter : IComparer<myDihPDBpair>
            {
                public int Compare(myDihPDBpair _p1, myDihPDBpair _p2)
                {
                    return _p1.myDihedralValue.CompareTo(_p2.myDihedralValue);
                }
            }
            public class myClusterGroupSorter : IComparer<List<int>>
            {
                public int Compare(List<int> _a1, List<int> _a2)
                {
                    int inta1 = _a1[0];
                    int inta2 = _a2[0];
                    return inta1.CompareTo(inta2);
                }
            }
            public class myLoopKeySorter : IComparer<string>
            {
                public int Compare(string _l1, string _l2)
                {
                    string loopType1 = _l1.Substring(0, 1);
                    string loopType2 = _l2.Substring(0, 1);
                    if (loopType1 != loopType2)
                    { return _l2.CompareTo(_l1); } // should rank L before H
                    string numeral1 = _l1.Substring(1, 1);
                    string numeral2 = _l2.Substring(1, 1);
                    if (numeral1 != numeral2)
                    { return (Convert.ToInt32(numeral1)).CompareTo(Convert.ToInt32(numeral2)); }
                    // next, find the length string and compare
                    numeral1 = _l1.Substring(3, _l1.IndexOf("-", 3) - 3);
                    numeral2 = _l2.Substring(3, _l2.IndexOf("-", 3) - 3);
                    if (numeral1 != numeral2)
                    { return (Convert.ToInt32(numeral1)).CompareTo(Convert.ToInt32(numeral2)); }
                    return _l1.CompareTo(_l2); // this handles cis-trans hash            
                }
            }
            public class myTurnDataObject
            {
                public myTurnDataObject()
                {
                    myTurnNameObj = (object)"empty";
                    myTurnTurnIDObj = (object)"notDef";
                    myClusterID = -1;
                }
                public myTurnDataObject(string _name, string _loopKey, string _turnID, int _clID, ArrayList _dih)
                {
                    myTurnNameObj = (object)_name;
                    myLoopKeyObj = (object)_loopKey;
                    myTurnTurnIDObj = (object)_turnID;
                    myClusterID = _clID;
                    myDihedrals = _dih;
                }
                // member variables
                public object myTurnNameObj = new object();
                public object myLoopKeyObj = new object();
                public object myTurnTurnIDObj = new object();
                public int myClusterID = new int();
                public ArrayList myDihedrals = new ArrayList();
            }
            #endregion
        
            
        }

       



        }

    

    

        
    


    

