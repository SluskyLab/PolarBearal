using System;
using System.Collections.Generic;
using System.Linq;
using betaBarrel.AtomParser;
using System.Collections;
using System.IO;

namespace betaBarrel
{
    class Program
    {
        //Change these locations and it will propagate to all other locations

        public static string parameterFile = @"Z:\Documents\PhD\SluskyLab\PolarBearal\par_hbond_1.txt"; //This parameter file holds donor/acceptors in every amino acid
        
        public static string polyDBDirectory = @"Z:\Documents\PhD\SluskyLab\PolyBarrelsDB\";
        public static string MacpolyDBDir = "/Users/meghan/Documents/PhD/SluskyLab/PolyBarrelsDB/";
        //public static string polyDBDirectory = @"Z:\Documents\PhD\SluskyLab\RenumbFrom1PDBs\";
        //public static string MacpolyDBDir = "/Users/meghan/Documents/PhD/SluskyLab/RenumbFrom1PDBs/";
        //public static string polyoutDirectory = @"Z:\Documents\PhD\SluskyLab\Output\PolyBarrels\";
        //public static string MacpolyDBDir = @"Z:\Dropbox (Slusky Lab)\MFranklin\StrucAlignTempCopy\v4_July\";
        //public static string polyDBDirectory = "/Volumes/G_Drive/StructureAlign/NirBenTal/Joanna_and_Meghan/v4_July/pdb_files/";
        //public static string polyoutDirectory = @"Z:\Dropbox (Slusky Lab)\MFranklin\StrucAlignTempCopy\v4_July\";
        public static string polyoutDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\";
        //public static string polyoutDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\RosettaLoops\";

        public static string monoDBDirectory = @"Z:\Documents\PhD\SluskyLab\MonoDB\";
        public static string MacmonoDBDir = "/Users/meghan/Documents/PhD/SluskyLab/MonoDB/";
        public static string monooutDirectory = @"Z:\Documents\PhD\SluskyLab\Output\monomers_v5\";
        //Remember to change lines to seqID for rel num or to ResNum to abs num in getLoopSeq function
        //public static string monooutDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\";
        //public static string monooutDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\RosettaLoops\";

        //public static string monoDBDirectory = @"Z:\Documents\PhD\v6_2018\OMBBNetwork\PDBFiles/";
        //public static string MacmonoDBDir = "/Users/meghan/Documents/PhD/v6_2018/OMBBNetwork/PDBFiles/";
        //public static string monooutDirectory = @"Z:\Documents\PhD\v6_2018\OMBBNetwork\";

        //public static string monoDBDirectory = @"Z:\Documents\PhD\SluskyLab\8-12Scaffolds\ForPBearal";
        //public static string MacmonoDBDir = "/Users/meghan/Documents/PhD/SluskyLab/8-12Scaffolds/ForPBearal";
        //public static string monooutDirectory = @"Z:\Documents\PhD\SluskyLab\8-12Scaffolds\PBearalData\";

        public static string solDBDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\LoopTwists\SolubleTwists\ThorntonDBMod\";
        public static string MacsolDBDir = "/Users/meghan/Documents/PhD/SluskyLab/BarrelChars/LoopTurns/LoopTwists/SolubleTwists/ThorntonDBMod/";
        public static string soloutDirectory = @"Z:\Documents\PhD\SluskyLab\BarrelChars\LoopTurns\LoopTwists\SolubleTwists\";

        public static Dictionary<string, AminoAcid> AADict = SharedFunctions.makeAADict();
        //The values in this dictionary are transcribed from CHARMM36 all-hydrogen topology file for proteins, May 2011
        public static Dictionary<Tuple<string, string>, double> partialChargesDict = new Dictionary<Tuple<string, string>, double>
                {
                    {new Tuple<string, string>("ARG", "HE"), 0.44}, 
                    {new Tuple<string, string>("ARG", "NE"), -0.70}, 
                    {new Tuple<string, string>("ARG", "NH1"), -0.80}, 
                    {new Tuple<string, string>("ARG", "HH11"), 0.46}, 
                    {new Tuple<string, string>("ARG", "HH12"), 0.46}, 
                    {new Tuple<string, string>("ARG", "NH2"), -0.80}, 
                    {new Tuple<string, string>("ARG", "HH21"), 0.46}, 
                    {new Tuple<string, string>("ARG", "HH22"), 0.46},

                    {new Tuple<string, string>("ASN", "ND2"), -0.62}, 
                    {new Tuple<string, string>("ASN", "HD21"), 0.32}, 
                    {new Tuple<string, string>("ASN", "HD22"), 0.40},
                    {new Tuple<string, string>("ASN", "CG"), 0.55}, 
                    {new Tuple<string, string>("ASN", "OD1"), -0.55},

                    {new Tuple<string, string>("ASP", "CG"), 0.62}, 
                    {new Tuple<string, string>("ASP", "OD1"), -0.76},
                    {new Tuple<string, string>("ASP", "OD2"), -0.76},

                    {new Tuple<string, string>("CYS", "HG"), 0.16}, 
                    {new Tuple<string, string>("CYS", "SG"), -0.23},

                    {new Tuple<string, string>("GLN", "NE2"), -0.62}, 
                    {new Tuple<string, string>("GLN", "HE21"), 0.32}, 
                    {new Tuple<string, string>("GLN", "HE22"), 0.30},
                    {new Tuple<string, string>("GLN", "CD"), 0.55}, 
                    {new Tuple<string, string>("GLN", "OE1"), -0.55},

                    {new Tuple<string, string>("GLU", "CD"), 0.62}, 
                    {new Tuple<string, string>("GLU", "OE1"), -0.76},
                    {new Tuple<string, string>("GLU", "OE2"), -0.76},

                    {new Tuple<string, string>("GLY", "O"), -0.51},

                    {new Tuple<string, string>("HIS", "NE2"), -0.36}, 
                    {new Tuple<string, string>("HIS", "HE2"), 0.32}, 
                    {new Tuple<string, string>("HIS", "ND1"), -0.70},

                    {new Tuple<string, string>("LYS", "NZ"), -0.30}, 
                    {new Tuple<string, string>("LYS", "HZ1"), 0.33},
                    {new Tuple<string, string>("LYS", "HZ2"), 0.33},
                    {new Tuple<string, string>("LYS", "HZ3"), 0.33},

                    {new Tuple<string, string>("SER", "HG1"), 0.43}, 
                    {new Tuple<string, string>("SER", "OG1"), -0.66},
 
                    {new Tuple<string, string>("THR", "HG1"), 0.43}, 
                    {new Tuple<string, string>("THR", "OG1"), -0.66}, 

                    {new Tuple<string, string>("TRP", "NE1"), -0.51}, 
                    {new Tuple<string, string>("TRP", "HE1"), 0.37},

                    {new Tuple<string, string>("TYR", "HH"), 0.43}, 
                    {new Tuple<string, string>("TYR", "OH"), -0.54}

                };

        static public void menu()
        {
            Console.WriteLine("1. Run PolarBearal");
            Console.WriteLine("2. Run Ellipse");
            Console.WriteLine("3. Run Ellipse Database");
            Console.WriteLine("4. Run Soluble Database");
            Console.WriteLine("5. Run single Protein");
            Console.WriteLine("6. Check out organisms");
            Console.WriteLine("8. Run PBarrel");
            Console.WriteLine("7. ");

        }

        static void Main(string[] args)
        {
            DateTime startTime = DateTime.Now;
            string choice = "";
            while (choice != "7")
            {
                menu();
                choice = Console.ReadLine();
                switch (choice)
                {
                    case "1":
                        PolarBearal.RunPolarBearal();
                        break;
                    case "2":
                        BarrelEllipse.RollTide();
                        break;
                    case "4":
                        startSolubleBarrel();
                        break;
                    case "5":
                        runProtein();
                        break;
                    case "6":
                        getOrganism newO = new getOrganism();
                        newO.runOrganism();
                        break;
                    case "8":
                        startPolyBarrel();
                        break;
                    default:
                        choice = "7";
                        break;
                }
            }


            Console.WriteLine("Program took from: \n start: {0} \n   end: {1}", startTime, DateTime.Now);

            return;
        }

        public static void runProtein()
        {
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();
            Dictionary<string, AminoAcid> AADict = SharedFunctions.makeAADict();

            Console.Write("\nPDBid: ");
            string input = Console.ReadLine();
            //string fileName = input;

            Program.runBetaBarrel(input, ref AADict, ref partialChargesDict);
            
        }

        public static Barrel runBetaBarrel(string fileName, ref Dictionary<string, AminoAcid> _aaDict, ref Dictionary<Tuple<string, string>, double> partialChargesDict)
            {
                //new peptide analysis
                //PeptideAnalysis myModule = new PeptideAnalysis();
                Directory.SetCurrentDirectory(monoDBDirectory);

                //insert name of xml file for doing things
                string pdbXMLfile = fileName + ".xml";
                string pdbName = fileName;
                AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();

                //Test to see if there is an xml file for this pdb ID name; if not, attempt to use the pdb file instead
                if (File.Exists(pdbXMLfile))
                {
                    AtomParser.XmlAtomParser myXMLparser = new AtomParser.XmlAtomParser();

                    {
                        Console.WriteLine("opened {0}", pdbXMLfile);
                        string emptyString = ""; // tells it to pass all atoms
                        myXMLparser.ParseXmlFileAndPassStorage(pdbXMLfile, emptyString, ref myAtomCat);
                    }
                }
                else
                {
                    string fileName2 = pdbName + ".pdb";
                    Console.WriteLine("opened {0}", fileName2);
                    myAtomCat = readPdbFile(fileName2, ref partialChargesDict);
                }
            pdbName = pdbName.Substring(0, 4).ToUpper();

                int chainNum = 0;
                //if (pdbName == "3RFZ" || pdbName == "3KVN") chainNum = 1;

                int stop = myAtomCat.ChainAtomList.Count();
                Console.Write("Protein Class {0}", chainNum);

                Protein newProt = new Protein(ref myAtomCat, chainNum, pdbName, monoDBDirectory);

                //Console.WriteLine("ChainNum {0}", chainNum);

                Console.Write("creating barrel class..");

                Barrel myBarrel = new Barrel(newProt.Chains[0], newProt);

                //if (myBarrel.Strands.Count < 8) return;
            // myBarrel.getTiltsByAA(ref _aaDict);
            //myBarrel.getTwist(ref _aaDict);


            //PolarBearal polarTest = new PolarBearal(ref myBarrel);

            return (myBarrel);
            }

        static public void startPolyBarrel()
        {
            //for making barrel
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();
            
            //PDBIDs that will be run
            string fileOfPDBs = polyDBDirectory + "PolyDBList_v4.txt"; //input file with list of polymeric xml files
            //string fileOfPDBs = polyoutDirectory + "OtherBarrelsList.txt";
    
            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    string fileLocation2 = polyoutDirectory + "AllBarrelChar.txt";
                    using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
                    {
                        string newLine = "PDB" + "\t\t" + "Total Chains" + "\t" + "Total Strands" + "\t" + "Length" + "\t" + "AvgLength" + "\t" + "MinLength" + "\t" + "MaxLength" + "\t" + "Radius" + "\t" + "Barrel Tilt";
                        file.WriteLine(newLine);
                    // Read and display lines from the file until the end of the file is reached.
                        while ((line = sr.ReadLine()) != null)
                        {
                            string[] splitLine = Array.FindAll<string>(((string)line).Split(
                                new char[] { ' ', '\t', ',' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                            string pdb = splitLine[0];
                            if (pdb != "IDs")
                            {
                                PBarrel myBarrel = runPBetaBarrel(pdb, ref AADict);
                                /*string char1 = myBarrel.PdbName;
                                string char2 = myBarrel.Axis.Length.ToString();
                                string char7 = myBarrel.StrandLength.Average().ToString();
                                string char8 = myBarrel.StrandLength.Min().ToString();
                                string char9 = myBarrel.StrandLength.Max().ToString();
                                string char3 = myBarrel.AvgRadius.ToString();
                                string char4 = myBarrel.protoBarrel.Count().ToString();
                                string char5 = myBarrel.Strands.Count().ToString();
                                string char6 = myBarrel.AvgTilt.ToString();
                                newLine = char1 + "\t" + char4 + "\t" + char5 + "\t" + char2 + "\t" + char7 + "\t" + char8 + "\t" + char9 + "\t" + char3 + "\t" + char6; 
                                file.WriteLine(newLine);*/
                            }
                        }
                    }
                    
                }
            }
            else
            {
                Console.WriteLine("could not open {0}", fileOfPDBs);
                Console.ReadLine();

            }
        }

        public static PBarrel runPBetaBarrel(string pdb, ref Dictionary<string, AminoAcid> _aaDict)
        {
            //Be sure that your chains are in alphabetical order in PDB structure; otherwise, alter the names so that they are or PBarrel will not know what position they belong in.
            //new peptide analysis
            //PeptideAnalysis myModule = new PeptideAnalysis();
            Directory.SetCurrentDirectory(polyDBDirectory);
            
            //insert name of xml file for doing things
            string pdbXMLfile = pdb + ".xml";
            string pdbName = pdb;
            AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
            
            //Test to see if there is an xml file for this pdb ID name; if not, attempt to use the pdb file instead
            if (File.Exists(pdbXMLfile))
            {
                AtomParser.XmlAtomParser myXMLparser = new AtomParser.XmlAtomParser();
                
                {
                    Console.WriteLine("opened {0}", pdbXMLfile);
                    string emptyString = ""; // tells it to pass all atoms
                    myXMLparser.ParseXmlFileAndPassStorage(pdbXMLfile, emptyString, ref myAtomCat);

                }
            }
            else
            {
                string fileName = pdb + ".pdb";
                Console.WriteLine("opened {0}", fileName);
                myAtomCat = readPdbFile(fileName, ref partialChargesDict);
            }
            
            //PolyProtein newProt = new PolyProtein(ref myAtomCat, 0, pdbName); // For considering only One chain at a time in polymeric beta barrel
            PolyProtein newProt = new PolyProtein(ref myAtomCat, pdbName); //For considering all chains

            Console.WriteLine("creating barrel class..");

            PBarrel myBarrel = new PBarrel(newProt, polyoutDirectory, polyDBDirectory);

            return (myBarrel);

        }

        static public void startSolubleBarrel()
        {
            //for making barrel
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();

            //PDBIDs that will be run
            string fileOfPDBs = solDBDirectory + "Thornton1989.txt"; //input file with list of pdb files

            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    string fileLocation2 = soloutDirectory + "AllSolubleChar.txt";
                    using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
                    {
                        string newLine = "PDB" + "\t\t" + "Total Chains" + "\t" + "Total Strands" + "\t" + "Length" + "\t" + "AvgLength" + "\t" + "MinLength" + "\t" + "MaxLength";
                        file.WriteLine(newLine);
                        // Read and display lines from the file until the end of the file is reached.
                        while ((line = sr.ReadLine()) != null)
                        {
                            string[] splitLine = Array.FindAll<string>(((string)line).Split(
                                new char[] { ' ', '\t', ',' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                            string pdb = splitLine[0];
                            if (pdb != "IDs")
                            {
                                SBarrel myBarrel = runSolBarrel(pdb, ref AADict, ref partialChargesDict);
                                /*string char1 = myBarrel.PdbName;
                                string char2 = myBarrel.Axis.Length.ToString();
                                string char7 = myBarrel.StrandLength.Average().ToString();
                                string char8 = myBarrel.StrandLength.Min().ToString();
                                string char9 = myBarrel.StrandLength.Max().ToString();
                                string char4 = myBarrel.protoBarrel.Count().ToString();
                                string char5 = myBarrel.Strands.Count().ToString();
                                newLine = char1 + "\t" + char4 + "\t" + char5 + "\t" + char2 + "\t" + char7 + "\t" + char8 + "\t" + char9 + "\t" + char3 + "\t" + char6;
                                file.WriteLine(newLine);*/
                            }
                        }
                    }

                }
            }
            else
            {
                Console.WriteLine("could not open {0}", fileOfPDBs);
                Console.ReadLine();

            }
        }

        public static SBarrel runSolBarrel(string pdb, ref Dictionary<string, AminoAcid> _aaDict, ref Dictionary<Tuple<string, string>, double> partialChargesDict)
        {
            Directory.SetCurrentDirectory(solDBDirectory);

            string fileName2 = pdb + ".pdb";
            Console.WriteLine("opened {0}", fileName2);
            AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
            myAtomCat = readPdbFile(fileName2, ref partialChargesDict);

            string pdbName = pdb.Substring(0, 4).ToUpper();

            int chainNum = 0;
            //if (pdbName == "3RFZ" || pdbName == "3KVN") chainNum = 1;

            int stop = myAtomCat.ChainAtomList.Count();
            Console.Write("Protein Class {0}", chainNum);

            Protein newProt = new Protein(ref myAtomCat, chainNum, pdbName, solDBDirectory);

            //Console.WriteLine("ChainNum {0}", chainNum);

            Console.Write("creating barrel class..");

            SBarrel myBarrel = new SBarrel(newProt.Chains[0], newProt);

            return (myBarrel);
        }
            //extracts atoms from pdb file
            public static AtomCategory readPdbFile(string pdbfilename, ref Dictionary<Tuple<string, string>, double> partialCharges)
            {
                AtomCategory myAtomCategory = new AtomCategory();
                
				ChainAtoms chainAtoms = new ChainAtoms();
                ArrayList atomList = new ArrayList();

                string preAsymId = "";
                string asymID = "";
				
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
                            
                            string[] spLine = Array.FindAll<string>(((string)line).Split(
                            new char[] { ' ', '\t', ',', ';' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                            List<string> splitLine = new List<string>(spLine);
                            
                            if (splitLine[0] == "ENDMDL") break;
                            
                            if (splitLine[0] == "ATOM")
                            {
                                AtomInfo myAtom = new AtomInfo();
                                asymID = line.Substring(21,1);
                                myAtom.altConfID = "";
                                myAtom.atomId = Convert.ToInt32(line.Substring(6,5).Trim());
                                myAtom.atomName = line.Substring(12,4).Trim();
                                myAtom.residue = line.Substring(17,3);
                                myAtom.authResidue = myAtom.residue;
                                myAtom.authSeqId = line.Substring(22,4).Trim();
                                if (Convert.ToInt32(myAtom.authSeqId) != previousSeqID) seqID = seqID + 1;
                                myAtom.seqId = seqID.ToString();
                                previousSeqID = Convert.ToInt32(myAtom.authSeqId);
                                
                                AtomParser.Coordinate coord = new Coordinate();
                                //Console.WriteLine("{0}, {1}, {2}", line.Substring(30, 8).Trim(), line.Substring(38, 8).Trim(), line.Substring(46,8).Trim());
                                coord = new Coordinate(Convert.ToDouble(line.Substring(30,8).Trim()), Convert.ToDouble(line.Substring(38,8).Trim()), Convert.ToDouble(line.Substring(46,8).Trim()));
                                myAtom.xyz = coord;
                                myAtom.atomType = line.Substring(76,2).Trim();
                                //This handful of conditions are occasionally different and throw off the check for hydrogen bonding, etc
                                if (myAtom.atomName == "HG") myAtom.atomName = "HG1";
                                if (myAtom.atomName == "OG") myAtom.atomName = "OG1";
                                if ((myAtom.atomName == "NH1" && myAtom.residue == "ARG")||(myAtom.atomName == "NZ" && myAtom.residue == "LYS")) myAtom.atomType = "N1+";
                                if ((myAtom.atomName == "OD2" && myAtom.residue == "ASP")||(myAtom.atomName == "OE2" && myAtom.residue == "GLU")) myAtom.atomType = "O1-";

                                myAtom.bFac = Convert.ToDouble(line.Substring(60,6).Trim());
                                myAtom.occupancy = 1.00;
                            
                                    //Compare chain ID. If this is now a new chain, write existing list, and create empty new ones.
                                    if (preAsymId != asymID && preAsymId != "" && atomList.Count > 0)
                                    {
                                        chainAtoms.AsymChain = preAsymId;
                                        chainAtoms.AuthAsymChain = preAsymId;
                                        chainAtoms.EntityID = "1"; // problem with int to string in new version
                                        chainAtoms.PolymerType = "-";

                                        AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                                        atomList.CopyTo(atomArray);
                                        chainAtoms.CartnAtoms = atomArray;
                                        myAtomCategory.AddChainAtoms(chainAtoms);
                                        atomList = new ArrayList();
                                        chainAtoms = new ChainAtoms();
                                    }
                                    //Add atom to current list of atoms if it is not a water molecule
                                    if (myAtom.residue.ToUpper() != "HOH")
                                    {
                                        if (myAtom.atomType == "H")
                                        {
                                            Tuple<string, string> key = new Tuple<string, string> (myAtom.residue, myAtom.atomName);
                                            if (partialCharges.ContainsKey(key) == true) atomList.Add(myAtom);
                                        }
                                        else atomList.Add(myAtom);
                                    }
                                    preAsymId = asymID;

                                
                            }
                        }
                        //Capture final chain
                        if (atomList.Count > 0)
                        { 
                            chainAtoms.AsymChain = preAsymId;
                            chainAtoms.AuthAsymChain = preAsymId;
                            chainAtoms.EntityID = "1";
                            chainAtoms.PolymerType = "-";
                            AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                            atomList.CopyTo(atomArray);
                            chainAtoms.CartnAtoms = atomArray;
                            myAtomCategory.AddChainAtoms(chainAtoms);
                            atomList = new ArrayList();
                            chainAtoms = new ChainAtoms();
                        }
				         
                    }
                }

                myAtomCategory.Resolution = 2.5;

                return myAtomCategory;
            }

        
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
                            new char[] { ' ', '\t', ',' }), delegate (string s)
                            {
                                return !String.IsNullOrEmpty(s);
                            });
                    string[] splitIntoWords2 = Array.FindAll<string>(_s2.Split(
                            new char[] { ' ', '\t', ',' }), delegate (string s)
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

    

