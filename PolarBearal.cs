

using System;
using System.Collections.Generic;
using System.IO;
using System.Windows.Media.Media3D;
using System.Linq;


namespace betaBarrelProgram
{
    public class PolarBearal
    {
        //polarBearal Variable 
        public static string path = Global.MONO_OUTPUT_DIR + "PolarBearal\\";
        double zone = 13;
        static public void menu()
        {
            Console.WriteLine("1. Create Protein Database");
            Console.WriteLine("2. Run Protein Database");
            Console.WriteLine("4. Quit");
        }

        static public void RunPolarBearal()
        {
            /*Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
            string file = @"\PolarBearalResults.txt";
            using (StreamWriter output = new System.IO.StreamWriter(path + file))
            { }

            file = @"\ExamineZ.txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(path + file))
            { }

            file = @"\PinNoutByProtein.txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(path + file))
            {
                output.Write("\n {0} \t {1} \t {2} \t {3} \t {4} \t {5}", "proteinID", "AAs", "IN", "Pin", "OUT", "Nout");
            }

            file = @"\DisplayAngles.txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(path + file))
            {
                output.Write("{0}\t{1}\t{2}", "aa", "angle", "inward facing");
				}*/

            StartPolarBearal();
        }

        static public void StartPolarBearal()
        {
            string choice = "";
            while (choice != "4")
            {
                menu();
                choice = Console.ReadLine();
                switch (choice)
                {
                    case "1":
                        CreateBetaBarrelProteinDatabase();
                        break;
                    case "2":
                        RunBetaBarrelProteinDatabase();
                        break;
                    default:
                        choice = "4";
                        break;
                    case "5":
                        CreateBetaBarrelProtein();
                        break;
                }
            }
        }

        static public void CreateBetaBarrelProtein()
        {
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();
            Dictionary<string, AminoAcid> AADict = SharedFunctions.makeAADict();

            Console.WriteLine("please give pdb:");
            string insert = Console.ReadLine();
            string fileName = insert + ".xml";

            Program.runBetaBarrel(fileName, ref AADict, ref Global.partialChargesDict);
        }


        static public void CreateBetaBarrelProteinDatabase()
        {
            Dictionary<string, int> pdbBeta = new Dictionary<string, int>();
            
            //string fileOfPDBs = @"Z:\Documents\PhD\SluskyLab\MonoDB\MonoDBAlts.txt"; //input file with list of xml files
            //string fileOfPDBs = @"Z:\Documents\PhD\SluskyLab\8-12Scaffolds\ScaffoldListPB.txt"; //input file with list of xml files
            string fileOfPDBs = @"Z:\Documents\PhD\SluskyLab\MonoDB\MonoDBList_v5_85.txt"; //input file with list of xml files
            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    string fileLocation2 = path + "AllBarrelChar.txt";
                    using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileLocation2))
                    {
                        string newLine = "PDB" + "\t\t" + "Total Strands" +"\t" + "Length" + "\t" + "AvgLength" + "\t" + "MinLength" + "\t" + "MaxLength" + "\t" + "Radius" + "\t" + "Barrel Tilt";
                        file.WriteLine(newLine);
                        // Read and display lines from the file until the end of the file is reached.
                        while ((line = sr.ReadLine()) != null)
                        {
                            string[] splitLine = line.Split(new char[] { ' ', '\t', ',' });
                            string pdb = splitLine[0];
                            //Console.Write(pdb);
                            if (pdb != "IDs")
                            {
                                string fileName = pdb;
                                //string fileName = pdb + ".pdb";
                                Barrel myBarrel = Program.runBetaBarrel(fileName, ref Program.AADict, ref Program.partialChargesDict);
                                string char1 = myBarrel.PdbName;
                                string char2 = myBarrel.Axis.Length.ToString();
                                string char7 = myBarrel.StrandLength.Average().ToString();
                                string char8 = myBarrel.StrandLength.Min().ToString();
                                string char9 = myBarrel.StrandLength.Max().ToString();
                                string char3 = myBarrel.AvgRadius.ToString();
                                string char4 = myBarrel.Strands.Count.ToString();
                                string char5 = myBarrel.AvgTilt.ToString();
                                string char6 = myBarrel.ShearNum.ToString();
                                string char10 = myBarrel.PrevTwists.Average().ToString();
                                newLine = char1 + "\t" + char4 + "\t" + char2 + "\t" + char7 + "\t" + char8 + "\t" + char9 + "\t" + char3 + "\t" + char5 + "\t" + char6 + "\t" + char10;
                                file.WriteLine(newLine);
                                //Console.WriteLine("Number of Proteins: {0} \t AAs: {1} \t Double Checked Directions: {2}", totalProteins, totalAAs, numDoubleChecks);
                            }
                        }
                    }
                }
                //Console.WriteLine("Number of Proteins: {0} \t AAs: {1} \t Double Checked Directions: {2}", totalProteins, totalAAs, numDoubleChecks);
            }
            else
            {
                Console.WriteLine("I am in {0}", System.IO.Directory.GetCurrentDirectory());
                Console.WriteLine("could not open {0}", fileOfPDBs);
                Console.ReadLine();

            }
        }

        static public void RunBetaBarrelProteinDatabase()
        {
            string fileOfPDBs = Global.MacmonoDBDir; //input file with list of monomeric xml files

            if (File.Exists(fileOfPDBs))
            {
                using (StreamReader sr = new StreamReader(fileOfPDBs))
                {
                    String line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        string[] splitLine = Array.FindAll<string>(((string)line).Split(
                            new char[] { ' ', '\t', ',' }), delegate (string s) { return !String.IsNullOrEmpty(s); });
                        string pdb = splitLine[0];
                        if (pdb != "IDs")
                        {
                            string fileName = pdb;
                            PolarBearal polarRetest = new PolarBearal(fileName);
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


        //Constructor 
        //uses a previously extracted barrel data to create a quick model of a barrel
        public PolarBearal(string PDBid)
        {
            largestAltSeq = 0;
            proteinID = PDBid;
            totalProteins++;

            string file = @"\betaBarrel_aaOnly\" + PDBid + ".txt";
            string[] FindChains = File.ReadAllLines(path + file);
            numChains = FindChains.Length;

            simpleBarrel = new aa[numChains][];
            Queue<aa> chain = new Queue<aa>();
            file = @"\betaBarrel\2" + PDBid + ".txt";
            string[] BarrelChains = File.ReadAllLines(path + file);

            int curChain = 0;
            char[] delPunc = { ':', ',' };
            foreach (string amino in BarrelChains)
            {
                string[] preaa = amino.Split(delPunc);

                if (preaa[0] != curChain.ToString())
                {
                    simpleBarrel[curChain] = chain.ToArray();
                    chain.Clear();
                    curChain++;
                }

                aa tempaa = new aa(preaa[1]);

                if (preaa[2] == "True")
                {
                    tempaa.Inward = true;
                }
                else
                {
                    tempaa.Inward = false;
                }

                tempaa.ResNum = Convert.ToInt16(preaa[3]);

                tempaa.height = Double.Parse(preaa[4]);

                tempaa.angle = Double.Parse(preaa[5]);
                //Console.Write(preaa[0] + ":" + preaa[1] + "," + preaa[2] + "," + preaa[3] + "," + preaa[4] + "\n");

                if (tempaa.height > -zone && tempaa.height < zone)
                {
                    chain.Enqueue(tempaa);
                    file = @"\DisplayAngles.txt";
                    using (StreamWriter output = File.AppendText(path + file))
                    {
                        output.Write("\n{0}\t{1}\t{2}", tempaa.m_aa_ID, tempaa.angle, tempaa.Inward);
                    }
                }
            }

            simpleBarrel[curChain] = chain.ToArray();
            seperatedBarrel = new Queue<aa[][]>();
            SeperateBarrel();

            printZ();
            PolarBearalCalculations();
            PinNoutByProtein();


        }

        //Constructor 
        //uses xml and the rest of the program to discover a barrel
        //extracts important data and add it to file for quicker extraction in future program runs
        public PolarBearal(ref Barrel myBarrel)
        {
            numChains = myBarrel.Strands.Count;
            proteinID = myBarrel.PdbName;
            largestAltSeq = 0;
            totalProteins++;

            //2D array to store data;[chain number] [amino acid number]
            simpleBarrel = new aa[numChains][];

            //create simple barrel
            for (int curChain = 0; curChain < numChains; curChain++)
            {
                simpleBarrel[curChain] = new aa[myBarrel.Strands[curChain].NumOfRes + 1];
                //add aa
                for (int cur_aa = 0; cur_aa <= myBarrel.Strands[curChain].NumOfRes; cur_aa++)
                {
                    //add aa member variables here
                    simpleBarrel[curChain][cur_aa] = new aa(myBarrel.Strands[curChain].Residues[cur_aa].OneLetCode);
                    simpleBarrel[curChain][cur_aa].Inward = myBarrel.Strands[curChain].Residues[cur_aa].Inward;
                    simpleBarrel[curChain][cur_aa].ResNum = myBarrel.Strands[curChain].Residues[cur_aa].ResNum;
                    simpleBarrel[curChain][cur_aa].X = myBarrel.Strands[curChain].Residues[cur_aa].BackboneCoords["CA"].X;
                    simpleBarrel[curChain][cur_aa].Y = myBarrel.Strands[curChain].Residues[cur_aa].BackboneCoords["CA"].Y;
                    simpleBarrel[curChain][cur_aa].height = myBarrel.Strands[curChain].Residues[cur_aa].Z;

                    simpleBarrel[curChain][cur_aa].angle = Vector3D.AngleBetween(myBarrel.Strands[curChain].Residues[cur_aa].BackboneCoords["CA"] - ((myBarrel.Strands[curChain].Residues[cur_aa].BackboneCoords["N"] + myBarrel.Strands[curChain].Residues[cur_aa].BackboneCoords["C"]) / 2), myBarrel.Axis);

                }
            }

            string file = @"\betaBarrel_aaOnly\" + proteinID + ".txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(path + file))
            {
                for (int curChain = 0; curChain < simpleBarrel.GetLength(0); curChain++)
                {
                    //add aa
                    for (int cur_aa = 0; cur_aa < simpleBarrel[curChain].Length; cur_aa++)
                    {
                        output.Write(" {0} ", simpleBarrel[curChain][cur_aa].m_aa_ID);
                    }
                    output.WriteLine();
                }
            }
            string _pdb;
            int _pdb_strands;
            string _res;
            int _res_num;
            int _res_strand;
            double _res_ca_x;
            double _res_ca_y;
            double _res_ca_z;
            bool _inward;
            double _angle;

            file = @"\betaBarrelStrands\" + proteinID + ".txt";
            using (System.IO.StreamWriter output = new System.IO.StreamWriter(path + file))
            {
                output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", "_pdb", "_pdb_strands", "_res", "_res_num", "_res_strand", "_res_ca_x", "_res_ca_y", "_res_ca_z", "_inward", "_angle");
                for (int curChain = 0; curChain < simpleBarrel.GetLength(0); curChain++)
                {
                    //add aa
                    for (int cur_aa = 0; cur_aa < simpleBarrel[curChain].Length; cur_aa++)
                    {
                        _pdb = proteinID;
                        _pdb_strands = simpleBarrel.GetLength(0);
                        _res = simpleBarrel[curChain][cur_aa].m_aa_ID;
                        _res_num = simpleBarrel[curChain][cur_aa].ResNum + 1;
                        _res_strand = curChain + 1;
                        _res_ca_x = simpleBarrel[curChain][cur_aa].X;
                        _res_ca_y = simpleBarrel[curChain][cur_aa].Y;
                        _res_ca_z = simpleBarrel[curChain][cur_aa].height;
                        _inward = simpleBarrel[curChain][cur_aa].Inward;
                        _angle = simpleBarrel[curChain][cur_aa].angle;

                        output.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", _pdb, _pdb_strands, _res, _res_num, _res_strand, _res_ca_x, _res_ca_y, _res_ca_z, _inward, _angle);
                    }
                }
            }
        }


        //breaks a strand into substrands
        public void SeperateBarrel()
        {
            Queue<aa[]> seperatedChain = new Queue<aa[]>();
            int subChainCounter = 0;
            Queue<aa> subChain = new Queue<aa>();
            aa[] chain;

            for (int curChain = 0; curChain < numChains; curChain++)
            {
                chain = simpleBarrel[curChain];

                //break apart chain
                for (int cur_aa = 0; cur_aa < chain.Length; cur_aa++)
                {
                    if (subChain.Count == 0)
                    {
                        subChain.Enqueue(chain[cur_aa]);
                    }
                    else if (!alternating_aa(chain[cur_aa - 1], chain[cur_aa]))
                    {
                        //check for an increased largest chain count
                        if (subChain.Count > largestAltSeq) largestAltSeq = subChain.Count;

                        aa[] subChain_aa = new aa[subChain.Count];
                        subChainCounter = 0;
                        foreach (aa number in subChain)
                        {
                            subChain_aa[subChainCounter++] = number;
                        }

                        subChain.Clear();
                        subChain.Enqueue(chain[cur_aa]);

                        seperatedChain.Enqueue(subChain_aa);
                    }
                    else
                    {
                        subChain.Enqueue(chain[cur_aa]);
                    }
                }
                //adds last chain if it ended with alternating sequances
                if (subChain.Count > 0)
                {
                    if (subChain.Count > largestAltSeq) largestAltSeq = subChain.Count;

                    aa[] subChain_aa = new aa[subChain.Count];
                    subChainCounter = 0;
                    foreach (aa number in subChain)
                    {
                        subChain_aa[subChainCounter++] = number;
                    }
                    subChain.Clear();
                    seperatedChain.Enqueue(subChain_aa);
                }


                aa[][] newChain = new aa[seperatedChain.Count][];
                for (int curSubChain = 0; 0 < seperatedChain.Count; curSubChain++)
                {
                    newChain[curSubChain] = seperatedChain.Dequeue();
                }

                seperatedBarrel.Enqueue(newChain);
            }
        }

        //determines if two given AA are PtoN or NtoP
        private bool alternating_aa(aa prevAA, aa curAA)
        {
            if (prevAA != null && curAA != null)
            {
                if ('P' == prevAA.mPolarity || prevAA.mPolarity == 'Q')
                {
                    if ('N' == curAA.mPolarity || curAA.mPolarity == 'Q') return (true);
                    else return (false);
                }
                else if ('N' == prevAA.mPolarity || prevAA.mPolarity == 'Q')
                {
                    if ('P' == curAA.mPolarity || curAA.mPolarity == 'Q') return (true);
                    else return (false);
                }
                else//no other way to be alternating
                {
                    return (false);
                }
            }
            return (false);
        }


        //Displays IN, Pin, OUT, Nout for each protein
        public void PinNoutByProtein()
        {
            int AAs = 0, Pin = 0, Nout = 0, IN = 0, OUT = 0;

            foreach (aa[] strand in this.simpleBarrel)
            {
                foreach (aa amino in strand)
                {
                    AAs++;

                    if (amino.Inward)
                    {
                        IN++;
                        if (amino.mPolarity == 'P' || amino.mPolarity == 'Q') Pin++;
                    }

                    if (!amino.Inward)
                    {
                        OUT++;
                        if (amino.mPolarity == 'N' || amino.mPolarity == 'Q') Nout++;
                    }
                }
            }

            string file = @"\PinNoutByProtein.txt";
            using (System.IO.StreamWriter output = File.AppendText(path + file))
            {
                output.Write("\n {0} \t {1} \t {2} \t {3} \t {4} \t {5}", proteinID, AAs, IN, Pin, OUT, Nout);
            }
        }

        //For normal calculations
        public void PolarBearalCalculations()
        {
            totalChains += simpleBarrel.Count();
            for (int i = 0; i < simpleBarrel.Length; i++)
            {
                strandSizes[simpleBarrel[i].Length]++;
            }

            foreach (aa[] strand in simpleBarrel)
            {
                foreach (aa amino in strand)
                {
                    _aa[IntOfeRes(amino.m_aa_ID)]++;
                    if (amino.Inward) aaInward[IntOfeRes(amino.m_aa_ID)]++;
                    else aaOutward[IntOfeRes(amino.m_aa_ID)]++;
                }
            }

            //numbers to make a graph of occurances by chain length
            int[] OcurrancesPerLengthOfAltSeq = new int[largestAltSeq + 1];

            aa[][] temp2D_aa;
            int tempLength;
            for (int curChain = 0; curChain < numChains; curChain++)
            {
                temp2D_aa = seperatedBarrel.Dequeue();
                for (int curSubChain = 0; curSubChain < temp2D_aa.Length; curSubChain++)
                {
                    tempLength = temp2D_aa[curSubChain].Length;
                    OcurrancesPerLengthOfAltSeq[tempLength]++;
                    seqCount[tempLength]++;
                }
                seperatedBarrel.Enqueue(temp2D_aa);
            }

            string file = @"\Polar\PolarBearalGraphOccurances.txt";
            using (StreamWriter sw = File.AppendText(path + file))
            {
                sw.Write("{0} \t", proteinID);
                foreach (int occurances in OcurrancesPerLengthOfAltSeq)
                {
                    sw.Write("\t {0}", occurances);
                }
                sw.Write("\n");
            }

            AlternatingFacing();
        }

        //converts amino acid 1 letter ID into associated interger
        public int IntOfeRes(string aa)
        {
            switch (aa)
            {
                case "L":
                    return (1);
                case "W":
                    return (2);
                case "P":
                    return (3);
                case "V":
                    return (4);
                case "I":
                    return (5);
                case "F":
                    return (6);
                case "Y":
                    return (7);
                case "H":
                    return (8);
                case "A":
                    return (9);
                case "M":
                    return (10);
                case "G":
                    return (11);
                case "T":
                    return (12);
                case "Q":
                    return (13);
                case "D":
                    return (14);
                case "N":
                    return (15);
                case "K":
                    return (16);
                case "S":
                    return (17);
                case "R":
                    return (18);
                case "E":
                    return (19);
                case "C":
                    return (20);
                default:
                    return (0);
            }
        }

        //use to examine specific values at different Zs
        public void printZ()
        {
            foreach (aa[] strand in simpleBarrel)
            {
                bool first = true;
                aa previousAA = strand[0];
                int strandNum = 0;
                foreach (aa amino in strand)
                {
                    if (!first)
                    {
                        if (amino.height > -zone && amino.height < zone && previousAA.height > -zone && previousAA.height < zone)
                        {
                            if (amino.mPolarity == 'P' || amino.mPolarity == 'Q') ZP[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            if (amino.mPolarity == 'N' || amino.mPolarity == 'Q') ZN[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;

                            if (amino.mPolarity == 'P' && previousAA.mPolarity == 'P')// || previousAA.mPolarity == 'Q' || amino.mPolarity == 'Q')
                            {
                                examineZP[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            }
                            else if (amino.mPolarity == 'N' && previousAA.mPolarity == 'N')// || previousAA.mPolarity == 'Q' || amino.mPolarity == 'Q')
                            {
                                examineZN[(int)Math.Floor(amino.height) + Convert.ToInt16(zone)]++;
                            }
                            else { }


                        }
                        strandNum++;
                    }
                    first = false;
                    previousAA = amino;
                }
            }
        }

        public void angleBreakDown()
        {
            //foreach()
        }
        //Examine Alternating Directionality
        void AlternatingFacing()
        {
            bool prevAAInward = true;
            double prevAAHeight = 0.0;
            char prevAAPol = 'P';
            int aaCount, seqCount, seqCount2;
            aa aa1_PreDev = new aa("A"), aa2_Dev = new aa("A"), aa3_PostDev = new aa("A");
            foreach (aa[] strand in this.simpleBarrel)
            {
                aaCount = 0;
                seqCount = 0;
                seqCount2 = 0;
                bool range = false;
                foreach (aa amino in strand)
                {
                    if (amino.height > -zone && amino.height < zone)
                    {
                        totalAAs++;
                        range = true;
                        if (aaCount > 0)
                        {
                            //totalAAs++;
                            if (seqCount == 0)
                            {
                                seqCount++;
                                //if (amino.Inward == false && amino.mPolarity == 'N' || amino.Inward == true && amino.mPolarity == 'P') seqCount++;
                            }
                            else
                            {
                                if (amino.Inward != prevAAInward) seqCount++;
                                else
                                {
                                    AlternatingInOutSeqLengths[seqCount]++;
                                    if (seqCount != 1)
                                    {
                                        if (amino.Inward && prevAAInward) aaDevInward[IntOfeRes(amino.m_aa_ID)]++;
                                        else if (!amino.Inward && !prevAAInward) aaDevOutward[IntOfeRes(amino.m_aa_ID)]++;
                                    }
                                    seqCount = 1;
                                }



                                if (seqCount2 == 0)
                                {
                                    seqCount2++;
                                    aa3_PostDev = amino;
                                }
                                else if (seqCount2 == 1)
                                {
                                    seqCount2++;
                                    aa2_Dev = aa3_PostDev;
                                    aa3_PostDev = amino;
                                }
                                else
                                {
                                    seqCount2++;
                                    aa1_PreDev = aa2_Dev;
                                    aa2_Dev = aa3_PostDev;
                                    aa3_PostDev = amino;

                                    if (aa1_PreDev.mPolarity == 'P' && aa2_Dev.mPolarity == 'P')
                                    {
                                        aap1[IntOfeRes(aa1_PreDev.m_aa_ID)]++;
                                        aap2[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                                        aan3[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                                    }

                                    if (aa1_PreDev.mPolarity == 'N' && aa2_Dev.mPolarity == 'N')
                                    {
                                        aan1[IntOfeRes(aa1_PreDev.m_aa_ID)]++;
                                        aan2[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                                        aap3[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                                    }
                                }
                            }
                        }
                        //else if (amino.Inward == false && amino.mPolarity == 'N' || amino.Inward == true && amino.mPolarity == 'P') seqCount++;
                        else seqCount++;
                    }


                    prevAAInward = amino.Inward;
                    prevAAHeight = amino.height;
                    prevAAPol = amino.mPolarity;
                    aaCount++;
                }
                if (aa3_PostDev.mPolarity == 'P' && aa2_Dev.mPolarity == 'P')
                {
                    aap1[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                    aap2[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                    aan3[21]++;
                }

                if (aa3_PostDev.mPolarity == 'N' && aa2_Dev.mPolarity == 'N')
                {
                    aan1[IntOfeRes(aa2_Dev.m_aa_ID)]++;
                    aan2[IntOfeRes(aa3_PostDev.m_aa_ID)]++;
                    aap3[21]++;
                }
                if (range) AlternatingInOutSeqLengths[seqCount]++;
            }
        }

        //member class aa
        class aa
        {
            //contructors
            public aa(string aa_ID)
            {
                m_aa_ID = aa_ID;
                mPolarity = 'U';
                assignPolarity();
                Inward = false;
                ResNum = -100;
                height = -100;
            }

            //assigns polarity from amino acid ID
            public void assignPolarity()
            {
                string ID = this.m_aa_ID;

                if (P.Contains(ID))
                {
                    mPolarity = 'P';
                }
                else if (N.Contains(ID))
                {
                    mPolarity = 'N';
                }
                else
                {
                    mPolarity = 'Q';
                }
            }

            //member variables
            public string m_aa_ID;//stores amino acid's 3 letter abreviation
            public char mPolarity;//should be P,N, or U(represent unassigned polarity}
            public bool Inward;//is true if amino acid is facing inward, false if facing outward
            public int ResNum;//residue number in protein
            public double X;//height along z axis
            public double Y;//height along z axis
            public double height;//height along z axis
            public double angle;
        }


        //member variables
        string proteinID;
        aa[][] simpleBarrel;
        Queue<aa[][]> seperatedBarrel;
        int numChains, largestAltSeq;
        public static int totalProteins = 0, totalChains = 0, totalAAs = 0;

        public static int[] seqCount = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        public static string P = "DEKNQRSTCH", N = "AFILVWYMP", Q = "G";

        public static int[] _aa = new int[21];
        public static int[] aaInward = new int[21];
        public static int[] aaOutward = new int[21];

        public static int[] aaDevInward = new int[21];
        public static int[] aaDevOutward = new int[21];

        public static int[] strandSizes = new int[50];
        public static int[] ZP = new int[150];
        public static int[] ZN = new int[150];
        public static int[] AlternatingInOutSeqLengths = new int[100];
        public static int[] examineZP = new int[125];
        public static int[] examineZN = new int[125];
        public static int numDoubleChecks = 0;

        //new data 10-3-15 Ryan
        public static int[] aap1 = new int[21];
        public static int[] aan1 = new int[21];
        public static int[] aap2 = new int[21];
        public static int[] aan2 = new int[21];
        public static int[] aap3 = new int[22];
        public static int[] aan3 = new int[22];
    }
}
