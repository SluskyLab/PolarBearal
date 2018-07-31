/*
**  File: PDBInfo.cs (previously getOrganism.cs)
**  Started: pre 6/01/2015
**  Contributors: Joanna Slusky, Ryan Feehan
**  Overview: Prints relevant PDB information
**
**  About: Opens and reads .pdb headers for relevant information about that pdb and outputs it to file
**
**  Last Edited: 7/21/2017 by Ryan Feehan
**  NEEDS UPDATE, pdb headers have changed since written
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.IO;

namespace betaBarrelProgram
{
    public class PDBInfo
    {
       
        public PDBInfo()
        {
            runPDBInfo();
        }

        static public void runPDBInfo()
        {

            getPDBInfo(Global.MACMONODBDIR, Global.MONO_DB_DIR, Global.MONO_OUTPUT_DIR + "PDBinfo.txt");
            getPDBInfo(Global.MACPOLYDBDIR, Global.POLY_DB_DIR, Global.POLY_OUTPUT_DIR + "PDBinfo.txt");
        }
        static public void getPDBInfo(string fileOfPDBs, string DB_dir, string Output_File)
        {
            Console.WriteLine("Starting at {0}", DateTime.Now);
            List<string> pdbList = new List<string>();

            using (System.IO.StreamWriter fileOut = new System.IO.StreamWriter(Output_File))
            {
                fileOut.WriteLine("pdb\torganism\tmolecule\tauthorOligomericState\tcompOligomericState\tmutated\tunp#\tunp");
            }

            if (File.Exists(fileOfPDBs))
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
                            pdbList.Add(pdb);
                    }
                }
            }
            else Console.WriteLine("pdb {0} does not exist!", fileOfPDBs);
            for (int pdbNum = 0; pdbNum < pdbList.Count; pdbNum++)
            {
                string pdbName = pdbList[pdbNum];
                string pdbFile = DB_dir + pdbName + ".pdb";
                /*
                // fetches .pdb from internet if needed
                if (File.Exists(pdbFile)==false)
                {
                    Console.WriteLine("downloading {0}", pdbName);
                    string url = "http://www.rcsb.org/pdb/files/" + pdbName + ".pdb";
                    
                    WebClient myWebClient = new WebClient();
                    myWebClient.DownloadFile(url, pdbFile);
                }
                */

                 if (File.Exists(pdbFile))
                {
                    string organismName1 = "0";
                    string organismName2 = "0";
                    List<string> molecule = new List<string>();
                    string authorOlig = "0";
                    string compOlig = "0";
                    string unpNum = "0";
                    string unp = "0";
                    string date = "0";

                    bool isMutated = false;

                    using (StreamReader sr = new StreamReader(pdbFile))
                    {
                        String line;
                        // Read and display lines from the file until the end of
                        // the file is reached.
                        int lineNum = 0;
                        while ((line = sr.ReadLine()) != null)
                        {
                          
                                string[] splitLine = Array.FindAll<string>(((string)line).Split(
                                new char[] { ' ', '\t', ',',';' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                                for (int wordNum = 0; wordNum < splitLine.Count(); wordNum++)
                                {
                                    if (lineNum == 0)
                                    {
                                        if ((splitLine[wordNum][0].ToString() == "0" || splitLine[wordNum][0].ToString() == "1" || splitLine[wordNum][0].ToString() == "2" || splitLine[wordNum][0].ToString() == "3") && date == "0") date = splitLine[wordNum];
                                    }
                                        if (splitLine[wordNum] == "ORGANISM_SCIENTIFIC:" && organismName1 == "0")
                                    {
                                        organismName1 = splitLine[wordNum + 1];
                                        organismName2 = splitLine[wordNum + 2];
                                    }
                                    if (splitLine[wordNum] == "MUTATION:" && splitLine[wordNum + 1] == "YES") isMutated = true;
                                    if (splitLine[wordNum] == "MOLECULE:" && molecule.Count == 0)
                                    {
                                        for (int wordNum1 = wordNum+1; wordNum1 < splitLine.Count(); wordNum1++)
                                        {
                                            molecule.Add(splitLine[wordNum1]);
                                        }
                                       
                                    }
                                    if (splitLine.Count()>6 && splitLine[wordNum] == "AUTHOR" && splitLine[wordNum+1] == "DETERMINED"&&splitLine[wordNum+2] == "BIOLOGICAL"&&splitLine[wordNum+3] == "UNIT:" &&authorOlig == "0")
                                    {
                                        authorOlig = splitLine[wordNum + 4];
                                    }
                                    if (splitLine.Count()>6&& splitLine[wordNum] == "SOFTWARE" && splitLine[wordNum+1] == "DETERMINED"&&splitLine[wordNum+2] == "QUATERNARY"&&splitLine[wordNum+3] == "STRUCTURE:" &&compOlig == "0")
                                    {
                                        compOlig = splitLine[wordNum + 4];
                                    }
                                    if (unpNum=="0" && splitLine[wordNum] == "UNP" && splitLine[wordNum+2].Contains("_"))
                                    {
                                        unpNum = splitLine[wordNum + 1];
                                        unp = splitLine[wordNum + 2];
                                    }
                                    
                                }
                            lineNum++;
                            
                        }
                    }
                    Console.WriteLine("processing {0}", pdbName);
                    using (System.IO.StreamWriter fileOut = new System.IO.StreamWriter(Output_File, true))
                    {
                        fileOut.Write("{0}\t{1}\t", pdbName, organismName1 + " " + organismName2);
                        for (int molNum = 0; molNum < molecule.Count; molNum++)
                        {
                            fileOut.Write("{0} ", molecule[molNum]);
                        }
                        fileOut.WriteLine("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}",authorOlig, compOlig, isMutated, unpNum, unp, date);

                    }
                }
                else Console.WriteLine("pdb {0} does not exist!", pdbName);
            }
          
        }
    }
}
