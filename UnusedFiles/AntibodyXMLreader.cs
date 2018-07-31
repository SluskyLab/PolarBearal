using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.IO;
using System.Xml;
using System.Data;

namespace betaBarrelProgram
{
    class AntibodyXMLreader
    {
        public AntibodyXMLreader()
        {
        }

        public AntibodyXMLreader(string _xmlFileWithPath)
        {
            Console.WriteLine("Before xml input, tables count = " + xmlDataSet.Tables.Count);
            ReadAntibodyXMLdb(_xmlFileWithPath);
            Console.WriteLine("After xml input, tables count = " + xmlDataSet.Tables.Count);
        }

        public void ReadAntibodyXMLdb(string _xmlFileWithPath)
        {
            xmlDataSet.ReadXml(_xmlFileWithPath);
        }

        public List<string> GetPDBids()
        {
            List<string> pdbIDs = new List<string>();
            foreach (DataRow myDataRow in xmlDataSet.Tables["antibody"].Rows)
            {
                pdbIDs.Add((string)myDataRow["pdb"]);
            }
            return pdbIDs;
        }

        public List<string> GetPDBidsWithoutDuplicates()
        {
            List<string> fullPDBlist = new List<string>();
            List<string> pdbsToSkip = new List<string>();
            List<string> pdbIDs = new List<string>();
            foreach (DataRow myDataRow in xmlDataSet.Tables["antibody"].Rows)
            {
                if (fullPDBlist.Contains((string)myDataRow["pdb"]))
                {
                    pdbsToSkip.Add((string)myDataRow["pdb"]);
                    Console.WriteLine("skipping {0} in pdblist: duplicated in master list", (string)myDataRow["pdb"]);
                }
                fullPDBlist.Add((string)myDataRow["pdb"]);
            }
            foreach (object pdbObj in fullPDBlist)
            {
                if (!pdbsToSkip.Contains((string)pdbObj))
                {
                    pdbIDs.Add((string)pdbObj);
                }
            }
            return pdbIDs;
        }

        public ArrayList GetResolutions()
        {
            ArrayList resolutionValues = new ArrayList();
            foreach (DataRow myDataRow in xmlDataSet.Tables["experimental_info"].Rows)
            {
                string resolutionString = (string)myDataRow["resolution"];
                if (resolutionString == "N/A")
                {
                    resolutionValues.Add((double)-1);
                }
                else
                {
                    resolutionString = resolutionString.Remove(resolutionString.Length - 1, 1);
                    resolutionValues.Add(Convert.ToDouble(resolutionString));
                }
            }
            return resolutionValues;
        }

        public Hashtable GetResolutionHashtable()
        {
            List<string> nonRedundantListOfPDBs = new List<string>();
            nonRedundantListOfPDBs = this.GetPDBidsWithoutDuplicates();
            return this.GetResolutionHashtable(nonRedundantListOfPDBs);
        }

        public Hashtable GetResolutionHashtable(List<string> _pdbsForList)
        {
            Hashtable resolutionHash = new Hashtable();
            List<string> fullPDBlist = this.GetPDBids();
            ArrayList resolutionValues = new ArrayList();
            foreach (DataRow myDataRow in xmlDataSet.Tables["experimental_info"].Rows)
            {
                string resolutionString = (string)myDataRow["resolution"];
                if (resolutionString == "N/A")
                {
                    resolutionValues.Add((double)-1);
                }
                else
                {
                    resolutionString = resolutionString.Remove(resolutionString.Length - 1, 1);
                    resolutionValues.Add(Convert.ToDouble(resolutionString));
                }
            }
            if (fullPDBlist.Count != resolutionValues.Count)
            {
                Console.WriteLine("fullPDBlist.Count {0} != resolutionValues.Count {1} in GetResolutionHashtable.",
                    fullPDBlist.Count, resolutionValues.Count);
                Console.ReadLine();
            }
            for (int fullPDBindex = 0; fullPDBindex < fullPDBlist.Count; fullPDBindex++)
            {
                if (_pdbsForList.Contains(fullPDBlist[fullPDBindex]))
                {
                    if (resolutionHash.ContainsKey(fullPDBlist[fullPDBindex]))
                    {
                        Console.WriteLine("resHash contains key {0} in GetResolutionHashtable.", (string)fullPDBlist[fullPDBindex]);
                        Console.ReadLine();
                    }
                    else
                    {
                        resolutionHash.Add(fullPDBlist[fullPDBindex], resolutionValues[fullPDBindex]);
                    }
                }
            }
            return resolutionHash;
        }
        /// <summary>
        /// returns max B-fac and resolution values for all pdbs in AB database without duplicate entries
        /// </summary>
        /// <param name="_myModule"></param>
        /// <param name="_pathToInputXMLPDBfiles"></param>
        /// <returns>Hashtable: keyed by pdbname, value: ArrayList:
        /// [0]: (double) resolution
        /// [1]: (Hashtable) keys: loopType values: (double) maxBfac</returns>
        public Hashtable GetQualityStatisticsHash(ref PeptideAnalysis _myModule, string _pathToInputXMLPDBfiles)
        {
            Hashtable qualityStatsHash = new Hashtable();
            List<string> nonredundantPDBlist = new List<string>();
            nonredundantPDBlist = GetPDBidsWithoutDuplicates();
            Hashtable resolutionHash = new Hashtable();
            Hashtable bfactorHash = new Hashtable();
            resolutionHash = GetResolutionHashtable(nonredundantPDBlist);
            //public Hashtable GetMaxBfactorHash(ref AntibodyXMLreader _xmlReader, string _pathToInputXMLPDBfiles, ArrayList _pdbList)
            bfactorHash = _myModule.GetMaxBfactorHash(this, _pathToInputXMLPDBfiles, nonredundantPDBlist);
            ICollection resKeys = resolutionHash.Keys;
            foreach (string pdbKey in resKeys)
            {
                if (bfactorHash.ContainsKey(pdbKey))
                {
                    ArrayList valuesForPDB = new ArrayList();
                    valuesForPDB.Add((double)resolutionHash[pdbKey]); // double
                    valuesForPDB.Add((Hashtable)bfactorHash[pdbKey]); // hash by loop type: values are doubles (bfac) 
                    qualityStatsHash.Add((string)pdbKey, valuesForPDB);
                }
                else
                {
                    Console.WriteLine("bfactorHash does not contain key {0} in GetQualityStatisticsHash", pdbKey);
                    Console.ReadLine();
                }
            }

            return qualityStatsHash;
        }


        public ArrayList GetLoopSequences(string _loopNumber)
        {
            ArrayList loopSequences = new ArrayList();

            foreach (DataRow myDataRow in xmlDataSet.Tables["cdr"].Rows)
            {
                if ((string)myDataRow["id"] == _loopNumber)
                {
                    loopSequences.Add((string)myDataRow["sequence"]);
                }
            }
            return loopSequences;
        }

        public void OpenPDBxmlFile(string _pdbXMLfile)
        {
            AtomParser.XmlAtomParser xmlAtomParser = 
                new AtomParser.XmlAtomParser();
            
        }

        public void ReportOnXmlDataSet()
        {
            Console.WriteLine(xmlDataSet.ToString()); // format?
            Console.WriteLine(xmlDataSet.Tables.ToString());
            Console.WriteLine(xmlDataSet.Tables.Count + "\n");
            for (int i = 0; i < xmlDataSet.Tables.Count; i++)
            {
                Console.Write(i + "\t");
                Console.Write(xmlDataSet.Tables[i]);
                Console.Write("\n");
            }
            Console.WriteLine();
            //int TTC = 1;
 
            /*for (int i = 0; i < xmlDataSet.Tables[TTC].Columns.Count; i++)
            {
                Console.Write(i + "\t");
                Console.Write(xmlDataSet.Tables[TTC].Columns[i] + "\n");
                Console.WriteLine("caption:\t" + xmlDataSet.Tables[TTC].Columns[i].Caption);
            }*/

            int setTheTableVal = 9;
            DataTable testDT = xmlDataSet.Tables[setTheTableVal];
            Console.WriteLine("\nName of DataTable in question:\n" + xmlDataSet.Tables[setTheTableVal]);

            int breakCounter = new int();
            breakCounter = 0;
            foreach (DataRow myDataRow in testDT.Rows)
            {
                //Console.WriteLine(breakCounter + "\t" + myDataRow.ToString());
                //Console.WriteLine("trying to print pdb = " + myDataRow["pdb"]);
                //string testThisString = (string)myDataRow["pdb"];
                //Console.WriteLine("testThisString = " + testThisString);
                object[] itemArrayDump = myDataRow.ItemArray;
                for (int j = 0; j < itemArrayDump.Length; j++)
                {
                    Console.WriteLine(j + "\t" + itemArrayDump[j]);
                }

                breakCounter++;
                if (breakCounter > 20) { break; }
                
            }
            int columnCounter = new int();
            columnCounter = 0;
            Console.WriteLine("\n\nNow for columns:\n");
            foreach (DataColumn myDataColumn in testDT.Columns)
            {
                Console.WriteLine(columnCounter + ":\t" + myDataColumn.Caption);
                columnCounter++;
            }
        }

        // member variables
        DataSet xmlDataSet = new DataSet();
    }
}
