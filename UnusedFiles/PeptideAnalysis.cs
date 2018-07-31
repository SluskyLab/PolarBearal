using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.IO;
using System.Xml;
using System.Data;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;
using System.Threading;

namespace betaBarrelProgram
{
    /// <summary>
    /// contains tools for analyzing antibody loop info, incl. dihedral info, sequence info, etc...
    /// </summary>
    class PeptideAnalysis
    {
        public PeptideAnalysis()
        {
        }

        public Hashtable ReadInScoringMatrix(string _txtFileWithPath)
        {
            Hashtable scoreTable = new Hashtable();
            StreamReader scoreReader = new StreamReader(_txtFileWithPath);
            string scoreline = null;
            List<string> linesOfScore = new List<string>();
            List<string> columnLabels = new List<string>();
            while ((scoreline = scoreReader.ReadLine()) != null)
            { linesOfScore.Add(scoreline); }
            scoreReader.Close();
            for (int lineCtr = 0; lineCtr < linesOfScore.Count; lineCtr++)
            {
                string[] oneLineSplit = Array.FindAll<string>(linesOfScore[lineCtr].Split(
                    new char[] { ' ', '\t', ',' }), delegate(string s) { return !String.IsNullOrEmpty(s); });
                if (lineCtr == 0)
                {
                    foreach (string aLabel in oneLineSplit)
                    { columnLabels.Add(aLabel); }
                }
                else
                {
                    Hashtable hashOfValuesForOneLine = new Hashtable();
                    for (int colCtr = 1; colCtr < oneLineSplit.Length; colCtr++)
                    { hashOfValuesForOneLine.Add(columnLabels[colCtr - 1], Convert.ToInt32(oneLineSplit[colCtr])); }
                    scoreTable.Add(oneLineSplit[0], hashOfValuesForOneLine);
                }
            }
            return scoreTable;
        }

        public double AlignmentScore(ref Hashtable _scoreMatrix, string _seq1, string _seq2)
        {
            double scoreToReturn = new double();
            return scoreToReturn;
        }

        public double CalculateTorsion(double[] _atom1, double[] _atom2, double[] _atom3, double[] _atom4)
        {
            double[] vec1 = new double[3] { _atom2[0] - _atom1[0], _atom2[1] - _atom1[1], _atom2[2] - _atom1[2] };
            double[] vec2 = new double[3] { _atom3[0] - _atom2[0], _atom3[1] - _atom2[1], _atom3[2] - _atom2[2] };
            double[] vec3 = new double[3] { _atom4[0] - _atom3[0], _atom4[1] - _atom3[1], _atom4[2] - _atom3[2] };

            double[] projVec = new double[3] { Math.Sqrt(DotProduct(vec2, vec2)) * vec1[0], Math.Sqrt(DotProduct(vec2, vec2)) * vec1[1],
                Math.Sqrt(DotProduct(vec2, vec2)) * vec1[2] };

            return (Math.Atan2((DotProduct(projVec, CrossProduct(vec2, vec3))), 
                DotProduct(CrossProduct(vec1, vec2), CrossProduct(vec2, vec3))))*180/Math.PI;
        }

        private double[] CrossProduct(double[] _vec1, double[] _vec2)
        {
            double[] crossProductValue = new double[3];
            crossProductValue[0] = _vec1[1] * _vec2[2] - _vec1[2] * _vec2[1];
            crossProductValue[1] = _vec1[2] * _vec2[0] - _vec1[0] * _vec2[2];
            crossProductValue[2] = _vec1[0] * _vec2[1] - _vec1[1] * _vec2[0];
            return crossProductValue;
        }

        private double DotProduct(double[] _vec1, double[] _vec2)
        {
            return (_vec1[0] * _vec2[0] + _vec1[1] * _vec2[1] + _vec1[2] * _vec2[2]);
        }

        public ArrayList CalculateAllLoopDihedrals(ref AntibodyXMLreader _xmlReader, string _pathToInputXMLPDBfiles)
        {
            ArrayList dataToReturn = new ArrayList();
            List<string> pdbIDs = new List<string>();
            Hashtable loopSequencesHash = new Hashtable();
            pdbIDs = _xmlReader.GetPDBids();
            string[] listOfLoops = new string[6] { "L1", "L2", "L3", "H1", "H2", "H3" };
            foreach (string loopLabel in listOfLoops)
            {
                loopSequencesHash.Add(loopLabel, _xmlReader.GetLoopSequences(loopLabel));
            }

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Count; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToInputXMLPDBfiles + (string)pdbIDs[pdbIndex] + ".xml";
                if (File.Exists(fileToOpen))
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in listOfLoops)
                    {
                        string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex];
                        int loopLength = new int();
                        if (loopSequence != "N/A")
                        {
                            loopLength = loopSequence.Length;
                        }
                        else
                        {
                            loopLength = -1;
                        }
                        ArrayList dihedralsOfOneLoop = new ArrayList();
                        ArrayList allDataFromCLD = new ArrayList();
                        ArrayList allDataForHash = new ArrayList();
                        allDataFromCLD = CalculateLoopDihedrals(ref myAtomCategory, loopSequence);
                        for (int seqIndex = 0; seqIndex < allDataFromCLD.Count; seqIndex++)
                        {
                            dihedralsOfOneLoop.Add((double[])((ArrayList)allDataFromCLD[seqIndex])[3]);
                        }
                        allDataForHash.Add(loopSequence);
                        allDataForHash.Add(loopLength);
                        allDataForHash.Add(dihedralsOfOneLoop);
                        loopInfoForOnePDB.Add(loopID, allDataForHash);
                    }
                    dataFromOnePDB.Add((string)pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    dataToReturn.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            return dataToReturn;
        }

        public ArrayList CalculateAllLoopDihedralsWithOmega(ref AntibodyXMLreader _xmlReader, string _pathToInputXMLPDBfiles)
        {
            ArrayList dataToReturn = new ArrayList();
            List<string> pdbIDs = new List<string>();
            Hashtable loopSequencesHash = new Hashtable();
            pdbIDs = _xmlReader.GetPDBids();
            string[] listOfLoops = new string[6] { "L1", "L2", "L3", "H1", "H2", "H3" };
            foreach (string loopLabel in listOfLoops)
            {
                loopSequencesHash.Add(loopLabel, _xmlReader.GetLoopSequences(loopLabel));
            }

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Count; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToInputXMLPDBfiles + (string)pdbIDs[pdbIndex] + ".xml";
                if (File.Exists(fileToOpen))
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in listOfLoops)
                    {
                        string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex];
                        int loopLength = new int();
                        if (loopSequence != "N/A")
                        {
                            loopLength = loopSequence.Length;
                        }
                        else
                        {
                            loopLength = -1;
                        }
                        ArrayList dihedralsOfOneLoop = new ArrayList();
                        ArrayList allDataFromCLD = new ArrayList();
                        ArrayList allDataForHash = new ArrayList();
                        allDataFromCLD = CalculateLoopDihedralsWithOmega(ref myAtomCategory, loopSequence);
                        for (int seqIndex = 0; seqIndex < allDataFromCLD.Count; seqIndex++)
                        {
                            dihedralsOfOneLoop.Add((double[])((ArrayList)allDataFromCLD[seqIndex])[3]);
                        }
                        allDataForHash.Add(loopSequence);
                        allDataForHash.Add(loopLength);
                        allDataForHash.Add(dihedralsOfOneLoop);
                        loopInfoForOnePDB.Add(loopID, allDataForHash);
                    }
                    dataFromOnePDB.Add((string)pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    dataToReturn.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            return dataToReturn;
        }

        public ArrayList CalculateAllLoopDihedralsWithOmega(string _pathToInputXMLPDBfiles, Hashtable _andrSeqHash, List<string> _loopIDs)
        {
            // changing to reflect addition of chain IDs 
            ArrayList dataToReturn = new ArrayList();
            Hashtable loopSequencesHash = new Hashtable();
            ICollection pdbAndChainKeys = _andrSeqHash.Keys;
            string[] pdbIDs = new string[pdbAndChainKeys.Count];
            int keyCounter = new int();
            keyCounter = 0;
            foreach (string pCK in pdbAndChainKeys)
            {
                pdbIDs[keyCounter] = pCK;
                keyCounter++;
            }

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Length; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToInputXMLPDBfiles + ((string)pdbIDs[pdbIndex]).Substring(0,4) + ".xml";
                if (File.Exists(fileToOpen))
                {
                    Hashtable hashForThisPDB = (Hashtable)_andrSeqHash[pdbIDs[pdbIndex]];
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    string chainID = null;
                    if (pdbIDs[pdbIndex].Length > 4)
                    {
                        chainID = pdbIDs[pdbIndex].Substring(4, 1);
                    }
                    else
                    {
                        chainID = "_";
                    }
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in _loopIDs)
                    {
                        string loopSequence = null;
                        //string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex]; // original
                        if (hashForThisPDB.ContainsKey(loopID))
                        {
                            loopSequence = (string)hashForThisPDB[loopID];
                        }
                        else
                        {
                            //loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex]; // orig
                            loopSequence = "?";
                        }
                        int loopLength = new int();
                        if ((loopSequence == "N/A") || (loopSequence == "?"))
                        {
                            loopLength = -1;
                        }
                        else
                        {
                            loopLength = loopSequence.Length;
                        }
                        ArrayList dihedralsOfOneLoop = new ArrayList();
                        ArrayList allDataFromCLD = new ArrayList();
                        ArrayList allDataForHash = new ArrayList();
                        allDataFromCLD = CalculateLoopDihedralsWithOmega(ref myAtomCategory, loopSequence, chainID);
                        for (int seqIndex = 0; seqIndex < allDataFromCLD.Count; seqIndex++)
                        {
                            dihedralsOfOneLoop.Add((double[])((ArrayList)allDataFromCLD[seqIndex])[3]);
                        }
                        allDataForHash.Add(loopSequence);
                        allDataForHash.Add(loopLength);
                        allDataForHash.Add(dihedralsOfOneLoop);
                        loopInfoForOnePDB.Add(loopID, allDataForHash);
                    }
                    dataFromOnePDB.Add((string)pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    dataToReturn.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            return dataToReturn;
        }

        public Hashtable GetMaxBfactorHash(AntibodyXMLreader _xmlReader, string _pathToInputXMLPDBfiles, List<string> _pdbList)
        {
            Hashtable maxBfacHash = new Hashtable();
            ArrayList allPDBdataInDB = new ArrayList();
            allPDBdataInDB = GetAllLoopMaxBfactors(ref _xmlReader, _pathToInputXMLPDBfiles);

            for (int allPDBindex = 0; allPDBindex < allPDBdataInDB.Count; allPDBindex++)
            {
                ArrayList dataOfOnePDB = (ArrayList)allPDBdataInDB[allPDBindex];
                if (_pdbList.Contains((string)dataOfOnePDB[0]))
                {
                    maxBfacHash.Add((string)dataOfOnePDB[0], (Hashtable)dataOfOnePDB[1]);
                }
            }
            return maxBfacHash;
        }

        public Hashtable GetMaxBfactorHash(string _pathToInputXMLPDBfiles, List<string> _pdbList, Hashtable _andrSeqHash, List<string> _loopIDs)
        {
            //_pdbList has no chain labels.  _andrSeqHash keys do
            Hashtable maxBfacHash = new Hashtable();
            ArrayList allPDBdataInDB = new ArrayList();
            allPDBdataInDB = GetAllLoopMaxBfactors(_pathToInputXMLPDBfiles, _andrSeqHash, _loopIDs);

            for (int allPDBindex = 0; allPDBindex < allPDBdataInDB.Count; allPDBindex++)
            {
                ArrayList dataOfOnePDB = (ArrayList)allPDBdataInDB[allPDBindex];
                if (_pdbList.Contains(((string)dataOfOnePDB[0]).Substring(0,4)))
                {
                    maxBfacHash.Add((string)dataOfOnePDB[0], (Hashtable)dataOfOnePDB[1]);
                }
            }
            return maxBfacHash;
        }

        public Hashtable GetResolutionHash(string _pathToPdbXmlFiles, List<string> _nrPDBs)
        {
            Hashtable resolutionHash = new Hashtable();
            foreach (string pdbID in _nrPDBs)
            {
                string fileToOpen = _pathToPdbXmlFiles + pdbID + ".xml";
                if (File.Exists(fileToOpen))
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    resolutionHash.Add(pdbID, myAtomCategory.Resolution);
                }
            }
            return resolutionHash;
        }

        public ArrayList GetAllLoopMaxBfactors(ref AntibodyXMLreader _xmlReader, string _pathToInputXMLPDBfiles)
            {
            ArrayList dataToReturn = new ArrayList();
            List<string> pdbIDs = new List<string>();
            Hashtable loopSequencesHash = new Hashtable();
            pdbIDs = _xmlReader.GetPDBids();
            string[] listOfLoops = new string[6] { "L1", "L2", "L3", "H1", "H2", "H3" };
            foreach (string loopLabel in listOfLoops)
            {
                loopSequencesHash.Add(loopLabel, _xmlReader.GetLoopSequences(loopLabel));
            }

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Count; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToInputXMLPDBfiles + (string)pdbIDs[pdbIndex] + ".xml";
                if (File.Exists(fileToOpen))
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in listOfLoops)
                    {
                        string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex];
                        int loopLength = new int();
                        if (loopSequence != "N/A")
                        {
                            loopLength = loopSequence.Length;
                        }
                        else
                        {
                            loopLength = -1;
                        }
                        double maxBfacOfOneLoop = new double();
                        maxBfacOfOneLoop = GetMaxBfacOfLoop(ref myAtomCategory, loopSequence);
                        loopInfoForOnePDB.Add(loopID, maxBfacOfOneLoop);
                    }
                    dataFromOnePDB.Add((string)pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    dataToReturn.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            return dataToReturn;
        }

        /// <summary>
        /// Returns a list of all detected chains, their author ChainIDs, and their sequences
        /// </summary>
        /// <param name="_pdbFileName">filename including path</param>
        /// <returns>ArrayList, each entry is ArrayList of [0]:AuthAsymId [1]:sequence [2]:bool:isJustTrace</returns>
        public ArrayList GetSequencesForAllChains(string _pdbFileName)
        {
            ArrayList seqsToReturn = new ArrayList();
            if (File.Exists(_pdbFileName))
            {
                AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                string atomToGrab = ""; // gets all atoms in file
                myXmlAtomParser.ParseXmlFileAndPassStorage(_pdbFileName, atomToGrab, ref myAtomCategory);
                AtomParser.ChainAtoms[] myChainAtomsArray = myAtomCategory.BackboneAtomList();
                AtomParser.AtomInfo[] currentChain; // this is how you refer to a single chain
                for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
                {
                    List<string> chainIDandSeqPair = new List<string>();
                    ArrayList chainIDseqAndTraceBool = new ArrayList();
                    chainIDseqAndTraceBool.Add(myChainAtomsArray[chainMarker].AuthAsymChain);
                    currentChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                    string theSequence = GiveSequenceOfChain(ref currentChain);
                    chainIDseqAndTraceBool.Add(theSequence);
                    // now: check if sequence is CA-only
                    bool alphaCarbonOnly = new bool();
                    alphaCarbonOnly = IsAlphaCarbonOnly(ref currentChain);
                    chainIDseqAndTraceBool.Add(alphaCarbonOnly);
                    seqsToReturn.Add(chainIDseqAndTraceBool);
                }
            }
            return seqsToReturn;
        }


        public ArrayList GetAllLoopMaxBfactors(string _pathToInputXMLPDBfiles, Hashtable _andrSeqHash, List<string> _loopIDs)
        {
            ArrayList dataToReturn = new ArrayList();
            List<string> pdbIDs = new List<string>();
            Hashtable loopSequencesHash = new Hashtable();
            ICollection pdbKeys = _andrSeqHash.Keys;
            foreach (string pdbID in pdbKeys)
            { pdbIDs.Add(pdbID); }
            //string[] listOfLoops = new string[6] { "L1", "L2", "L3", "H1", "H2", "H3" };

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Count; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToInputXMLPDBfiles + pdbIDs[pdbIndex].Substring(0,4) + ".xml";
                Hashtable hashForThisPDB = null;
                if (_andrSeqHash.ContainsKey(pdbIDs[pdbIndex]))
                { hashForThisPDB = (Hashtable)_andrSeqHash[pdbIDs[pdbIndex]]; }
                if (File.Exists(fileToOpen) && hashForThisPDB != null)
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    string chainID = pdbIDs[pdbIndex].Substring(4, 1);
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in _loopIDs)
                    {
                        //string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex]; // old version
                        string loopSequence = null;
                        if (hashForThisPDB.ContainsKey(loopID))
                        {
                            loopSequence = (string)hashForThisPDB[loopID];
                        }
                        else
                        {
                            loopSequence = "?";
                        }
                        int loopLength = new int();
                        if ((loopSequence == "N/A") || (loopSequence == "?"))
                        {
                            loopLength = -1;
                        }
                        else
                        {
                            loopLength = loopSequence.Length;
                        }
                        double maxBfacOfOneLoop = new double();
                        maxBfacOfOneLoop = GetMaxBfacOfLoop(ref myAtomCategory, loopSequence, chainID);
                        loopInfoForOnePDB.Add(loopID, maxBfacOfOneLoop);
                    }
                    dataFromOnePDB.Add(pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    dataToReturn.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            return dataToReturn;
        }

        public Hashtable GetQualityStatisticsHash(string _pathToInputXMLPDBfiles, Hashtable _andrSeqHash, Hashtable _confEHash, List<string> _loopIDs)
        {
            Hashtable qualityStatsHash = new Hashtable();
            List<string> nonredundantPDBlist = new List<string>();
            ICollection seqPDBKeys = _andrSeqHash.Keys; // no chainIDs!
            foreach (string thePDBandChainID in seqPDBKeys)
            {
                string pdbID = thePDBandChainID.Substring(0, 4);
                if (!nonredundantPDBlist.Contains(pdbID))
                { nonredundantPDBlist.Add(pdbID); }
            }

            Hashtable bfactorHash = new Hashtable();
            Hashtable resolutionHash = new Hashtable();

            bfactorHash = GetMaxBfactorHash(_pathToInputXMLPDBfiles, nonredundantPDBlist, _andrSeqHash, _loopIDs);
            resolutionHash = GetResolutionHash(_pathToInputXMLPDBfiles, nonredundantPDBlist);


            ICollection bfacKeys = bfactorHash.Keys; // does have chain IDs!
            foreach (string pdbKey in bfacKeys)
            {
                if (_confEHash.ContainsKey(pdbKey) && resolutionHash.ContainsKey(pdbKey.Substring(0, 4)))
                {
                    ArrayList valuesForPDB = new ArrayList();
                    valuesForPDB.Add((double)resolutionHash[pdbKey.Substring(0, 4)]); // double
                    valuesForPDB.Add((Hashtable)bfactorHash[pdbKey]); // hash by loop type: values are doubles (bfac)
                    valuesForPDB.Add((Hashtable)_confEHash[pdbKey]);
                    qualityStatsHash.Add((string)pdbKey, valuesForPDB);
                }
                else
                {
                    if (!resolutionHash.ContainsKey(pdbKey.Substring(0, 4)))
                    { Console.WriteLine("resolutionHash does not contain key {0} in GetQualityStatisticsHash", pdbKey.Substring(0, 4)); }
                    if (!_confEHash.ContainsKey(pdbKey))
                    { Console.WriteLine("_confEHash does not contain key {0} in GetQualityStatisticsHash", pdbKey); }
                }
            }

            return qualityStatsHash;
        }

        public ArrayList GetAllLoopBfacs(ref AntibodyXMLreader _xmlReader, string _pathToXMLfiles)
        {
            ArrayList allBfacsByPDBandLoop = new ArrayList();
            // start of paste here
            List<string> pdbIDs = new List<string>();
            Hashtable loopSequencesHash = new Hashtable();
            pdbIDs = _xmlReader.GetPDBids();
            string[] listOfLoops = new string[6] { "L1", "L2", "L3", "H1", "H2", "H3" };
            foreach (string loopLabel in listOfLoops)
            {
                loopSequencesHash.Add(loopLabel, _xmlReader.GetLoopSequences(loopLabel));
            }

            for (int pdbIndex = 0; pdbIndex < pdbIDs.Count; pdbIndex++)
            {
                ArrayList filesThatDoNotExist = new ArrayList();
                string fileToOpen = _pathToXMLfiles + (string)pdbIDs[pdbIndex] + ".xml";
                if (File.Exists(fileToOpen))
                {
                    AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                    AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                    string atomToGrab = ""; // gets all atoms in file
                    Hashtable loopInfoForOnePDB = new Hashtable();
                    ArrayList dataFromOnePDB = new ArrayList();
                    myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                    foreach (string loopID in listOfLoops)
                    {
                        string loopSequence = (string)((ArrayList)loopSequencesHash[loopID])[pdbIndex];
                        int loopLength = new int();
                        if (loopSequence != "N/A")
                        {
                            loopLength = loopSequence.Length;
                        }
                        else
                        {
                            loopLength = -1;
                        }
                        ArrayList bfactorsOfOneLoop = new ArrayList();
                        bfactorsOfOneLoop = GetBackboneBfactorsOfOneLoop(ref myAtomCategory, loopSequence);
                        loopInfoForOnePDB.Add(loopID, bfactorsOfOneLoop);
                    }
                    dataFromOnePDB.Add((string)pdbIDs[pdbIndex]);
                    dataFromOnePDB.Add(loopInfoForOnePDB);
                    allBfacsByPDBandLoop.Add(dataFromOnePDB);
                }
                else
                {
                    filesThatDoNotExist.Add(fileToOpen); // can do something with this information later
                }

            }
            // end of paste here

            return allBfacsByPDBandLoop;
        }

        public double GetMaxBfacOfLoop(ref AtomParser.AtomCategory _atomCategory, string _loopSequence)
        {
            double maxValueOfLoop = new double();
            maxValueOfLoop = -1;

            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();

            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                maxValueOfLoop = (double)999;
                return maxValueOfLoop;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get Bfactors & compare to maxValue so far
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] BfactorsOfOneResidue = GetBackboneBfactorsOfOneResidue(ref chainWithLoop, indexOfAtomInBackbone);
                foreach (double atomicBfactor in BfactorsOfOneResidue)
                {
                    if (atomicBfactor > maxValueOfLoop)
                    {
                        maxValueOfLoop = atomicBfactor;
                    }
                }
            }

            return maxValueOfLoop;
        }

        public double GetMaxBfacOfLoop(ref AtomParser.AtomCategory _atomCategory, string _loopSequence, string _chainID)
        {
            double maxValueOfLoop = new double();
            maxValueOfLoop = -1;

            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();
            indexOfLoopPosition = -1;

            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                if (myChainAtomsArray[chainMarker].AuthAsymChain != _chainID)
                { continue; }
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                maxValueOfLoop = (double)999;
                return maxValueOfLoop;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms(); // out of bounds here
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get Bfactors & compare to maxValue so far
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] BfactorsOfOneResidue = GetBackboneBfactorsOfOneResidue(ref chainWithLoop, indexOfAtomInBackbone);
                foreach (double atomicBfactor in BfactorsOfOneResidue)
                {
                    if (atomicBfactor > maxValueOfLoop)
                    {
                        maxValueOfLoop = atomicBfactor;
                    }
                }
            }

            return maxValueOfLoop;
        }

        public ArrayList GetBackboneBfactorsOfOneLoop(ref AtomParser.AtomCategory _atomCategory, string _loopSequence)
        {
            ArrayList BfactorsOfLoopByResidue = new ArrayList();

            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();

            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                return BfactorsOfLoopByResidue;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get Bfactors & compare to maxValue so far
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] BfactorsOfOneResidue = GetBackboneBfactorsOfOneResidue(ref chainWithLoop, indexOfAtomInBackbone);
                BfactorsOfLoopByResidue.Add(BfactorsOfOneResidue);
            }

            return BfactorsOfLoopByResidue;
        }

        /// <summary>
        /// calculates dihedral angles for given loop and pdbIndex (from GetPDBids() in AntibodyXMLreader)
        /// </summary>
        /// <param name="_atomCategory">ref to AtomCategory</param>
        /// <param name="_loopSequence">Sequence of loop(string)</param>
        /// <returns>ArrayList of an ArrayList for each position, which consists of:
        /// [0]: position in loop(int)
        /// [1]: seqID of that position(string)
        /// [2]: residueType of that position(string)
        /// [3]: phi,psi of that position (double[2])
        /// </returns>
        public ArrayList CalculateLoopDihedrals(ref AtomParser.AtomCategory _atomCategory, string _loopSequence)
        {

            ArrayList dihedrals = new ArrayList();
            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            //int[] indexOfSeqString = new int[2];
            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();
            

            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                ArrayList returnValueForMissingLoop = new ArrayList();
                int positionInLoop = new int();
                positionInLoop = -1;
                string seqID = "-1";
                string resType = "X";
                double[] phiPsi = new double[2] { 999, 999 };
                returnValueForMissingLoop.Add(positionInLoop);
                returnValueForMissingLoop.Add(seqID);
                returnValueForMissingLoop.Add(resType);
                returnValueForMissingLoop.Add(phiPsi);
                dihedrals.Add(returnValueForMissingLoop);
                return dihedrals;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);
            //indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[indexOfLoopPosition]);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // do this tomorrow: for each currentPosition value (which is the position within the loop sequence):
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get position to calcphipsi
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] onePhiPsi = new double[2];
                onePhiPsi = CalculatePhiAndPsi(ref chainWithLoop, indexOfAtomInBackbone);
                dataPointForOnePosition.Add(currentPosition);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].seqId);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].residue);
                dataPointForOnePosition.Add(onePhiPsi);
                dihedrals.Add(dataPointForOnePosition);
            }
            return dihedrals;
        }


        /// <summary>
        /// calculates dihedral angles for given loop and pdbIndex (from GetPDBids() in AntibodyXMLreader)
        /// </summary>
        /// <param name="_atomCategory">ref to AtomCategory</param>
        /// <param name="_loopSequence">Sequence of loop(string)</param>
        /// <returns>ArrayList of an ArrayList for each position, which consists of:
        /// [0]: position in loop(int)
        /// [1]: seqID of that position(string)
        /// [2]: residueType of that position(string)
        /// [3]: phi,psi,omega of that position (double[3])
        /// </returns>
        public ArrayList CalculateLoopDihedralsWithOmega(ref AtomParser.AtomCategory _atomCategory, string _loopSequence)
        {

            ArrayList dihedrals = new ArrayList();
            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            //int[] indexOfSeqString = new int[2];
            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();


            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                ArrayList returnValueForMissingLoop = new ArrayList();
                int positionInLoop = new int();
                positionInLoop = -1;
                string seqID = "-1";
                string resType = "X";
                double[] phiPsiOmega = new double[3] { 999, 999, 999 };
                returnValueForMissingLoop.Add(positionInLoop);
                returnValueForMissingLoop.Add(seqID);
                returnValueForMissingLoop.Add(resType);
                returnValueForMissingLoop.Add(phiPsiOmega);
                dihedrals.Add(returnValueForMissingLoop);
                return dihedrals;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);
            //indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[indexOfLoopPosition]);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // do this tomorrow: for each currentPosition value (which is the position within the loop sequence):
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get position to calcphipsi
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] onePhiPsiOmega = new double[3];
                onePhiPsiOmega = CalculatePhiPsiAndOmega(ref chainWithLoop, indexOfAtomInBackbone);
                dataPointForOnePosition.Add(currentPosition);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].seqId);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].residue);
                dataPointForOnePosition.Add(onePhiPsiOmega);
                dihedrals.Add(dataPointForOnePosition);
            }
            return dihedrals;
        }

        public ArrayList CalculateLoopDihedralsWithOmega(ref AtomParser.AtomCategory _atomCategory, string _loopSequence, string _chainID)
        {

            ArrayList dihedrals = new ArrayList();
            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            //int[] indexOfSeqString = new int[2];
            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();
            indexOfLoopPosition = -1;


            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                string authAsymCh = myChainAtomsArray[chainMarker].AuthAsymChain;
                if (authAsymCh != _chainID)
                { continue; }
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    chainLabel = chainMarker;
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                ArrayList returnValueForMissingLoop = new ArrayList();
                int positionInLoop = new int();
                positionInLoop = -1;
                string seqID = "-1";
                string resType = "X";
                double[] phiPsiOmega = new double[3] { 999, 999, 999 };
                returnValueForMissingLoop.Add(positionInLoop);
                returnValueForMissingLoop.Add(seqID);
                returnValueForMissingLoop.Add(resType);
                returnValueForMissingLoop.Add(phiPsiOmega);
                dihedrals.Add(returnValueForMissingLoop);
                return dihedrals;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);
            //indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[indexOfLoopPosition]);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                // do this tomorrow: for each currentPosition value (which is the position within the loop sequence):
                // 1.) get the seqID of that position from listOfCAseqIds[indexOfLoopPosition + currentPosition]
                // 2.) find the ALL BACKBONE index position via FindAtomBySeqID()
                // 3.) use that index to get position to calcphipsi
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(indexOfLoopPosition + currentPosition)]);
                double[] onePhiPsiOmega = new double[3];
                onePhiPsiOmega = CalculatePhiPsiAndOmega(ref chainWithLoop, indexOfAtomInBackbone);
                dataPointForOnePosition.Add(currentPosition);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].seqId);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].residue);
                dataPointForOnePosition.Add(onePhiPsiOmega);
                dihedrals.Add(dataPointForOnePosition);
            }
            return dihedrals;
        }

        public ArrayList CalculateLoopDihedralsWithOmega(ref AtomParser.AtomCategory _atomCategory, string _loopSequence, string _chainID, string _FASTAentryOfFirst)
        {
            ArrayList dihedrals = new ArrayList();
            // first, find loop sequence
            // then, get position of that sequence into chain object
            // finally, get torsions and sequence info, seqID info

            //int[] indexOfSeqString = new int[2];
            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCategory.BackboneAtomList();
            int chainLabel = new int();
            int indexOfLoopPosition = new int();
            indexOfLoopPosition = -1;


            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                string authAsymCh = myChainAtomsArray[chainMarker].AuthAsymChain;
                if (authAsymCh != _chainID)
                { continue; }
                AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                indexOfLoopPosition = FindThisSequence(ref tempChain, _loopSequence);
                if (indexOfLoopPosition > -1)
                {
                    int chainCtr = new int();
                    chainCtr = indexOfLoopPosition;
                    while (tempChain[indexOfLoopPosition].seqId != _FASTAentryOfFirst)
                    {
                        if (chainCtr + _loopSequence.Length > tempChain.Length) { indexOfLoopPosition = -1; break; }
                        AtomParser.AtomInfo[] smallerTempChain = new AtomParser.AtomInfo[tempChain.Length - indexOfLoopPosition - 1]; // new chain starts one position after start of found previous sequence
                        for (int atomCtr = chainCtr + 1; atomCtr < tempChain.Length; atomCtr++)
                        { smallerTempChain[atomCtr - chainCtr - 1] = tempChain[atomCtr]; }
                        int nextPosition = FindThisSequence(ref smallerTempChain, _loopSequence);
                        if (nextPosition < 0)
                        { indexOfLoopPosition = -1; break; }
                        chainCtr += nextPosition + 1;
                        indexOfLoopPosition = chainCtr;
                    }
                    if (indexOfLoopPosition > 0)
                    {
                        chainLabel = chainMarker;
                    }
                    break;  // found loop sequence, so exit with information
                }
            }

            if (indexOfLoopPosition == -1) // didn't find sequence
            {
                ArrayList returnValueForMissingLoop = new ArrayList();
                int positionInLoop = new int();
                positionInLoop = -1;
                string seqID = "-1";
                string resType = "X";
                double[] phiPsiOmega = new double[3] { 999, 999, 999 };
                returnValueForMissingLoop.Add(positionInLoop);
                returnValueForMissingLoop.Add(seqID);
                returnValueForMissingLoop.Add(resType);
                returnValueForMissingLoop.Add(phiPsiOmega);
                dihedrals.Add(returnValueForMissingLoop);
                return dihedrals;
            }

            AtomParser.AtomInfo[] chainWithLoop = myChainAtomsArray[chainLabel].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref chainWithLoop);
            int startingLoopPosition = listOfCAseqIds.IndexOf(_FASTAentryOfFirst);
            //indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[indexOfLoopPosition]);

            for (int currentPosition = 0; currentPosition < _loopSequence.Length; currentPosition++)
            {
                int indexOfAtomInBackbone = new int();
                ArrayList dataPointForOnePosition = new ArrayList();
                indexOfAtomInBackbone = FindAtomBySeqID(ref chainWithLoop, (string)listOfCAseqIds[(startingLoopPosition + currentPosition)]);
                double[] onePhiPsiOmega = new double[3];
                onePhiPsiOmega = CalculatePhiPsiAndOmega(ref chainWithLoop, indexOfAtomInBackbone);
                dataPointForOnePosition.Add(currentPosition);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].seqId);
                dataPointForOnePosition.Add(chainWithLoop[indexOfAtomInBackbone].residue);
                dataPointForOnePosition.Add(onePhiPsiOmega);
                dihedrals.Add(dataPointForOnePosition);
            }
            return dihedrals;
        }


        /// <summary>
        /// convert three letter aa code to one letter code
        /// </summary>
        /// <param name="threeLetters">AA three letter code</param>
        /// <returns></returns>
        public string threeToOne(string threeLetters)
        {
            int threeIndex = three2OneTable.IndexOf(threeLetters);
            if (threeIndex == -1) // not found
                return "X";
            else
                return three2OneTable.Substring(threeIndex + 5, 1);

        }

        public string GiveSequenceOfChain(ref AtomParser.AtomInfo[] _chain)
        {
            string chainSequence = "";
            foreach (AtomParser.AtomInfo _atom in _chain)
            {
                if (_atom.atomName == "CA" && (_atom.altConfID == "" || _atom.altConfID == "A"))
                {
                    chainSequence += threeToOne(_atom.residue);
                }
            }
            return chainSequence;
        }

        public bool IsAlphaCarbonOnly(ref AtomParser.AtomInfo[] _chain)
        {
            bool justTrace = new bool();
            justTrace = true;
            foreach (AtomParser.AtomInfo _atom in _chain)
            {
                if (_atom.atomName == "N")
                { justTrace = false; }
            }
            return justTrace;
        }

        /// <summary>
        /// Gets Sequences Flanking the supplied sequence from a pdb XML
        /// </summary>
        /// <param name="_pdbID">pdbID with chainID</param>
        /// <param name="_pathToPDBXMLs">pathToPDB-xml with frontslash</param>
        /// <param name="_NtermRes">number of res to supply for Nterm flank</param>
        /// <param name="_CtermRes">number of res to supply for Cterm flank</param>
        /// <param name="_loopSeq">Sequence for which flanks are needed</param>
        /// <returns>list-string of [0]:ntermflank [1]:querySeq [2]:C-term flank (X for no amino acid)</returns>
        public List<string> GetFlankingSequences(string _pdbID, string _pathToPDBXMLs, int _NtermRes, int _CtermRes, string _loopSeq)
        {
            List<string> sequencesToReturn = new List<string>();
            string NtermFlank = null;
            string CtermFlank = null;
            string fileToOpen = _pathToPDBXMLs + _pdbID.Substring(0, 4) + ".xml";
            if (File.Exists(fileToOpen))
            {
                AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
                AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
                string chainID = null;
                string atomToGrab = ""; // gets all atoms in the file
                if (_pdbID.Length > 4)
                {
                    chainID = _pdbID.Substring(4, 1);
                }
                else
                {
                    chainID = "_";
                }
                myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
                AtomParser.ChainAtoms[] myChainAtomsArray = myAtomCategory.BackboneAtomList();
                int indexOfLoopPosition = -1;

                for (int chainCtr = 0; chainCtr < myChainAtomsArray.Length; chainCtr++)
                {
                    string authAsynCh = myChainAtomsArray[chainCtr].authAsymChain;
                    if (authAsynCh == chainID)
                    {
                        AtomParser.AtomInfo[] tempChain = myChainAtomsArray[chainCtr].BackboneAtoms();
                        string theWholeChainSeq = GiveSequenceOfChain(ref tempChain);
                        indexOfLoopPosition = theWholeChainSeq.IndexOf(_loopSeq);
                        if (indexOfLoopPosition != -1)
                        {
                            // 01234567Loop01234567
                            // 01234567890123456789
                            // first, get the N-term flank
                            if (indexOfLoopPosition > (_NtermRes - 1)) // long enough that there's a full 8 residues before loop
                            {
                                NtermFlank = theWholeChainSeq.Substring(indexOfLoopPosition - _NtermRes, _NtermRes);
                                sequencesToReturn.Add(NtermFlank);
                            }
                            else
                            {
                                // add an "X" for each position missing at the start
                                NtermFlank = "";
                                for (int xCtr = 0; xCtr < (_NtermRes - indexOfLoopPosition); xCtr++)
                                { NtermFlank += "X"; }
                                if (indexOfLoopPosition < 0)
                                { NtermFlank += theWholeChainSeq.Substring(0, indexOfLoopPosition); }
                                sequencesToReturn.Add(NtermFlank);
                            }
                            sequencesToReturn.Add(_loopSeq);
                            // now, get the C-term flank
                            int testThisValue = new int();
                            testThisValue = theWholeChainSeq.Length;
                            int testThisNext = new int();
                            testThisNext = _loopSeq.Length;
                            if (theWholeChainSeq.Length > (indexOfLoopPosition + _loopSeq.Length + _CtermRes))
                            {
                                CtermFlank = theWholeChainSeq.Substring((indexOfLoopPosition + _loopSeq.Length), _CtermRes);
                                sequencesToReturn.Add(CtermFlank);
                            }
                            else // loop too close to end of chain: add X characters
                            {
                                CtermFlank = "";
                                if (theWholeChainSeq.Length > (indexOfLoopPosition + _loopSeq.Length)) // at least one residue after loop
                                { CtermFlank = theWholeChainSeq.Substring((indexOfLoopPosition + _loopSeq.Length)); }
                                for (int xCtr = 0; xCtr < (indexOfLoopPosition + _loopSeq.Length + _CtermRes - theWholeChainSeq.Length); xCtr++)
                                { CtermFlank += "X"; }
                                sequencesToReturn.Add(CtermFlank);
                            }                           
                        }
                        else // didn't find loop sequence
                        {
                            NtermFlank = "";
                            CtermFlank = "";
                            for (int xCtr = 0; xCtr < _NtermRes; xCtr++)
                            { NtermFlank += "X"; }
                            for (int xCtr = 0; xCtr < _CtermRes; xCtr++)
                            { CtermFlank += "X"; }
                            sequencesToReturn.Add(NtermFlank);
                            sequencesToReturn.Add(_loopSeq);
                            sequencesToReturn.Add(CtermFlank);
                        }
                    }
                }
            }            
            return sequencesToReturn;
        }

        public ArrayList GiveSequenceAndResidueIDsOfChain(string _pathToPdbXmlFiles, string _pdbIDWithChain)
        {
            // ************** working here ******************
            ArrayList arrayToReturn = new ArrayList(); // [0]: string (sequence) [1] ArrayList (identifiers)
            string pdbWithoutChainID = _pdbIDWithChain.Substring(0, 4);
            string chainID = _pdbIDWithChain.Substring(4, 1);
            string fileToOpen = _pathToPdbXmlFiles + pdbWithoutChainID + ".xml";
            AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
            AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
            string atomToGrab = "";
            myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);
            AtomParser.ChainAtoms[] myChainAtomsArray = myAtomCategory.BackboneAtomList();
            AtomParser.AtomInfo[] correctChain;
            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                if (myChainAtomsArray[chainMarker].AuthAsymChain != chainID)
                { continue; }
                //AtomParser.AtomInfo[] correctChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                correctChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                string sequence = GiveSequenceOfChain(ref correctChain);
                arrayToReturn.Add(sequence);
                ArrayList identifiers = new ArrayList();
                identifiers = MakeCAseqIDarray(ref correctChain);
                arrayToReturn.Add(identifiers);
                break;
            }
            return arrayToReturn;
        }

        public int FindThisSequence(ref AtomParser.AtomInfo[] _chain, string _sequence)
        {
            string chainSequence = GiveSequenceOfChain(ref _chain);
            return chainSequence.IndexOf(_sequence);
        }

        public double EndToEndDist(string _pdbAndCh, string _seq, string _pathToAntibodyXMLlib)
        {
            double theDistance = new double();
            string chainID = _pdbAndCh.Substring(4, 1);
            string pdbXMLIDwithPath = _pathToAntibodyXMLlib + _pdbAndCh.Substring(0, 4) + ".xml";
            string fetchAtomType = "";
            AtomParser.AtomCategory theAtomCat = new AtomParser.AtomCategory();
            AtomParser.XmlAtomParser myXmlAtomP = new AtomParser.XmlAtomParser();
            myXmlAtomP.ParseXmlFileAndPassStorage(pdbXMLIDwithPath, fetchAtomType, ref theAtomCat);

            AtomParser.ChainAtoms[] myChainAtomsArray = theAtomCat.ChainAtomList;
            int correctChainNumber = new int();
            correctChainNumber = 9999;
            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                if (myChainAtomsArray[chainMarker].AuthAsymChain != chainID)
                { continue; }
                correctChainNumber = chainMarker;
                break;
            }
            AtomParser.AtomInfo[] correctChain = myChainAtomsArray[correctChainNumber].BackboneAtoms();
            ArrayList listOfCAseqIds = new ArrayList();
            listOfCAseqIds = MakeCAseqIDarray(ref correctChain);
            int indexOfLoopPosition = FindThisSequence(ref correctChain, _seq);
            int indexOfFirstAtom = FindAtomBySeqID(ref correctChain, (string)listOfCAseqIds[indexOfLoopPosition]);
            int indexOfLastAtom = FindAtomBySeqID(ref correctChain, (string)listOfCAseqIds[indexOfLoopPosition + _seq.Length - 1]);
            double[] atom1 = new double[3];
            double[] atom2 = new double[3];
            atom1[0] = correctChain[indexOfFirstAtom].xyz.X;
            atom1[1] = correctChain[indexOfFirstAtom].xyz.Y;
            atom1[2] = correctChain[indexOfFirstAtom].xyz.Z;

            atom2[0] = correctChain[indexOfLastAtom].xyz.X;
            atom2[1] = correctChain[indexOfLastAtom].xyz.Y;
            atom2[2] = correctChain[indexOfLastAtom].xyz.Z;
            theDistance = Math.Sqrt(Math.Pow((atom1[0] - atom2[0]), 2) + Math.Pow((atom1[1] - atom2[1]), 2) + Math.Pow((atom1[2] - atom2[2]), 2));

            return theDistance;
        }

        public AtomParser.AtomInfo[] FindCorrectChain(ref AtomParser.AtomCategory _atomCat, string _authChainID)
        {
            //AtomParser.ChainAtoms[] myChainAtomsArray = _atomCat.BackboneAtomList(); // for backbone atoms only
            AtomParser.ChainAtoms[] myChainAtomsArray = _atomCat.ChainAtomList;
            AtomParser.AtomInfo[] correctChain;
            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                if (myChainAtomsArray[chainMarker].AuthAsymChain != _authChainID)
                { continue; }
                //correctChain = myChainAtomsArray[chainMarker].BackboneAtoms(); // returns backbone atoms only
                correctChain = myChainAtomsArray[chainMarker].CartnAtoms;
                return correctChain;
            }
            AtomParser.AtomInfo[] dummyChain = new AtomParser.AtomInfo[0];
            return dummyChain;
        }

        /// <summary>
        /// returns residue identifiers (author-generated) for start and end of input sequence
        /// </summary>
        /// <param name="_pathToXmls">path to pdb-xml files</param>
        /// <param name="_pdbAndChainID">pdbandAuthorChainID</param>
        /// <param name="_seq">loopSequence</param>
        /// <returns>strings of residueNum/IDs for start and end residues of loop</returns>
        public ArrayList FindStartAndEndResiNumbers(string _pathToXmls, string _pdbAndChainID, string _seq)
        {
            ArrayList resiNumValuesToReturn = new ArrayList();
            string badValueOne = "-1";
            string badValueTwo = "-1";
            resiNumValuesToReturn.Add(badValueOne);
            resiNumValuesToReturn.Add(badValueTwo);
            string pdbIDwithoutChain = _pdbAndChainID.Substring(0, 4);
            string chainID = _pdbAndChainID.Substring(4, 1);
            string fileToOpen = _pathToXmls + pdbIDwithoutChain + ".xml";
            AtomParser.AtomCategory myAtomCategory = new AtomParser.AtomCategory();
            AtomParser.XmlAtomParser myXmlAtomParser = new AtomParser.XmlAtomParser();
            string atomToGrab = "";
            myXmlAtomParser.ParseXmlFileAndPassStorage(fileToOpen, atomToGrab, ref myAtomCategory);

            //authSeqID is the value we're searching
            AtomParser.ChainAtoms[] myChainAtomsArray = myAtomCategory.BackboneAtomList();
            AtomParser.AtomInfo[] correctChain;
            for (int chainMarker = 0; chainMarker < myChainAtomsArray.Length; chainMarker++)
            {
                if (myChainAtomsArray[chainMarker].AuthAsymChain != chainID)
                { continue; }
                correctChain = myChainAtomsArray[chainMarker].BackboneAtoms();
                string sequence = GiveSequenceOfChain(ref correctChain);
                ArrayList identifiers = MakeCAseqIDarray(ref correctChain);
                // return an error code if the sequence length is not equal to the number of identifiers
                int seqLength = sequence.Length;
                // now, find search sequence in sequence string, then get start and end indices
                int startIndex = new int();
                startIndex = sequence.IndexOf(_seq);
                int endIndex = new int();
                endIndex = startIndex + _seq.Length - 1;
                // check for errors here
                bool thereIsAnError = new bool();
                thereIsAnError = false;
                if (sequence.Length != identifiers.Count)
                { thereIsAnError = true; }
                if (startIndex < 0)
                { thereIsAnError = true; }
                if (sequence.Length < (startIndex + _seq.Length))
                { thereIsAnError = true; }
                if (!thereIsAnError)
                {
                    resiNumValuesToReturn[0] = (string)identifiers[startIndex];
                    resiNumValuesToReturn[1] = (string)identifiers[endIndex];
                    break;
                }
            }
            return resiNumValuesToReturn;
        }

        /// <summary>
        /// finds int[] (position,chain#) of sequence in protein
        /// </summary>
        /// <param name="_atomCat">Protein container</param>
        /// <param name="_seq">sequence string</param>
        /// <returns>int[2] of (position, chain#)</returns>
        public int[] FindSequenceInProtein(ref AtomParser.AtomCategory _atomCat, string _seq)
        {
            int[] indexPosition = new int[2] { -1, -1 };
            int chainCounter = 0;
            foreach (AtomParser.ChainAtoms _chain in _atomCat.ChainAtomList)
            {
                AtomParser.AtomInfo[] thisChain = _chain.CartnAtoms;
                indexPosition[0] = GiveSequenceOfChain(ref thisChain).IndexOf(_seq);

                if (indexPosition[0] > -1)
                {
                    indexPosition[1] = chainCounter;
                    return indexPosition;
                }
                chainCounter++;
            }
            return indexPosition;
        }

        public ArrayList MakeCAseqIDarray(ref AtomParser.AtomInfo[] _chain)
        {
            ArrayList atomSeqAndIDArray = new ArrayList();
            foreach (AtomParser.AtomInfo _atom in _chain)
            {
                if (_atom.atomName == "CA" && (_atom.altConfID == "" || _atom.altConfID == "A"))
                {
                    atomSeqAndIDArray.Add(_atom.seqId);
                }
            }
            
            return atomSeqAndIDArray;
        }
        /// <summary>
        /// Given seqID, returns position of atom in AtomInfo[] array
        /// </summary>
        /// <param name="_chain">AtomParser.AtomInfo[] array, by ref</param>
        /// <param name="_seqID">residue # in string form</param>
        /// <returns>index in _chain of appropriate atom, or -1 if none found</returns>
        public int FindAtomBySeqID(ref AtomParser.AtomInfo[] _chain, string _seqID)
        {
            for (int atomIndex = 0; atomIndex < _chain.Length; atomIndex++)
            {
                if ((_chain[atomIndex].seqId == _seqID) && (_chain[atomIndex].atomName == "CA") && (_chain[atomIndex].altConfID == "" || _chain[atomIndex].altConfID == "A"))
                {
                    return atomIndex;
                }
            }
            return -1;
        }

        public int LocateAtomByID(ref AtomParser.AtomInfo[] _chain, int _atomID)
        {
            for (int indexValue = 0; indexValue < _chain.Length; indexValue++)
            {
                if (_chain[indexValue].atomId == _atomID && (_chain[indexValue].altConfID == "" || _chain[indexValue].altConfID == "A"))
                {
                    return indexValue;
                }
            }
            Console.WriteLine("Didn't find atom {0} in LocateAtomByID, should abort now!\n", _atomID);
            Console.ReadLine();
            return -1;
        }

        public List<AtomParser.AtomInfo> LocateAtomsOfResidue(ref AtomParser.AtomInfo[] _chain, string _fastaResID)
        {
            List<AtomParser.AtomInfo> atomsInRes = new List<AtomParser.AtomInfo>();
            foreach (AtomParser.AtomInfo anAtom in _chain)
            {
                if (anAtom.seqId == _fastaResID)
                { atomsInRes.Add(anAtom); }
            }

            return atomsInRes;
        }

        public ArrayList FindAllTurnsInThisCDR(string _pdbIDchainAndLoopID, string _pathToXMLfiles, ArrayList _cdrInfo, double _turnDistThreshold)
        {
            ArrayList emptyInfoToReturn = new ArrayList();
            ArrayList ALwithinEmpty = new ArrayList();
            string errorMsg = "Error calculating turns for this CDR";
            ALwithinEmpty.Add(errorMsg);
            emptyInfoToReturn.Add(ALwithinEmpty);
            ArrayList turnInfoToReturn = new ArrayList();
            //first, open the file. find the chain.  find each residue by fasta entry
            AtomParser.XmlAtomParser myxmlparser = new AtomParser.XmlAtomParser();
            AtomParser.AtomCategory myAtomCat = new AtomParser.AtomCategory();
            string xmlFile = _pathToXMLfiles + _pdbIDchainAndLoopID.Substring(0, 4).ToLower() + ".xml";
            myxmlparser.ParseXmlFileAndPassStorage(xmlFile, "", ref myAtomCat);
            AtomParser.AtomInfo[] correctChain = FindCorrectChain(ref myAtomCat, _pdbIDchainAndLoopID.Substring(4, 1));
            List<List<AtomParser.AtomInfo>> residuesOfCDR = new List<List<AtomParser.AtomInfo>>();
            string seq = (string)_cdrInfo[0];
            List<string> fastaValues = (List<string>)_cdrInfo[1];
            if (fastaValues.Count != seq.Length)
            { return emptyInfoToReturn; } // just in case there is a mismatch in CDR lengths (header seq versus res-by-res entries)
            List<int> alphaCarbonIndices = new List<int>();
            for (int resctr = 0; resctr < fastaValues.Count; resctr++)
            { 
                residuesOfCDR.Add(LocateAtomsOfResidue(ref correctChain, fastaValues[resctr]));
                if (residuesOfCDR[residuesOfCDR.Count - 1].Count == 0)
                { return emptyInfoToReturn; } // residue was not found, so quitting
                if ((threeToOne(residuesOfCDR[residuesOfCDR.Count - 1][0].residue) != seq.Substring(resctr, 1)) && (threeToOne(residuesOfCDR[residuesOfCDR.Count - 1][0].residue) != "X"))
                { return emptyInfoToReturn; } // this is to ensure that the sequence is correct

                for (int atmCtr = 0; atmCtr < residuesOfCDR[resctr].Count; atmCtr++)
                {
                    if (residuesOfCDR[resctr][atmCtr].atomName == "CA" && (residuesOfCDR[resctr][atmCtr].altConfID == "" || residuesOfCDR[resctr][atmCtr].altConfID == "A"))
                    {
                        for (int chainAtmCtr = 0; chainAtmCtr < correctChain.Length; chainAtmCtr++)
                        {
                            if (correctChain[chainAtmCtr] == residuesOfCDR[resctr][atmCtr])
                            { alphaCarbonIndices.Add(chainAtmCtr); break; }
                        }
                    }
                }
            }
            if (alphaCarbonIndices.Count != seq.Length) { return emptyInfoToReturn; } // check to make sure 

            // now: calculate all the dihedrals, find the conformation, find all the turns, give info in correct format
            List<double[]> allDihedralsForCDR = new List<double[]>();
            for (int resctr = 0; resctr < fastaValues.Count; resctr++)
            {
                allDihedralsForCDR.Add(CalculatePhiPsiAndOmega(ref correctChain, alphaCarbonIndices[resctr]));
                if (resctr > 2)
                {
                    if ((correctChain[alphaCarbonIndices[resctr - 3]] - correctChain[alphaCarbonIndices[resctr]]) <= _turnDistThreshold)
                    {
                        ArrayList oneTurnInfo = new ArrayList();
                        ArrayList dihedralsOfTurn = new ArrayList(); // needs to be in this format to calculate conf hash
                        for (int dihResCtr = resctr - 3; dihResCtr < resctr + 1; dihResCtr++)
                        {
                            dihedralsOfTurn.Add(allDihedralsForCDR[dihResCtr][0]);
                            dihedralsOfTurn.Add(allDihedralsForCDR[dihResCtr][1]);
                            dihedralsOfTurn.Add(allDihedralsForCDR[dihResCtr][2]);
                        }
                        oneTurnInfo.Add(fastaValues[resctr - 3]); // fastaID of res1
                        oneTurnInfo.Add(seq.Substring(resctr - 3, 4)); // sequence of turn
                        oneTurnInfo.Add(correctChain[alphaCarbonIndices[resctr - 3]] - correctChain[alphaCarbonIndices[resctr]]);
                        oneTurnInfo.Add(dihedralsOfTurn);
                        oneTurnInfo.Add(resctr - (int)2); // one-based value of start of turn
                        turnInfoToReturn.Add(oneTurnInfo);
                    }
                }
            }
            // turn info: fastaID of res2, seq, CA 1 to 4 dist, dihedrals, confhash            
            
            return turnInfoToReturn;
        }

        /// <summary>
        /// calc backbone dihedrals given AtomInfo[] ref of backbone atoms and CA position in the array 
        /// </summary>
        /// <param name="_chain">AtomInfo[] of backbone atoms</param>
        /// <param name="_atomIndex">index of alpha carbon of residue of interest in ALL BACKBONE array, NOT JUST CA!</param>
        /// <returns>double[2] of phi, psi(returns 999 if can't calculate value)</returns>
        public double[] CalculatePhiAndPsi(ref AtomParser.AtomInfo[] _chain, int _atomIndex)
        {
            double[] phiAndPsi = new double[2] { 999, 999 };

            if (_chain[_atomIndex].atomName != "CA")
            {
                return phiAndPsi;
            }

            int[] atomIndices = new int[5] { -1, -1, _atomIndex, -1, -1 };

            //order: 
            // [0] C' (index-1)
            // [1] N
            // [2] CA
            // [3] C'
            // [4] N (index+1)

            // now, identify all the atoms necessary
            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            // this would be prettier by using residue index #s but "seqID" is a string, not an int
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    break; // this means that the previous residue was hit before atoms could be found, i.e. missing atoms
                }

                if (_chain[backIndex].atomName == "C" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[0] = backIndex;
                    break; // found the C-atom, can exit now
                }
            }

            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[backIndex].atomName == "N" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[1] = backIndex;
                    break; // found the N-atom, can exit now
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "C" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[3] = forwardIndex;
                    break;
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "N" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[4] = forwardIndex;
                    break;
                }
            }

            if ((atomIndices[1] > -1) && (atomIndices[3] > -1))
            {
                double[] atom1 = new double[3];
                double[] atom2 = new double[3];
                double[] atom3 = new double[3];

                atom1[0] = _chain[atomIndices[1]].xyz.X;
                atom1[1] = _chain[atomIndices[1]].xyz.Y;
                atom1[2] = _chain[atomIndices[1]].xyz.Z;

                atom2[0] = _chain[atomIndices[2]].xyz.X;
                atom2[1] = _chain[atomIndices[2]].xyz.Y;
                atom2[2] = _chain[atomIndices[2]].xyz.Z;

                atom3[0] = _chain[atomIndices[3]].xyz.X;
                atom3[1] = _chain[atomIndices[3]].xyz.Y;
                atom3[2] = _chain[atomIndices[3]].xyz.Z;

                if ((atomIndices[0] > -1))
                {
                    double[] atom0 = new double[3];
                    atom0[0] = _chain[atomIndices[0]].xyz.X;
                    atom0[1] = _chain[atomIndices[0]].xyz.Y;
                    atom0[2] = _chain[atomIndices[0]].xyz.Z;
                    phiAndPsi[0] = CalculateTorsion(atom0, atom1, atom2, atom3);
                }
                if ((atomIndices[4] > -1))
                {
                    double[] atom4 = new double[3];
                    atom4[0] = _chain[atomIndices[4]].xyz.X;
                    atom4[1] = _chain[atomIndices[4]].xyz.Y;
                    atom4[2] = _chain[atomIndices[4]].xyz.Z;
                    phiAndPsi[1] = CalculateTorsion(atom1, atom2, atom3, atom4);
                }
            }
            return phiAndPsi;
        }

        double[] CalculatePhiPsiAndOmega(ref AtomParser.AtomInfo[] _chain, int _atomIndex)
        {
            double[] phiPsiAndOmega = new double[3] { 999, 999, 999 };

            if (_chain[_atomIndex].atomName != "CA")
            {
                return phiPsiAndOmega;
            }

            int[] atomIndices = new int[6] { -1, -1, -1, _atomIndex, -1, -1 };

            //order: 
            // [0] CA (index-1)
            // [1] C' (index-1)
            // [2] N
            // [3] CA
            // [4] C'
            // [5] N (index+1)

            // now, identify all the atoms necessary

            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            // this would be prettier by using residue index #s but "seqID" is a string, not an int
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[0] = backIndex;
                    break; // found CA for previous atom
                }
            }

            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    break; // this means that the previous residue was hit before atoms could be found, i.e. missing atoms
                }

                if (_chain[backIndex].atomName == "C" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[1] = backIndex;
                    break; // found the C-atom, can exit now
                }
            }

            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[backIndex].atomName == "N" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[2] = backIndex;
                    break; // found the N-atom, can exit now
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "C" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[4] = forwardIndex;
                    break;
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "N" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[5] = forwardIndex;
                    break;
                }
            }

            if ((atomIndices[1] > -1) && (atomIndices[3] > -1))
            {
                double[] atom1 = new double[3];
                double[] atom2 = new double[3];
                double[] atom3 = new double[3];

                atom1[0] = _chain[atomIndices[1]].xyz.X;
                atom1[1] = _chain[atomIndices[1]].xyz.Y;
                atom1[2] = _chain[atomIndices[1]].xyz.Z;

                atom2[0] = _chain[atomIndices[2]].xyz.X;
                atom2[1] = _chain[atomIndices[2]].xyz.Y;
                atom2[2] = _chain[atomIndices[2]].xyz.Z;

                atom3[0] = _chain[atomIndices[3]].xyz.X;
                atom3[1] = _chain[atomIndices[3]].xyz.Y;
                atom3[2] = _chain[atomIndices[3]].xyz.Z;

                if ((atomIndices[0] > -1))
                {
                    double[] atom0 = new double[3];
                    atom0[0] = _chain[atomIndices[0]].xyz.X;
                    atom0[1] = _chain[atomIndices[0]].xyz.Y;
                    atom0[2] = _chain[atomIndices[0]].xyz.Z;
                    phiPsiAndOmega[2] = CalculateTorsion(atom0, atom1, atom2, atom3);
                }
                if ((atomIndices[4] > -1) && (atomIndices[5]  > -1))
                {
                    double[] atom4 = new double[3];
                    double[] atom5 = new double[3];
                    atom4[0] = _chain[atomIndices[4]].xyz.X;
                    atom4[1] = _chain[atomIndices[4]].xyz.Y;
                    atom4[2] = _chain[atomIndices[4]].xyz.Z;
                    atom5[0] = _chain[atomIndices[5]].xyz.X;
                    atom5[1] = _chain[atomIndices[5]].xyz.Y;
                    atom5[2] = _chain[atomIndices[5]].xyz.Z;
                    phiPsiAndOmega[0] = CalculateTorsion(atom1, atom2, atom3, atom4);
                    phiPsiAndOmega[1] = CalculateTorsion(atom2, atom3, atom4, atom5);
                }
            }
            return phiPsiAndOmega;
        }

        public double[] GetBackboneBfactorsOfOneResidue(ref AtomParser.AtomInfo[] _chain, int _atomIndex)
        {
            double[] Bfactors = new double[4] { (double)999, (double)999, (double)999, (double)999 };
            // && (_atom.altConfID == "" || _atom.altConfID == "A")
            if (_chain[_atomIndex].atomName != "CA")
            {
                return Bfactors;
            }

            int[] atomIndices = new int[4] { -1, _atomIndex, -1, -1 };

            //order: 
            // [0] N
            // [1] CA
            // [2] C'
            // [3] O

            // now, identify all the atoms necessary

            for (int backIndex = _atomIndex - 1; backIndex > -1; backIndex--)
            {
                if (_chain[backIndex].atomName == "CA" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[backIndex].atomName == "N" && (_chain[backIndex].altConfID == "" || _chain[backIndex].altConfID == "A"))
                {
                    atomIndices[0] = backIndex;
                    break; // found the N-atom, can exit now
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "C" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[2] = forwardIndex;
                    break;
                }
            }

            for (int forwardIndex = _atomIndex + 1; forwardIndex < _chain.Length; forwardIndex++)
            {
                if (_chain[forwardIndex].atomName == "CA" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    break;
                }

                if (_chain[forwardIndex].atomName == "O" && (_chain[forwardIndex].altConfID == "" || _chain[forwardIndex].altConfID == "A"))
                {
                    atomIndices[3] = forwardIndex;
                    break;
                }
            }

            for (int atomCounter = 0; atomCounter < atomIndices.Length; atomCounter++)
            {
                if (atomIndices[atomCounter] > -1)
                {
                    Bfactors[atomCounter] = _chain[atomIndices[atomCounter]].bFac;
                }
            }

            return Bfactors;
        }

        // K-methods (K-means and K-medoids) start here

        public double CalcKDistance(double _angle1, double _angle2)
        {
            return (2 * (1 - System.Math.Cos((_angle1 - _angle2) * Math.PI / 180)));
        }

        public double AngleDifference(double _angle1, double _angle2)
        {
            return (Math.Abs((_angle1 - _angle2 + (double)180) % (double)360 - (double)180));
        }

        /// <summary>
        /// returns K-method Loop Distance D (double)
        /// </summary>
        /// <param name="_dihedralsFromLoop1">ArrayList of dihedrals (doubles) from loop 1</param>
        /// <param name="_dihedralsFromLoop2">ArrayList of dihedrals (doubles) from loop 2</param>
        /// <returns>loop distance as double (k-methods)</returns>
        public double CalculateLoopDistance(ArrayList _dihedralsFromLoop1, ArrayList _dihedralsFromLoop2)
        {
            if (_dihedralsFromLoop1.Count != _dihedralsFromLoop2.Count)
            {
                Console.WriteLine("unequal Loop Lengths in PeptideAnalysis.CalculateLoopDistance");
                Console.ReadLine();
                return 10000;
            }
            double loopDistance = new double();
            loopDistance = 0;
            for (int loopIndex = 0; loopIndex < _dihedralsFromLoop1.Count; loopIndex++)
            {
                loopDistance += CalcKDistance(((double)_dihedralsFromLoop1[loopIndex]), ((double)_dihedralsFromLoop2[loopIndex]));
            }
            return loopDistance;
        }

        public double CalculateLoopDistanceNoOmega(PeptideAnalysis.myLoopDataObject _l1, PeptideAnalysis.myLoopDataObject _l2)
        {
            if (_l1.myDihValues.Count != _l2.myDihValues.Count)
            {
                Console.WriteLine("unequal Loop Lengths in PeptideAnalysis.CalculateLoopDistance");
                Console.ReadLine();
                return 10000;
            }
            double loopDistance = new double();
            loopDistance = 0;
            for (int loopIndex = 0; loopIndex < _l1.myDihValues.Count; loopIndex = loopIndex + 3)
            {
                loopDistance += CalcKDistance((double)_l1.myDihValues[loopIndex], (double)_l2.myDihValues[loopIndex]);
                loopDistance += CalcKDistance((double)_l1.myDihValues[loopIndex + 1], (double)_l2.myDihValues[loopIndex + 1]);
            }
            return loopDistance;
        }

        /// <summary>
        /// Calculates the average angle (K-methods
        /// </summary>
        /// <param name="_angles"></param>
        /// <returns>double average angle (degrees)</returns>
        public double CalcAverageAngle(ArrayList _angles)
        {
            if (_angles.Count == 0)
            {
                Console.WriteLine("trying to calc avg angle with zero angles in set, error!"); return 55555;
            }
            double cosSum = new double();
            double sinSum = new double();
            cosSum = 0;
            sinSum = 0;
            foreach (double oneAngle in _angles)
            {
                cosSum += (1 / _angles.Count) * Math.Cos(oneAngle * Math.PI / 180);
                sinSum += (1 / _angles.Count) * Math.Sin(oneAngle * Math.PI / 180);
            }
            return Math.Atan2(sinSum, cosSum)*180/Math.PI;
        }
        /// <summary>
        /// Calculates the average loop conformation given an input arraylist of loop dihedrals(K-methods)
        /// </summary>
        /// <param name="_dihedralsOneLoopAllPDBs">
        /// input ArrayList: each entry is an ArrayList (one per pdbID) containing:
        /// * [0] pdbname
        /// * [1] ArrayList of dihedrals for one loop: each entry is a double[3] of phi,psi,omega</param>
        /// <returns>ArrayList of doubles: average phi,psi,omega for each loop position(so phi0,psi0,omega0,phi1,...</returns>
        ArrayList CalculateAverageLoopConformation(ArrayList _dihedralsByPDB)
        {
            ArrayList averageLoopConf = new ArrayList();
            ArrayList dihedralsByLoopPosition = new ArrayList();
            string[] dihedralType = { "phi", "psi", "omega" };
            // for each loop position and each dihedral type, need an arraylist of angles
            for (int loopPosition = 0; loopPosition < ((ArrayList)((ArrayList)dihedralsByLoopPosition[0])[1]).Count; loopPosition++)
            {
                ArrayList listForLoopPosition = new ArrayList();
                for (int dih = 0; dih < dihedralType.Length; dih++)
                {
                    ArrayList listForLoopPositionAndDihedralType = new ArrayList();
                    for (int pdbIndex = 0; pdbIndex < _dihedralsByPDB.Count; pdbIndex++)
                    {
                        listForLoopPositionAndDihedralType.Add(((double[])(((ArrayList)((ArrayList)dihedralsByLoopPosition[pdbIndex])[1])[loopPosition]))[dih]);
                    }
                    listForLoopPosition.Add(listForLoopPositionAndDihedralType);
                }
                dihedralsByLoopPosition.Add(listForLoopPosition);
            }
            for (int loopPosition = 0; loopPosition < dihedralsByLoopPosition.Count; loopPosition++)
            {
                for (int dih = 0; dih < dihedralType.Length; dih++)
                {
                    averageLoopConf.Add(CalcAverageAngle((ArrayList)((ArrayList)dihedralsByLoopPosition[loopPosition])[dih]));
                }
            }
            
            return averageLoopConf;
        }

        /// <summary>
        /// K-methods: outputs distance matrix
        /// </summary>
        /// <param name="_dihedralsByPDB">ArrayList of ArrayLists(one per PDB)
        /// each PDB ArrayList contains:
        /// [0]: pdbname
        /// [1]: ArrayList of dihedrals(phi0,psi0,omega0,phi1,psi1,omega1,...)</param>
        /// <returns>distanceMatrix[,]</returns>
        public double[,] GenerateDistanceMatrix(ref List<myLoopDataObject> _loopObjsByPDB)
        {
            double[,] distanceMatrix = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];
            for (int pdbIndex = 0; pdbIndex < _loopObjsByPDB.Count; pdbIndex++)
            {
                for (int secondPdbIndex = pdbIndex; secondPdbIndex < _loopObjsByPDB.Count; secondPdbIndex++)
                {
                    if (pdbIndex == secondPdbIndex)
                    {
                        distanceMatrix[pdbIndex, secondPdbIndex] = 0;
                    }
                    else
                    {
                        // get two arrayLists of dihedrals and pass to CalculateLoopDistance
                        //ArrayList dihedralsOfFirstPDB = ((ArrayList)((ArrayList)(_dihedralsByPDB[pdbIndex]))[1]);
                        //ArrayList dihedralsOfSecondPDB = ((ArrayList)((ArrayList)(_dihedralsByPDB[secondPdbIndex]))[1]);
                        ArrayList dihedralsOfFirstPDB = (_loopObjsByPDB[pdbIndex]).myDihValues;
                        ArrayList dihedralsOfSecondPDB = (_loopObjsByPDB[secondPdbIndex]).myDihValues;
                        distanceMatrix[pdbIndex, secondPdbIndex] = CalculateLoopDistance(dihedralsOfFirstPDB, dihedralsOfSecondPDB);
                        distanceMatrix[secondPdbIndex, pdbIndex] = distanceMatrix[pdbIndex, secondPdbIndex];
                    }
                }
            }
            return distanceMatrix;
        }

        /// <summary>
        /// Clusters by K-medoids method
        /// </summary>
        /// <param name="_numberOfClusters">number of clusters in sim</param>
        /// <param name="_dihedralsByPDB">ArrayList of ArrayLists (one per PDB)
        /// each PDB arraylist contains:
        /// [0]: pdbname
        /// [1]: ArrayList of dihedrals(phi0,psi0,omega0,phi1,psi1,omega1,...)</param>
        /// <param name="_scatterResults">passes scatter results back to main program or calling function</param>
        /// <returns>int[] of cluster membership by PDB index</returns>
        public int[] ClusterByKMedoids(int _numberOfClusters, List<myLoopDataObject> _dihedralsByPDB, int _maxNumberOfIterations, bool _withOutliers, ref double _scatterResults, ref ArrayList _centerOfEachCluster)
        {
            int[] bestClusterData = new int[_dihedralsByPDB.Count];
            ArrayList dihedralsOfFirstPDB = new ArrayList();
            dihedralsOfFirstPDB = _dihedralsByPDB[0].myDihValues;
            int loopLength = new int();
            loopLength = (int)(dihedralsOfFirstPDB.Count/3);
            // setup: calculate matrix of distances
            // need: *list of pdbs, dihedrals
            // to create distance matrix, need to pass ArrayList of ArrayLists (one per dihedral-position)[loop0][loop1]
            double[,] distanceMatrix = new double[_dihedralsByPDB.Count, _dihedralsByPDB.Count];
            distanceMatrix = GenerateDistanceMatrix(ref _dihedralsByPDB);

            Random random = new Random();
            double withinClusterScatter = new double();
            int[] bestClusterCenters = new int[_numberOfClusters];

            for (int iterationCounter = 0; iterationCounter < _maxNumberOfIterations; iterationCounter++)
            {
                int[] clusterData = new int[_dihedralsByPDB.Count];
                int[] clusterCenters = new int[_numberOfClusters];

                for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                {
                    clusterData[pdbIndex] = -1;
                }

                for (int clusterIndex = 0; clusterIndex < clusterCenters.Length; clusterIndex++)
                {
                    clusterCenters[clusterIndex] = -1;
                }

                // assign centers randomly
                
                for (int clusterIndex = 0; clusterIndex < clusterCenters.Length; clusterIndex++)
                {
                    clusterCenters[clusterIndex] = -1; // stays at -1 if cluster has no members
                    while (clusterCenters[clusterIndex] == -1)
                    {
                        int potentialCenter = new int();
                        potentialCenter = random.Next(0, clusterData.Length);
                        bool isCenterAlreadyAssigned = false;
                        for (int checkClusterIndex = 0; checkClusterIndex < clusterCenters.Length; checkClusterIndex++)
                        {
                            if (clusterCenters[checkClusterIndex] == potentialCenter)
                            {
                                isCenterAlreadyAssigned = true;
                            }
                        }
                        if (!isCenterAlreadyAssigned)
                        {
                            clusterCenters[clusterIndex] = potentialCenter;
                        }
                    }
                }
                

                //k-medoids++ D-squared weighting
                /*
                clusterCenters[0] = random.Next(0, clusterData.Length);
                for (int clusterIndex = 1; clusterIndex < clusterCenters.Length; clusterIndex++)
                {
                    // calculate the weights
                    double normalizationFactor = new double();
                    normalizationFactor = 0.000;
                    double[] newCenterProbabilities = new double[clusterData.Length];
                    double[] runningTotal = new double[clusterData.Length];
                    for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                    {
                        newCenterProbabilities[pdbIndex] = distanceMatrix[pdbIndex, clusterCenters[0]];
                        for (int clstrInd2 = 1; clstrInd2 < clusterCenters.Length; clstrInd2++)
                        {
                            if (clusterCenters[clstrInd2] > -1)
                            {
                                if (distanceMatrix[pdbIndex, clusterCenters[clstrInd2]] < newCenterProbabilities[pdbIndex])
                                {
                                    newCenterProbabilities[pdbIndex] = distanceMatrix[pdbIndex, clusterCenters[clstrInd2]];
                                }
                            }
                        }
                        normalizationFactor += newCenterProbabilities[pdbIndex];
                        runningTotal[pdbIndex] = normalizationFactor;
                    }
                    // generate a random number*normalization factor, look for it in set
                    double randomValForNextCenter = new double();
                    randomValForNextCenter = normalizationFactor * random.NextDouble();
                    for (int pdbIndex = 1; pdbIndex < clusterData.Length; pdbIndex++)
                    {
                        if ((randomValForNextCenter > runningTotal[pdbIndex - 1]) &&
                            (randomValForNextCenter >= runningTotal[pdbIndex]))
                        {
                            clusterCenters[clusterIndex] = pdbIndex;
                        }
                    }
                }*/

                // matrix of distance sums needed anyway

                double[] distanceSums = new double[_dihedralsByPDB.Count];

                // iterate until clusters don't change:
                bool haveClustersChanged = new bool();
                haveClustersChanged = true;
                // in loop: *calc distance sums, *find minimum distance sum & assign center, *reassign cluster IDs
                while (haveClustersChanged)
                {
                    int[] previousIterationClusterData = new int[_dihedralsByPDB.Count];
                    previousIterationClusterData = clusterData;

                    // re-assign cluster assignments if necessary
                    haveClustersChanged = false;
                    for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                    {
                        for (int clusterCounter = 0; clusterCounter < _numberOfClusters; clusterCounter++)
                        {
                            if (clusterData[pdbIndex] == -1) // move into first empty cluster, if any.  Otherwise, move into current cluster
                            {
                                clusterData[pdbIndex] = clusterCounter;
                            }
                            else
                            {
                                if (clusterCenters[clusterCounter] > -1)
                                {
                                    if (distanceMatrix[pdbIndex, clusterCenters[clusterData[pdbIndex]]] > distanceMatrix[pdbIndex, clusterCenters[clusterCounter]])
                                    {
                                        clusterData[pdbIndex] = clusterCounter;
                                    }
                                }
                            }
                        }
                    }

                    // calculate distance sums
                    for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                    {
                        distanceSums[pdbIndex] = 0;
                        for (int secondPdb = 0; secondPdb < clusterData.Length; secondPdb++)
                        {
                            if (clusterData[pdbIndex] == clusterData[secondPdb])
                            {
                                distanceSums[pdbIndex] += distanceMatrix[pdbIndex, secondPdb];
                            }
                        }
                    }
                    // reassign center of each cluster at current medoid
                    for (int clusterIndex = 0; clusterIndex < clusterCenters.Length; clusterIndex++)
                    {
                        for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                        {
                            if (clusterData[pdbIndex] == clusterIndex)
                            {
                                if (clusterCenters[clusterIndex] > -1)
                                {
                                    if (distanceSums[clusterCenters[clusterIndex]] > distanceSums[pdbIndex])
                                    {
                                        clusterCenters[clusterIndex] = pdbIndex;
                                    }
                                }
                                else
                                {
                                    clusterCenters[clusterIndex] = pdbIndex;
                                }
                            }
                        }
                    }

                    // remove outliers from their clusters

                    // new outlier method: exclude those more than 90 degrees away from any phi,psi,omega value of center value
                    if (_withOutliers)
                    {
                        double degreeThreshold = new double();
                        degreeThreshold = 45;
                        for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                        {
                            int memberCluster = new int();
                            memberCluster = clusterData[pdbIndex];
                            if (memberCluster > -1)
                            {
                                if (FindOutlierByDegrees(ref _dihedralsByPDB, pdbIndex, clusterCenters[memberCluster], degreeThreshold))
                                {
                                    clusterData[pdbIndex] = -1;
                                }
                            }
                        }
                    }
                    

                    if (previousIterationClusterData == clusterData)
                    {
                        haveClustersChanged = false;
                    }
                    else
                    {
                        haveClustersChanged = true;
                    }
                }
                // calculate within-cluster scatter; if better than best previous, replace value
                double currentClusterScatter = new double();
                currentClusterScatter = 0;
                for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
                {
                    if (clusterData[pdbIndex] > -1)
                    {
                        currentClusterScatter += distanceMatrix[pdbIndex, clusterCenters[clusterData[pdbIndex]]];
                    }
                }
                if (iterationCounter == 0)
                {
                    bestClusterData = clusterData;
                    bestClusterCenters = clusterCenters;
                    withinClusterScatter = currentClusterScatter;
                }
                else
                {
                    if (currentClusterScatter < withinClusterScatter)
                    {
                        bestClusterData = clusterData;
                        bestClusterCenters = clusterCenters;
                        withinClusterScatter = currentClusterScatter;
                    }
                }
            }
            // prepare the data for return
            for (int clusterIndex = 0; clusterIndex < bestClusterCenters.Length; clusterIndex++)
            {
                // pass back indices instead
                // makes life much easier when formatting changes
                _centerOfEachCluster.Add(bestClusterCenters[clusterIndex]);
            }
            _scatterResults = withinClusterScatter;
            return bestClusterData;
        }

        bool FindOutlierByDegrees(ref List<myLoopDataObject> _dihedralsByPdb, int _pdbIndex, int _centerIndex, double _degreeThreshold)
        {
            ArrayList dihedralsOfPdb = _dihedralsByPdb[_pdbIndex].myDihValues;
            ArrayList dihedralsOfCenter = _dihedralsByPdb[_centerIndex].myDihValues;
            for (int dihIndex = 0; dihIndex < dihedralsOfPdb.Count; dihIndex++)
            {
                double maxValue = new double();
                double minValue = new double();
                maxValue = (double)dihedralsOfCenter[dihIndex] + _degreeThreshold;
                minValue = (double)dihedralsOfCenter[dihIndex] - _degreeThreshold;
                if (Math.Pow(((double)dihedralsOfPdb[dihIndex] - (double)dihedralsOfCenter[dihIndex]), 2) > Math.Pow(_degreeThreshold,2))
                {
                    if (!((maxValue > 180) && ((double)dihedralsOfPdb[dihIndex] < maxValue - 360)))
                    {
                        if (!((minValue < -180) && ((double)dihedralsOfPdb[dihIndex] > minValue + 360)))
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        public List<double> TransformACoordWithTheseusInput(ref List<double> _position, ref List<List<double>> _rotMatrix, ref List<double> _transVector)
        {
            List<double> returnValues = new List<double>();
            List<double> workingVal = new List<double>();
            workingVal.Add(_position[0]);
            workingVal.Add(_position[1]);
            workingVal.Add(_position[2]);
            workingVal[0] -= _transVector[0];
            workingVal[1] -= _transVector[1];
            workingVal[2] -= _transVector[2];
            returnValues.Add(workingVal[0] * _rotMatrix[0][0] + workingVal[1] * _rotMatrix[1][0] + workingVal[2] * _rotMatrix[2][0]);
            returnValues.Add(workingVal[0] * _rotMatrix[0][1] + workingVal[1] * _rotMatrix[1][1] + workingVal[2] * _rotMatrix[2][1]);
            returnValues.Add(workingVal[0] * _rotMatrix[0][2] + workingVal[1] * _rotMatrix[1][2] + workingVal[2] * _rotMatrix[2][2]);
            return returnValues;
        }

        public double[] CalculateAngleStatistics(ArrayList _setOfAngles)
        {
            double[] statsToReturn = new double[4]; // [0]: average value
                                                       // [1]: circular variance
                                                       // [2]: circular standard deviation
                                                       // [3]: dispersion
            if(_setOfAngles.Count==0)
            {
                Console.WriteLine("Tried to calculate avg value of empty set.  Exiting PeptideAnalysis.CalculateAngleStatistics");
                Console.ReadLine();
                statsToReturn[0] = 999;
                return statsToReturn;
            }
            
            double cValue = new double();
            double sValue = new double();
            double rValue = new double();
            double dispConstant = new double();

            cValue = 0;
            sValue = 0;
            dispConstant = 0;

            for (int angleIndex = 0; angleIndex < _setOfAngles.Count; angleIndex++)
            {
                cValue += Math.Cos((double)_setOfAngles[angleIndex] * Math.PI / (double)180);
                sValue += Math.Sin((double)_setOfAngles[angleIndex] * Math.PI / (double)180);
            }

            cValue = cValue / ((double)_setOfAngles.Count);
            sValue = sValue / ((double)_setOfAngles.Count);

            rValue = Math.Sqrt(Math.Pow(cValue, (double)2) + Math.Pow(sValue, (double)2));

            statsToReturn[0] = Math.Atan2(sValue, cValue) * (double)180 / Math.PI;
            statsToReturn[1] = (double)1 - rValue;
            statsToReturn[2] = Math.Sqrt((double)-2 * Math.Log(rValue));

            for (int angleIndex = 0; angleIndex < _setOfAngles.Count; angleIndex++)
            {
                dispConstant += Math.Cos(((double)_setOfAngles[angleIndex] - statsToReturn[0]) * Math.PI / (double)180);
            }
            dispConstant = dispConstant / ((double)_setOfAngles.Count);
            statsToReturn[3] = Math.Acos(dispConstant) * (double)180 / Math.PI;
            
            return statsToReturn;
        }

        public double CalculateTransFraction(ArrayList _setOfOmegaAngles)
        {
            int transCount = new int();
            transCount = 0;
            for (int angleIndex = 0; angleIndex < _setOfOmegaAngles.Count; angleIndex++)
            {
                if (((double)_setOfOmegaAngles[angleIndex] > (double)90) || ((double)_setOfOmegaAngles[angleIndex] < (double)-90))
                {
                    transCount++;
                }
            }
            return (double)transCount / (double)_setOfOmegaAngles.Count;
        }

        /// <summary>
        /// Clusters by affinity method
        /// </summary>
        /// <param name="_numberOfClusters">number of clusters in sim: passed back</param>
        /// <param name="_dihedralsByPDB">ArrayList of ArrayLists (one per PDB)
        /// each PDB arraylist contains:
        /// [0]: pdbname
        /// [1]: ArrayList of dihedrals(phi0,psi0,omega0,phi1,psi1,omega1,...)</param>
        /// [2]: sequence (for data collection later)
        /// <param name="_scatterResults">passes scatter results back to main program or calling function</param>
        /// <returns>int[] of cluster membership by PDB index</returns>
        public void ClusterByAffinity(ref int _numberOfClusters, ref List<myLoopDataObject> _loopObjsByPDB, int _maxIterationCount, ref double _scatterResults, ref ArrayList _centerOfEachCluster, ref ArrayList _centerForEachPoint, ref ArrayList _LinesToWrite, double _selfSimFactor, string _progFileName)
        {
            bool useSimulatedAnnealing = new bool();
            useSimulatedAnnealing = false;
            Console.WriteLine("Started clustering at {0}", DateTime.Now);
            int[] bestClusterData = new int[_loopObjsByPDB.Count];
            int[] centersForEachPoint = new int[_loopObjsByPDB.Count];
            int[] centersForEachPointOnPreviousIteration = new int[_loopObjsByPDB.Count];
            //int numberOfIterationsWithConstantCenterAssignments = new int();
            //numberOfIterationsWithConstantCenterAssignments = 0;
            //int numberOfIterationsWithConstantTotalClusters = new int();
            //numberOfIterationsWithConstantTotalClusters = 0;
            int totalNumberOfClusters = new int();
            totalNumberOfClusters = 0;
            int previousTotalNumberOfClusters = new int();
            string doNotWriteFlag = "donotwrite";
            bool writeProgFiles = new bool();
            writeProgFiles = true;
            int lastExecutedIteration = new int();
            lastExecutedIteration = 0;
            if (_progFileName == doNotWriteFlag)
            {
                writeProgFiles = false;
            }
            ArrayList linesToWrite = new ArrayList();
            //bool needToBreak = new bool();
            //needToBreak = false;
            bool reportOnIteration = new bool();
            //bool firstReportOfConstCenters = new bool();
            //firstReportOfConstCenters = false;
            //bool firstReportOfConstClusters = new bool();
            //firstReportOfConstClusters = false;
            
            // setup: calculate matrix of distances
            // need: *list of pdbs, dihedrals
            // to create distance matrix, need to pass ArrayList of ArrayLists (one per dihedral-position)[loop0][loop1]
            //double[,] distanceMatrix = new double[_dihedralsByPDB.Count, _dihedralsByPDB.Count];
            //double[,] similarityMatrix = new double[_dihedralsByPDB.Count, _dihedralsByPDB.Count];
            
            //testing with number line here!!
            
            double[,] distanceMatrix = new double[_loopObjsByPDB.Count,_loopObjsByPDB.Count];
            double[,] similarityMatrix = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];

            distanceMatrix = GenerateDistanceMatrix(ref _loopObjsByPDB);
            
            
            // ************** opening up a temporary write file (pdbs,centers)
            ArrayList progFileLinesToWrite = new ArrayList();

            //*************** end of write file block
            
            
            //distanceMatrix = GenerateDistanceMatrix(_dihedralsByPDB);

            // now, calculate the similarities
            // self-similarity can be an adjustable parameter
            double selfSimilarityConstant = new double();
            int countForAverage = new int();
            countForAverage = 0;
            selfSimilarityConstant = 0;
            for (int simIndex = 0; simIndex < _loopObjsByPDB.Count; simIndex++)
            {
                for (int simIndex2 = 0; simIndex2 < _loopObjsByPDB.Count; simIndex2++)
                {
                    similarityMatrix[simIndex, simIndex2] = -1 * distanceMatrix[simIndex, simIndex2];
                    if (simIndex != simIndex2)
                    {
                        selfSimilarityConstant += similarityMatrix[simIndex, simIndex2];
                        countForAverage++;
                    }
                }
            }
            selfSimilarityConstant = _selfSimFactor*selfSimilarityConstant / countForAverage;
            for (int simIndex = 0; simIndex < _loopObjsByPDB.Count; simIndex++)
            {
                similarityMatrix[simIndex, simIndex] = selfSimilarityConstant;
            }

            // more setup: initialize availabilities, set update parameter lambda, declare responsibility matrix

            double[,] availabilities = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];
            double[,] responsibilities = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];
            double updateLambda = new double();
            updateLambda = 0.5;
            for (int avIndex = 0; avIndex < _loopObjsByPDB.Count; avIndex++)
            {
                for (int avIndex2 = 0; avIndex2 < _loopObjsByPDB.Count; avIndex2++)
                {
                    availabilities[avIndex, avIndex2] = 0;
                    responsibilities[avIndex, avIndex2] = 0;
                }
            }

            updateLambda = 0.5;
            
            //double withinClusterScatter = new double();
            bool haveSeenIdenticalCenterAssignments = new bool();
            haveSeenIdenticalCenterAssignments = false;
            bool reportNextClusters = new bool();
            reportNextClusters = false;
            int numberOfIterationsWithConstantClusterAssignments = new int();
            numberOfIterationsWithConstantClusterAssignments = 0;
            int[] previousClusterAssignments = new int[_loopObjsByPDB.Count];
            for (int loopIndex = 0; loopIndex < _loopObjsByPDB.Count; loopIndex++)
            { previousClusterAssignments[loopIndex] = 0; }

            for (int iterationCounter = 0; iterationCounter < _maxIterationCount; iterationCounter++)
            {
                centersForEachPoint.CopyTo(centersForEachPointOnPreviousIteration, 0); // deep copy?
                previousTotalNumberOfClusters = totalNumberOfClusters;
                if (iterationCounter % 20 == 0)
                { reportOnIteration = true; }
                else
                { reportOnIteration = false; }

                // relaxing parameter
                // doesn't do anything with clustering of AB conformations, because they cluster before 250 anyways
                // (and there the max iteration count is 500)
                if ((useSimulatedAnnealing) && ((iterationCounter > _maxIterationCount / 2) || haveSeenIdenticalCenterAssignments))
                {
                    updateLambda = updateLambda * (double)0.95;
                }

                //step 1: update responsibilities
                if (reportOnIteration)
                { Console.WriteLine("Started algorithm loop {0} at {1}", iterationCounter, DateTime.Now); }
                // find the value of k' (!=k) that maxs the quantity (a(i,k') + s(i,k')) for fixed i
                
                double maxValueForResp = new double();
                for (int respCounter = 0; respCounter < _loopObjsByPDB.Count; respCounter++) // loop over i
                {
                    for (int cenCounter = 0; cenCounter < _loopObjsByPDB.Count; cenCounter++)
                    {
                        bool hasMaxValueBeenSet = new bool();
                        hasMaxValueBeenSet = false;

                        for (int findingMaxInd = 0; findingMaxInd < _loopObjsByPDB.Count; findingMaxInd++)
                        {
                            if (findingMaxInd == cenCounter)
                            {
                                continue;
                            }
                            double currentQuant = new double();
                            currentQuant = similarityMatrix[respCounter, findingMaxInd];
                            if (respCounter != cenCounter)
                            {
                                currentQuant += availabilities[respCounter, findingMaxInd];
                            }
                            if (hasMaxValueBeenSet)
                            {
                                if (currentQuant > maxValueForResp)
                                {
                                    maxValueForResp = currentQuant;
                                }
                            }
                            else
                            {
                                maxValueForResp = currentQuant;
                                hasMaxValueBeenSet = true;
                            }
                        }
                        responsibilities[respCounter, cenCounter] = (1 - updateLambda) * responsibilities[respCounter, cenCounter] + 
                            updateLambda * (similarityMatrix[respCounter, cenCounter] - maxValueForResp); // check this if using annealing: lambda versus 1-lambda backwards?
                    }
                }

                // step 2: update availabilities
                if (reportOnIteration)
                { Console.WriteLine("Finished responsibilities, starting availabilities at {0}", DateTime.Now); }
                for (int avCtr = 0; avCtr < _loopObjsByPDB.Count; avCtr++)
                {
                    for (int cenCtr = 0; cenCtr < _loopObjsByPDB.Count; cenCtr++)
                    {
                        if (avCtr == cenCtr)
                        {
                            double accumValueForSelfAv = new double();
                            accumValueForSelfAv = 0;
                            for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                            {
                                if (pointCtr == cenCtr)
                                {
                                    continue;
                                }
                                if (responsibilities[pointCtr, cenCtr] < 0)
                                {
                                    accumValueForSelfAv += responsibilities[pointCtr, cenCtr];
                                }
                            }
                            availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr] +
                                updateLambda * accumValueForSelfAv; 
                        }
                        else
                        {
                            double accumulateQuantity = new double();
                            accumulateQuantity = responsibilities[cenCtr, cenCtr];
                            for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                            {
                                if ((pointCtr == avCtr) || (pointCtr == cenCtr))
                                {
                                    continue;
                                }
                                if (responsibilities[pointCtr, cenCtr] > 0)
                                {
                                    accumulateQuantity += responsibilities[pointCtr, cenCtr];
                                }
                            }
                            if (accumulateQuantity > 0)
                            {
                                availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr];
                            }
                            else
                            {
                                availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr] + 
                                    updateLambda * accumulateQuantity;
                            }
                        }
                    }
                }

                // Step 3: compute centers
                // debugging variables here
                // define clustering after this step?

                if (reportOnIteration)
                { Console.WriteLine("Finished availabilities, computing centers at {0}", DateTime.Now); }
                for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                {
                    double maxValForSum = new double();
                    double testQuantity = new double();
                    int exemplarForThisPoint = new int();
                    exemplarForThisPoint = 0;
                    maxValForSum = availabilities[pointCtr, 0] + responsibilities[pointCtr, 0];
                    for (int cenCtr = 1; cenCtr < _loopObjsByPDB.Count; cenCtr++)
                    {
                        testQuantity = availabilities[pointCtr, cenCtr] + responsibilities[pointCtr, cenCtr];
                        if (testQuantity > maxValForSum)
                        {
                            maxValForSum = testQuantity;
                            exemplarForThisPoint = cenCtr;
                        }
                    }
                    
                    centersForEachPoint[pointCtr] = exemplarForThisPoint;
                }

                // see if the center/exemplar assignments are the same as the previous iteration


                
                //code block: convergence by constant cluster ID
                int[] thisIterClusterData = new int[_loopObjsByPDB.Count];
                bool thisIterClusteringFailed = new bool();
                int[] thisIterCenters = ClusterWithSorting(ref centersForEachPoint, ref thisIterClusterData, ref thisIterClusteringFailed);
                if (thisIterClusteringFailed && writeProgFiles)
                {
                    string oopsstring = "attempt at clustering in loop failed in affinity clustering for progfile " + _progFileName + " at iteration " + iterationCounter.ToString();
                    linesToWrite.Add(oopsstring);
                    progFileLinesToWrite.Add(oopsstring);
                }
                if (thisIterClusteringFailed) { Console.WriteLine("clustering failed in affinity clustering. error!"); }
                if ((iterationCounter < 200) || (iterationCounter % 1000 == 0) || (reportNextClusters))
                {
                    int debugNumberClusters = new int();
                    debugNumberClusters = thisIterCenters.Length;
                    string debugLine = iterationCounter.ToString() + "\t" + debugNumberClusters.ToString();
                    progFileLinesToWrite.Add(debugLine);
                    if (iterationCounter > 200)
                    {
                        if ((iterationCounter - 10) % 1000 == 0) { reportNextClusters = false; }
                        else { reportNextClusters = true; }
                    }
                }
                bool testedAsEqual = new bool();
                testedAsEqual = true;
                for (int dIndex = 0; dIndex < thisIterClusterData.Length; dIndex++)
                {
                    if (thisIterClusterData[dIndex] != previousClusterAssignments[dIndex])
                    { testedAsEqual = false; break; }
                }

                if (testedAsEqual)
                {
                    numberOfIterationsWithConstantClusterAssignments++;
                    if (numberOfIterationsWithConstantClusterAssignments > 3)
                    {
                        if (iterationCounter > 10)
                        { lastExecutedIteration = iterationCounter; break; }
                    }
                }
                else
                { numberOfIterationsWithConstantClusterAssignments = 0; previousClusterAssignments = thisIterClusterData; }

                // converge by constant centers
                /*
                areTheCentersConstant = true;
                for (int thePdbCtr = 0; thePdbCtr < _loopObjsByPDB.Count; thePdbCtr++)
                {
                    if (centersForEachPoint[thePdbCtr] != centersForEachPointOnPreviousIteration[thePdbCtr])
                    {
                        areTheCentersConstant = false;
                    }
                }
                if (areTheCentersConstant)
                {
                    Console.WriteLine("centers constant from previous iteration at it " + iterationCounter);
                    numberOfIterationsWithConstantCenterAssignments++;
                    haveSeenIdenticalCenterAssignments = true;
                    if (numberOfIterationsWithConstantCenterAssignments > 4)
                    { lastExecutedIteration = iterationCounter; break; }
                }
                else
                {
                    numberOfIterationsWithConstantCenterAssignments = 0;
                    if (centersForEachPoint == centersForEachPointOnPreviousIteration)
                    {
                        Console.WriteLine("centers erroneously equal at it = {0}", iterationCounter);
                        if (writeToProgressFile && writeProgFiles)
                        {
                            string progLineToWrite = "centers erroneously equal at it = " + iterationCounter.ToString();
                            progFileLinesToWrite.Add(progLineToWrite);
                        }
                    }
                }

                if ((numberOfIterationsWithConstantCenterAssignments > 0) && writeProgFiles)
                {
                    string progLineToWrite = "centers constant for " + numberOfIterationsWithConstantCenterAssignments.ToString() +
                        " iterations at iteration " + iterationCounter.ToString();
                    //progFileLinesToWrite.Add(progLineToWrite); // 1234567890 add this back!
                }
                */
                // end of converge by constant centers
                // output for FS
                
                if ((iterationCounter == _maxIterationCount - 1) && writeProgFiles)
                {
                    string progLine = "Hit max iterations at it " + iterationCounter.ToString() +
                        ", did not converge according to test, exiting.";
                    progFileLinesToWrite.Add(progLine);
                    lastExecutedIteration = iterationCounter;
                }
            }

            // get _centerForEachPoint ready to return

            foreach (int intCenter in centersForEachPoint)
            {
                _centerForEachPoint.Add(intCenter);
            }
            // construct clusters from exemplar data
            Console.WriteLine("Starting cluster construction at {0}", DateTime.Now);

            bool finalClusteringFailed = new bool();
            int[] centersForEachCluster = ClusterWithSorting(ref centersForEachPoint, ref bestClusterData, ref finalClusteringFailed);
            if (finalClusteringFailed)
            {
                string oopsstring = "attempt at clustering in wrapup failed in affinity clustering for progfile " + _progFileName;
                linesToWrite.Add(oopsstring);
                if (writeProgFiles)
                { progFileLinesToWrite.Add(oopsstring); }
            }
            if (writeProgFiles)
            {
                string concatFilename = @"F:\andreasAB\abClusteringNR\progFiles\concatProg.log";
                StreamWriter progFileW = new StreamWriter(concatFilename, true);
                progFileW.WriteLine("for run {0}", _progFileName);
                foreach (object lineObj in progFileLinesToWrite)
                {
                    progFileW.WriteLine((string)lineObj);
                }
                progFileW.WriteLine("Finished at iteration {0}", lastExecutedIteration);
                progFileW.Close();
                Console.WriteLine("{0} written.", _progFileName);
            }

            foreach (int clCenter in centersForEachCluster)
            {
                _centerOfEachCluster.Add(clCenter);
            }

            _numberOfClusters = centersForEachCluster.Length;

            // get the algorithm ready for output
           
            // compute scatter
            //_centerOfEachCluster
            _scatterResults = 0;
            for (int pointCtr = 0; pointCtr < bestClusterData.Length; pointCtr++)
            {
                if ((bestClusterData[pointCtr] < _centerOfEachCluster.Count) && (bestClusterData[pointCtr] > -1))
                {
                    _scatterResults += distanceMatrix[pointCtr, (int)_centerOfEachCluster[bestClusterData[pointCtr]]];
                }
            }
            _LinesToWrite = linesToWrite;
            Console.WriteLine("about to return at {0}", DateTime.Now);
            
            // sort clusters by length: fix bestClusterData (int[]) and _centerOfEachCluster (ArrayList)

            //return bestClusterData;
            for (int pdbIndex = 0; pdbIndex < _loopObjsByPDB.Count; pdbIndex++)
            {
                _loopObjsByPDB[pdbIndex].myClusterID = bestClusterData[pdbIndex];
            }
            return;
        }
        /// <summary>
        /// Allows for direct input of distance matrix
        /// </summary>
        /// <param name="_loopObjsByPDB">dummy loop objects: only name and clusterID are imporant</param>
        /// <param name="_distanceMatrix">type double[,]: must be same order as loop data object</param>
        /// <param name="_scatterResults">scatter value</param>
        /// <param name="_maxIterationCount">max # iterations</param>
        /// <param name="_selfSimFactor">selfSim scaling factor</param>
        public void ClusterByAffinityWithDistanceMatrixInput(ref List<myLoopDataObject> _loopObjsByPDB, ref double[,] _distanceMatrix, ref double _scatterResults, int _maxIterationCount, double _selfSimFactor)
        {
            int _numberOfClusters = new int();
            ArrayList _centerOfEachCluster = new ArrayList();
            ArrayList _centerForEachPoint = new ArrayList();
            ArrayList _LinesToWrite = new ArrayList();
            string _progFileName = "NOWRITE";



            bool useSimulatedAnnealing = new bool();
            useSimulatedAnnealing = false;
            Console.WriteLine("Started clustering at {0}", DateTime.Now);
            int[] bestClusterData = new int[_loopObjsByPDB.Count];
            int[] centersForEachPoint = new int[_loopObjsByPDB.Count];
            int[] centersForEachPointOnPreviousIteration = new int[_loopObjsByPDB.Count];
            //int numberOfIterationsWithConstantCenterAssignments = new int();
            //numberOfIterationsWithConstantCenterAssignments = 0;
            //int numberOfIterationsWithConstantTotalClusters = new int();
            //numberOfIterationsWithConstantTotalClusters = 0;
            int totalNumberOfClusters = new int();
            totalNumberOfClusters = 0;
            int previousTotalNumberOfClusters = new int();
            string doNotWriteFlag = "donotwrite";
            bool writeProgFiles = new bool();
            writeProgFiles = true;
            int lastExecutedIteration = new int();
            lastExecutedIteration = 0;
            if (_progFileName == doNotWriteFlag)
            {
                writeProgFiles = false;
            }
            ArrayList linesToWrite = new ArrayList();
            //bool needToBreak = new bool();
            //needToBreak = false;
            bool reportOnIteration = new bool();
            //bool firstReportOfConstCenters = new bool();
            //firstReportOfConstCenters = false;
            //bool firstReportOfConstClusters = new bool();
            //firstReportOfConstClusters = false;

            // setup: calculate matrix of distances
            // need: *list of pdbs, dihedrals
            // to create distance matrix, need to pass ArrayList of ArrayLists (one per dihedral-position)[loop0][loop1]
            //double[,] distanceMatrix = new double[_dihedralsByPDB.Count, _dihedralsByPDB.Count];
            //double[,] similarityMatrix = new double[_dihedralsByPDB.Count, _dihedralsByPDB.Count];

            //testing with number line here!!

            double[,] distanceMatrix = _distanceMatrix;
            double[,] similarityMatrix = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];

            // ************** opening up a temporary write file (pdbs,centers)
            ArrayList progFileLinesToWrite = new ArrayList();

            //*************** end of write file block


            //distanceMatrix = GenerateDistanceMatrix(_dihedralsByPDB);

            // now, calculate the similarities
            // self-similarity can be an adjustable parameter
            double selfSimilarityConstant = new double();
            int countForAverage = new int();
            countForAverage = 0;
            selfSimilarityConstant = 0;
            for (int simIndex = 0; simIndex < _loopObjsByPDB.Count; simIndex++)
            {
                for (int simIndex2 = 0; simIndex2 < _loopObjsByPDB.Count; simIndex2++)
                {
                    similarityMatrix[simIndex, simIndex2] = -1 * distanceMatrix[simIndex, simIndex2];
                    if (simIndex != simIndex2)
                    {
                        selfSimilarityConstant += similarityMatrix[simIndex, simIndex2];
                        countForAverage++;
                    }
                }
            }
            selfSimilarityConstant = _selfSimFactor * selfSimilarityConstant / countForAverage;
            for (int simIndex = 0; simIndex < _loopObjsByPDB.Count; simIndex++)
            {
                similarityMatrix[simIndex, simIndex] = selfSimilarityConstant;
            }

            // more setup: initialize availabilities, set update parameter lambda, declare responsibility matrix

            double[,] availabilities = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];
            double[,] responsibilities = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];
            double updateLambda = new double();
            updateLambda = 0.5;
            for (int avIndex = 0; avIndex < _loopObjsByPDB.Count; avIndex++)
            {
                for (int avIndex2 = 0; avIndex2 < _loopObjsByPDB.Count; avIndex2++)
                {
                    availabilities[avIndex, avIndex2] = 0;
                    responsibilities[avIndex, avIndex2] = 0;
                }
            }

            updateLambda = 0.5;

            //double withinClusterScatter = new double();
            bool haveSeenIdenticalCenterAssignments = new bool();
            haveSeenIdenticalCenterAssignments = false;
            bool reportNextClusters = new bool();
            reportNextClusters = false;
            int numberOfIterationsWithConstantClusterAssignments = new int();
            numberOfIterationsWithConstantClusterAssignments = 0;
            int[] previousClusterAssignments = new int[_loopObjsByPDB.Count];
            for (int loopIndex = 0; loopIndex < _loopObjsByPDB.Count; loopIndex++)
            { previousClusterAssignments[loopIndex] = 0; }

            for (int iterationCounter = 0; iterationCounter < _maxIterationCount; iterationCounter++)
            {
                centersForEachPoint.CopyTo(centersForEachPointOnPreviousIteration, 0); // deep copy?
                previousTotalNumberOfClusters = totalNumberOfClusters;
                if (iterationCounter % 20 == 0)
                { reportOnIteration = true; }
                else
                { reportOnIteration = false; }

                // relaxing parameter
                // doesn't do anything with clustering of AB conformations, because they cluster before 250 anyways
                // (and there the max iteration count is 500)
                if ((useSimulatedAnnealing) && ((iterationCounter > _maxIterationCount / 2) || haveSeenIdenticalCenterAssignments))
                {
                    updateLambda = updateLambda * (double)0.95;
                }

                //step 1: update responsibilities
                if (reportOnIteration)
                { Console.WriteLine("Started algorithm loop {0} at {1}", iterationCounter, DateTime.Now); }
                // find the value of k' (!=k) that maxs the quantity (a(i,k') + s(i,k')) for fixed i

                double maxValueForResp = new double();
                for (int respCounter = 0; respCounter < _loopObjsByPDB.Count; respCounter++) // loop over i
                {
                    for (int cenCounter = 0; cenCounter < _loopObjsByPDB.Count; cenCounter++)
                    {
                        bool hasMaxValueBeenSet = new bool();
                        hasMaxValueBeenSet = false;

                        for (int findingMaxInd = 0; findingMaxInd < _loopObjsByPDB.Count; findingMaxInd++)
                        {
                            if (findingMaxInd == cenCounter)
                            {
                                continue;
                            }
                            double currentQuant = new double();
                            currentQuant = similarityMatrix[respCounter, findingMaxInd];
                            if (respCounter != cenCounter)
                            {
                                currentQuant += availabilities[respCounter, findingMaxInd];
                            }
                            if (hasMaxValueBeenSet)
                            {
                                if (currentQuant > maxValueForResp)
                                {
                                    maxValueForResp = currentQuant;
                                }
                            }
                            else
                            {
                                maxValueForResp = currentQuant;
                                hasMaxValueBeenSet = true;
                            }
                        }
                        responsibilities[respCounter, cenCounter] = (1 - updateLambda) * responsibilities[respCounter, cenCounter] +
                            updateLambda * (similarityMatrix[respCounter, cenCounter] - maxValueForResp); // check this if using annealing: lambda versus 1-lambda backwards?
                    }
                }

                // step 2: update availabilities
                if (reportOnIteration)
                { Console.WriteLine("Finished responsibilities, starting availabilities at {0}", DateTime.Now); }
                for (int avCtr = 0; avCtr < _loopObjsByPDB.Count; avCtr++)
                {
                    for (int cenCtr = 0; cenCtr < _loopObjsByPDB.Count; cenCtr++)
                    {
                        if (avCtr == cenCtr)
                        {
                            double accumValueForSelfAv = new double();
                            accumValueForSelfAv = 0;
                            for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                            {
                                if (pointCtr == cenCtr)
                                {
                                    continue;
                                }
                                if (responsibilities[pointCtr, cenCtr] < 0)
                                {
                                    accumValueForSelfAv += responsibilities[pointCtr, cenCtr];
                                }
                            }
                            availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr] +
                                updateLambda * accumValueForSelfAv;
                        }
                        else
                        {
                            double accumulateQuantity = new double();
                            accumulateQuantity = responsibilities[cenCtr, cenCtr];
                            for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                            {
                                if ((pointCtr == avCtr) || (pointCtr == cenCtr))
                                {
                                    continue;
                                }
                                if (responsibilities[pointCtr, cenCtr] > 0)
                                {
                                    accumulateQuantity += responsibilities[pointCtr, cenCtr];
                                }
                            }
                            if (accumulateQuantity > 0)
                            {
                                availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr];
                            }
                            else
                            {
                                availabilities[avCtr, cenCtr] = (1 - updateLambda) * availabilities[avCtr, cenCtr] +
                                    updateLambda * accumulateQuantity;
                            }
                        }
                    }
                }

                // Step 3: compute centers
                // debugging variables here
                // define clustering after this step?

                if (reportOnIteration)
                { Console.WriteLine("Finished availabilities, computing centers at {0}", DateTime.Now); }
                for (int pointCtr = 0; pointCtr < _loopObjsByPDB.Count; pointCtr++)
                {
                    double maxValForSum = new double();
                    double testQuantity = new double();
                    int exemplarForThisPoint = new int();
                    exemplarForThisPoint = 0;
                    maxValForSum = availabilities[pointCtr, 0] + responsibilities[pointCtr, 0];
                    for (int cenCtr = 1; cenCtr < _loopObjsByPDB.Count; cenCtr++)
                    {
                        testQuantity = availabilities[pointCtr, cenCtr] + responsibilities[pointCtr, cenCtr];
                        if (testQuantity > maxValForSum)
                        {
                            maxValForSum = testQuantity;
                            exemplarForThisPoint = cenCtr;
                        }
                    }

                    centersForEachPoint[pointCtr] = exemplarForThisPoint;
                }

                // see if the center/exemplar assignments are the same as the previous iteration



                //code block: convergence by constant cluster ID
                int[] thisIterClusterData = new int[_loopObjsByPDB.Count];
                bool thisIterClusteringFailed = new bool();
                int[] thisIterCenters = ClusterWithSorting(ref centersForEachPoint, ref thisIterClusterData, ref thisIterClusteringFailed);
                if (thisIterClusteringFailed && writeProgFiles)
                {
                    string oopsstring = "attempt at clustering in loop failed in affinity clustering for progfile " + _progFileName + " at iteration " + iterationCounter.ToString();
                    linesToWrite.Add(oopsstring);
                    progFileLinesToWrite.Add(oopsstring);
                }
                if (thisIterClusteringFailed) { Console.WriteLine("clustering failed in affinity clustering. error!"); }
                if ((iterationCounter < 200) || (iterationCounter % 1000 == 0) || (reportNextClusters))
                {
                    int debugNumberClusters = new int();
                    debugNumberClusters = thisIterCenters.Length;
                    string debugLine = iterationCounter.ToString() + "\t" + debugNumberClusters.ToString();
                    progFileLinesToWrite.Add(debugLine);
                    if (iterationCounter > 200)
                    {
                        if ((iterationCounter - 10) % 1000 == 0) { reportNextClusters = false; }
                        else { reportNextClusters = true; }
                    }
                }
                bool testedAsEqual = new bool();
                testedAsEqual = true;
                for (int dIndex = 0; dIndex < thisIterClusterData.Length; dIndex++)
                {
                    if (thisIterClusterData[dIndex] != previousClusterAssignments[dIndex])
                    { testedAsEqual = false; break; }
                }

                if (testedAsEqual)
                {
                    numberOfIterationsWithConstantClusterAssignments++;
                    if (numberOfIterationsWithConstantClusterAssignments > 3)
                    {
                        if (iterationCounter > 10)
                        { lastExecutedIteration = iterationCounter; break; }
                    }
                }
                else
                { numberOfIterationsWithConstantClusterAssignments = 0; previousClusterAssignments = thisIterClusterData; }

                // converge by constant centers
                /*
                areTheCentersConstant = true;
                for (int thePdbCtr = 0; thePdbCtr < _loopObjsByPDB.Count; thePdbCtr++)
                {
                    if (centersForEachPoint[thePdbCtr] != centersForEachPointOnPreviousIteration[thePdbCtr])
                    {
                        areTheCentersConstant = false;
                    }
                }
                if (areTheCentersConstant)
                {
                    Console.WriteLine("centers constant from previous iteration at it " + iterationCounter);
                    numberOfIterationsWithConstantCenterAssignments++;
                    haveSeenIdenticalCenterAssignments = true;
                    if (numberOfIterationsWithConstantCenterAssignments > 4)
                    { lastExecutedIteration = iterationCounter; break; }
                }
                else
                {
                    numberOfIterationsWithConstantCenterAssignments = 0;
                    if (centersForEachPoint == centersForEachPointOnPreviousIteration)
                    {
                        Console.WriteLine("centers erroneously equal at it = {0}", iterationCounter);
                        if (writeToProgressFile && writeProgFiles)
                        {
                            string progLineToWrite = "centers erroneously equal at it = " + iterationCounter.ToString();
                            progFileLinesToWrite.Add(progLineToWrite);
                        }
                    }
                }

                if ((numberOfIterationsWithConstantCenterAssignments > 0) && writeProgFiles)
                {
                    string progLineToWrite = "centers constant for " + numberOfIterationsWithConstantCenterAssignments.ToString() +
                        " iterations at iteration " + iterationCounter.ToString();
                    //progFileLinesToWrite.Add(progLineToWrite); // 1234567890 add this back!
                }
                */
                // end of converge by constant centers
                // output for FS

                if ((iterationCounter == _maxIterationCount - 1) && writeProgFiles)
                {
                    string progLine = "Hit max iterations at it " + iterationCounter.ToString() +
                        ", did not converge according to test, exiting.";
                    progFileLinesToWrite.Add(progLine);
                    lastExecutedIteration = iterationCounter;
                }
            }

            // get _centerForEachPoint ready to return

            foreach (int intCenter in centersForEachPoint)
            {
                _centerForEachPoint.Add(intCenter);
            }
            // construct clusters from exemplar data
            Console.WriteLine("Starting cluster construction at {0}", DateTime.Now);

            bool finalClusteringFailed = new bool();
            int[] centersForEachCluster = ClusterWithSorting(ref centersForEachPoint, ref bestClusterData, ref finalClusteringFailed);
            if (finalClusteringFailed)
            {
                string oopsstring = "attempt at clustering in wrapup failed in affinity clustering for progfile " + _progFileName;
                linesToWrite.Add(oopsstring);
                if (writeProgFiles)
                { progFileLinesToWrite.Add(oopsstring); }
            }
            if (writeProgFiles)
            {
                string concatFilename = @"F:\andreasAB\abClusteringNR\progFiles\concatProg.log";
                StreamWriter progFileW = new StreamWriter(concatFilename, true);
                progFileW.WriteLine("for run {0}", _progFileName);
                foreach (object lineObj in progFileLinesToWrite)
                {
                    progFileW.WriteLine((string)lineObj);
                }
                progFileW.WriteLine("Finished at iteration {0}", lastExecutedIteration);
                progFileW.Close();
                Console.WriteLine("{0} written.", _progFileName);
            }

            foreach (int clCenter in centersForEachCluster)
            {
                _centerOfEachCluster.Add(clCenter);
            }

            _numberOfClusters = centersForEachCluster.Length;

            // get the algorithm ready for output

            // compute scatter
            //_centerOfEachCluster
            _scatterResults = 0;
            for (int pointCtr = 0; pointCtr < bestClusterData.Length; pointCtr++)
            {
                if ((bestClusterData[pointCtr] < _centerOfEachCluster.Count) && (bestClusterData[pointCtr] > -1))
                {
                    _scatterResults += distanceMatrix[pointCtr, (int)_centerOfEachCluster[bestClusterData[pointCtr]]];
                }
            }
            _LinesToWrite = linesToWrite;
            Console.WriteLine("about to return at {0}", DateTime.Now);

            // sort clusters by length: fix bestClusterData (int[]) and _centerOfEachCluster (ArrayList)

            //return bestClusterData;
            for (int pdbIndex = 0; pdbIndex < _loopObjsByPDB.Count; pdbIndex++)
            {
                _loopObjsByPDB[pdbIndex].myClusterID = bestClusterData[pdbIndex];
            }
            return;
        }

        public ArrayList MakingCluster(ref int[] _centersByPoint, ref ArrayList _pointsByCenter, 
            ref bool[] _pointsAdded, int _pointToAdd, bool _wholeCluster)
        {
            ArrayList pointsInThisCycle = new ArrayList();
            if (!_pointsAdded[_pointToAdd])
            {
                _pointsAdded[_pointToAdd] = true;
                pointsInThisCycle.Add(_pointToAdd);
                ArrayList pointsFromExemplar = new ArrayList();
                pointsFromExemplar = MakingCluster(ref _centersByPoint, ref _pointsByCenter,
                    ref _pointsAdded, _centersByPoint[_pointToAdd], false);
                foreach (object pObject in pointsFromExemplar)
                {
                    if (!pointsInThisCycle.Contains(pObject))
                    {
                        pointsInThisCycle.Add(pObject);
                    }
                }
                ArrayList pointsAsExemplar = new ArrayList();
                pointsAsExemplar = ((ArrayList)_pointsByCenter[_pointToAdd]);
                for (int exPointIndex = 0; exPointIndex < pointsAsExemplar.Count; exPointIndex++)
                {
                    ArrayList pointsFromThisRound = new ArrayList();
                    pointsFromThisRound = MakingCluster(ref _centersByPoint, ref _pointsByCenter, 
                        ref _pointsAdded, (int)pointsAsExemplar[exPointIndex], false);
                    foreach (object pObject in pointsFromThisRound)
                    {
                        if (!pointsInThisCycle.Contains(pObject))
                        {
                            pointsInThisCycle.Add(pObject);
                        }
                    }
                }

            }
            if (_wholeCluster)
            {
                int[] pointsToSort = new int[pointsInThisCycle.Count];
                for (int pointIndex = 0; pointIndex < pointsInThisCycle.Count; pointIndex++)
                {
                    pointsToSort[pointIndex] = (int)pointsInThisCycle[pointIndex];
                }
                Array.Sort(pointsToSort);
                for (int pointIndex = 0; pointIndex < pointsInThisCycle.Count; pointIndex++)
                {
                    pointsInThisCycle[pointIndex] = pointsToSort[pointIndex];
                }
            }
            return pointsInThisCycle;
        }

        /// <summary>
        /// Clusters by new method.
        /// </summary>
        /// <param name="_centersByPoint"></param>
        /// <param name="_rawClusterData">cluster membership (as integer) indexed by point.  Value of -100:unassigned</param>
        /// <returns>int[] of centers, one per cluster</returns>
        public int[] ClusterByNewMethod(ref int[] _centersByPoint, ref int[] _rawClusterData, ref bool _clusteringFailed)
        {
            ArrayList returnVariable = new ArrayList();
            ArrayList groupedIntoClusters = new ArrayList();
            List<clusterContainer> containerList = new List<clusterContainer>();
            containerList = FromCenterListToClusters(ref _centersByPoint, ref _clusteringFailed);
            int[] centersByCluster = new int[containerList.Count];
            for (int clusterIndex = 0; clusterIndex < containerList.Count; clusterIndex++)
            {
                int[] infoForOneCluster = new int[containerList[clusterIndex].memberPDBs.Count];
                for (int memberIndex = 0; memberIndex < containerList[clusterIndex].memberPDBs.Count; memberIndex++)
                {
                    infoForOneCluster[memberIndex] = (int)containerList[clusterIndex].memberPDBs[memberIndex];
                }
                groupedIntoClusters.Add(infoForOneCluster);
                centersByCluster[clusterIndex] = containerList[clusterIndex].centerPDBValue;
            }

            // ****************
            int[] clusterData = new int[_centersByPoint.Length];
            ArrayList checkCenterOfEachCluster = new ArrayList();
            for (int pdbIndex = 0; pdbIndex < clusterData.Length; pdbIndex++)
            {
                _rawClusterData[pdbIndex] = -100; // to look for unassigned values
            }

            for (int listIndex = 0; listIndex < containerList.Count; listIndex++)
            {
                for (int memberIndex = 0; memberIndex < containerList[listIndex].memberPDBs.Count; memberIndex++)
                {
                    if (((int)containerList[listIndex].memberPDBs[memberIndex] > -1) &&
                        ((int)containerList[listIndex].memberPDBs[memberIndex] < _rawClusterData.Length))
                    {
                        _rawClusterData[(int)containerList[listIndex].memberPDBs[memberIndex]] = listIndex;
                    }
                }
            }

            // ****************
            return centersByCluster;
        }

        public int[] ClusterWithSorting(ref int[] _centersByPoint, ref int[] _rawClusterData, ref bool _clusteringFailed)
        {
            int[] centersByCluster = ClusterByNewMethod(ref _centersByPoint, ref _rawClusterData, ref _clusteringFailed);
            // for each cluster, sort the membership.  Sort by cluster length, then by top value
            // assumes no non-assigned points
            // number of clusters == centersByCluster.Count
            List<List<int>> clustersToSort = new List<List<int>>();
            for (int clIndex = 0; clIndex < centersByCluster.Length; clIndex++)
            {
                List<int> oneCluster = new List<int>();
                clustersToSort.Add(oneCluster);
            }
            for (int pIndex = 0; pIndex < _rawClusterData.Length; pIndex++)
            {
                clustersToSort[_rawClusterData[pIndex]].Add(pIndex);
            }
            List<int> firstIndexForKeepingTrack = new List<int>(); // this records the first index for each cluster pre-sorting
            for (int clIndex = 0; clIndex < centersByCluster.Length; clIndex++)
            { firstIndexForKeepingTrack.Add(clustersToSort[clIndex][0]); }

            clustersToSort.Sort(new sortClusters());

            // now: re-assign _rawClusterData and centersByCluster.  _centersByPoint is unchanged

            int[] newCenterVector = new int[centersByCluster.Length];
            for (int clIndex = 0; clIndex < centersByCluster.Length; clIndex++)
            {
                // find clustersToSort[clIndex][0] in firstIndex record.  that was the original clusterIndex
                newCenterVector[clIndex] = centersByCluster[firstIndexForKeepingTrack.IndexOf(clustersToSort[clIndex][0])];
            }
            for (int clIndex = 0; clIndex < centersByCluster.Length; clIndex++)
            {
                centersByCluster[clIndex] = newCenterVector[clIndex];
                for (int ctIndex = 0; ctIndex < (clustersToSort[clIndex]).Count; ctIndex++)
                {
                    _rawClusterData[clustersToSort[clIndex][ctIndex]] = clIndex;
                }
            }
            return centersByCluster;
        }

        private List<clusterContainer> FromCenterListToClusters(ref int[] _centersByPoint, ref bool _clusteringFailed)
        {
            //ArrayList clusters = new ArrayList();
            
            ArrayList pointsByCenter = new ArrayList();
            List<clusterContainer> clusters = new List<clusterContainer>();


            for (int pointIndex = 0; pointIndex < _centersByPoint.Length; pointIndex++)
            {
                List<int> pointsList = new List<int>();
                pointsByCenter.Add(pointsList);
            }

            for (int pointIndex = 0; pointIndex < _centersByPoint.Length; pointIndex++)
            {
                if (_centersByPoint[pointIndex] > -1 && _centersByPoint[pointIndex] < _centersByPoint.Length)
                {
                    ((List<int>)pointsByCenter[_centersByPoint[pointIndex]]).Add(pointIndex);
                }
            }
            int clusterContainerCount = new int();
            clusterContainerCount = 0;

            for (int pointIndex = 0; pointIndex < pointsByCenter.Count; pointIndex++)
            {
                List<int> currentList = (List<int>)pointsByCenter[pointIndex];
                if (currentList.Count > 0)
                {
                    if (!currentList.Contains(pointIndex))
                    { currentList.Add(pointIndex); }
                    int[] centerAndItsCenter = new int[2];
                    centerAndItsCenter[0] = pointIndex;
                    centerAndItsCenter[1] = _centersByPoint[pointIndex];
                    clusterContainer myClusterContainer = new clusterContainer(currentList, centerAndItsCenter, clusterContainerCount);
                    clusters.Add(myClusterContainer);
                    clusterContainerCount++;
                }
            }

            // sort by length

            clusters.Sort(new comparingClusters());
            clusters.Reverse();

            // merge
            // if a cluster contains the center of another cluster, merge the two

            // added merge code starts here

            // merge code ends here
            bool continueTheMerging = new bool();
            continueTheMerging = true;
            while (continueTheMerging)
            {
                continueTheMerging = false;
                for (int clusterIndex1 = 0; clusterIndex1 < clusters.Count; clusterIndex1++)
                {
                    for (int clusterIndex2 = clusterIndex1 + 1; clusterIndex2 < clusters.Count; clusterIndex2++)
                    {
                        // add the new merging criterion here
                        List<int[]> candidates1 = clusters[clusterIndex1].mergeCandidatesWithCenters;
                        List<int[]> candidates2 = clusters[clusterIndex2].mergeCandidatesWithCenters;
                        if (candidates1.Count > 0)
                        {
                            foreach (int[] cPair in candidates1)
                            {
                                if (clusters[clusterIndex2].memberPDBs.Contains(cPair[1]))
                                {
                                    continueTheMerging = true;
                                    break;
                                }
                            }
                        }
                        if (!continueTheMerging)
                        {
                            if (candidates2.Count > 0)
                            {
                                foreach (int[] cPair in candidates2)
                                {
                                    if (clusters[clusterIndex1].memberPDBs.Contains(cPair[1]))
                                    {
                                        continueTheMerging = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (continueTheMerging)
                        {
                            clusters[clusterIndex1].mergeWithCluster(clusters[clusterIndex2]);
                            clusters.RemoveAt(clusterIndex2);
                        }
                        if (continueTheMerging) // if a merge has occured, restart the algorithm
                        { break; }
                    }
                    if (continueTheMerging)
                    { break; }
                }
            }

            // resort

            clusters.Sort(new comparingClusters());
            clusters.Reverse();

            // now test!!

            // test that for each point, only one cluster contains the point
            // also, test that for each point, its center is in the cluster

            _clusteringFailed = false;

            for (int cli1 = 0; cli1 < clusters.Count; cli1++)
            {
                List<int> members1 = clusters[cli1].memberPDBs;
                foreach (int aMember in members1)
                {
                    if (aMember < 0 || aMember >= _centersByPoint.Length)
                    { _clusteringFailed = true; break; }
                    if (!members1.Contains(_centersByPoint[aMember]))
                    { _clusteringFailed = true; break; }
                }
                
                for (int cli2 = cli1 + 1; cli2 < clusters.Count; cli2++)
                {
                    if (_clusteringFailed)
                    { break; }
                    List<int> members2 = clusters[cli2].memberPDBs;
                    foreach (int aMember in members1)
                    {
                        if (members2.Contains(aMember))
                        { _clusteringFailed = true; }
                    }
                }
                if (_clusteringFailed)
                { break; }
            }

            return clusters;
        }
        public void GunzipXMLfiles(List<string> _pdbList, string _sourceDir, string _outputDir)
        {
            foreach (string pdbFile in _pdbList)
            {
                System.Diagnostics.Process gunzipProc = new System.Diagnostics.Process();
                gunzipProc.EnableRaisingEvents = false;
                gunzipProc.StartInfo.FileName = "\"C:\\Program Files (x86)\\FCCC\\BioDownloader\\Auxiliary\\gz-decompress-minigzip.bat\"";
                string argsForGunzip = "\"";
                argsForGunzip += _sourceDir + "\\";
                argsForGunzip += pdbFile;
                argsForGunzip += ".xml.gz\" \"" + _outputDir + "\\\"";
                gunzipProc.StartInfo.Arguments = argsForGunzip;
                gunzipProc.Start();
                gunzipProc.WaitForExit();
            }
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
        // data in the output: 
        // [0]: pdbname
        // [1]: quality stats (res, bfac, confE)
        // [2]: sequence
        // [3]: cluster membership (int)
        // [4]: dihedral values
        [Serializable]
        public class myLoopDataObject : ICloneable
        {
            public myLoopDataObject()
            {
                myPdbIDobj = null;
                myResolution = 999;
                myBfac = 999;
                myConfE = 999;
                mySeqObj = null;
                myClusterID = -999;
                myDihValues = null;
                myLoopIDobj = null;
            }
            public myLoopDataObject(string _pdbID, string _loopID, string _seq, int _clusterID, double _res, double _bfac, double _confE, ArrayList _dihValues)
            {
                myPdbIDobj = (object)_pdbID;
                myLoopIDobj = (object)_loopID;
                mySeqObj = (object)_seq;
                myClusterID = _clusterID;
                myResolution = _res;
                myBfac = _bfac;
                myConfE = _confE;
                myDihValues = _dihValues;
            }
            public string PrintMyLoopDataObject()
            {
                string returnString = null;
                returnString = (string)myPdbIDobj + "," + myClusterID + "," + (string)mySeqObj + "," + myResolution + "," + myBfac + "," + myConfE;
                foreach (object dih in myDihValues)
                {
                    double theValue = (double)dih;
                    string valToPrint = theValue.ToString("#0.00");
                    returnString += "," + valToPrint;
                }
                return returnString;
            }
            public string PrintMyLoopDataObjectShortVersion()
            {
                string returnString = null;
                returnString = (string)myPdbIDobj + "," + myClusterID + "," + (string)mySeqObj + "," + myResolution + "," + myBfac + "," + myConfE;
                return returnString;
            }
            // member variables
            public object myPdbIDobj = new object();
            public double myResolution = new double();
            public double myBfac = new double();
            public double myConfE = new double();
            public object mySeqObj = new object();
            public int myClusterID = new int();
            public ArrayList myDihValues = new ArrayList();
            public object myLoopIDobj = new object();


            public object Clone() // deep copy
            {
                MemoryStream ms = new MemoryStream();
                BinaryFormatter loopBF = new BinaryFormatter();
                loopBF.Serialize(ms, this);
                ms.Position = 0;
                myLoopDataObject newLoopObj = (myLoopDataObject)loopBF.Deserialize(ms);
                ms.Close();
                return newLoopObj;
            }


        }
        public class myLoopSorter : IComparer<myLoopDataObject>
        {
            public int Compare(myLoopDataObject _loop1, myLoopDataObject _loop2)
            {
                return (_loop1.myClusterID == _loop2.myClusterID) ?
                    ((string)_loop1.myPdbIDobj).CompareTo((string)_loop2.myPdbIDobj) : _loop1.myClusterID.CompareTo(_loop2.myClusterID);
            }
        }
        public class myClusteredLoopsSorter : IComparer<List<myLoopDataObject>>
        {
            public int Compare(List<myLoopDataObject> _cl1, List<myLoopDataObject> _cl2)
            {
                return (_cl1.Count == _cl2.Count) ? ((string)_cl1[0].myPdbIDobj).CompareTo((string)_cl2[0].myPdbIDobj) : _cl2.Count.CompareTo(_cl1.Count); // ranks larger cluster before smaller
            }
        }

        private class clusterContainer
        {
            public clusterContainer(List<int> _memberList, int[] _centerValueAndItsCenter, List<int> _clusterLabel)
            {
                memberPDBs = _memberList;
                numberOfMembers = _memberList.Count;
                centerPDBValue = _centerValueAndItsCenter[0];
                mergeCandidatesWithCenters.Add(_centerValueAndItsCenter);
                foreach (int oneLabel in _clusterLabel)
                { clusterLabel.Add(oneLabel); }
            }

            public clusterContainer(List<int> _memberList, int[] _centerValueAndItsCenter, int _clusterLabel)
            {
                memberPDBs = _memberList;
                numberOfMembers = _memberList.Count;
                centerPDBValue = _centerValueAndItsCenter[0];
                mergeCandidatesWithCenters.Add(_centerValueAndItsCenter);
                clusterLabel.Add(_clusterLabel);
            }

            public clusterContainer(List<int> _memberList, List<int[]> _mergeCandidates, int _centerValue, List<int> _clusterLabel)
            {
                memberPDBs = _memberList;
                mergeCandidatesWithCenters = _mergeCandidates;
                numberOfMembers = _memberList.Count;
                centerPDBValue = _centerValue;
                foreach (int oneLabel in _clusterLabel)
                { clusterLabel.Add(oneLabel); }
            }

            public clusterContainer(List<int> _memberList, List<int[]> _mergeCandidates, int _centerValue, int _clusterLabel)
            {
                memberPDBs = _memberList;
                mergeCandidatesWithCenters = _mergeCandidates;
                numberOfMembers = _memberList.Count;
                centerPDBValue = _centerValue;
                clusterLabel.Add(_clusterLabel);
            }

            public clusterContainer()
            {
                numberOfMembers = 0;
            }

            public void mergeWithCluster(clusterContainer _otherCluster)
            {
                foreach (int _memberPDB in _otherCluster.memberPDBs)
                {
                    if (!memberPDBs.Contains(_memberPDB))
                    {
                        memberPDBs.Add(_memberPDB);
                    }
                }
                if (_otherCluster.numberOfMembers > numberOfMembers)
                {
                    centerPDBValue = _otherCluster.centerPDBValue;
                }
                foreach (int aLabel in _otherCluster.clusterLabel)
                { clusterLabel.Add(aLabel); }
                numberOfMembers = memberPDBs.Count;
                List<int[]> newMergeList = new List<int[]>();
                foreach (int[] memberCenterPair in mergeCandidatesWithCenters)
                {
                    int centerValue = memberCenterPair[1];
                    //if (!(_otherCluster.memberPDBs.Contains(centerValue) && memberPDBs.Contains(centerValue)))
                    { newMergeList.Add(memberCenterPair); }
                }
                foreach (int[] memberCenterPair in _otherCluster.mergeCandidatesWithCenters)
                {
                    int centerValue = memberCenterPair[1];
                    //if (!(_otherCluster.memberPDBs.Contains(centerValue) && memberPDBs.Contains(centerValue)))
                    {
                        if (!newMergeList.Contains(memberCenterPair))
                        { newMergeList.Add(memberCenterPair); }
                    }
                }
                mergeCandidatesWithCenters = newMergeList;
            }

            public List<int> memberPDBs = new List<int>();
            public List<int[]> mergeCandidatesWithCenters = new List<int[]>(); // List of int[2]s: [0]: merge candidate [1]: its center
            public int numberOfMembers = new int();
            public int centerPDBValue = new int();
            public List<int> clusterLabel = new List<int>();
        }

        private class comparingClusters : IComparer<clusterContainer>
        {
            public int Compare(clusterContainer _cl1, clusterContainer _cl2)
            {
                return _cl1.numberOfMembers - _cl2.numberOfMembers;
            }
        }
        private class sortClusters : IComparer<List<int>>
        {
            public int Compare(List<int> _cluster1, List<int> _cluster2)
            {
                if (_cluster1.Count != _cluster2.Count)
                {
                    return _cluster2.Count - _cluster1.Count; // should rank the longer cluster first
                }
                return _cluster1[0] - _cluster2[0];
            }
        }
    }
}