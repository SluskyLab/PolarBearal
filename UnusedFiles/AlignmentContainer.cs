using System;
using System.Text;
using System.Collections;

namespace betaBarrelProgram
{
    /// <summary>
    /// each entry in the hash table corresponds to a position in the multiple alignment
    /// each key is the name of a given position within that alignment
    /// each value is a hash of protein identifiers, for which each value is an array of strings:
    ///     key: positionID
    ///     Value: AlignmentPosition object
    ///     note: contains iterators, etc.
    /// </summary>
    class AlignmentContainer
    {
        // constructors and destructors
        public AlignmentContainer()
        {
            sourceOfAlignment = "unknown";
            numberOfPositions = 0;
        }

        //public member functions
        public void AddNewPosition(string _position, string _previousPosition, string _nextPosition)
        // note: this function doesn't take into account insertion
        // requires manual resolution of _prevK, _nextK
        {
            AlignmentPosition newPosition = new AlignmentPosition(_previousPosition,_nextPosition);
            alnHash.Add(_position, newPosition);
            numberOfPositions++;
            return;
        }

        public void AddNewPosition(string _position) // overload to add position to end of alignment
        {
            // find last position
            string currentLastPosition = GetLastPosition();
            string endKey = "end";
            AlignmentPosition oldLastPosition = new AlignmentPosition();
            oldLastPosition = (AlignmentPosition)alnHash[currentLastPosition];
            oldLastPosition.SetNextAlignmentPosition(_position);
            AlignmentPosition newPosition = new AlignmentPosition();
            newPosition.SetPreviousAlignmentPosition(currentLastPosition);
            newPosition.SetNextAlignmentPosition(endKey);
            alnHash.Add(_position, newPosition);
            numberOfPositions++;
            return;
        }

        public void AddResidueToPosition(string _position, string _protID, string[] _resInfo)
        {
            //AlignmentPosition tmpRefHolder = new AlignmentPosition(); // making reference to get to object
            //tmpRefHolder = (AlignmentPosition)alnHash[_position]; // assigning reference to object
            //tmpRefHolder.addResidueToAlignmentPosition(_protID, _resInfo); // test this to see if it works // nope, doesn't work
            if(!alnHash.ContainsKey(_position))
            {
                Console.WriteLine("reporting that alnHash does not contain key " + _position);
                return;
            }

            ((AlignmentPosition)alnHash[_position]).AddResidueToAlignmentPosition(_protID, _resInfo);

            if (!proteinsInAlignment.Contains(_protID))
            {
                proteinsInAlignment.Add(_protID);
            }
            return;
        }

        public bool InsertPosition(string _position, string _previousPosition, string _nextPosition)
            // inserts alignment position object in between noted positions
        {
            if (alnHash.ContainsKey(_previousPosition) && alnHash.ContainsKey(_nextPosition))
            {
                if (GetNextPositionKey(_previousPosition) == _nextPosition)
                {
                    AlignmentPosition newPosition = new AlignmentPosition(_previousPosition, _nextPosition);
                    alnHash.Add(_position, newPosition);
                    AlignmentPosition prevP = new AlignmentPosition();
                    AlignmentPosition nextP = new AlignmentPosition();
                    prevP = (AlignmentPosition)alnHash[_previousPosition];
                    nextP = (AlignmentPosition)alnHash[_nextPosition];
                    prevP.SetNextAlignmentPosition(_position);  // need to test this! Does this change underlying obj?
                    nextP.SetPreviousAlignmentPosition(_position); // ditto
                    numberOfPositions++;
                    return true;
                }
                else
                { return false; }
            }
            else
            { return false; }
        }

        public bool DeletePosition(string _position) // deletes position, moves Next and Prev keys and returns true if position exists
        {
            if (alnHash.ContainsKey(_position))
            {
                string previousPosition = GetPreviousPositionKey(_position);
                string nextPosition = GetNextPositionKey(_position);
                if (previousPosition != "start")
                {
                    AlignmentPosition prevP = new AlignmentPosition();
                    prevP = (AlignmentPosition)alnHash[previousPosition];
                    prevP.SetNextAlignmentPosition(nextPosition);
                }

                if (nextPosition != "end")
                {
                    AlignmentPosition nextP = new AlignmentPosition();
                    nextP = (AlignmentPosition)alnHash[nextPosition];
                    nextP.SetPreviousAlignmentPosition(previousPosition);
                }
                Hashtable infoFromDeletedPosition = ((AlignmentPosition)alnHash[_position]).GetResidueInfo(); // need to get position info before it's deleted
                alnHash.Remove(_position);
                // clean up protein list here: if this is the last position containing said protein,
                // remove it from the list

                string[] allKeys = new string[alnHash.Count];
                alnHash.Keys.CopyTo(allKeys, 0);

                string[] listOfProteinsInDeletedPosition = new string[infoFromDeletedPosition.Count];
                infoFromDeletedPosition.Keys.CopyTo(listOfProteinsInDeletedPosition, 0);
                bool[] proteinStillPresent = new bool[listOfProteinsInDeletedPosition.Length];

                for (int deletedPositionProteinIndex = 0; deletedPositionProteinIndex < listOfProteinsInDeletedPosition.Length; deletedPositionProteinIndex++)
                {
                    proteinStillPresent[deletedPositionProteinIndex] = false;
                    for (int alnKeyIndex = 0; alnKeyIndex < allKeys.Length; alnKeyIndex++)
                    {
                        if (((AlignmentPosition)alnHash[allKeys[alnKeyIndex]]).IsProteinPresent(listOfProteinsInDeletedPosition[deletedPositionProteinIndex]))
                        {
                            proteinStillPresent[deletedPositionProteinIndex] = true;
                            break;
                        }
                    }
                }
                for (int deletedPositionProteinIndex = 0; deletedPositionProteinIndex < listOfProteinsInDeletedPosition.Length; deletedPositionProteinIndex++)
                {
                    if (!proteinStillPresent[deletedPositionProteinIndex])
                    {
                        if(proteinsInAlignment.Contains(listOfProteinsInDeletedPosition[deletedPositionProteinIndex]))
                        {
                            proteinsInAlignment.RemoveAt(proteinsInAlignment.IndexOf(listOfProteinsInDeletedPosition[deletedPositionProteinIndex]));
                        }
                    }
                } // this loop gets rid of proteins that are no longer in alignment after deletion of position
                numberOfPositions--;
                return true;
            }
            else
            { return false; }
        }

        // member functions for reading in different types of alignments

        public bool ReadInPOSAalignment(ArrayList _RawLinesFromFile)
            // note: this algorithm assumes that there aren't multiple lines per protein
            // if there are (like ClustalW format), then this input format must be changed
        {
            Hashtable sequencesAndPositions = new Hashtable();

            ArrayList ProteinLabels = new ArrayList();
            ArrayList POSAfragLabels = new ArrayList(); // ArrayList of POSAfragment labels for keys later
            ArrayList POSAfragsByLine = new ArrayList(); // put data in here before sorting and putting in next AL
            ArrayList POSAfragsByProtein = new ArrayList(); // arraylist of ArrayLists of alignment string frags
          
            char commentSignal = '#'; // requires single quotes for type char

            sourceOfAlignment = "POSA";

            Console.Write("direct input from file, by line, separated by lines of asterisks\n\n\n");
            for (int x = 0; x < _RawLinesFromFile.Count; x++)
            {
                Console.WriteLine(x + ": " + _RawLinesFromFile[x]);
                Console.WriteLine("*************************************");
            }
            Console.ReadLine();

            for(int lineIndex = 0; lineIndex < _RawLinesFromFile.Count; lineIndex++)
            {
                string[] words = Array.FindAll<string>(((string)_RawLinesFromFile[lineIndex]).Split(
                    new char[] { ' ', '\t' }), delegate(string s)
                {
                    return !String.IsNullOrEmpty(s);
                }); // (return type here should be string?)


                string alignmentLine = String.Join(" ", words); // string for alignment with " " delim

                if (((string)_RawLinesFromFile[lineIndex])[0] == commentSignal)
                { // is comment line: alignment are labels
                    for(int k = 1; k < words.Length; k++ )
                    { POSAfragLabels.Add(words[k]); } // these will become the labels for the position hash keys
                }
                else
                { // is alignment line: contains protein names and alignment names
                  // this method is messy, but I need to separate out the PDB name (which can have numbers)
                        // from the POSA alignment (where {number} are all treated as delimiters to be skipped)
                    string[] arrayToGetFirstWord = Array.FindAll<string>(alignmentLine.Split(
                        new char[] { ' ', '\t', '{', '}' }),
                        delegate(string t) { return !String.IsNullOrEmpty(t); });
                    string[] lessArray = new string[arrayToGetFirstWord.Length - 1];
                    for (int wordArrayIndex = 1; wordArrayIndex < arrayToGetFirstWord.Length; wordArrayIndex++)
                    {
                        lessArray[wordArrayIndex - 1] = arrayToGetFirstWord[wordArrayIndex];
                    }
                    string holdingString = String.Join("", lessArray);
                    string[] alignLineWords = Array.FindAll<string>(holdingString.Split(
                        new char[] { ' ', '\t', '{', '}', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' }),
                        delegate(string t) { return !String.IsNullOrEmpty(t); });
                    if (alignLineWords.Length >= 1)
                    {
                        ArrayList alignmentFragsOfLine = new ArrayList();
                        alignmentFragsOfLine.Add(arrayToGetFirstWord[0]); // adding protein ID
                        for (int j = 0; j < alignLineWords.Length; j++) //  [0] is not proteinID
                        {
                            alignmentFragsOfLine.Add(alignLineWords[j]);
                        }
                        POSAfragsByLine.Add(alignmentFragsOfLine);
                        if(!ProteinLabels.Contains((string)alignmentFragsOfLine[0]))
                        {
                            ProteinLabels.Add((string)alignmentFragsOfLine[0]);
                        }
                    }
                }
            }
            // now, to sort POSAfragsByLine into appropriate lines
            for (int i = 0; i < ProteinLabels.Count; i++)
            {
                ArrayList fragsOfOneProtein = new ArrayList();
                for (int j = 0; j < POSAfragsByLine.Count; j++)
                {
                    if (((ArrayList)POSAfragsByLine[j])[0] == ProteinLabels[i]) // this can allow for multiple lines per protein
                    {
                        for (int k = 1; k < ((ArrayList)POSAfragsByLine[j]).Count; k++)
                        {
                            fragsOfOneProtein.Add(((ArrayList)POSAfragsByLine[j])[k]);
                        }
                    }
                }
                POSAfragsByProtein.Add(fragsOfOneProtein);
            }

            // REMINDER: what's going on so far
            // at this point: POSAfragsByProtein contain ArrayLists of sequence fragments, in S1, S2, .. POSA order
            // ArrayLists in POSAfragsByProtein are in by protein in order of ArrayList ProteinLabels

            // POSAfragLabels has S1, S2, ... labels
            // Now: to make the AlignmentContainer Object for this alignment

            // in order for this next section to work, all proteins must have the same # of POSA fragments (and same # of res in each fragment)

            string previousPositionString = "start";

            if (((ArrayList)POSAfragsByProtein[0]).Count != POSAfragLabels.Count)
            { 
                Console.WriteLine("Number of protein fragment labels in header not consistent with number of fragments in protein sequences in POSA input, Warning!"); 
                Console.ReadLine();
                return false; // bails out of the subroutine with error
            }


            //POSAfragsByProtein[protein][fragment]

            
            for (int byFragment = 0; byFragment < POSAfragLabels.Count; byFragment++) // one loop per POSA label (S1, S2, ... )
            {
                for (int byPosition = 0; byPosition < ((string)((ArrayList)POSAfragsByProtein[0])[byFragment]).Length; byPosition++) // one loop per POSA position ** problem here
                {
                    string positionString = (string)POSAfragLabels[byFragment];
                    positionString += "_";
                    positionString += (byPosition + 1).ToString();
                    // now, set previous and next position strings
                    string nextPositionString;
                    if (byPosition == ((string)((ArrayList)POSAfragsByProtein[0])[byFragment]).Length - 1)
                    {
                        if (byFragment == POSAfragLabels.Count - 1) // end of the alignment
                        {
                            nextPositionString = "end";
                        }
                        else // end of the fragment, but not on the last fragment
                        {
                            nextPositionString = (string)POSAfragLabels[byFragment + 1];
                            nextPositionString += "_1";
                        }
                    }
                    else // not at end of fragment
                    {
                        nextPositionString = (string)POSAfragLabels[byFragment];
                        nextPositionString += "_";
                        nextPositionString += (byPosition + 2).ToString();
                    }
                    // now that previous and next positions are set, add the position
                    AddNewPosition(positionString, previousPositionString, nextPositionString);
                    previousPositionString = positionString; //ready for the next position
                    for (int byProtein = 0; byProtein < POSAfragsByProtein.Count; byProtein++)
                    {
                        char charOfResidueToAdd = ((string)((ArrayList)POSAfragsByProtein[byProtein])[byFragment])[byPosition];
                        char[] tempCharArray = new char[1];
                        tempCharArray[0] = charOfResidueToAdd;
                        string residueToAdd = new string(tempCharArray);
                        string chainOfResidue = "Z"; // these are dummy values because the POSA
                        string stringOfResiduePosition = "9999"; // alignment has no info about them

                        string[] resInfo = new string[3];
                        resInfo[0] = residueToAdd;
                        resInfo[1] = chainOfResidue;
                        resInfo[2] = stringOfResiduePosition;
                        AddResidueToPosition(positionString, (string)ProteinLabels[byProtein], resInfo);
                    }
                }
            }
            

            return true;
        }

        public bool ReadInClustalWAlignment(ArrayList _RawLinesFromFile)
        {
            ArrayList proteinLabels = new ArrayList();
            ArrayList proteinSequences = new ArrayList(); // strings of sequence by protein
            Hashtable proteinsAndSequences = new Hashtable(); // keys->proteinLabels, values -> sequences (as StringBuilders, so append works)

            sourceOfAlignment = "ClustalW";

            for (int lineIndex = 0; lineIndex < _RawLinesFromFile.Count; lineIndex++)
            {
                if (((string)_RawLinesFromFile[lineIndex]).Length != 0) // skipping garbage lines
                {
                    if (!(((string)_RawLinesFromFile[lineIndex])[0] == ' ' || ((string)_RawLinesFromFile[lineIndex])[0] == '\t' || ((string)_RawLinesFromFile[lineIndex])[0] == '\n')) // all conditions for skipping the line
                    {
                        string[] words = Array.FindAll<string>(((string)_RawLinesFromFile[lineIndex]).Split(
                        new char[] { ' ', '\t' }), delegate(string s)
                        {
                            return !String.IsNullOrEmpty(s);
                        });
                        if (!(words[0] == "CLUSTAL") && words.Length > 1)
                        {
                            //body of loop here
                            if (!(proteinLabels.Contains(words[0])))
                            {
                                proteinLabels.Add(words[0]); // I want to keep it in order of ClustalW file
                                StringBuilder convertSequenceToSB = new StringBuilder(words[1]);
                                proteinsAndSequences.Add(words[0], convertSequenceToSB); // handle error here too? nah
                            }
                            else
                            {
                                if (proteinsAndSequences.ContainsKey(words[0]))
                                {
                                    ((StringBuilder)proteinsAndSequences[words[0]]).Append(words[1]);
                                }
                                else
                                {
                                    Console.WriteLine("Error in ClustalW file input.  Hash keys not consistent with list of protein labels from file");
                                    Console.WriteLine("Aborting CW file read");
                                    Console.ReadLine();
                                    return false;
                                }
                            }
                        }
                    }
                }
            }

            // now, to check the StringBuilder values... if they aren't all the same size, there's a problem

            int alignmentLength = new int();
            string[] keysToCWAlnHash = new string[proteinsAndSequences.Count];
            proteinsAndSequences.Keys.CopyTo(keysToCWAlnHash, 0);
            alignmentLength = ((StringBuilder)proteinsAndSequences[keysToCWAlnHash[0]]).Length;

            for (int i = 1; i < keysToCWAlnHash.Length; i++) // skip [0] ... that's where alignmentLength comes from
            {
                if (alignmentLength != ((StringBuilder)proteinsAndSequences[keysToCWAlnHash[i]]).Length)
                {
                    Console.WriteLine("Alignment sequence length (with gaps) for " + keysToCWAlnHash[i] + " is not equal to reference alignment length for " + keysToCWAlnHash[0]);
                    Console.Write("Protein tags in Hashtable and their associated alignment sequence lengths:\n\n");
                    for (int hashIndex = 0; hashIndex < keysToCWAlnHash.Length; hashIndex++)
                    {
                        Console.WriteLine(keysToCWAlnHash[hashIndex] + "\t" + ((StringBuilder)proteinsAndSequences[keysToCWAlnHash[i]]).Length);
                    }
                    Console.WriteLine("About to bail from ClustalW-reading subroutine");
                    Console.ReadLine();
                    return false;
                }
            }

            string previousPositionLabel = "start";


            for (int positionIndex = 0; positionIndex < ((StringBuilder)proteinsAndSequences[keysToCWAlnHash[0]]).Length; positionIndex++)
            {
                
                string positionLabel = "CW_" + (positionIndex + 1).ToString();
                string nextPositionLabel;
                if(positionIndex != (((StringBuilder)proteinsAndSequences[keysToCWAlnHash[0]]).Length - 1))
                {
                    nextPositionLabel = "CW_" + (positionIndex + 2).ToString();
                }
                else
                {
                    nextPositionLabel = "end";
                }

                AddNewPosition(positionLabel, previousPositionLabel, nextPositionLabel);
                for (int protIndex = 0; protIndex < proteinLabels.Count; protIndex++)
                {
                    char charOfResidueToAdd = ((StringBuilder)proteinsAndSequences[proteinLabels[protIndex]]).ToString()[positionIndex];
                    char[] tempCharArray = new char[1];
                    tempCharArray[0] = charOfResidueToAdd;
                    string[] currentResInfo = new string[3];
                    currentResInfo[1] = "Z";
                    currentResInfo[2] = "9999";
                    string residueToAdd = new string(tempCharArray);
                    currentResInfo[0] = residueToAdd;
                    AddResidueToPosition(positionLabel, (string)proteinLabels[protIndex], currentResInfo);
                }
                previousPositionLabel = positionLabel;
            }
            return true;
        }

        // accessors

        public string GetPreviousPositionKey(string _position)
        {
            AlignmentPosition tmpHolder = new AlignmentPosition();
            tmpHolder = (AlignmentPosition)alnHash[_position];
            return tmpHolder.GetPreviousAlignmentPosition();
        }

        public string GetNextPositionKey(string _position)
        {
            AlignmentPosition tmpHolder = new AlignmentPosition();
            tmpHolder = (AlignmentPosition)this.alnHash[_position];
            return tmpHolder.GetNextAlignmentPosition();
        }


        public string GetFirstPosition()
        {
            string[] allKeys = new string[alnHash.Count];
            alnHash.Keys.CopyTo(allKeys, 0);
            string lookForThis = "start";
            string currentPosition = allKeys[0];
            string previousPosition = "foobar";
            while (true)
            {
                previousPosition = GetPreviousPositionKey(currentPosition);
                if (previousPosition == lookForThis)
                {
                    break;
                }
                else
                {
                    currentPosition = previousPosition;
                }
            }
            return currentPosition;
        }

        public string GetLastPosition()
        {
            string[] allKeys = new string[alnHash.Count];
            alnHash.Keys.CopyTo(allKeys, 0);
            string lookForThis = "end";
            string currentPosition = allKeys[allKeys.Length - 1];
            string nextPosition = "foobar";
            while (true)
            {
                nextPosition = GetNextPositionKey(currentPosition);
                if (nextPosition == lookForThis)
                {
                    break;
                }
                else
                {
                    currentPosition = nextPosition;
                }
            }
            return currentPosition;
        }

        public string[] GetAllKeys() // gets all keys (in their order in the hash, not alignment order
        {
            string[] allKeys = new string[alnHash.Count];
            alnHash.Keys.CopyTo(allKeys, 0);
            return allKeys;
        }

        public string[] GetAllKeysInOrder() // gets all keys in alignment order
        {
            string currentPosition = GetFirstPosition();
            string[] allKeysInOrder = new string[alnHash.Count];
            int keyCounter = 0;
            while (currentPosition != "end")
            {
                allKeysInOrder[keyCounter] = currentPosition;
                currentPosition = GetNextPositionKey(currentPosition);
                keyCounter++;
            }
            return allKeysInOrder;
        }

        public bool IsProteinPresentAtPosition(string _position, string _protID)
        {
            return ((AlignmentPosition)alnHash[_position]).IsProteinPresent(_protID);
        }

        public string GetResidue(string _position, string _protID)
        {
            return ((AlignmentPosition)alnHash[_position]).GetResidueFromProtein(_protID);
        }

        public void PrintAlignmentToConsole()
        {
            string[] allKeysInOrder = GetAllKeysInOrder();
            ArrayList totalAlignment = new ArrayList();
            for (int i = 0; i < proteinsInAlignment.Count; i++) // sets up the proteinsInAlignment array (one entry per protein)
            {
                StringBuilder proteinLine = new StringBuilder();
                totalAlignment.Add(proteinLine);
            }

            for (int keyIndex = 0; keyIndex < allKeysInOrder.Length; keyIndex++)
            {
                for (int protIndex = 0; protIndex < proteinsInAlignment.Count; protIndex++)
                {
                    if (IsProteinPresentAtPosition(allKeysInOrder[keyIndex], (string)proteinsInAlignment[protIndex]))
                    {
                        ((StringBuilder)totalAlignment[protIndex]).Append(GetResidue(allKeysInOrder[keyIndex], (string)proteinsInAlignment[protIndex]));
                    }
                    else
                    {
                        ((StringBuilder)totalAlignment[protIndex]).Append("-");
                    }
                }
            }

            Console.Write("\nAlignment:\n\n");
            for (int lineIndex = 0; lineIndex < proteinsInAlignment.Count; lineIndex++)
            {
                Console.Write((StringBuilder)totalAlignment[lineIndex]);
                Console.Write("\t" + (string)proteinsInAlignment[lineIndex] + "\n");
            }
            Console.WriteLine("Source of alignment = " + sourceOfAlignment);
            return;
        }

        public ArrayList GetListOfProteinsInAlignment()
        {
            return proteinsInAlignment;
        }

        public void ReportEachResidueFromEachPosition() // this is present for debugging
        {
            string[] allKeysInOrder = GetAllKeysInOrder();
            for (int keyIndex = 0; keyIndex < allKeysInOrder.Length; keyIndex++)
            {
                Console.WriteLine("key: " + allKeysInOrder[keyIndex]);
                Hashtable tmpResInfo = new Hashtable();
                tmpResInfo = ((AlignmentPosition)alnHash[allKeysInOrder[keyIndex]]).GetResidueInfo();
                
                string[] resKeys = new string[tmpResInfo.Count];
                tmpResInfo.Keys.CopyTo(resKeys, 0);
                for (int resKeyIndex = 0; resKeyIndex < resKeys.Length; resKeyIndex++)
                {
                    Console.WriteLine(resKeys[resKeyIndex] + "\t" + ((string[])tmpResInfo[resKeys[resKeyIndex]])[0]);
                }
                
            }

            return;
        }

        // member variables
        private Hashtable alnHash = new Hashtable (); // key: positionID, value:alignmentPosition
        private ArrayList proteinsInAlignment = new ArrayList(); // list of proteins included in alignment
        private int numberOfPositions = new int();
        private string sourceOfAlignment;
    }
}
