using System.Collections;

namespace betaBarrelProgram
{
    /// <summary>
    /// a container object for each position in the alignment
    /// consists of:
    /// * a Hashtable of (key) ProteinID and (Value) string[]:
    ///     [0]: resType
    ///     [1]: chain ID
    ///     [2]: position# (saved as string)
    /// * string of key to previous position in alignment
    /// * string of key to next position in alignment
    /// </summary>
    class AlignmentPosition
    {
        //constructors
        public AlignmentPosition()
        {
            previousKey = "null";
            nextKey = "null";
        }

        public AlignmentPosition(string _previousKey, string _nextKey)
        {
            previousKey = _previousKey;
            nextKey = _nextKey;
        }

        public AlignmentPosition(string _previousKey, string _nextKey, string _protID, string[] _resInfo)
        {
            previousKey = _previousKey;
            nextKey = _nextKey;
            positionHash.Add(_protID, _resInfo);
        }

        // public member methods
        public void SetPreviousAlignmentPosition(string _previousPositionKey)
        {
            previousKey = _previousPositionKey;
            return;
        }

        public void SetNextAlignmentPosition(string _nextPositionKey)
        {
            nextKey = _nextPositionKey;
            return;
        }

        public void AddResidueToAlignmentPosition(string _protID, string[] _resInfo)
        {
            positionHash.Add(_protID, _resInfo);
            return;
        }

        // accessors
        public string GetPreviousAlignmentPosition()
        {
            return previousKey;
        }

        public string GetNextAlignmentPosition()
        {
            return nextKey;
        }

        public Hashtable GetResidueInfo()
        {
            return positionHash;
        }

        public string GetResidueFromProtein(string _protID)
        {
            string[] resInfo = new string[3];
            resInfo = (string[])positionHash[_protID];
            return resInfo[0]; // residue position
        }

        public bool IsProteinPresent(string _protID)
        {
            return positionHash.ContainsKey(_protID);
        }

        // member variables
        private Hashtable positionHash = new Hashtable ();
        private string previousKey;
        private string nextKey;
    }
}
