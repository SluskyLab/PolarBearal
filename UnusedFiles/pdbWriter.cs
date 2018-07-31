using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.IO;
using System.Xml;
using System.Data;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;

namespace betaBarrelProgram.AtomParser
{
    class pdbWriter
    {
        public string WriteFromAtomCategory(string _pdbId, string _chainId, AtomCategory _atomCategoryToPrint, string _pathToOutput)
        {
            ChainAtoms[] myChainAtomsArray = _atomCategoryToPrint.ChainAtomList;
            
            string emptyString = "dakj";
            return emptyString;
        }
        #region write atoms to a file
        /// <summary>
        /// write list of atoms into a file
        /// </summary>
        /// <param name="pdbId">pdb entry</param>
        /// <param name="chain">a list of atoms</param>
        /// <param name="fileNum">an integer number for file sequence id</param>
        /// <param name="remark">description of the file, can be empty</param>
        public string WriteAtoms(string _pdbId, string _chainId, AtomInfo[] _chain, string _filePath)
        {
            string fileName = _pdbId;
            fileName += ".pdb";
            if (_chainId == "")
            {
                _chainId = "A";
            }
            FileStream fileStream = new FileStream(Path.Combine(_filePath, fileName), FileMode.Create, FileAccess.Write);
            StreamWriter fileWriter = new StreamWriter(fileStream);
            string header = "HEADER    " + _pdbId + "                    " + DateTime.Now;
            fileWriter.WriteLine(header);

            try
            {
                string line = "";
                int atomCount = 1;
                foreach (AtomInfo atom in _chain)
                {
                    line = "ATOM  ";
                    string atomIdStr = atomCount.ToString();
                    line += atomIdStr.PadLeft(5, ' ');
                    line += " ";
                    string atomName = atom.atomName;
                    if (atomName != "" && atom.atomType != "H" && atomName.Length < 4)
                    {
                        atomName = " " + atomName;
                    }
                    line += atomName.PadRight(4, ' ');
                    if (atom.altConfID == "")
                    { line += " "; }
                    else
                    { line += atom.altConfID; }
                    line += atom.residue;
                    line += " ";
                    line += _chainId;
                    line += atom.seqId.PadLeft(4, ' ');
                    line += "    ";
                    line += FormatDoubleString(atom.xyz.X, 4, 3);
                    line += FormatDoubleString(atom.xyz.Y, 4, 3);
                    line += FormatDoubleString(atom.xyz.Z, 4, 3);
                    // next line was: line += "  1.00"; (for dummy occupancy)
                    // replace with:
                    line += atom.occupancy.ToString("##0.00").PadLeft(6);
                    // end of occupancy replace block

                    // next line was: line += "  0.00"; (for dummy bfactor)
                    // replace with:
                    //string bFacString = _inputAtom.bFac.ToString("##0.00");
                    line += atom.bFac.ToString("##0.00").PadLeft(6);
                    // end of bfactor replace block
                    line += "    ";
                    line += atom.atomType;
                    fileWriter.WriteLine(line);
                    atomCount++;
                }
                fileWriter.WriteLine("END");
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
            return Path.Combine(_filePath, fileName);
        }
        #endregion

        /// <summary>
        /// WriteSingleLine from old line and new coordinates
        /// </summary>
        /// <param name="_oldLine">(string)original pdb line</param>
        /// <param name="_coords">(List-double)new coordinates</param>
        /// <returns></returns>
        public string WriteSingleLine(string _oldLine, List<double> _coords)
        {
            string line = "";
            line += _oldLine.Substring(0, 30);
            line += FormatDoubleString(_coords[0], 4, 3);
            line += FormatDoubleString(_coords[1], 4, 3);
            line += FormatDoubleString(_coords[2], 4, 3);
            line += _oldLine.Substring(54);
            return line;
        }

        public string WriteSingleLine(AtomInfo _inputAtom, int _atomCounter, string _chainId)
        {
            string line = "";

            line = "ATOM  ";
            string atomIdStr = _atomCounter.ToString();
            line += atomIdStr.PadLeft(5, ' ');
            line += " ";
            string atomName = _inputAtom.atomName;
            if (atomName != "" && _inputAtom.atomType != "H" && atomName.Length < 4)
            {
                atomName = " " + atomName;
            }
            line += atomName.PadRight(4, ' ');
            if (_inputAtom.altConfID == "")
            { line += " "; }
            else
            { line += _inputAtom.altConfID; }
            line += _inputAtom.residue;
            line += " ";
            line += _chainId;
            line += _inputAtom.seqId.PadLeft(4, ' ');
            line += "    ";
            line += FormatDoubleString(_inputAtom.xyz.X, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Y, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Z, 4, 3);
            // next line was: line += "  1.00"; (for dummy occupancy)
            // replace with:
            line += _inputAtom.occupancy.ToString("##0.00").PadLeft(6);
            // end of occupancy replace block

            // next line was: line += "  0.00"; (for dummy bfactor)
            // replace with:
            //string bFacString = _inputAtom.bFac.ToString("##0.00");
            line += _inputAtom.bFac.ToString("##0.00").PadLeft(6);
            // end of bfactor replace block
            line += "    ";
            line += _inputAtom.atomType;
            return line;

        }
        /// <summary>
        /// This is the version that renumbers residues, insertion codes, and chain assignments
        /// Feed it info, outputs a pdb line for one atom
        /// </summary>
        /// <param name="_inputAtom"></param>
        /// <param name="_atomCounter"></param>
        /// <param name="_chainID"></param>
        /// <param name="_newResi"></param>
        /// <param name="_newChain"></param>
        /// <returns></returns>
        public string WriteSingleLine(AtomInfo _inputAtom, int _atomCounter, string _newResi, string _newChain)
        {
            string line = "";

            line = "ATOM  ";
            string atomIdStr = _atomCounter.ToString();
            line += atomIdStr.PadLeft(5, ' ');
            line += " ";
            string atomName = _inputAtom.atomName;
            if (atomName != "" && _inputAtom.atomType != "H" && atomName.Length < 4)
            {
                atomName = " " + atomName;
            }
            line += atomName.PadRight(4, ' ');
            // next line was: line += " ";
            // replace by:
            if (_inputAtom.altConfID == "")
            { line += " "; }
            else
            { line += _inputAtom.altConfID; }
            // end of replace block
            line += _inputAtom.residue;
            line += " ";
            line += _newChain;
            //line += _inputAtom.seqId.PadLeft(4, ' ');
            line += _newResi.PadLeft(4, ' '); // replaces this line: line += _inputAtom.seqId.PadLeft(4, ' ');
            line += "    ";
            line += FormatDoubleString(_inputAtom.xyz.X, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Y, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Z, 4, 3);
            // next line was: line += "  1.00"; (for dummy occupancy)
            // replace with:
            line += _inputAtom.occupancy.ToString("##0.00").PadLeft(6);
            // end of occupancy replace block

            // next line was: line += "  0.00"; (for dummy bfactor)
            // replace with:
            //string bFacString = _inputAtom.bFac.ToString("##0.00");
            line += _inputAtom.bFac.ToString("##0.00").PadLeft(6);
            // end of bfactor replace block
            line += "    ";
            line += _inputAtom.atomType;
            return line;
        }

        /// <summary>
        /// Writes a pdbline given an atom as input (for use for those residues with insertion codes)
        /// </summary>
        /// <param name="_inputAtom"></param>
        /// <param name="_atomCounter"></param>
        /// <param name="_newResi"></param>
        /// <param name="_newChain"></param>
        /// <param name="_insertionCode"></param>
        /// <returns></returns>
        public string WriteSingleLine(AtomInfo _inputAtom, int _atomCounter, string _newResi, string _newChain, string _insertionCode)
        {
            string line = "";

            line = "ATOM  ";
            string atomIdStr = _atomCounter.ToString();
            line += atomIdStr.PadLeft(5, ' ');
            line += " ";
            string atomName = _inputAtom.atomName;
            if (atomName != "" && _inputAtom.atomType != "H" && atomName.Length < 4)
            {
                atomName = " " + atomName;
            }
            line += atomName.PadRight(4, ' ');
            // next line was: line += " ";
            // replace by:
            if (_inputAtom.altConfID == "")
            { line += " "; }
            else
            { line += _inputAtom.altConfID; }
            // end of replace block
            line += _inputAtom.residue;
            line += " ";
            line += _newChain;
            //line += _inputAtom.seqId.PadLeft(4, ' ');
            line += _newResi.PadLeft(4, ' '); // replaces this line: line += _inputAtom.seqId.PadLeft(4, ' ');
            line += _insertionCode;
            line += "   ";
            line += FormatDoubleString(_inputAtom.xyz.X, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Y, 4, 3);
            line += FormatDoubleString(_inputAtom.xyz.Z, 4, 3);
            // next line was: line += "  1.00"; (for dummy occupancy)
            // replace with:
            line += _inputAtom.occupancy.ToString("##0.00").PadLeft(6);
            // end of occupancy replace block

            // next line was: line += "  0.00"; (for dummy bfactor)
            // replace with:
            //string bFacString = _inputAtom.bFac.ToString("##0.00");
            line += _inputAtom.bFac.ToString("##0.00").PadLeft(6);
            // end of bfactor replace block
            line += "    "; // old version from QiFang
            //line += "           "; // version consistent with pdb
            line += _inputAtom.atomType;
            return line;
        }

        public void WriteWholeProtein(ref AtomParser.AtomCategory _atomCat, string _fileWithPath)
        {
            StreamWriter protWriter = new StreamWriter(_fileWithPath, false);
            foreach (AtomParser.ChainAtoms thisChainAtoms in _atomCat.ChainAtomList)
            {
                for (int atomIndex = 0; atomIndex < thisChainAtoms.CartnAtoms.Length; atomIndex++)
                {
                    protWriter.WriteLine(WriteSingleLine(thisChainAtoms.CartnAtoms[atomIndex],
                        (atomIndex + 1), thisChainAtoms.AsymChain));
                }
            }
            protWriter.Close();
            return;
        }

        public void WriteProteinFragment(ref AtomParser.AtomCategory _atomCat, int _fragStart,
            int _fragEnd, int _chainIndex, string _chainID, string _fileWithPath)
        {
            // this is complicated because the indices refer to alpha-carbon-only, while this function writes
            // whole thing
            AtomParser.ChainAtoms theChainAtoms = new ChainAtoms();
            AtomParser.ChainAtoms onlyCAChainAtoms = new ChainAtoms();
            theChainAtoms = _atomCat.ChainAtomList[_chainIndex];
            onlyCAChainAtoms = _atomCat.CalphaAtomList()[_chainIndex];
            string firstCAauthSeqID = onlyCAChainAtoms.CartnAtoms[_fragStart].authSeqId;
            string lastCAauthSeqID = onlyCAChainAtoms.CartnAtoms[_fragEnd].authSeqId;
            StreamWriter protWriter = new StreamWriter(_fileWithPath, false);
            int atomCounter = new int();
            atomCounter = 1;
            bool writingInFragment = new bool();
            bool hitLastCA = new bool();
            writingInFragment = false;
            hitLastCA = false;
            for (int atomIndex = 0; atomIndex < theChainAtoms.CartnAtoms.Length; atomIndex++)
            {
                if (theChainAtoms.CartnAtoms[atomIndex].authSeqId == firstCAauthSeqID)
                {
                    writingInFragment = true;
                }
                if (writingInFragment)
                {
                    if (!hitLastCA)
                    {
                        if (theChainAtoms.CartnAtoms[atomIndex].authSeqId == lastCAauthSeqID)
                        {
                            hitLastCA = true;
                        }
                    }
                    else
                    {
                        if (theChainAtoms.CartnAtoms[atomIndex].authSeqId != lastCAauthSeqID) // hit the residue after the last residue in fragment
                        {
                            break;
                        }
                    }
                    protWriter.WriteLine(WriteSingleLine(theChainAtoms.CartnAtoms[atomIndex],
                            atomCounter, _chainID));
                }
            }
            protWriter.Close();
            return;
        }

        public void WriteCAOnlyFragment(ref AtomParser.AtomCategory _atomCat, int _fragStart,
            int _fragEnd, int _chainIndex, string _chainID, string _fileWithPath)
        {
            AtomParser.ChainAtoms theChainAtoms = new ChainAtoms();
            theChainAtoms = (_atomCat.CalphaAtomList())[_chainIndex];
            StreamWriter protWriter = new StreamWriter(_fileWithPath, false);
            int atomCounter = new int();
            atomCounter = 1;
            for (int atomIndex = _fragStart; atomIndex < _fragEnd; atomIndex++)
            {
                protWriter.WriteLine(WriteSingleLine(theChainAtoms.CartnAtoms[atomIndex],
                    atomCounter, _chainID));
            }
            protWriter.Close();
            return;
        }


        #region format strings
        /// <summary>
        /// format a double into a string 
        /// (8.3) (1234.123)
        /// </summary>
        /// <param name="val"></param>
        /// <returns></returns>
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
        #endregion
    }
}
