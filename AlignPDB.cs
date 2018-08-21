/*
**  File: AlignPDB.cs
**  Started: 7/31/2017
**  Contributors: Ryan Feehan
**  Overview:  Make GIF for 4rjw to 5dl6 16 strand to 18 strand transformation
**
**  About:  
**  Last Edited: 
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
    public class AlignPDB
    {
        public static string ALIGN_OUTPUT = @"\\psf\Home\Desktop\SluskyLab\strandChange_16_18\";

        /*
        ** Pre:     
        ** Post:      
        ** About:    
        */
        public AlignPDB()
        {

            //get 16 strand beta barrel
            first_protein = null;
            first_barrel = null;
            SharedFunctions.runBetaBarrel_RYAN("4rjw", ref first_protein, ref first_barrel);

            //get 18 strand beta barrel
            last_protein = null;
            last_barrel = null;
            SharedFunctions.runBetaBarrel_RYAN("5dl6", ref last_protein, ref last_barrel);

            Align();


            WritePDB("first", first_protein);
        }

        /* align proteins on first strand of the transformation */
        public void Align()
        {
            WritePDB("initial_first", first_protein);
            WritePDB("initial_last", last_protein);
            //TODO: ALIGN FIRST AND LAST PDBS
            int start_first = 0;
            int start_last = 0;
            double change_x = 0;
            double change_y = 0;
            double change_z = 0;
            /*  PDB 4rjw Strand 1 residues 199-210 
            **  PDB 5dl6 Strand 1 residues 189-200 ?
            */
            //find start residues
            foreach(Res res in first_protein.Chains[0].Residues)
            {
                if (res.SeqID == 199) break;
                start_first++;
            }
            foreach (Res res in last_protein.Chains[0].Residues)
            {
                if (res.SeqID == 189) break;
                start_last++;
            }
            //find average change in backbone coords
            for (int i = start_first, j = start_last; i <= start_first + 11; i++, j++)
            {
                change_x += first_protein.Chains[0].Residues[i].BackboneCoords["CA"].X - last_protein.Chains[0].Residues[i].BackboneCoords["CA"].X;
                change_y += first_protein.Chains[0].Residues[i].BackboneCoords["CA"].Y - last_protein.Chains[0].Residues[i].BackboneCoords["CA"].Y;
                change_z += first_protein.Chains[0].Residues[i].BackboneCoords["CA"].Z - last_protein.Chains[0].Residues[i].BackboneCoords["CA"].Z;
            }
            change_x = change_x / 11;
            change_y = change_y / 11;
            change_z = change_z / 11;

            Vector3D change = new Vector3D(change_x, change_y, change_z);
            //adjust last protein
            for(int i = 0; i < last_protein.Chains[0].Residues.Count; i++)
            {
                for (int j = 0; j < last_protein.Chains[0].Residues[i].Atoms.Count; j++)
                {
                    last_protein.Chains[0].Residues[i].Atoms[j].Coords -= change;
                }
            }

            WritePDB("aligned_last", last_protein);
        }

        public void CalculateChange()
        {
            //TODO: Calculate Res Change from first to last pdb
        }
        public void FramePDBs()
        {
            //TODO: Make 8 Frame PDBs
        }

        /*
        ** About:    
        */
        public string WritePDB(string _pdbId, Protein _protein)
        {
            string fileName = _pdbId;
            fileName += ".pdb";
           
            FileStream fileStream = new FileStream(Path.Combine(ALIGN_OUTPUT, fileName), FileMode.Create, FileAccess.Write);
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
            return Path.Combine(ALIGN_OUTPUT, fileName);
        }

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

        Protein first_protein;
        Barrel first_barrel;
        Protein last_protein;
        Barrel last_barrel;

    }

}