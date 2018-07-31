using System;
using System.Data;
using System.IO;
using System.Collections;
using System.Xml;

namespace betaBarrelProgram.AtomParser
{ 
	//Summary description for PdbXmlAtomParser.
	public class XmlAtomParser
	{
		public XmlAtomParser()
		{

		}

		//parse one XML file
        public void ParseXmlFile(string thisXmlFile, string atomType)
		{
			int xmlIndex = thisXmlFile.LastIndexOf ("\\");
			// 4 character for Pdb entry ID
			string pdbId = thisXmlFile.Substring (xmlIndex + 1, 4);

			// <PDBx:struct_biol_genCategory>
			// the category to generate the biological units
			// BuGenCategory buGenCat = new BuGenCategory ();

			// <PDBx:atom_siteCategory>
			// the coordinates of atoms, 
			AtomCategory atomCat = new AtomCategory ();

			// <PDBx:cellCategory>
			// the cystal1 record
			// CrystalInfo crystalInfo = new CrystalInfo ();
			
			XmlDocument xmlDoc = new XmlDocument();
			try
			{
				xmlDoc.Load (thisXmlFile);
				// Create an XmlNamespaceManager for resolving namespaces.
				XmlNamespaceManager nsManager = new XmlNamespaceManager(xmlDoc.NameTable);
				nsManager.AddNamespace("PDBx", "http://deposit.pdb.org/pdbML/pdbx.xsd");
				
				Hashtable entityPolyTypeHash = new Hashtable ();

				// if there are protein chains, no, return
				///// parse polymer type of an entity
                /*
				bool hasProtein = false;
				ParseEntityPolyCategory(ref xmlDoc, entityPolyTypeHash, ref nsManager, out hasProtein);
				if (! hasProtein)
				{
					entryCrystal = null;
					return ;
				}	
                */

				///////////////
				// parse atom_sitescategory <PDBx:fract_transf_matrix11>
				//ParseFractTransfMatrix(ref xmlDoc, ref entryCrystal.scaleCat, ref nsManager);
				
				//////////////
				// Parse Cryst1 record
				//ParseCryst1 (ref xmlDoc, ref entryCrystal.crystal1, ref nsManager);

				///////////////
				// Parse PDBx:struct_biol_genCategory
				//ParseBuGenCategory(ref xmlDoc, ref entryCrystal.buGenCat, ref nsManager);

				///////////////
				// Parse atom 
				//ParseAtoms(ref xmlDoc, ref entryCrystal.atomCat, ref nsManager, entityPolyTypeHash, atomType);
                ParseAtoms(ref xmlDoc, ref atomCat, ref nsManager, entityPolyTypeHash, atomType);
				
			}
			catch (Exception ex)
			{
				throw new Exception (string.Format ("Parse {0} Errors: {1}", pdbId, ex.Message));
			}
		}

/// <summary>
/// added this function myself .... so that atomCat can be referenced elsewhere
/// </summary>
/// <param name="thisXmlFile"></param>
/// <param name="atomType"></param>
/// <param name="PassedAtomCat"></param>
        public void ParseXmlFileAndPassStorage(string thisXmlFile, string atomType, ref AtomCategory PassedAtomCat)
		{
			int xmlIndex = thisXmlFile.LastIndexOf ("\\");
			// 4 character for Pdb entry ID
			string pdbId = thisXmlFile.Substring (xmlIndex + 1, 4);
			
			XmlDocument xmlDoc = new XmlDocument();
			try
			{
				xmlDoc.Load (thisXmlFile);
				// Create an XmlNamespaceManager for resolving namespaces.
                XmlNamespaceManager nsManager = new XmlNamespaceManager(xmlDoc.NameTable);

                string xmlNameSpace = xmlDoc.DocumentElement.Attributes["xmlns:PDBx"].InnerText;

                nsManager.AddNamespace("PDBx", xmlNameSpace);

                Hashtable entityPolyTypeHash = new Hashtable();

				// Parse atom 
				//ParseAtoms(ref xmlDoc, ref entryCrystal.atomCat, ref nsManager, entityPolyTypeHash, atomType);
                ParseAtoms(ref xmlDoc, ref PassedAtomCat, ref nsManager, entityPolyTypeHash, atomType);
				
			}
			catch (Exception ex)
			{
				throw new Exception (string.Format ("Parse {0} Errors: {1}", pdbId, ex.Message));
			}
		}


        #region parse atom new format
        /// <summary>

        /// parse the coordinate of C alphas 

        /// </summary>

        /// <param name="xmlDoc"></param>

        /// <param name="calphaInfoHash"></param>

        /// <param name="nsManager"></param>
        //private void ParseAtomsOld(ref XmlDocument xmlDoc, ref AtomCategory atomCat, ref XmlNamespaceManager nsManager,
       //     string retrievingAtomType)

        private void ParseAtoms(ref XmlDocument xmlDoc, ref AtomCategory atomCat, ref XmlNamespaceManager nsManager,
        Hashtable entityPolyTypeHash, string retrievingAtomType)
        {

            XmlNodeList atomNodeList = xmlDoc.DocumentElement.SelectNodes("descendant::PDBx:atom_siteCategory/PDBx:atom_site", nsManager);

            int atomId = 0;
            string asymId = "";
            string preAsymId = "";
            string preAuthAsymId = "";
            string authAsymId = "";
            string preEntityId = "";
            string preResidue = "";
            string entityId = "";
            string residue = "";
            string authResidue = "";
            string seqId = "";
            string authSeqId = "";
            double cartnX = 0.0;
            double cartnY = 0.0;
            double cartnZ = 0.0;
            string atomType = "-";
            string atomName = "-";
            int modelNum = 1;
            int heterResidueNum = 0;

            double Bfac = 0.0; // added this value here (BN)
            double res = -1; // added this value here (BN)
            double occ = 0.0; // added this value here (BN) for occupancy

            // string polymerType = "";

            // find the resolution here-- added (BN)
            XmlNodeList resolutionNodes = xmlDoc.DocumentElement.SelectNodes("descendant::PDBx:refineCategory/PDBx:refine", nsManager);
            
			foreach (XmlNode resNode in resolutionNodes) // this should only run once
            {
                XmlNodeList childNodes = resNode.ChildNodes;
                foreach (XmlNode refineNode in childNodes)
                {
                    if (refineNode.Name.ToLower().IndexOf("pdbx:ls_d_res_high") > -1)
                    {
                        if (refineNode.InnerText.ToString() != "")
                        {
                            string resString = refineNode.InnerText.ToString();
                            res = Convert.ToDouble(resString);
                            atomCat.Resolution = res;
                        }
                        break;
                    }
                }
            }
            string altConfID = ""; // added this value
            // xml tag called "B_iso_or_equiv"

            ChainAtoms chainAtoms = new ChainAtoms();
            ArrayList atomList = new ArrayList();
            bool isAtomNeeded = false;
			//Run through each atom
            foreach (XmlNode atomNode in atomNodeList)
            {
                isAtomNeeded = false;
                atomId = System.Convert.ToInt32(atomNode.Attributes[0].InnerText.ToString());
                XmlNodeList atomInfoNodeList = atomNode.ChildNodes;
                
				//Go through each "child node" and define all properties for a single atom
				foreach (XmlNode atomInfoNode in atomInfoNodeList)
                {
                    if (atomInfoNode.Name.ToLower() == "pdbx:type_symbol")
                    {
                        atomType = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_atom_id")
                    {
                        atomName = atomInfoNode.InnerText;
						//I think the if/elif sections will just pull out residues depending on the input of retrieving atom type, which is currently set to null
                        //This means that it always returns true
                        if (retrievingAtomType == "CA" || retrievingAtomType == "CB")
                        {
                            if (atomInfoNode.InnerText.ToUpper() != retrievingAtomType)
                            {
                                isAtomNeeded = false;
                                break;
                            }
                            else
                            {
                                isAtomNeeded = true;
                                continue;
                            }
                        }

                        else if (retrievingAtomType == "CA_CB") 
                        {
                            if (atomInfoNode.InnerText.ToUpper() != "CA" &&
                            atomInfoNode.InnerText.ToUpper() != "CB")
                            {
                                isAtomNeeded = false;
                                break;
                            }
                            else
                            {
                                isAtomNeeded = true;
                                continue;
                            }
                        }
                        else
                        {
                            isAtomNeeded = true;
                        }
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:pdbx_pdb_model_num")
                    {
                        modelNum = Convert.ToInt16(atomInfoNode.InnerText);
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_comp_id")
                    {
                        residue = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:auth_comp_id")
                    {
                        authResidue = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_asym_id")
                    {
                        asymId = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower().IndexOf("pdbx:auth_asym_id") > -1)
                    {
                        authAsymId = atomInfoNode.InnerText.ToString();
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_entity_id")
                    {
                        entityId = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_seq_id")
                    {
                        seqId = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:auth_seq_id")
                    {
                        authSeqId = atomInfoNode.InnerText;
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:cartn_x")
                    {
                        cartnX = System.Convert.ToDouble(atomInfoNode.InnerText);
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:cartn_y")
                    {
                        cartnY = System.Convert.ToDouble(atomInfoNode.InnerText);
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:cartn_z")
                    {
                        cartnZ = System.Convert.ToDouble(atomInfoNode.InnerText);
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:b_iso_or_equiv")
                    {
                        Bfac = System.Convert.ToDouble(atomInfoNode.InnerText.ToString());
                        continue; //added: check
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:label_alt_id")
                    {
                        altConfID = atomInfoNode.InnerText.ToString();
                        continue;
                    }
                    if (atomInfoNode.Name.ToLower() == "pdbx:occupancy")
                    {
                        occ = System.Convert.ToDouble(atomInfoNode.InnerText.ToString());
                        continue;
                    }
                }
                
				//If this has crossed into the next chain, write previous residues to chain and start chain list over
				if (preAsymId != asymId && preAsymId != "" && atomList.Count > 0)
                {
                    chainAtoms.AsymChain = preAsymId;
                    chainAtoms.AuthAsymChain = preAuthAsymId;
                    chainAtoms.EntityID = preEntityId; // problem with int to string in new version
                    if (entityPolyTypeHash.ContainsKey(preEntityId))
                    {
                        chainAtoms.PolymerType = entityPolyTypeHash[preEntityId].ToString();
                    }
                    else
                    {
                        chainAtoms.PolymerType = "-";
                    }
                    AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                    atomList.CopyTo(atomArray);
                    chainAtoms.CartnAtoms = atomArray;
                    atomCat.AddChainAtoms(chainAtoms);
                    atomList = new ArrayList();
                    chainAtoms = new ChainAtoms();
                    heterResidueNum = 0;
                    preResidue = "";
                }

                if (modelNum > 1) // only pick up the model with model number 1
                {
                    break;
                }
                if (isAtomNeeded && residue.ToUpper() != "HOH")
                {
                    if (seqId == "")
                    {
                        if (preResidue != residue)
                        {
                            heterResidueNum++;
                        }
                        seqId = heterResidueNum.ToString();
                    }
                    AtomInfo atomInfo = new AtomInfo();
                    atomInfo.atomId = atomId;
                    atomInfo.atomType = atomType;
                    atomInfo.atomName = atomName;
                    atomInfo.seqId = seqId;
                    atomInfo.authSeqId = authSeqId;
                    atomInfo.residue = residue;
                    atomInfo.authResidue = authResidue;
                    atomInfo.xyz.X = cartnX;
                    atomInfo.xyz.Y = cartnY;
                    atomInfo.xyz.Z = cartnZ;
                    atomInfo.bFac = Bfac;
                    atomInfo.altConfID = altConfID;
                    atomInfo.occupancy = occ;
                    /* if (entityPolyTypeHash.ContainsKey(entityId))
                    {
                    polymerType = entityPolyTypeHash[entityId].ToString();
                    }
                    else
                    {
                    polymerType = "-";
                    }
                    atomCat.AddAtom(entityId, asymId, polymerType, atomInfo);*/
                    atomList.Add(atomInfo);
                }
                preAsymId = asymId;
                preAuthAsymId = authAsymId;
                preEntityId = entityId;
                preResidue = residue;
            }
            // add the last one
            if (atomList.Count > 0)
            {
                chainAtoms.AsymChain = asymId;
                chainAtoms.AuthAsymChain = authAsymId;
                chainAtoms.EntityID = entityId;
                if (entityPolyTypeHash.ContainsKey(entityId))
                {
                    chainAtoms.PolymerType = entityPolyTypeHash[entityId].ToString();
                }
                else
                {
                    chainAtoms.PolymerType = "-";
                }
                AtomInfo[] atomArray = new AtomInfo[atomList.Count];
                atomList.CopyTo(atomArray);
                chainAtoms.CartnAtoms = atomArray;
                atomCat.AddChainAtoms(chainAtoms);
            }
        }

        #endregion
	}
}
