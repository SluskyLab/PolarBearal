using System;
using System.Collections;

namespace betaBarrelProgram.AtomParser
{
	/// <summary>
	/// A protein chain
	/// </summary>
	public class ChainAtoms
	{
        internal string authAsymChain = "";
		internal string asymChain = "";
		internal string entityId = "";
		internal string polymerType = "";
		// store the list of atoms for this asymmetric chain
		public ArrayList cartnAtoms = null;
		public ChainAtoms()
		{
			cartnAtoms = new ArrayList ();
		}

		public ChainAtoms(ChainAtoms extChain)
		{
            this.authAsymChain = extChain.AuthAsymChain;
			this.asymChain = extChain.AsymChain;
			this.entityId = extChain.entityId;
			this.polymerType = extChain.PolymerType;
			foreach (AtomInfo atom in extChain.CartnAtoms)
			{
				this.cartnAtoms.Add (atom.Clone ());
			}
		}

		#region properties
		// asymmetric chain id
		public string AsymChain
		{
			get
			{
				return asymChain;
			}
			set
			{
				asymChain = value;
			}
		}

		// entity id
		public string EntityID
		{
			get 
			{
				return entityId;
			}
			set
			{
				entityId = value;
			}
		}
		public string PolymerType
		{
			get
			{
				return polymerType;
			}
			set
			{
				polymerType = value;
			}
		}
        public string AuthAsymChain
        {
            get
            {
                return authAsymChain;
            }
            set
            {
                authAsymChain = value;
            }
        }

		// the list of atoms for this asymmetric chain
		public AtomInfo [] CartnAtoms
		{
			get
			{
				AtomInfo[] atomInfoList = new AtomInfo[cartnAtoms.Count];
				cartnAtoms.CopyTo (atomInfoList);
				return atomInfoList;
			}
			set
			{
				if (value == null)
				{
					return;
				}
				AtomInfo[] atomInfoList = (AtomInfo[]) value;
				cartnAtoms.Clear ();
				cartnAtoms.AddRange (atomInfoList);
			}
		}
		#endregion

		#region atom types
		/// <summary>
		/// c-alpha atoms
		/// </summary>
		public AtomInfo[] CalphaAtoms ()
		{
			ArrayList calphaList = new ArrayList ();
			foreach (AtomInfo atom in cartnAtoms)
			{
				if (atom.atomName.ToUpper () == "CA")
				{
					calphaList.Add (new AtomInfo(atom));
				}
			}
			AtomInfo[] calphaArray = new AtomInfo [calphaList.Count];
			calphaList.CopyTo (calphaArray);
			return calphaArray;
		}

		/// <summary>
		/// c-beta atoms
		/// </summary>
		public AtomInfo[] CbetaAtoms ()
		{
			ArrayList cbetaList = new ArrayList ();
			foreach (AtomInfo atom in cartnAtoms)
			{
				if (atom.atomName.ToUpper () == "CB")
				{
					cbetaList.Add (new AtomInfo(atom));
				}
			}
			AtomInfo[] cbetaArray = new AtomInfo [cbetaList.Count];
			cbetaList.CopyTo (cbetaArray);
			return cbetaArray;
		}

		/// <summary>
		/// c-beta and C-alpha atoms
		/// </summary>
		public AtomInfo[] CalphaCbetaAtoms ()
		{
			ArrayList calphaCbetaList = new ArrayList ();
			foreach (AtomInfo atom in cartnAtoms)
			{
				if (atom.atomName.ToUpper () == "CB" || atom.atomName.ToUpper () == "CA")
				{
					calphaCbetaList.Add (new AtomInfo(atom));
				}
			}
			// for some PDB entries, there is only Calpha
			/*	if (calphaCbetaList.Count == 0)
				{
					foreach (AtomInfo atom in cartnAtoms)
					{
						if (atom.atomName.ToUpper () == "CA")
						{
							calphaCbetaList.Add (new AtomInfo(atom));
						}
					}
				}*/
			AtomInfo[] calphaCbetaArray = new AtomInfo [calphaCbetaList.Count];
			calphaCbetaList.CopyTo (calphaCbetaArray);
			return calphaCbetaArray;
		}


		/// <summary>
		/// get backbone atoms
		/// </summary>
		public AtomInfo[] BackboneAtoms ()
		{			
			ArrayList backboneList = new ArrayList ();
			foreach (AtomInfo atom in cartnAtoms)
			{
				if (atom.atomName == "CA" || atom.atomName == "N" ||
					atom.atomName == "C" || atom.atomName == "O" || atom.atomName == "OXT")
				{
					backboneList.Add (new AtomInfo(atom));
				}
			}

			AtomInfo[] backboneArray = new AtomInfo [backboneList.Count];
			backboneList.CopyTo (backboneArray);
			return backboneArray;			
		}
		#endregion

		/// <summary>
		/// add an atom to the list
		/// </summary>
		/// <param name="atom"></param>
		/// <returns></returns>
		public int AddAtom(AtomInfo atom)
		{
			return cartnAtoms.Add (atom);
		}
		
		/// <summary>
		/// number of residues based on number of C-alpha
		/// </summary>
		/// <returns></returns>
		public int GetNumOfResidues ()
		{
			return this.CalphaAtoms ().Length;
		}


		#region sort by atom ID
		/// <summary>
		/// sort chains by atom ID
		/// </summary>
		public void SortByAtomId ()
		{
			if (cartnAtoms == null)
			{
				return ;
			}
			AtomInfo[] chain = this.CartnAtoms;
			Hashtable seqHash = new Hashtable ();
			foreach (AtomInfo atom in chain)
			{
				// sould have unique atom id
				seqHash.Add (atom.atomId, atom);
			}
			ArrayList atomIdList = new ArrayList (seqHash.Keys);
			atomIdList.Sort ();
			AtomInfo[] sortedChain = new AtomInfo [chain.Length];
			int count = 0;
			foreach (object atomId in atomIdList)
			{
				sortedChain[count] = (AtomInfo)seqHash[atomId];
				count ++;
			}
			this.CartnAtoms = sortedChain;
		}
		#endregion
	}
}
