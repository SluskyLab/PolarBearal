using System;
using System.Collections;

namespace betaBarrelProgram.AtomParser
{
	/// <summary>
	/// category for atoms of each asymmetric chains
	/// </summary>
	public class AtomCategory
	{
		public ArrayList chainAtomsList = new ArrayList ();
        internal double resolution = new double();

		public ChainAtoms[] ChainAtomList
		{
			get
			{
				ChainAtoms [] chainList = new ChainAtoms [chainAtomsList.Count];
				chainAtomsList.CopyTo (chainList);
				return chainList;
			}
			set
			{
				if (value == null) return;
				chainAtomsList.Clear ();
				chainAtomsList.AddRange ((ChainAtoms[])value);
			}
		}

        public double Resolution
        {
            get
            {
                return resolution;
            }
            set
            {
                resolution = value;
            }
        }

		/// <summary>
		/// chains and their Calpha atoms
		/// </summary>
		/// <returns></returns>
		public ChainAtoms[] CalphaAtomList ()
		{
			ChainAtoms[] calphaChainList = new ChainAtoms [chainAtomsList.Count];
			for (int i = 0; i < chainAtomsList.Count; i ++)
			{
				ChainAtoms chainCalphaAtoms = new ChainAtoms ();
				chainCalphaAtoms.CartnAtoms = ((ChainAtoms)chainAtomsList[i]).CalphaAtoms();
				chainCalphaAtoms.AsymChain = ((ChainAtoms)chainAtomsList[i]).AsymChain;
                chainCalphaAtoms.AuthAsymChain = ((ChainAtoms)chainAtomsList[i]).AuthAsymChain;
				chainCalphaAtoms.EntityID = ((ChainAtoms)chainAtomsList[i]).EntityID;
				chainCalphaAtoms.PolymerType = ((ChainAtoms)chainAtomsList[i]).PolymerType;
				calphaChainList[i] = chainCalphaAtoms;
			}
			return calphaChainList;
		}

		/// <summary>
		/// chains and their Cbeta atoms
		/// </summary>
		/// <returns></returns>
		public ChainAtoms[] CbetaAtomList ()
		{
			ChainAtoms[] cbetaChainList = new ChainAtoms [chainAtomsList.Count];
			for (int i = 0; i < chainAtomsList.Count; i ++)
			{
				ChainAtoms chainCbetaAtoms = new ChainAtoms ();
				chainCbetaAtoms.CartnAtoms = ((ChainAtoms)chainAtomsList[i]).CbetaAtoms ();
				chainCbetaAtoms.AsymChain = ((ChainAtoms)chainAtomsList[i]).AsymChain;
                chainCbetaAtoms.AuthAsymChain = ((ChainAtoms)ChainAtomList[i]).AuthAsymChain;
				chainCbetaAtoms.EntityID = ((ChainAtoms)chainAtomsList[i]).EntityID;
				chainCbetaAtoms.PolymerType = ((ChainAtoms)chainAtomsList[i]).PolymerType;
				cbetaChainList[i] = chainCbetaAtoms;
			}
			return cbetaChainList;
		}

		/// <summary>
		/// chains and their Calpha and Cbeta atoms
		/// </summary>
		/// <returns></returns>
		public ChainAtoms[] CalphaCbetaAtomList ()
		{
			ChainAtoms[] calphaCbetaChainList = new ChainAtoms [chainAtomsList.Count];
			for (int i = 0; i < chainAtomsList.Count; i ++)
			{
				ChainAtoms chainCalphaCbetaAtoms = new ChainAtoms ();
				chainCalphaCbetaAtoms.CartnAtoms = ((ChainAtoms)chainAtomsList[i]).CalphaCbetaAtoms ();
				chainCalphaCbetaAtoms.AsymChain = ((ChainAtoms)chainAtomsList[i]).AsymChain;
                chainCalphaCbetaAtoms.AuthAsymChain = ((ChainAtoms)chainAtomsList[i]).AuthAsymChain;
				chainCalphaCbetaAtoms.EntityID = ((ChainAtoms)chainAtomsList[i]).EntityID;
				chainCalphaCbetaAtoms.PolymerType = ((ChainAtoms)chainAtomsList[i]).PolymerType;
				calphaCbetaChainList[i] = chainCalphaCbetaAtoms;
			}
			return calphaCbetaChainList;
		}

		/// <summary>
		/// chains and their backbone atoms
		/// </summary>
		public ChainAtoms[] BackboneAtomList ()
		{
			ChainAtoms[] backboneChainList = new ChainAtoms [chainAtomsList.Count];
			for (int i = 0; i < chainAtomsList.Count; i ++)
			{
				ChainAtoms backboneChainAtoms = new ChainAtoms ();
				backboneChainAtoms.CartnAtoms = ((ChainAtoms)chainAtomsList[i]).BackboneAtoms ();
				backboneChainAtoms.AsymChain = ((ChainAtoms)chainAtomsList[i]).AsymChain;
                backboneChainAtoms.AuthAsymChain = ((ChainAtoms)chainAtomsList[i]).AuthAsymChain;
				backboneChainAtoms.EntityID = ((ChainAtoms)chainAtomsList[i]).EntityID;
				backboneChainAtoms.PolymerType = ((ChainAtoms)chainAtomsList[i]).PolymerType;
				backboneChainList[i] = backboneChainAtoms;
			}
			return backboneChainList;			
		}

		/// <summary>
		/// add asymmetry chain
		/// </summary>
		/// <param name="chainAtoms"></param>
		public void AddChainAtoms(ChainAtoms chainAtoms)
		{
			chainAtomsList.Add (chainAtoms);
		}

        
	}
}

