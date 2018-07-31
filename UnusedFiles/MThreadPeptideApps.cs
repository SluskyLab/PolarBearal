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
    class MThreadPeptideApps
    {

        public void ClusterByAffinityWithMultiThreading(ref int _numberOfClusters, ref List<PeptideAnalysis.myLoopDataObject> _loopObjsByPDB, int _maxIterationCount, ref ArrayList _centerOfEachCluster, ref ArrayList _centerForEachPoint, ref ArrayList _LinesToWrite, double _selfSimFactor, string _progFileName, ref PeptideAnalysis _myModule)
        {
            bool reportOnEachIteration = new bool();
            reportOnEachIteration = true;
            bool useSimulatedAnnealing = new bool();
            useSimulatedAnnealing = false;
            Console.WriteLine("Started clustering at {0}", DateTime.Now);
            int[] bestClusterData = new int[_loopObjsByPDB.Count];
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

            double[,] distanceMatrix = new double[_loopObjsByPDB.Count, _loopObjsByPDB.Count];

            distanceMatrix = _myModule.GenerateDistanceMatrix(ref _loopObjsByPDB);

            //generateMultiThreadInfoHere
            int numberOfRowsPerThread = _loopObjsByPDB.Count / totalNumberOfThreads; // forces answer to be rounded-off int

            //end multithreadInfo


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
                    localSimilarities[simIndex, simIndex2] = -1 * distanceMatrix[simIndex, simIndex2];
                    if (simIndex != simIndex2)
                    {
                        selfSimilarityConstant += localSimilarities[simIndex, simIndex2];
                        countForAverage++;
                    }
                }
            }
            selfSimilarityConstant = _selfSimFactor * selfSimilarityConstant / countForAverage;
            for (int simIndex = 0; simIndex < _loopObjsByPDB.Count; simIndex++)
            {
                localSimilarities[simIndex, simIndex] = selfSimilarityConstant;
            }

            // more setup: initialize availabilities, set update parameter lambda, declare responsibility matrix

            double updateLambda = new double();
            updateLambda = 0.5;
            for (int avIndex = 0; avIndex < _loopObjsByPDB.Count; avIndex++)
            {
                for (int avIndex2 = 0; avIndex2 < _loopObjsByPDB.Count; avIndex2++)
                {
                    localAvailabilities[avIndex, avIndex2] = 0;
                    localResponsibilities[avIndex, avIndex2] = 0;
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
                previousTotalNumberOfClusters = totalNumberOfClusters;
                if (iterationCounter % 20 == 0)
                { reportOnIteration = true; }
                else
                { reportOnIteration = false; }
                if (reportOnEachIteration)
                { reportOnIteration = true; }

                // relaxing parameter
                // doesn't do anything with clustering of AB conformations, because they cluster before 250 anyways
                // (and there the max iteration count is 500)
                if ((useSimulatedAnnealing) && ((iterationCounter > _maxIterationCount / 2) || haveSeenIdenticalCenterAssignments))
                {
                    updateLambda = updateLambda * (double)0.95;
                }


                if (reportOnIteration)
                { Console.WriteLine("Started algorithm loop {0} at {1}", iterationCounter, DateTime.Now); }

                //step 1: update responsibilities
                // find the value of k' (!=k) that maxs the quantity (a(i,k') + s(i,k')) for fixed i
                List<Thread> threadsForResp = new List<Thread>();
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    Thread aTaskThread = new Thread(this.CalculateResponsibilitiesMT);
                    threadsForResp.Add(aTaskThread);
                }
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    ArrayList paramsListForPass = new ArrayList();
                    paramsListForPass.Add((numberOfRowsPerThread * thrCtr));
                    if (thrCtr < totalNumberOfThreads - 1)
                    { paramsListForPass.Add((numberOfRowsPerThread * (thrCtr + 1)) - 1); }
                    else
                    { paramsListForPass.Add(_loopObjsByPDB.Count - 1); }
                    paramsListForPass.Add(updateLambda);
                    threadsForResp[thrCtr].Start((object)paramsListForPass);
                }

                // now, wait until the threads are done

                thrCountdown.Wait();
                thrCountdown.AddCount(totalNumberOfThreads);

                // step 2: update availabilities
                if (reportOnIteration)
                { Console.WriteLine("Finished responsibilities, starting availabilities at {0}", DateTime.Now); }
                
                List<Thread> threadsForAv = new List<Thread>();
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    Thread aTaskThread = new Thread(this.CalculateAvailabilitiesMT);
                    threadsForAv.Add(aTaskThread);
                }
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    ArrayList paramsListForPass = new ArrayList();
                    paramsListForPass.Add((numberOfRowsPerThread * thrCtr));
                    if (thrCtr < totalNumberOfThreads - 1)
                    { paramsListForPass.Add((numberOfRowsPerThread * (thrCtr + 1)) - 1); }
                    else
                    { paramsListForPass.Add(_loopObjsByPDB.Count - 1); }
                    paramsListForPass.Add(updateLambda);
                    threadsForAv[thrCtr].Start((object)paramsListForPass);
                }
                // then, wait until the threads finish

                thrCountdown.Wait();
                thrCountdown.AddCount(totalNumberOfThreads);
                
                // Step 3: compute centers
                // debugging variables here
                // define clustering after this step?

                if (reportOnIteration)
                { Console.WriteLine("Finished availabilities, computing centers at {0}", DateTime.Now); }

                List<Thread> threadsForCtr = new List<Thread>();
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    Thread aTaskThread = new Thread(this.ComputeCentersMT);
                    threadsForCtr.Add(aTaskThread);
                }
                for (int thrCtr = 0; thrCtr < totalNumberOfThreads; thrCtr++)
                {
                    ArrayList paramsListForPass = new ArrayList();
                    paramsListForPass.Add((numberOfRowsPerThread * thrCtr));
                    if (thrCtr < totalNumberOfThreads - 1)
                    { paramsListForPass.Add((numberOfRowsPerThread * (thrCtr + 1)) - 1); }
                    else
                    { paramsListForPass.Add(_loopObjsByPDB.Count - 1); }
                    threadsForCtr[thrCtr].Start((object)paramsListForPass);
                }

                thrCountdown.Wait();
                thrCountdown.AddCount(totalNumberOfThreads);

                // see if the center/exemplar assignments are the same as the previous iteration
                //bool writeToProgressFile = new bool();


                if (reportOnIteration)
                { Console.WriteLine("Finished computing centers, calculating clusters at {0}", DateTime.Now); }

                //code block: convergence by constant cluster ID
                int[] thisIterClusterData = new int[_loopObjsByPDB.Count];
                bool thisIterClusteringFailed = new bool();
                // don't need to multithread this-- not O(N-squared)?

                int[] thisIterCenters = _myModule.ClusterWithSorting(ref localCenterAssignments, ref thisIterClusterData, ref thisIterClusteringFailed);

                // end of Clustering MT block
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

                // output for FS

                if (reportOnIteration)
                { Console.WriteLine("Finished iteration {0} at {1}", iterationCounter.ToString(), DateTime.Now); }

                if ((iterationCounter == _maxIterationCount - 1) && writeProgFiles)
                {
                    string progLine = "Hit max iterations at it " + iterationCounter.ToString() +
                        ", did not converge according to test, exiting.";
                    progFileLinesToWrite.Add(progLine);
                    lastExecutedIteration = iterationCounter;
                }
            }

            // get _centerForEachPoint ready to return

            foreach (int intCenter in localCenterAssignments)
            {
                _centerForEachPoint.Add(intCenter);
            }
            // construct clusters from exemplar data
            Console.WriteLine("Starting cluster construction at {0}", DateTime.Now);

            bool finalClusteringFailed = new bool();
            int[] centersForEachCluster = _myModule.ClusterWithSorting(ref localCenterAssignments, ref bestClusterData, ref finalClusteringFailed);
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
        public void CalculateResponsibilitiesMT(object _paramsObject)
        {
            ArrayList paramsValues = (ArrayList)_paramsObject;
            int startingRow = (int)paramsValues[0];
            int endingRow = (int)paramsValues[1];
            double updateLambda = (double)paramsValues[2];
            double maxValueForResp = new double();
            for (int respCounter = startingRow; respCounter < endingRow + 1; respCounter++) // loop over i
            {
                for (int cenCounter = 0; cenCounter < totalNumberOfStructures; cenCounter++)
                {
                    bool hasMaxValueBeenSet = new bool();
                    hasMaxValueBeenSet = false;

                    for (int findingMaxInd = 0; findingMaxInd < totalNumberOfStructures; findingMaxInd++)
                    {
                        if (findingMaxInd == cenCounter)
                        {
                            continue;
                        }
                        double currentQuant = new double();
                        currentQuant = localSimilarities[respCounter, findingMaxInd];
                        if (respCounter != cenCounter)
                        {
                            currentQuant += localAvailabilities[respCounter, findingMaxInd];
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
                    localResponsibilities[respCounter, cenCounter] = (1 - updateLambda) * localResponsibilities[respCounter, cenCounter] +
                        updateLambda * (localSimilarities[respCounter, cenCounter] - maxValueForResp); // check this if using annealing: lambda versus 1-lambda backwards?
                }
            }
            thrCountdown.Signal();
            return;
        }

        public void CalculateAvailabilitiesMT(object _paramsObject)
        {
            ArrayList paramsFromMain = (ArrayList)_paramsObject;
            int firstRow = (int)paramsFromMain[0];
            int lastRow = (int)paramsFromMain[1];
            double updateLambda = (double)paramsFromMain[2];
            for (int avCtr = firstRow; avCtr < lastRow + 1; avCtr++)
            {
                for (int cenCtr = 0; cenCtr < totalNumberOfStructures; cenCtr++)
                {
                    if (avCtr == cenCtr)
                    {
                        double accumValueForSelfAv = new double();
                        accumValueForSelfAv = 0;
                        for (int pointCtr = 0; pointCtr < totalNumberOfStructures; pointCtr++)
                        {
                            if (pointCtr == cenCtr)
                            {
                                continue;
                            }
                            if (localResponsibilities[pointCtr, cenCtr] < 0)
                            {
                                accumValueForSelfAv += localResponsibilities[pointCtr, cenCtr];
                            }
                        }
                        localAvailabilities[avCtr, cenCtr] = (1 - updateLambda) * localAvailabilities[avCtr, cenCtr] +
                            updateLambda * accumValueForSelfAv;
                    }
                    else
                    {
                        double accumulateQuantity = new double();
                        accumulateQuantity = localResponsibilities[cenCtr, cenCtr];
                        for (int pointCtr = 0; pointCtr < totalNumberOfStructures; pointCtr++)
                        {
                            if ((pointCtr == avCtr) || (pointCtr == cenCtr))
                            {
                                continue;
                            }
                            if (localResponsibilities[pointCtr, cenCtr] > 0)
                            {
                                accumulateQuantity += localResponsibilities[pointCtr, cenCtr];
                            }
                        }
                        if (accumulateQuantity > 0)
                        {
                            localAvailabilities[avCtr, cenCtr] = (1 - updateLambda) * localAvailabilities[avCtr, cenCtr];
                        }
                        else
                        {
                            localAvailabilities[avCtr, cenCtr] = (1 - updateLambda) * localAvailabilities[avCtr, cenCtr] +
                                updateLambda * accumulateQuantity;
                        }
                    }
                }
            }
            thrCountdown.Signal();
            return;
        }

        public void ComputeCentersMT(object _paramsViaObject)
        {
            ArrayList paramsFromMain = (ArrayList)_paramsViaObject;
            int firstRow = (int)paramsFromMain[0];
            int lastRow = (int)paramsFromMain[1];
            for (int pointCtr = firstRow; pointCtr < lastRow + 1; pointCtr++)
            {
                double maxValForSum = new double();
                double testQuantity = new double();
                int exemplarForThisPoint = new int();
                exemplarForThisPoint = 0;
                maxValForSum = localAvailabilities[pointCtr, 0] + localResponsibilities[pointCtr, 0];
                for (int cenCtr = 1; cenCtr < totalNumberOfStructures; cenCtr++)
                {
                    testQuantity = localAvailabilities[pointCtr, cenCtr] + localResponsibilities[pointCtr, cenCtr];
                    if (testQuantity > maxValForSum)
                    {
                        maxValForSum = testQuantity;
                        exemplarForThisPoint = cenCtr;
                    }
                }
                localCenterAssignments[pointCtr] = exemplarForThisPoint;
            }

            thrCountdown.Signal();
            return;
        }

        // constructor?
        public MThreadPeptideApps(int _numberOfStructures, int _numberOfThreads)
        {
            totalNumberOfStructures = _numberOfStructures;
            totalNumberOfThreads = _numberOfThreads;
            localAvailabilities = new double[_numberOfStructures, _numberOfStructures];
            localResponsibilities = new double[_numberOfStructures, _numberOfStructures];
            localSimilarities = new double[_numberOfStructures, _numberOfStructures];
            localCenterAssignments = new int[_numberOfStructures];
            thrCountdown = new MyCountdown(_numberOfThreads);
        }

        // list of variables
        static double[,] localAvailabilities;
        static double[,] localResponsibilities;
        static double[,] localSimilarities;
        static int[] localCenterAssignments;
        int totalNumberOfStructures;
        int totalNumberOfThreads;
        MyCountdown thrCountdown;

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
