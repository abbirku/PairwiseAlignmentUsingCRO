using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class PopulationInitialization
    {
        int numOfColumns;
        MultipleSequenceInformation msi;
        List<MoleculeRepresentation> molArr;
        Random rand;
        int popSize;

        /// <summary>
        /// This constructor needs rand and msi information when the class is instantiated 
        /// </summary>
        /// <param name="rand"></param>
        /// <param name="msi"></param>
        public PopulationInitialization(Random rand, MultipleSequenceInformation msi, int popSize)
        {
            this.rand = rand;
            this.msi = msi;
            molArr = null;
            numOfColumns = 0;
            this.popSize = popSize;
        }

        /// <summary>
        /// Get Set Section of this Class
        /// This section return all the molecule information as an array
        /// and set number of columns as it is needed to generate random array
        /// </summary>
        /// <returns></returns>
        public List<MoleculeRepresentation> getPopulation()
        {
            return this.molArr;
        }

        public void setNumOfColumns(int columns)
        {
            this.numOfColumns = columns;
        }

        public List<int> getRandomArray()
        {
            int r = 0;
            int x = 0;

            bool[] is_used = new bool[numOfColumns];
            for (int i = 0; i < numOfColumns; i++)
            {
                is_used[i] = false;
            }
            List<int> permutationVector = new List<int>();

            while (x < numOfColumns)
            {
                r = rand.Next(0, numOfColumns);
                if (is_used[r])
                {
                    r = x;
                }
                if (!is_used[r])
                {
                    permutationVector.Add(r);
                    is_used[r] = true;
                    x++;
                }
            }

            return permutationVector;
        }

        public void createSolutions(double initialKE)
        {
            FitnessFunction fitFun = new FitnessFunction();

            OnWallInCol tempOnWall = new OnWallInCol(rand);
            Decomposition tempDecom = new Decomposition(rand);
            IntMolIneCol tempIntMol = new IntMolIneCol(rand);
            Synthesis tempSyn = new Synthesis(rand);

            //MoleculeRepresentation[] molReArr = new MoleculeRepresentation[popSize];
            List<MoleculeRepresentation> molReArr = new List<MoleculeRepresentation>();
            for (int i = 0; i < popSize; i++)
            {
                MoleculeRepresentation molRe = new MoleculeRepresentation(msi);
                molRe.id = i;
                molReArr.Add(molRe);
                molReArr[i].createMoleculeMatrix();
            }
            SingleSequenceInformation[] singSeInfoArr = msi.getTheArray();

            for (int i = 0; i < popSize; i++)
            {
                //molReArr[i].createMoleculeMatrix();
                numOfColumns = molReArr[i].getNumOfColumns();
                char[,] arr = molReArr[i].getMoleculeMatrix();
                for (int j = 0; j < msi.getNumOfSequences(); j++)
                {
                    List<int> vec = getRandomArray();
                    List<int> tempVec = vec.GetRange(0, singSeInfoArr[j].getTheSequenceLength());
                    tempVec.Sort();
                    for (int k = 0; k < singSeInfoArr[j].getTheSequenceLength(); k++)
                    {
                        arr[j, tempVec[k]] = msi.getTheCharacter(j, k);//arr is the omega
                    }
                }
                molReArr[i].setMoleculeMatrix(arr);
                int PE = fitFun.alignmentScore(arr, msi.getNumOfSequences(), molReArr[i].getNumOfColumns());
                molReArr[i].setMolPE(PE);
                molReArr[i].setMolKE(initialKE);
                molReArr[i].setNumHit(0);
                molReArr[i].setMoleculeMinStructure(arr);
                molReArr[i].setMinPE(PE);
                molReArr[i].setMinHit(0);
                
                //testing on wall ineffective collision
                /*
                Console.WriteLine("Molecule i = "+i+" has "+molReArr[i].getNumOfSequences()+" sequences.");
                char[,] tempChar = tempOnWall.onWallInEffCollision(molReArr[i]);
                for(int x=0;x<molReArr[i].getNumOfSequences();x++){
                    for (int y= 0; y<molReArr[i].getNumOfColumns();y++ )
                    {
                        Console.Write(tempChar[x,y]);
                    }
                    Console.WriteLine();
                }
                Console.WriteLine();
                */

                //Testing decomposition
                //tempDecom.decomposition(molReArr[i],msi);


                //testing Inter molecular ineffective collision
                /*tempIntMol.intMolIneCollision1(molReArr[i]);
                tempIntMol.intMolIneCollision2(molReArr[i]);*/

                

            }

            //tempSyn.synthesis(molReArr[0], molReArr[1], msi);

            molArr = molReArr;
        }

    }
}
