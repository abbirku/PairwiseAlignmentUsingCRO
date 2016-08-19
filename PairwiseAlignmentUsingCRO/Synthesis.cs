using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class Synthesis
    {
        Random rand;
        int numOfSequences;
        int numOfColumns;
        char[,] molArr1;
        char[,] molArr2;

        public Synthesis(Random rand)
        {
            this.rand = rand;
            numOfSequences = 0;
            numOfColumns = 0;
            molArr1 = null;
        }

        void printArr(char[,] arr, int numOfSquences, int numOfColumns)
        {
            for (int i = 0; i < numOfSequences; i++)
            {
                for (int j = 0; j < numOfColumns; j++)
                {
                    Console.Write(arr[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public MoleculeRepresentation synthesis(MoleculeRepresentation mol1, MoleculeRepresentation mol2, MultipleSequenceInformation msi)
        {
            MoleculeRepresentation[] molReArr = new MoleculeRepresentation[2];
            for (int i = 0; i < 2; i++)
            {
                molReArr[i] = new MoleculeRepresentation(msi);
                molReArr[i].createMoleculeMatrix();
            }

            this.numOfSequences = mol1.getNumOfSequences();
            this.numOfColumns = mol1.getNumOfColumns();
            this.molArr1 = mol1.getMoleculeMatrix();
            this.molArr2 = mol2.getMoleculeMatrix();

            char[,] wp1 = molReArr[0].getMoleculeMatrix();//womega prime 1
            char[,] wp2 = molReArr[1].getMoleculeMatrix();//womega prime 2

            double d = numOfSequences / 2.0;
            int bottomSection = (int)Math.Ceiling(d);
            int upperSection = numOfSequences - bottomSection;

            for (int i = 0; i < upperSection; i++)
            {
                for (int j = 0; j < numOfColumns; j++)
                {
                    wp1[i, j] = molArr1[i, j];
                    wp2[i, j] = molArr2[i, j];
                }
            }

            for (int i = upperSection; i < numOfSequences; i++)
            {
                for (int j = 0; j < numOfColumns; j++)
                {
                    wp1[i, j] = molArr2[i, j];
                    wp2[i, j] = molArr1[i, j];
                }
            }

            molReArr[0].setMoleculeMatrix(wp1);
            molReArr[1].setMoleculeMatrix(wp2);

            /*printArr(molArr1,numOfSequences,numOfColumns);
            printArr(molArr2, numOfSequences, numOfColumns);
            printArr(wp1, numOfSequences, numOfColumns);
            printArr(wp2, numOfSequences, numOfColumns);*/

            d = rand.NextDouble();

            if (d <= 0.5)
            {
                return molReArr[0];
            }
            else
            {
                return molReArr[1];
            }
        }
    }
}
