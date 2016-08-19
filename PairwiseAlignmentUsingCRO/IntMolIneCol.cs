using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class IntMolIneCol
    {

        Random rand;
        char[,] molArr1;
        char[,] molArr2;
        int numOfSequences;
        int numOfColumns;

        public IntMolIneCol(Random rand)
        {
            this.rand = rand;
            molArr1 = null;
            molArr2 = null;
            numOfColumns = 0;
            numOfSequences = 0;
        }

        public int randSequence()
        {
            return rand.Next(0, numOfSequences);
        }

        int randGap1(int seqInd)
        {
            int r = 0;
            while (true)
            {
                r = rand.Next(0, numOfColumns);
                if (molArr1[seqInd, r] == '-')
                {
                    return r;
                }
            }
        }

        int randGap2(int seqInd)
        {
            int r = 0;
            while (true)
            {
                r = rand.Next(0, numOfColumns);
                if (molArr2[seqInd, r] == '-')
                {
                    return r;
                }
            }
        }

        int nearChar1(int seqInd, int spInd)
        {
            //checking right side
            for (int i = spInd; i < numOfColumns; i++)
            {
                if (molArr1[seqInd, i] != '-')
                {
                    return i;
                }
            }

            //checking left side
            for (int j = spInd; j > -1; j--)
            {
                if (molArr1[seqInd, j] != '-')
                {
                    return j;
                }
            }

            return spInd;
        }

        int nearChar2(int seqInd, int spInd)
        {
            //checking right side
            for (int i = spInd; i < numOfColumns; i++)
            {
                if (molArr2[seqInd, i] != '-')
                {
                    return i;
                }
            }

            //checking left side
            for (int j = spInd; j > -1; j--)
            {
                if (molArr2[seqInd, j] != '-')
                {
                    return j;
                }
            }

            return spInd;
        }

        public char[,] intMolIneCollision1(MoleculeRepresentation mol1)
        {
            this.molArr1 = mol1.getMoleculeMatrix();
            this.numOfSequences = mol1.getNumOfSequences();
            this.numOfColumns = mol1.getNumOfColumns();

            int rSeq1 = randSequence();
            int rGap1 = randGap1(rSeq1);
            int nChar1 = nearChar1(rSeq1, rGap1);

            char temp1 = molArr1[rSeq1, rGap1];
            molArr1[rSeq1, rGap1] = molArr1[rSeq1, nChar1];
            molArr1[rSeq1, nChar1] = temp1;

            return molArr1;
        }

        public char[,] intMolIneCollision2(MoleculeRepresentation mol2) {
            this.molArr2 = mol2.getMoleculeMatrix();
            this.numOfSequences = mol2.getNumOfSequences();
            this.numOfColumns = mol2.getNumOfColumns();

            int rSeq2 = randSequence();
            int rGap2 = randGap2(rSeq2);
            int nChar2 = nearChar2(rSeq2, rGap2);

            char temp2 = molArr2[rSeq2, rGap2];
            molArr2[rSeq2, rGap2] = molArr2[rSeq2, nChar2];
            molArr2[rSeq2, nChar2] = temp2;

            return molArr2;
        }


    }
}
