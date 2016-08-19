using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class OnWallInCol
    {
        char[,] molArr;
        Random rand;
        int numOfSequences;
        int numOfColumns;

        public OnWallInCol(Random rand)
        {
            this.molArr = null;
            this.rand = rand;
            this.numOfSequences = 0;
            this.numOfColumns = 0;
        }

        public char[,] getChangedMolMatrix()
        {
            return molArr;
        }

        public int randSequence()
        {
            return rand.Next(0, numOfSequences);
        }

        int randGap(int seqInd)
        {
            int r = 0;
            while (true)
            {
                r = rand.Next(0, numOfColumns);
                if (molArr[seqInd, r] == '-')
                {
                    return r;
                }
            }
        }

        int nearChar(int seqInd, int spInd)
        {
            //checking right side
            for (int i = spInd; i < numOfColumns; i++)
            {
                if (molArr[seqInd, i] != '-')
                {
                    return i;
                }
            }

            //checking left side
            for (int j = spInd; j > -1; j--)
            {
                if (molArr[seqInd, j] != '-')
                {
                    return j;
                }
            }

            return spInd;
        }

        public char[,] onWallInEffCollision(MoleculeRepresentation mol)
        {
            this.molArr = mol.getMoleculeMatrix();
            this.numOfSequences = mol.getNumOfSequences();
            this.numOfColumns = mol.getNumOfColumns();

            int rSeq = randSequence();
            int rGap = randGap(rSeq);
            int nChar = nearChar(rSeq, rGap);

            char temp = molArr[rSeq, rGap];
            molArr[rSeq, rGap] = molArr[rSeq, nChar];
            molArr[rSeq, nChar] = temp;

            //mol.setMoleculeMatrix(molArr);

            return molArr;
        }

    }
}
