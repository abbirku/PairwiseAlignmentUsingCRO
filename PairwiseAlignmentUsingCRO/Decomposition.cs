using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class Decomposition
    {
        Random rand;
        int numOfSequences;
        int numOfColumns;
        char[,] molArr;

        public Decomposition(Random rand)
        {
            this.rand = rand;
            numOfSequences = 0;
            numOfColumns = 0;
            molArr = null;
        }

        int randSequence(bool[] is_used)
        {
            int r = 0;
            while (true)
            {
                r = rand.Next(0, numOfSequences);
                if (!is_used[r])
                {
                    return r;
                }
            }
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

        int randGap(int seqInd,bool[] is_used) {
            int r = 0;
            int x = 100;
            int random_space = -1;

            while (x > 0)
            {
                r = rand.Next(0,numOfColumns);
                if (!is_used[r])
                {
                    if (molArr[seqInd, r] == '-')
                    {
                        random_space = r;
                        is_used[r] = true;
                        break;
                    }
                }
                x--;
            }
            return random_space;
        }

        int randAlph(int seqInd)
        {
            int r = 0;
            while (true)
            {
                r = rand.Next(0, numOfColumns);
                if (molArr[seqInd, r] != '-')
                {
                    return r;
                }
            }
        }

        public MoleculeRepresentation[] decomposition(MoleculeRepresentation mol, MultipleSequenceInformation msi)
        {
            MoleculeRepresentation[] tempMolArr = new MoleculeRepresentation[2];
            for (int i = 0; i<2;i++ )
            {
                tempMolArr[i] = new MoleculeRepresentation(msi);
                tempMolArr[i].createMoleculeMatrix();
            }

            this.numOfSequences = mol.getNumOfSequences();
            this.numOfColumns = mol.getNumOfColumns();
            this.molArr = mol.getMoleculeMatrix();

            char[,] molArr1 = tempMolArr[0].getMoleculeMatrix();
            char[,] molArr2 = tempMolArr[1].getMoleculeMatrix();
            

            bool[] is_used = new bool[mol.getNumOfColumns()];
            for (int i = 0; i<mol.getNumOfColumns();i++ )
            {
                is_used[i] = false;
            }

            int rSeq = randSequence(is_used);
            int rSpa = randGap(rSeq);
            int rAlph = randAlph(rSeq);

            //For molecule one -> change space operator
            //First copying the whole array
            for (int i = 0; i < mol.getNumOfSequences(); i++)
            {
                for (int j = 0; j < mol.getNumOfColumns(); j++)
                {
                    molArr1[i, j] = molArr[i, j];
                }
            }

            //changing a random row
            if (rSpa < rAlph)
            {
                int tempMolCol = rSpa;
                for (int i = rSpa + 1; i <= rAlph; i++)
                {
                    molArr1[rSeq, tempMolCol] = molArr[rSeq, i];
                    tempMolCol++;
                }
                molArr1[rSeq, rAlph] = '-';
            }
            else
            {
                int tempMolCol = rAlph + 1;
                for (int i = rAlph; i <= (rSpa - 1); i++)
                {
                    molArr1[rSeq, tempMolCol] = molArr[rSeq, i];
                    tempMolCol++;
                }
                molArr1[rSeq, rAlph] = '-';
            }
            tempMolArr[0].setMoleculeMatrix(molArr1);

            //printing
            /*for (int i = 0; i < mol.getNumOfSequences(); i++)
            {
                for (int j = 0; j < mol.getNumOfColumns(); j++)
                {
                    Console.Write(molArr1[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();*/

            //For molecule two -> merge space operator
            rSeq = randSequence(is_used);
            bool[] is_used2 = new bool[mol.getNumOfColumns()];
            for (int i = 0; i < mol.getNumOfColumns(); i++)
            {
                is_used2[i] = false;
            }
            int random_space_1 = randGap(rSeq,is_used2);
            int random_space_2 = randGap(rSeq,is_used2);

            //copying the whole matrix
            //Console.WriteLine("Original---->");
            for (int i = 0; i < mol.getNumOfSequences(); i++)
            {
                for (int j = 0; j < mol.getNumOfColumns(); j++)
                {
                    molArr2[i, j] = molArr[i, j];
                }
            }

            if (random_space_1 > random_space_2)
            {
                int temp = random_space_2;
                random_space_2 = random_space_1;
                random_space_1 = temp;
            }

            if (random_space_1 != -1 && random_space_2 != -1)
            {
                int mol_col = random_space_1 + 2;
                if (mol_col <= (mol.getNumOfColumns() - 1))
                {
                    for (int i = random_space_1 + 1; i <= (random_space_2 - 1); i++)
                    {
                        molArr2[rSeq, mol_col] = molArr[rSeq, i];
                        mol_col++;
                    }
                    molArr2[rSeq, random_space_1 + 1] = '-';
                }
            }
            tempMolArr[1].setMoleculeMatrix(molArr2);

            /*
            Console.WriteLine("Changed---->");
            for (int i = 0; i < mol.getNumOfSequences(); i++)
            {
                for (int j = 0; j < mol.getNumOfColumns(); j++)
                {
                    Console.Write(molArr2[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();*/

            return tempMolArr;
        }
    }
}
