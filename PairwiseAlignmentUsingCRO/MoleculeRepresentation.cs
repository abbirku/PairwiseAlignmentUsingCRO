using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class MoleculeRepresentation
    {
        int k_sequences;//number of sequence
        int maxLenOf_k_Sequences;//Maximum length
        int numOfColumns_w;//number of column
        MultipleSequenceInformation msi;

        char[,] moleculeMatrix;
        int molPE;
        double molKE;
        int numHit;
        char[,] minStruct;
        int minPE;
        int minHit;
        public int id;

        /*public MoleculeRepresentation()
        {
            this.k_sequences = 0;
            this.maxLenOf_k_Sequences = 0;
            this.numOfColumns_w = 0;
            this.msi = null;

            this.moleculeMatrix = null;
            this.molPE = 0;
            this.molKE = 0.0;
            this.numHit = 0;
            this.minStruct = null;
            this.minPE = 0;
            this.minHit = 0;
        }*/

        public MoleculeRepresentation(MultipleSequenceInformation msi)
        {
            this.k_sequences = msi.getNumOfSequences();
            this.maxLenOf_k_Sequences = 0;
            this.numOfColumns_w = 0;
            this.msi = msi;

            this.moleculeMatrix = null;
            this.molPE = 0;
            this.molKE = 0.0;
            this.numHit = 0;
            this.minStruct = null;
            this.minPE = 0;
            this.minHit = 0;
            this.id = 0;
        }

        //Get Set section of each field
        public void set_K_Sequences(int k_sequences)
        {
            this.k_sequences = k_sequences;
        }

        public int getNumOfSequences()
        {
            return this.k_sequences;
        }

        public void setMSI(MultipleSequenceInformation m)
        {
            this.msi = m;
        }

        public int getNumOfColumns()
        {
            return this.numOfColumns_w;
        }

        public void setMoleculeMatrix(char[,] chArr)
        {
            this.moleculeMatrix = chArr;
        }

        public char[,] getMoleculeMatrix()
        {
            return this.moleculeMatrix;
        }

        public void setMolPE(int PE)
        {
            this.molPE = PE;
        }

        public int getMolPE()
        {
            return this.molPE;
        }

        public void setMolKE(double KE)
        {
            this.molKE = KE;
        }

        public double getMolKE() {
            return this.molKE;
        }

        public void setNumHit(int hit)
        {
            this.numHit = hit;
        }

        public int getNumHit()
        {
            return this.numHit;
        }

        public void setMoleculeMinStructure(char[,] structure)
        {
            this.minStruct = structure;
        }

        public char[,] getMoleculeMinStructure()
        {
            return this.minStruct;
        }

        public void setMinPE(int minPE)
        {
            this.minPE = minPE;
        }

        public int getMinPE()
        {
            return this.minPE;
        }

        public void setMinHit(int minHit)
        {
            this.minHit = minHit;
        }

        public int getMinHit()
        {
            return this.minHit;
        }

        //Operation Section of this class

        //Calculates the column by multiplying 1.2 with max sequence length
        public void calculateColumn()
        {
            int maxLen = 0;
            SingleSequenceInformation[] ssiArrOb = msi.getTheArray();
            for (int i = 0; i < k_sequences; i++)
            {
                if (ssiArrOb[i].getTheSequenceLength() > maxLen)
                {
                    maxLen = ssiArrOb[i].getTheSequenceLength();
                }
            }
            maxLenOf_k_Sequences = maxLen;
            numOfColumns_w = (int)Math.Ceiling(1.2 * maxLenOf_k_Sequences);
            //Console.WriteLine(numOfColumns_w);
        }

        //Create Molecule matrix with blank space
        public void createMoleculeMatrix()
        {
            calculateColumn();
            moleculeMatrix = new char[k_sequences, numOfColumns_w];
            for (int i = 0; i < k_sequences; i++)
            {
                for (int j = 0; j < numOfColumns_w; j++)
                {
                    moleculeMatrix[i, j] = '-';
                }
            }
        }

        public void createMoleculeMatrix(int row, int column)
        {
            moleculeMatrix = new char[row, column];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < column; j++)
                {
                    moleculeMatrix[i, j] = '-';
                }
            }
            numOfColumns_w = column;
            k_sequences = row;
            /*parentMatrix = new char*[row];
            for (int i = 0; i < row; i++)
            {
                parentMatrix[i] = new char[column];
                for (int j = 0; j < numOfColumns_w; j++)
                {
                    parentMatrix[i][j] = '-';
                }
            }
            numOfColumns_w = column;*/
        }
    }
}
