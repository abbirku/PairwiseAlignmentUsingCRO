using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class MultipleSequenceInformation
    {
        SingleSequenceInformation[] ssiArr; //Single Sequence Information class
        int numOfSequences; //number of sequences

        public MultipleSequenceInformation() //constructor
        {
            this.ssiArr = null;
            this.numOfSequences = 0;
        }

        public void setNumOfSequences(int numOfSequences)
        {
            this.numOfSequences = numOfSequences;
        }

        public int getNumOfSequences() //send the number of sequences to the others
        {
            return this.numOfSequences;
        }

        public char getTheCharacter(int i, int j) //get the specific character
        {
            string str = this.ssiArr[i].getTheSequence();
            return str[j];
        }

        public SingleSequenceInformation[] getTheArray() //send the single sequence information array to others
        {
            return this.ssiArr;
        }

        public void createAMultipleSequence(string[] strArr, int num)//create multiple sequences
        {
            this.ssiArr = new SingleSequenceInformation[num]; //creating a Single sequence information array

            for (int i = 0; i < num; i++) //creating object for each of the array index. Note: Without it we can't do object operation on it.
            {   //This will give null exception
                ssiArr[i] = new SingleSequenceInformation();
            }

            for (int i = 0; i < num; i++) //putting information on each of the array object element
            {
                ssiArr[i].setTheSequence(strArr[i]);
                ssiArr[i].setTheSequenceLength(strArr[i].Length);
            }
        }

        public void collectAndCreateMultipleSequences(string fileName)// collect and create multiple sequences
        {
            GetInput gi = new GetInput();
            gi.getInput(fileName);
            this.createAMultipleSequence(gi.getSequences(), gi.getNumOfSequences());
            setNumOfSequences(gi.getNumOfSequences());
        }
    }
}
