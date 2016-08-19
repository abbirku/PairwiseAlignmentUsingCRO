using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class GetInput
    {
        int numOfSequences;
        string[] sequences;

        public GetInput()
        {
            numOfSequences = 0;
        }

        public int getNumOfSequences()
        {
            return numOfSequences;
        }
        public string[] getSequences()
        {
            return sequences;
        }

        public void getInput(string fileName)
        {

            try
            {
                sequences = System.IO.File.ReadAllLines(fileName);
                numOfSequences = sequences.Length;
            }
            catch (Exception exp)
            {
                Console.WriteLine(exp.ToString());
            }

        }
    }
}
