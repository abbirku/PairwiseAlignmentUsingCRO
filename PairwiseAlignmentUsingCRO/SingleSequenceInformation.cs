using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class SingleSequenceInformation
    {
        string theSequence;
        int sequenceLength;

        public SingleSequenceInformation()
        {
            this.theSequence = "";
            this.sequenceLength = 0;
        }

        public void setTheSequence(string sequence)
        {
            this.theSequence = sequence;
        }

        public string getTheSequence()
        {
            return this.theSequence;
        }

        public void setTheSequenceLength(int length)
        {
            this.sequenceLength = length;
        }

        public int getTheSequenceLength()
        {
            return this.sequenceLength;
        }
    }
}
