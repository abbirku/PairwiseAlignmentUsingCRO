using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class FitnessFunction
    {
        int[,] blastMatrix;
        Dictionary<char, int> dic;

        public FitnessFunction() { 
            blastMatrix = new int[,]{{5,-4,-4,-4,-10},
                                     {-4,5,-4,-4,-10},
                                     {-4,-4,5,-4,-10},
                                     {-4,-4,-4,5,-10},
                                     {-10,-10,-10,-10,0}};

            dic = new Dictionary<char, int>();
            dic.Add('A',0);
            dic.Add('T',1);
            dic.Add('C', 2);
            dic.Add('G', 3);
            dic.Add('-', 4);
        }

        int score(char c1,char c2) { 
            int x = dic[c1];
            int y = dic[c2];

            return blastMatrix[x,y];
        }

        int pairScore(char[,] alignmentArray,int numOfColumns,int sequence_i,int sequence_j) {
            int sum = 0;
            char c1='-', c2='-';
            for (int p = 0; p<numOfColumns;p++ )
            {
                c1 = alignmentArray[sequence_i, p];
                c2 = alignmentArray[sequence_j, p];
                sum += score(c1,c2);
            }
            return sum;
        }

        int generateAndCalculatePairScore(char[,] alignmentArray,int numOfSequences,int numOfColumns) {
            int sum = 0;
            for (int i = 0; i<numOfSequences-1;i++ )
            {
                for (int j = i+1; j<numOfSequences;j++ )
                {
                    sum += pairScore(alignmentArray,numOfColumns,i,j);
                }
            }
            return sum;
        }

        public int alignmentScore(char[,] alignmentArray,int numOfSequences,int numOfColumns) {
           int score = generateAndCalculatePairScore(alignmentArray,numOfSequences,numOfColumns)*(-1);
           return score;
        }
    }
}
