using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
                Random rand = new Random();
                FitnessFunction ft = new FitnessFunction();

                int popSize = 1000;
                double InitialKE = 1000.0;//KE data type changed
                double KElossRate = 0.2;
                double MoleColl = 0.0001;
                int decomThresh = 500;
                int synThresh = 5;
                int buffer = 100;
                int numOfIteration = 999999;
                string lines;

                // Write the string to a file.
                System.IO.StreamWriter file = new System.IO.StreamWriter("F:\\test.txt");

                
                /*Console.Write("Population Size = ");
                popSize = Convert.ToInt32(Console.ReadLine());
                Console.Write("KELossRate = ");
                KElossRate = Convert.ToDouble(Console.ReadLine());
                Console.Write("Initial KE = ");
                InitialKE = Convert.ToDouble(Console.ReadLine());
                Console.Write("MoleColl = ");
                MoleColl = Convert.ToDouble(Console.ReadLine());
                Console.Write("Alpha = ");
                decomThresh = Convert.ToInt32(Console.ReadLine());
                Console.Write("Beta = ");
                synThresh = Convert.ToInt32(Console.ReadLine());
                Console.Write("Buffer = ");
                buffer = Convert.ToInt32(Console.ReadLine());
                Console.Write("Number of iteration = ");
                numOfIteration = Convert.ToInt32(Console.ReadLine());*/

                for (int x = 0; x<20;x++ )
                {
                    lines = "Iteration = " + x;
                    file.WriteLine(lines);
                    //Testing GetInput Class
                    GetInput giOb = new GetInput();
                    giOb.getInput("input.txt");

                    //Testing Multiple Sequence Information
                    MultipleSequenceInformation msiOb = new MultipleSequenceInformation();
                    msiOb.collectAndCreateMultipleSequences("input.txt");

                    CRO_Algorithm croAlgo = new CRO_Algorithm(rand, popSize, KElossRate, InitialKE, MoleColl, decomThresh, synThresh, buffer, numOfIteration, msiOb);
                    croAlgo.run();
                    List<MoleculeRepresentation> population = croAlgo.getPopulation();
                    MoleculeRepresentation temp = null;
                    int minPE = population[0].getMinPE();
                    temp = population[0];
                    Console.WriteLine("Popsize = " + population.Count);

                    for (int i = 1; i < population.Count; i++)
                    {
                        if (population[i].getMinPE() < minPE)
                        {
                            minPE = population[i].getMinPE();
                            temp = population[i];
                        }
                    }


                    char[,] ftArr = temp.getMoleculeMinStructure();
                    Console.WriteLine("Alignment score = -" + ft.alignmentScore(ftArr, temp.getNumOfSequences(), temp.getNumOfColumns()));
                    string str = "Alignment score = -" + ft.alignmentScore(ftArr, temp.getNumOfSequences(), temp.getNumOfColumns());
                    file.WriteLine(str);



                    if (temp != null)
                    {
                        /*lines = "";
                        lines = "Temp is not null.";
                        file.WriteLine(lines);*/
                        lines = "";
                        lines = "Minimul PE = " + temp.getMinPE() + "";
                        file.WriteLine(lines);
                        char[,] theResult = temp.getMoleculeMinStructure();

                        for (int i = 0; i < temp.getNumOfSequences(); i++)
                        {
                            lines = "";
                            for (int j = 0; j < temp.getNumOfColumns(); j++)
                            {
                                lines += theResult[i, j];
                            }
                            file.WriteLine(lines);
                        }
                    }
                    else
                    {
                        lines = "";
                        lines = "Temp is null";
                        file.WriteLine(lines);
                    }
                    lines = "\n";
                    file.WriteLine(lines);
                }

                file.Close();
            }
            catch (Exception exp)
            {
                Console.WriteLine(exp.Message.ToString());
            }
            Console.ReadLine();

        }
    }
}
