using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PairwiseAlignmentUsingCRO
{
    class CRO_Algorithm
    {
        int popSize;
        double KElossRate;
        double initialKE;//KE data type changed
        double moleColl;
        int decomThresh;
        int synThresh;
        double buffer;
        Random rand;
        MoleculeRepresentation globalMol;
        int numOfIteration;
        int globalMinPE;
        int minPE;
        int num_of_iter;
        MultipleSequenceInformation msi;


        FitnessFunction fitOb;
        OnWallInCol onWallOb;
        Decomposition decomOb;
        IntMolIneCol intMolColOb;
        Synthesis synthesisOb;
        List<MoleculeRepresentation> population;

        public CRO_Algorithm(Random rand, int popSize, double KElossRate, double initialKE, double moleColl, int decomThresh, int synThresh, double buffer, int numOfIteration, MultipleSequenceInformation m)
        {
            this.rand = rand;
            this.popSize = popSize;
            this.KElossRate = KElossRate;
            this.initialKE = initialKE;//KE data type changed
            this.moleColl = moleColl;
            this.decomThresh = decomThresh;
            this.synThresh = synThresh;
            this.buffer = buffer;
            this.globalMol = new MoleculeRepresentation(m);
            this.numOfIteration = numOfIteration;
            this.globalMinPE = 99999999;
            this.minPE = 99999999;
            this.num_of_iter = 0;
            this.msi = m;


            this.fitOb = new FitnessFunction();
            onWallOb = new OnWallInCol(rand);
            decomOb = new Decomposition(rand);
            intMolColOb = new IntMolIneCol(rand);
            synthesisOb = new Synthesis(rand);

            population = null;
        }

        int randMole()
        {
            return rand.Next(0, population.Count);
        }

        public List<MoleculeRepresentation> getPopulation()
        {
            return population;
        }

        public int getNumOfIteration() {
            return this.num_of_iter;
        }

        void onWallIneffectiveCollision(MoleculeRepresentation selMol)
        {
            int PE, PEomegaP;
            double KE, KEomegaP;

            char[,] omegaP = onWallOb.onWallInEffCollision(selMol);
            PEomegaP = fitOb.alignmentScore(omegaP, selMol.getNumOfSequences(), selMol.getNumOfColumns());
            int hit = selMol.getNumHit();
            hit += 1;
            selMol.setNumHit(hit);

            PE = selMol.getMolPE();
            KE = selMol.getMolKE();

            if (PE + KE >= PEomegaP)
            {
                KEomegaP = ((PE - PEomegaP + KE) * KElossRate);
                buffer = buffer + (PE - PEomegaP + KE) * (1 - KElossRate);
                selMol.setMoleculeMatrix(omegaP);
                selMol.setMolPE(PEomegaP);
                selMol.setMolKE(KEomegaP);
                if (selMol.getMolPE() < selMol.getMinPE())
                {
                    selMol.setMoleculeMinStructure(selMol.getMoleculeMatrix());
                    selMol.setMinPE(selMol.getMolPE());
                    selMol.setMinHit(selMol.getNumHit());
                    //added
                    this.minPE = selMol.getMinPE();
                }
            }
        }



        void decompositionWork(MoleculeRepresentation selMol, int molIndex, List<MoleculeRepresentation> population)
        {
            Console.WriteLine("Decomposed");
            //variable declaration
            int PEomegaP1, PEomegaP2, PE;
            double del1, del2, del3, KE;
            double KEwp1, KEwp2;
            double Edec;
            int hit;

            PE = selMol.getMolPE();
            KE = selMol.getMolKE();

            //creating two decomposed molecule
            MoleculeRepresentation[] chanMols;
            chanMols = decomOb.decomposition(selMol, msi);

            //getting the decomposed structure
            char[,] omegaP1, omegaP2;
            omegaP1 = chanMols[0].getMoleculeMatrix();
            omegaP2 = chanMols[1].getMoleculeMatrix();

            //calculating their fitness
            PEomegaP1 = fitOb.alignmentScore(omegaP1, selMol.getNumOfSequences(), selMol.getNumOfColumns());
            PEomegaP2 = fitOb.alignmentScore(omegaP2, selMol.getNumOfSequences(), selMol.getNumOfColumns());

            if (PE + KE >= PEomegaP1 + PEomegaP2)
            {
                Edec = PE + KE - (PEomegaP1 + PEomegaP2);
                del3 = rand.NextDouble();
                KEwp1 = (Edec * del3);
                KEwp2 = (Edec * (1 - del3));

                chanMols[0].setMolPE(PEomegaP1);
                chanMols[1].setMolPE(PEomegaP2);
                chanMols[0].setMolKE(KEwp1);
                chanMols[1].setMolKE(KEwp2);

                chanMols[0].setMoleculeMinStructure(omegaP1);
                chanMols[1].setMoleculeMinStructure(omegaP2);

                chanMols[0].setMinPE(chanMols[0].getMolPE());
                chanMols[1].setMinPE(chanMols[1].getMolPE());

                population.RemoveAt(molIndex);
                population.Add(chanMols[0]);
                population.Add(chanMols[1]);

                //added
                if (minPE > chanMols[0].getMinPE())
                {
                    minPE = chanMols[0].getMinPE();
                }

                if (minPE > chanMols[1].getMinPE())
                {
                    minPE = chanMols[1].getMinPE();
                }

            }
            else
            {
                del1 = rand.NextDouble();
                del2 = rand.NextDouble();

                Edec = (PE + KE + (del1 * del2 * buffer) - (PEomegaP1 + PEomegaP2));
                if (Edec >= 0)
                {
                    buffer = buffer * (1 - del1 * del2);


                    del3 = rand.NextDouble();
                    KEwp1 = (Edec * del3);
                    KEwp2 = (Edec * (1 - del3));

                    chanMols[0].setMolPE(PEomegaP1);
                    chanMols[1].setMolPE(PEomegaP2);
                    chanMols[0].setMolKE(KEwp1);
                    chanMols[1].setMolKE(KEwp2);

                    chanMols[0].setMoleculeMinStructure(omegaP1);
                    chanMols[1].setMoleculeMinStructure(omegaP2);

                    chanMols[0].setMinPE(chanMols[0].getMolPE());
                    chanMols[1].setMinPE(chanMols[1].getMolPE());

                    /*if(!(molIndex>population.Count-1)){
                        
                    }*/

                    population.RemoveAt(molIndex);

                    population.Add(chanMols[0]);
                    population.Add(chanMols[1]);

                    //added
                    if (minPE > chanMols[0].getMinPE())
                    {
                        minPE = chanMols[0].getMinPE();
                    }

                    if (minPE > chanMols[1].getMinPE())
                    {
                        minPE = chanMols[1].getMinPE();
                    }


                }
                else
                {
                    hit = selMol.getNumHit();
                    hit += 1;
                    selMol.setNumHit(hit);
                }
            }
        }

        void interMolecularIneffectiveCollision(MoleculeRepresentation selMol1, MoleculeRepresentation selMol2)
        {
            char[,] wp1, wp2;
            int PEwp1, PEwp2, hit, PEw1, PEw2;
            double Einter, KEw1, KEw2, del4, KEwp1, KEwp2;

            PEw1 = selMol1.getMolPE();
            PEw2 = selMol2.getMolPE();
            KEw1 = selMol1.getMolKE();
            KEw2 = selMol2.getMolKE();

            //after collision structure
            wp1 = intMolColOb.intMolIneCollision1(selMol1);
            wp2 = intMolColOb.intMolIneCollision2(selMol2);

            //calculating their firness
            PEwp1 = fitOb.alignmentScore(wp1, selMol1.getNumOfSequences(), selMol1.getNumOfColumns());
            PEwp2 = fitOb.alignmentScore(wp2, selMol2.getNumOfSequences(), selMol2.getNumOfColumns());

            //assigning number of hit
            hit = selMol1.getNumHit();
            hit += 1;
            selMol1.setNumHit(hit);

            hit = selMol2.getNumHit();
            hit += 1;
            selMol2.setNumHit(hit);

            //calculating intermolecular posibility
            Einter = (PEw1 + PEw2 + KEw1 + KEw2) - (PEwp1 + PEwp2);
            if (Einter >= 0)
            {
                del4 = rand.NextDouble();

                //calculating KE for changed molecule
                KEwp1 = Einter * del4;
                KEwp2 = Einter * (1 - del4);

                selMol1.setMoleculeMatrix(wp1);
                selMol2.setMoleculeMatrix(wp2);

                selMol1.setMolPE(PEwp1);
                selMol2.setMolPE(PEwp2);

                selMol1.setMolKE(KEwp1);
                selMol2.setMolKE(KEwp2);

                if (selMol1.getMolPE() < selMol1.getMinPE())
                {
                    selMol1.setMoleculeMinStructure(selMol1.getMoleculeMatrix());
                    selMol1.setMinPE(selMol1.getMolPE());
                    selMol1.setMinHit(selMol1.getNumHit());
                    //added
                    this.minPE = selMol1.getMinPE();
                }
                if (selMol2.getMolPE() < selMol2.getMinPE())
                {
                    selMol2.setMoleculeMinStructure(selMol2.getMoleculeMatrix());
                    selMol2.setMinPE(selMol2.getMolPE());
                    selMol2.setMinHit(selMol2.getNumHit());
                    //added
                    this.minPE = selMol2.getMinPE();
                }
            }
        }

        void synthesisWork(MoleculeRepresentation selMol1, MoleculeRepresentation selMol2, int molIndex1, int molIndex2, List<MoleculeRepresentation> population)
        {
            Console.WriteLine("Synthesis");
            MoleculeRepresentation synMol;
            char[,] wp;
            int PEwp, PEw1, PEw2, hit;
            double KEwp, KEw1, KEw2;

            PEw1 = selMol1.getMolPE();
            PEw2 = selMol2.getMolPE();
            KEw1 = selMol1.getMolKE();
            KEw2 = selMol2.getMolKE();

            synMol = synthesisOb.synthesis(selMol1, selMol2, msi);
            wp = synMol.getMoleculeMatrix();

            PEwp = fitOb.alignmentScore(wp, synMol.getNumOfSequences(), synMol.getNumOfColumns());

            if ((PEw1 + PEw2 + KEw1 + KEw2) >= PEwp)
            {
                KEwp = (PEw1 + PEw2 + KEw1 + KEw2) - PEwp;
                synMol.setMolPE(PEwp);
                synMol.setMolKE(KEwp);
                synMol.setMoleculeMinStructure(wp);
                synMol.setMinPE(PEwp);

                int id1 = selMol1.id;
                int id2 = selMol2.id;

                for (int i = 0; i < population.Count; i++)
                {
                    if (population[i].id == id1)
                    {
                        population.RemoveAt(i);
                        break;
                    }
                }

                for (int i = 0; i < population.Count; i++)
                {
                    if (population[i].id == id2)
                    {
                        population.RemoveAt(i);
                        break;
                    }
                }


                /*MoleculeRepresentation[] remItem = new[] { selMol1, selMol2 };
                population.RemoveAll(x => remItem.Contains(x));*/

                population.Add(synMol);

                if (minPE > synMol.getMinPE())
                {
                    minPE = synMol.getMinPE();
                }
            }
            else
            {
                hit = selMol1.getNumHit();
                hit += 1;
                selMol1.setNumHit(hit);

                hit = selMol2.getNumHit();
                hit += 1;
                selMol2.setNumHit(hit);
            }
        }

        public void run()
        {
            //List<MoleculeRepresentation> population = new List<MoleculeRepresentation>();
            int iter = 0;
            double t;
            int molIndex;
            int anotherMolIndex;
            int numOfHit, minNumHit;
            int stop_criteria = 0;

            PopulationInitialization popInit = new PopulationInitialization(rand, msi, popSize);
            popInit.createSolutions(initialKE);
            this.population = popInit.getPopulation();
            int i = 0;

            var watch = System.Diagnostics.Stopwatch.StartNew();

            while (iter <= numOfIteration)
            {
                t = rand.NextDouble();
                if (t > moleColl)//Uni-molicular
                {
                    molIndex = randMole();
                    numOfHit = population[molIndex].getNumHit();
                    minNumHit = population[molIndex].getMinHit();
                    if ((numOfHit - minNumHit) > decomThresh)//decomposition
                    {
                        decompositionWork(population[molIndex], molIndex, population);
                        iter += 2;
                    }
                    else//on wall
                    {
                        onWallIneffectiveCollision(population[molIndex]);
                        iter += 1;
                    }
                }
                else//inter-molecular
                {
                    molIndex = randMole();
                    anotherMolIndex = randMole();
                    if (molIndex == anotherMolIndex)
                    {
                        anotherMolIndex = randMole();
                    }

                    //synthesis
                    if ((population[molIndex].getMolKE() <= synThresh) && (population[anotherMolIndex].getMolKE() <= synThresh))
                    {
                        synthesisWork(population[molIndex], population[anotherMolIndex], molIndex, anotherMolIndex, population);
                        iter += 2;
                    }
                    //inter nolecular ineffective collision
                    else
                    {
                        interMolecularIneffectiveCollision(population[molIndex], population[anotherMolIndex]);
                        iter++;
                    }
                }
                if (i == 10000)
                {
                    //Console.WriteLine("Now at = " + iter);
                    i = 0;
                }
                i++;

                if (minPE != globalMinPE)
                {
                    if((minPE < globalMinPE)){
                        globalMinPE = minPE;
                    }
                    stop_criteria = 0;
                }
                else
                {
                    stop_criteria += 1;
                }

                if(stop_criteria==50){
                    break;
                }
                this.num_of_iter += 1;

            }

            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;
            Console.WriteLine("Taken time = " + elapsedMs / 1000 + " seconds");

            //this.population = population;

        }

    }
}
