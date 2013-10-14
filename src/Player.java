import java.util.Arrays;
import java.util.Random;

class Player {
    private static final double MINIMUM_HIT_PROB                          = 0.6;
    private static final double MINIMUM_SEQUENCEPROB_FOR_BIRDS_OWN_MATRIX = 1;
    private static       int    NR_STATES                                 = 5;
    private static       int    SKIP_CALCULATING_SHOOT_UNTIL              = 50;
    private static       int    DO_NOT_SHOOT_UNTIL                        = 70;
    private static       int    MAX_NR_STEPS                              = 100;
    private static final int    MAXIMUM_SAVED_MATRICES                    = 10;


    private int lastRound = -1;
    private int currentBirdStep;

    double[][][] currentBirdsTransitionMatrices;
    double[][][] currentBirdsEmissionMatrices;
    double[][]   currentBirdsInitialStatePDVectors;
    double[]     lastStepGuess;
    boolean[][]  stepCorrectGuesses;
    private boolean guessedSpecies = false;

    double[][][][] savedBirdsTransitionMatrics;
    double[][][][] savedBirdsEmissionMatrices;
    private int[] currentGuessedSpecies;

    //    private double[] lGuessPercentage;


    //Information variables
    int numberOfCorrectGuesses  = 0;
    int numberOfGuesses         = 0;
    int numberOfShots           = 0;
    int numberOfHits            = 0;
    int numberOfBirdsStillAlive = 0;
    int totalNumberOfBirds      = 0;
    int[] birdDiedOrIsAtStep;

    // /constructor

    // /There is no data in the beginning, so not much should be done here.
    public Player() {
    }

    /**
     * Shoot!
     * <p/>
     * This is the function where you start your work.
     * <p/>
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each birds contains all past actions.
     * <p/>
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue   time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
        if (pState.getRound() != lastRound) {
            totalNumberOfBirds += pState.getNumBirds();
            numberOfBirdsStillAlive += pState.getNumBirds();
            lastRound = pState.getRound();
            currentBirdStep = 0;
            initialize(pState.getNumBirds());
            System.err.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!New round: " + pState.getRound());
        } else {
            currentBirdStep++;
        }
        if (currentBirdStep < SKIP_CALCULATING_SHOOT_UNTIL) {
            return cDontShoot;
        }
        int birdsAlive = 0;

        for (int i = 0; i < pState.getNumBirds(); i++) {
            Bird bird = pState.getBird(i);
            if (bird.isAlive()) {
                birdsAlive++;
                reCalculateHMM(i, bird, 50);
                birdDiedOrIsAtStep[i] = currentBirdStep;
            }
        }

        if (currentBirdStep == DO_NOT_SHOOT_UNTIL) {
            computeGuess(pState);
        }
        Action action;
        action = getBestShootingAction(pState);
        if (currentBirdStep < DO_NOT_SHOOT_UNTIL || pState.getRound() == 0) {
            action = cDontShoot;
        }


        //        System.out.println("!!!!TEST!!!!");
        // This line choose not to shoot
        if (pDue.remainingMs() < 1) {
            throw new RuntimeException("Time is out! Time left: " + pDue.remainingMs() + " ms");
        }
//        System.err.println("Shooting time left: " + pDue.remainingMs());
        if (!action.equals(cDontShoot)) {

            double emissionProb = calculateEmissionSequenceProb(currentBirdsTransitionMatrices[action.getBirdNumber()], currentBirdsEmissionMatrices[action.getBirdNumber()], currentBirdsInitialStatePDVectors[action.getBirdNumber()], currentBirdStep + 1, pState.getBird(action.getBirdNumber()));
            if (emissionProb < 100000) {
                action = cDontShoot;
            } else {
                numberOfShots++;
                System.err.println("Birds alive: " + birdsAlive);
                //                System.err.println("EmissionProb:" + emissionProb);
            }
        }
        return action;

        // This line would predict that bird 0 will move right and shoot at it
        //		return new Action(0, Constants.MOVE_RIGHT);
    }

    private void initialize(int numBirds) {
        currentGuessedSpecies = new int[numBirds];
        if (savedBirdsEmissionMatrices == null) {
            savedBirdsTransitionMatrics = new double[Constants.COUNT_SPECIES][MAXIMUM_SAVED_MATRICES][][];
            savedBirdsEmissionMatrices = new double[Constants.COUNT_SPECIES][MAXIMUM_SAVED_MATRICES][][];
        }

        birdDiedOrIsAtStep = new int[numBirds];
        lastStepGuess = new double[numBirds];
        stepCorrectGuesses = new boolean[numBirds][MAX_NR_STEPS];
        currentBirdsTransitionMatrices = new double[numBirds][NR_STATES][NR_STATES];
        currentBirdsEmissionMatrices = new double[numBirds][NR_STATES][Constants.COUNT_MOVE];
        currentBirdsInitialStatePDVectors = new double[numBirds][NR_STATES];

        for (int birdI = 0; birdI < numBirds; birdI++) {
            randomizeAndNormalizeMatrix(currentBirdsTransitionMatrices[birdI]);
            randomizeAndNormalizeMatrix(currentBirdsEmissionMatrices[birdI]); //TODO: Not randomize this one, be smart somehow
            randomizeAndNormalizeVector(currentBirdsInitialStatePDVectors[birdI]);
        }
    }

    private void randomizeAndNormalizeMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            randomizeAndNormalizeVector(matrix[i]);
        }
    }

    private void randomizeAndNormalizeVector(double[] vector) {
        Random r = new Random();
        double totalPercentage = 0;
        double tempValue;
        for (int i = 0; i < vector.length; i++) {
            tempValue = r.nextDouble();
            vector[i] = tempValue;
            totalPercentage += tempValue;
        }
        normalizeVector(vector, totalPercentage);
    }

    private void normalizeVector(double[] vectorToNormalize, double totalValue) {
        for (int i = 0; i < vectorToNormalize.length; i++) {
            vectorToNormalize[i] /= totalValue;
        }
    }

    private void reCalculateHMM(int birdIndex, Bird bird, int maximumIterations) {
        double[][] currentTM = currentBirdsTransitionMatrices[birdIndex];
        double[][] currentEM = currentBirdsEmissionMatrices[birdIndex];
        double[] currentV = currentBirdsInitialStatePDVectors[birdIndex];

        updateMatrices(currentTM, currentEM, currentV, bird, currentBirdStep + 1, maximumIterations);
        double[][] tempTM = new double[NR_STATES][NR_STATES];
        double[][] tempEM = new double[NR_STATES][Constants.COUNT_MOVE];
        double[] tempV = new double[NR_STATES];

        randomizeAndNormalizeMatrix(tempTM);
        randomizeAndNormalizeMatrix(tempEM);
        randomizeAndNormalizeVector(tempV);
        updateMatrices(tempTM, tempEM, tempV, bird, currentBirdStep + 1, maximumIterations);

        double emissionProb1 = calculateEmissionSequenceProb(currentTM, currentEM, currentV, birdDiedOrIsAtStep[birdIndex] + 1, bird);
        double emissionProb2 = calculateEmissionSequenceProb(tempTM, tempEM, tempV, birdDiedOrIsAtStep[birdIndex] + 1, bird);
        if (emissionProb2 > emissionProb1) {
            currentBirdsTransitionMatrices[birdIndex] = tempTM;
            currentBirdsEmissionMatrices[birdIndex] = tempEM;
            currentBirdsInitialStatePDVectors[birdIndex] = tempV;
        }
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     * <p/>
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue   time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        System.err.println(pDue.remainingMs());
        printRoundInfo();
        /*
         * Here you should write your clever algorithms to guess the species of
		 * each bird. This skeleton makes no guesses, better safe than sorry!
		 */
        computeFinalGuess(pState);

        printVector(currentGuessedSpecies);

        printSequenceLikeleyHoodForOwnMatrix(pState);

        return currentGuessedSpecies;
    }

    private void computeGuess(GameState pState) {
        currentGuessedSpecies = new int[pState.getNumBirds()];
        Arrays.fill(currentGuessedSpecies, Constants.SPECIES_UNKNOWN);

        for (int i = 0; i < pState.getNumBirds(); i++) {
            Bird bird = pState.getBird(i);
            double bestEmissionProb = 0;
            int specie = 0;
            for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
                double bestEmissionForS = 0;

                for (int t = 0; t < MAXIMUM_SAVED_MATRICES; t++) {
                    double[][] tm = savedBirdsTransitionMatrics[s][t];
                    double[][] em = savedBirdsEmissionMatrices[s][t];
                    if (em == null) {
                        continue;
                    }
                    double[] initialStatePDVector = new double[NR_STATES];
                    initialStatePDVector[0] = 1;
                    for (int j = 0; j < NR_STATES; j++) {


                        double emissionProb = calculateEmissionSequenceProb(tm, em, initialStatePDVector, birdDiedOrIsAtStep[i] + 1, bird);
                        if (emissionProb > bestEmissionForS) {
                            bestEmissionForS = emissionProb;
                        }

                        if (j != NR_STATES - 1) {
                            initialStatePDVector[j] = 0;
                            initialStatePDVector[j + 1] = 1;
                        }
                    }
                    if (bestEmissionForS > bestEmissionProb) {
                        bestEmissionProb = bestEmissionForS;
                        specie = s;
                    }
                }
                System.err.println(" (" + i + "-" + s + ") " + bestEmissionForS);
            }

            if (bestEmissionProb > 1000000000000000.0) {
                currentGuessedSpecies[i] = specie;
            }
        }


        //            for (int i = 0; i < pSpecies.length; i++) {
        //                System.err.print(" (" + i + ") " + pSpecies[i]);
        //            }
        //            System.err.println();


    }


    private void computeFinalGuess(GameState pState) {
        currentGuessedSpecies = new int[pState.getNumBirds()];
        Arrays.fill(currentGuessedSpecies, Constants.SPECIES_UNKNOWN);

        if (!guessedSpecies) {
            //JUst guess
            for (int i = 0; i < pState.getNumBirds(); ++i) {
                currentGuessedSpecies[i] = Constants.SPECIES_SKYLARK;
            }
            guessedSpecies = true;
        } else {

            for (int i = 0; i < pState.getNumBirds(); i++) {
                Bird bird = pState.getBird(i);
                double bestEmissionProb = 0;
                int specie = 0;
                for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
                    double bestEmissionForS = 0;

                    for (int t = 0; t < MAXIMUM_SAVED_MATRICES; t++) {
                        double[][] tm = savedBirdsTransitionMatrics[s][t];
                        double[][] em = savedBirdsEmissionMatrices[s][t];
                        if (em == null) {
                            continue;
                        }
                        double[] initialStatePDVector = new double[NR_STATES];
                        initialStatePDVector[0] = 1;
                        for (int j = 0; j < NR_STATES; j++) {


                            double emissionProb = calculateEmissionSequenceProb(tm, em, initialStatePDVector, birdDiedOrIsAtStep[i] + 1, bird);
                            if (emissionProb > bestEmissionForS) {
                                bestEmissionForS = emissionProb;
                            }

                            if (j != NR_STATES - 1) {
                                initialStatePDVector[j] = 0;
                                initialStatePDVector[j + 1] = 1;
                            }
                        }
                        if (bestEmissionForS > bestEmissionProb) {
                            bestEmissionProb = bestEmissionForS;
                            specie = s;
                        }
                    }
                    System.err.println(" (" + i + "-" + s + ") " + bestEmissionForS);
                }

                if (bestEmissionProb > 0) {
                    currentGuessedSpecies[i] = specie;
                } else {
                    currentGuessedSpecies[i] = Constants.SPECIES_UNKNOWN;
                    for (int k = 0; k < savedBirdsTransitionMatrics.length; k++) {
                        if (savedBirdsTransitionMatrics[k][0] == null) {
                            currentGuessedSpecies[i] = k;
                            break;
                        }
                    }
                }
            }


            //            for (int i = 0; i < pSpecies.length; i++) {
            //                System.err.print(" (" + i + ") " + pSpecies[i]);
            //            }
            //            System.err.println();

        }
    }

    private void printSequenceLikeleyHoodForOwnMatrix(GameState pState) {
        System.err.print("Bird own emission likelyhood");
        for (int i = 0; i < pState.getNumBirds(); i++) {
            double sequence = calculateEmissionSequenceProb(i, pState.getBird(i));
            System.err.println("Bird: " + i + ": emissionLikelyHood " + sequence);
        }
    }

    private void printVector(int[] o) {
        for (Object t : o) {
            System.err.print(t + " ");
        }
        System.err.println();
    }

    private void printVector(double[] o) {
        for (Object t : o) {
            System.err.print(t + " ");
        }
        System.err.println();
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird  the bird you hit
     * @param pDue   time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD: " + pBird + "!!!");
        numberOfHits++;
        numberOfBirdsStillAlive--;
    }

    private void printRoundInfo() {
        System.err.println("Shots: " + numberOfShots);
        System.err.println("Hits: " + numberOfHits);
        System.err.println("Total birds: " + totalNumberOfBirds);
        System.err.println("Birds still alive: " + numberOfBirdsStillAlive);
        System.err.println("Total hit percentage: " + (double) numberOfHits / numberOfShots);

        System.err.println("Birds alive percentage: " + (double) numberOfBirdsStillAlive / totalNumberOfBirds);
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState   the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue     time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
        for (int i = 0; i < pState.getNumBirds(); i++) {
            int specie = pSpecies[i];
            if (specie == Constants.SPECIES_BLACK_STORK) {
                System.err.println("-------------BLAAAAAAAAAAAAAAAACK STOOOOOOOOOOOOOOOOOOOOOORK---------------");
            }
            if (specie != -1) {
                for (int t = 0; t < MAXIMUM_SAVED_MATRICES; t++) {
                    if (savedBirdsTransitionMatrics[specie][t] == null) {
                        savedBirdsEmissionMatrices[specie][t] = currentBirdsEmissionMatrices[i];
                        savedBirdsTransitionMatrics[specie][t] = currentBirdsTransitionMatrices[i];
                        break;
                    }
                }
            }
        }
        //                Bird testBird = pState.getBird(0);
        //
        //                for (int i = 0; i < pState.getNumBirds(); i++) {
        //                    double[] initialStatePDVector = new double[pState.getNumBirds()];
        //                    initialStatePDVector[0] = 1;
        //                    double bestEmissionProb = 0;
        //                    for (int j = 0; j < currentBirdsTransitionMatrices.length; j++) {
        //                        double emissionProb = calculateEmissionSequenceProb(currentBirdsTransitionMatrices[i], currentBirdsEmissionMatrices[i], initialStatePDVector, birdDiedOrIsAtStep[0] + 1, testBird);
        //                        if (emissionProb > bestEmissionProb) {
        //                            bestEmissionProb = emissionProb;
        //                        }
        //                        if (j != currentBirdsTransitionMatrices.length - 1) {
        //                            initialStatePDVector[j] = 0;
        //                            initialStatePDVector[j + 1] = 1;
        //                        }
        //                    }
        //                    System.err.println(" (" + i + ") " + bestEmissionProb);
        //                }
        for (int i = 0; i < pSpecies.length; i++) {
            System.err.print(" (" + i + ") " + pSpecies[i]);
        }
        System.err.println();


        for (int g = 0; g < currentGuessedSpecies.length; g++) {
            if (currentGuessedSpecies[g] != Constants.SPECIES_UNKNOWN) {
                numberOfGuesses++;
                if (pSpecies[g] == currentGuessedSpecies[g]) {
                    numberOfCorrectGuesses++;
                }
            }
        }
        System.err.println("Guesses made: " + numberOfGuesses);
        System.err.println("Correct guesses made: " + numberOfCorrectGuesses);
        System.err.println("Percentage correct guessses: " + (double) numberOfCorrectGuesses / numberOfGuesses);
    }

    public static final Action cDontShoot = new Action(-1, -1);

    public double updateMatrices(double[][] transitionMatrix, double[][] emissionMatrix, double[] initialStatePDVector, Bird bird, int emissionSize, int maximumIterations) {
        double lastMaxProb = -1;
        double currentMaxProb = 0;
        double[][] forwardProb;
        double[][] backwardsProb;
        double[][] forwardBackwardAverageProb;
        double[][][] forwardBackwardAverageProb2;

        for (int iteration = 0; iteration < maximumIterations; iteration++) {
            forwardProb = new double[transitionMatrix.length][emissionSize];
            double[][] forwardProbWOSmooth = new double[transitionMatrix.length][emissionSize];
            double tempTotalProb = 0;
            for (int i = 0; i < transitionMatrix.length; i++) {
                double currentProb = initialStatePDVector[i] * emissionMatrix[i][bird.getObservation(0)];
                forwardProb[i][0] = currentProb;
                forwardProbWOSmooth[i][0] = currentProb;
                tempTotalProb += currentProb;
            }
            for (int i = 0; i < transitionMatrix.length; i++) {
                forwardProb[i][0] = forwardProb[i][0] / tempTotalProb;
            }

            for (int m = 1; m < emissionSize; m++) {
                int e = bird.getObservation(m);
                double tempTotalProb2 = 0;
                for (int i = 0; i < transitionMatrix.length; i++) {
                    double probForTransitionI = 0;
                    double probForTransitionIWOSmooth = 0;
                    for (int j = 0; j < transitionMatrix.length; j++) {
                        probForTransitionI += forwardProb[j][m - 1] * transitionMatrix[j][i];
                        probForTransitionIWOSmooth += forwardProbWOSmooth[j][m - 1] * transitionMatrix[j][i];
                    }
                    double probForTransitionIEmissionE = probForTransitionI * emissionMatrix[i][e];
                    double probForTransitionIEmissionEWOSmooth = probForTransitionIWOSmooth * emissionMatrix[i][e];
                    forwardProb[i][m] = probForTransitionIEmissionE;
                    forwardProbWOSmooth[i][m] = probForTransitionIEmissionEWOSmooth;
                    tempTotalProb2 += probForTransitionIEmissionE;
                }
                for (int i = 0; i < transitionMatrix.length; i++) {
                    forwardProb[i][m] = forwardProb[i][m] / tempTotalProb2;
                }
            }
            currentMaxProb = 0;
            for (int i = 0; i < transitionMatrix.length; i++) {
                if (forwardProbWOSmooth[i][forwardProbWOSmooth.length - 1] > currentMaxProb) {
                    currentMaxProb = forwardProbWOSmooth[i][forwardProb[0].length - 1];
                }
            }

            if (lastMaxProb != -1) {
                if (currentMaxProb < lastMaxProb) {
                    break;
                }
            }
            lastMaxProb = currentMaxProb;


            backwardsProb = new double[transitionMatrix.length][emissionSize];
            for (int i = 0; i < transitionMatrix.length; i++) {
                backwardsProb[i][backwardsProb[0].length - 1] = 1;
            }

            for (int m = emissionSize - 1; m > 0; m--) {
                int e = bird.getObservation(m);
                double totalProb = 0;
                for (int i = 0; i < transitionMatrix.length; i++) {
                    double probForTransitionI = 0;
                    for (int j = 0; j < transitionMatrix.length; j++) {
                        probForTransitionI += transitionMatrix[i][j] * emissionMatrix[j][e] * backwardsProb[j][m];
                    }
                    totalProb += probForTransitionI;
                    backwardsProb[i][m - 1] = probForTransitionI;
                }
                for (int i = 0; i < transitionMatrix.length; i++) {
                    backwardsProb[i][m - 1] = backwardsProb[i][m - 1] / totalProb;
                }
            }

            forwardBackwardAverageProb = new double[transitionMatrix.length][emissionSize];
            for (int m = 0; m < emissionSize; m++) {
                for (int i = 0; i < transitionMatrix.length; i++) {
                    double multiplication1 = forwardProb[i][m] * backwardsProb[i][m];
                    double multiplication2 = 0;
                    for (int j = 0; j < transitionMatrix.length; j++) {
                        multiplication2 += forwardProb[j][m] * backwardsProb[j][m];
                    }
                    forwardBackwardAverageProb[i][m] = multiplication1 / multiplication2;
                }
            }

            forwardBackwardAverageProb2 = new double[transitionMatrix.length][transitionMatrix.length][emissionSize];
            for (int m = 0; m < emissionSize - 1; m++) {
                int e = bird.getObservation(m + 1);
                double multiplication2 = 0;
                for (int i2 = 0; i2 < transitionMatrix.length; i2++) {
                    for (int j2 = 0; j2 < transitionMatrix.length; j2++) {
                        multiplication2 += forwardProb[i2][m] * transitionMatrix[i2][j2] * backwardsProb[j2][m + 1] * emissionMatrix[j2][e];
                    }
                }
                for (int i = 0; i < transitionMatrix.length; i++) {
                    double multiplication1 = 0;
                    for (int j = 0; j < transitionMatrix.length; j++) {
                        multiplication1 = forwardProb[i][m] * transitionMatrix[i][j] * backwardsProb[j][m + 1] * emissionMatrix[j][e];
                        forwardBackwardAverageProb2[i][j][m] = multiplication1 / multiplication2;
                    }
                }
            }


            for (int i = 0; i < transitionMatrix.length; i++) {
                initialStatePDVector[i] = forwardBackwardAverageProb[i][0];
            }

            for (int i = 0; i < transitionMatrix.length; i++) {
                for (int j = 0; j < transitionMatrix.length; j++) {
                    double multiplication1 = 0;
                    double multiplication2 = 0;
                    for (int m = 0; m < emissionSize - 1; m++) {
                        multiplication1 += forwardBackwardAverageProb2[i][j][m];
                        multiplication2 += forwardBackwardAverageProb[i][m];
                    }
                    transitionMatrix[i][j] = multiplication1 / multiplication2;
                }
            }

            for (int i = 0; i < transitionMatrix.length; i++) {
                for (int e = 0; e < emissionMatrix[0].length; e++) {
                    double multiplication1 = 0;
                    double multiplication2 = 0;
                    for (int m = 0; m < emissionSize - 1; m++) {
                        if (bird.getObservation(m) == e) {
                            multiplication1 += forwardBackwardAverageProb[i][m];
                        }
                        multiplication2 += forwardBackwardAverageProb[i][m];
                    }
                    emissionMatrix[i][e] = multiplication1 / multiplication2;
                }
            }
        }
        return currentMaxProb;
    }

    private Action getBestShootingAction(GameState pState) {
        Action returnAction = cDontShoot;
        int stepType = 0;
        int birdNumber = 0;
        double percentage = 0;
        for (int birdI = 0; birdI < pState.getNumBirds(); birdI++) {
            Bird bird = pState.getBird(birdI);

            if (bird.isAlive()) {
                double[] nextEmissionVector = getNextEmission(currentBirdsTransitionMatrices[birdI], currentBirdsEmissionMatrices[birdI], bird.getLastObservation());
                double totalProb = 0;
                for (int m = 0; m < nextEmissionVector.length; m++) {
                    double probOfStep = nextEmissionVector[m];
                    totalProb += probOfStep;
                    if (probOfStep > percentage) {
                        if (!isBirdUnknownOrBlack(birdI)) {

                            stepType = m;
                            birdNumber = birdI;
                            percentage = probOfStep;
                        }
                    }
                }


                if (totalProb < 0.99 || Double.isNaN(totalProb)) {
                    //                    System.err.println("TotalProb of bird: " + birdI + "is: " + totalProb);
                    randomizeAndNormalizeMatrix(currentBirdsTransitionMatrices[birdI]);
                    randomizeAndNormalizeMatrix(currentBirdsEmissionMatrices[birdI]);
                    randomizeAndNormalizeVector(currentBirdsInitialStatePDVectors[birdI]);
                    continue;
                }
                double emissionProb = calculateEmissionSequenceProb(currentBirdsTransitionMatrices[birdI], currentBirdsEmissionMatrices[birdI], currentBirdsInitialStatePDVectors[birdI], currentBirdStep + 1, bird);
                if (emissionProb < MINIMUM_SEQUENCEPROB_FOR_BIRDS_OWN_MATRIX) {
                    //                    System.err.println("EmissionProb: " + emissionProb);
                    randomizeAndNormalizeMatrix(currentBirdsTransitionMatrices[birdI]);
                    randomizeAndNormalizeMatrix(currentBirdsEmissionMatrices[birdI]);
                    randomizeAndNormalizeVector(currentBirdsInitialStatePDVectors[birdI]);
                    continue;
                }


                //                }
            }
        }
//        System.err.println("Percentage of hit: " + percentage);
        if (percentage > MINIMUM_HIT_PROB) {
            returnAction = new Action(birdNumber, stepType);
        }
        return returnAction;
    }

    private boolean isBirdUnknownOrBlack(int birdI) {
        if (currentGuessedSpecies[birdI] == Constants.SPECIES_UNKNOWN || currentGuessedSpecies[birdI] == Constants.SPECIES_BLACK_STORK) {
            return true;
        }
        return false;
    }

    private double[] getNextEmission(double[][] transitionMatrix, double[][] emissionMatrix, int lastEmission) {
        double[] lastStateVector = getLastStateVector(emissionMatrix, lastEmission);


        double[] nextStateVector = new double[transitionMatrix.length];
        for (int transitionRow = 0; transitionRow < transitionMatrix.length; transitionRow++) {
            for (int transitionCol = 0; transitionCol < transitionMatrix.length; transitionCol++) {
                nextStateVector[transitionCol] += lastStateVector[transitionRow] * transitionMatrix[transitionRow][transitionCol];
            }
        }

        double[] returnEmissionVector = new double[emissionMatrix[0].length];
        for (int emissionIndex = 0; emissionIndex < returnEmissionVector.length; emissionIndex++) {
            for (int stateProbIndex = 0; stateProbIndex < nextStateVector.length; stateProbIndex++) {
                returnEmissionVector[emissionIndex] += nextStateVector[stateProbIndex] * emissionMatrix[stateProbIndex][emissionIndex];
            }
        }
        return returnEmissionVector;
    }

    private double[] getLastStateVector(double[][] emissionMatrix, int lastEmission) {
        double[] returnVector = new double[emissionMatrix.length];
        double totalValue = 0;
        for (int i = 0; i < emissionMatrix.length; i++) {
            double value = emissionMatrix[i][lastEmission];
            returnVector[i] = value;
            totalValue += value;
        }
        normalizeVector(returnVector, totalValue);
        return returnVector;
    }

    private double calculateEmissionSequenceProb(int birdIndex, Bird bird) {
        return calculateEmissionSequenceProb(currentBirdsTransitionMatrices[birdIndex], currentBirdsEmissionMatrices[birdIndex], currentBirdsInitialStatePDVectors[birdIndex], birdDiedOrIsAtStep[birdIndex] + 1, bird);
    }

    private double calculateEmissionSequenceProb(double[][] transitionMatrix, double[][] emissionMatrix, double[] initialStatePDVector, int emissionSize, Bird bird) {

        double[][] probOfEmissionSequenceAtState = new double[transitionMatrix.length][emissionSize];


        for (int i = 0; i < transitionMatrix.length; i++) {

            double currentProb = initialStatePDVector[i] * emissionMatrix[i][bird.getObservation(0)];
            probOfEmissionSequenceAtState[i][0] = currentProb * Constants.COUNT_SPECIES;
        }


        for (int m = 1; m < emissionSize; m++) {
            int e = bird.getObservation(m);
            for (int i = 0; i < transitionMatrix.length; i++) {


                double probForTransitionIWOSmooth = 0;

                for (int j = 0; j < transitionMatrix.length; j++) {
                    probForTransitionIWOSmooth += probOfEmissionSequenceAtState[j][m - 1] * transitionMatrix[j][i];
                }


                if (e == -1) {
                    throw new RuntimeException();
                }
                double probForTransitionIEmissionEWOSmooth = probForTransitionIWOSmooth * emissionMatrix[i][e] * Constants.COUNT_SPECIES;
                probOfEmissionSequenceAtState[i][m] = probForTransitionIEmissionEWOSmooth;
            }
        }
        double emissionSequenceProb = 0;
        for (int i = 0; i < transitionMatrix.length; i++) {
            if (probOfEmissionSequenceAtState[i][probOfEmissionSequenceAtState[0].length - 1] > emissionSequenceProb) {

                emissionSequenceProb = probOfEmissionSequenceAtState[i][probOfEmissionSequenceAtState[0].length - 1];
            }
        }
        return emissionSequenceProb;
    }
}
