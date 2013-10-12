import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * @author mathiaslindblom
 */
public class Main {
    private static double[][] transisionMatrix;
    private static double[][] emissionMatrix;
    private static double[]   initialStatePDVector; //Says matrix in assignment but should/could be vector;
    private static int[]      emissionSequence;
    private static double totalProb = 0;
    private static int[]  mostLikelyPath;
    private static double probOfMostLikelyPath;

    private static double[][] forwardProb;
    private static double[][] backwardsProb;
    private static double     currentMaxProb;
    private static double     lastMaxProb;


    public static void main(String[] args) throws IOException {


        ArrayList<String[]> lines = new ArrayList<String[]>();

        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

        String line;
        while ((line = in.readLine()) != null) {
            if (line.contains(";")) {
                break;
            }

            lines.add(line.split(" "));
        }

        for (int n = 0; n < 4; n++) {
            String[] stringArray = lines.get(n);
            if (n == 3) {
                int size = Integer.parseInt(stringArray[0]);
                emissionSequence = new int[size];
                for (int i = 0; i < size; i++) {
                    emissionSequence[i] = Integer.parseInt(stringArray[1 + i]);
                }
                break;
            }
            int row = Integer.parseInt(stringArray[0]);
            int col = Integer.parseInt(stringArray[1]);
            double[][] matrix = new double[row][col];

            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    matrix[i][j] = Double.parseDouble(stringArray[col * i + j + 2]);
                }
            }
            if (n == 0) { transisionMatrix = matrix; }
            if (n == 1) { emissionMatrix = matrix; }
            if (n == 2) {
                if (matrix.length > 1) {
                    throw new RuntimeException();
                }
                initialStatePDVector = matrix[0];
            }
        }


        lastMaxProb = -1;

        for(int iteration = 0; iteration < 50; iteration++) {
            forwardProb = new double[transisionMatrix.length][emissionSequence.length];
            double[][] forwardProbWOSmooth = new double[transisionMatrix.length][emissionSequence.length];
            double tempTotalProb = 0;
            for (int i = 0; i < transisionMatrix.length; i++) {
                double currentProb = initialStatePDVector[i] * emissionMatrix[i][emissionSequence[0]];
                forwardProb[i][0] = currentProb;
                forwardProbWOSmooth[i][0] = currentProb;
                tempTotalProb += currentProb;
            }
            for (int i = 0; i < transisionMatrix.length; i++) {
                forwardProb[i][0] = forwardProb[i][0] / tempTotalProb;
            }


            mostLikelyPath = new int[emissionSequence.length];
            for (int m = 1; m < emissionSequence.length; m++) {
                int e = emissionSequence[m];
                double tempTotalProb2 = 0;
                for (int i = 0; i < transisionMatrix.length; i++) {
                    double probForTransitionI = 0;
                    double probForTransitionIWOSmooth = 0;
                    for (int j = 0; j < transisionMatrix.length; j++) {
                        probForTransitionI += forwardProb[j][m - 1] * transisionMatrix[j][i];
                        probForTransitionIWOSmooth += forwardProbWOSmooth[j][m-1] * transisionMatrix[j][i];
                    }
                    double probForTransitionIEmissionE = probForTransitionI * emissionMatrix[i][e];
                    double probForTransitionIEmissionEWOSmooth = probForTransitionIWOSmooth * emissionMatrix[i][e];
                    forwardProb[i][m] = probForTransitionIEmissionE;
                    forwardProbWOSmooth[i][m] = probForTransitionIEmissionEWOSmooth;
                    tempTotalProb2 +=probForTransitionIEmissionE;
                }
                for (int i = 0; i < transisionMatrix.length; i++) {
                    forwardProb[i][m] = forwardProb[i][m] / tempTotalProb2;
                }
            }
            currentMaxProb = 0;
            for (int i = 0; i < transisionMatrix.length; i++) {
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


            backwardsProb = new double[transisionMatrix.length][emissionSequence.length];
            for (int i = 0; i < transisionMatrix.length; i++) {
                backwardsProb[i][backwardsProb[0].length - 1] = 1;
            }

            for (int m = emissionSequence.length - 1; m > 0; m--) {
                int e = emissionSequence[m];
                double totalProb = 0;
                for (int i = 0; i < transisionMatrix.length; i++) {
                    double probForTransitionI = 0;
                    for (int j = 0; j < transisionMatrix.length; j++) {
                        probForTransitionI += transisionMatrix[i][j] * emissionMatrix[j][e] * backwardsProb[j][m];
                    }
                    totalProb += probForTransitionI;
                    backwardsProb[i][m - 1] = probForTransitionI;
                }
                for (int i = 0; i < transisionMatrix.length; i++) {
                    backwardsProb[i][m - 1] = backwardsProb[i][m - 1] / totalProb;
                }
            }

            double[][] forwardBackwardAverageProb = new double[transisionMatrix.length][emissionSequence.length];
            for (int m = 0; m < emissionSequence.length; m++) {
                for (int i = 0; i < transisionMatrix.length; i++) {
                    double multiplication1 = forwardProb[i][m] * backwardsProb[i][m];
                    double multiplication2 = 0;
                    for (int j = 0; j < transisionMatrix.length; j++) {
                        multiplication2 += forwardProb[j][m] * backwardsProb[j][m];
                    }
                    forwardBackwardAverageProb[i][m] = multiplication1 / multiplication2;
                }
            }

            double[][][] forwardBackwardAverageProb2 = new double[transisionMatrix.length][transisionMatrix.length][emissionSequence.length];
            for (int m = 0; m < emissionSequence.length - 1; m++) {
                int e = emissionSequence[m + 1];
                double multiplication2 = 0;
                for (int i2 = 0; i2 < transisionMatrix.length; i2++) {
                    for (int j2 = 0; j2 < transisionMatrix.length; j2++) {
                        multiplication2 += forwardProb[i2][m] * transisionMatrix[i2][j2] * backwardsProb[j2][m + 1] * emissionMatrix[j2][e];
                    }
                }
                for (int i = 0; i < transisionMatrix.length; i++) {
                    double multiplication1 = 0;
                    for (int j = 0; j < transisionMatrix.length; j++) {
                        multiplication1 = forwardProb[i][m] * transisionMatrix[i][j] * backwardsProb[j][m + 1] * emissionMatrix[j][e];
                        forwardBackwardAverageProb2[i][j][m] = multiplication1 / multiplication2;
                    }
                }
            }


            for (int i = 0; i < transisionMatrix.length; i++) {
                initialStatePDVector[i] = forwardBackwardAverageProb[i][0];
            }

            for (int i = 0; i < transisionMatrix.length; i++) {
                for (int j = 0; j < transisionMatrix.length; j++) {
                    double multiplication1 = 0;
                    double multiplication2 = 0;
                    for (int m = 0; m < emissionSequence.length - 1; m++) {
                        multiplication1 += forwardBackwardAverageProb2[i][j][m];
                        multiplication2 += forwardBackwardAverageProb[i][m];
                    }
                    transisionMatrix[i][j] = multiplication1 / multiplication2;
                }
            }

            for (int i = 0; i < transisionMatrix.length; i++) {
                for (int e = 0; e < emissionMatrix[0].length; e++) {
                    double multiplication1 = 0;
                    double multiplication2 = 0;
                    for (int m = 0; m < emissionSequence.length - 1; m++) {
                        if (emissionSequence[m] == e) {
                            multiplication1 += forwardBackwardAverageProb[i][m];
                        }
                        multiplication2 += forwardBackwardAverageProb[i][m];
                    }
                    emissionMatrix[i][e] = multiplication1 / multiplication2;
                }
            }

        }

        System.out.print(transisionMatrix.length + " " + transisionMatrix.length);
        for(int i = 0; i < transisionMatrix.length; i++){
            for(int j = 0; j < transisionMatrix.length; j++){
                System.out.print(" " + Math.round(transisionMatrix[i][j]*1000000.0)/1000000.0);
            }
        }
        System.out.println();

        System.out.print(emissionMatrix.length + " " + emissionMatrix[0].length);
        for(int i = 0; i < emissionMatrix.length; i++){
            for(int m = 0; m < emissionMatrix[0].length; m++){
                System.out.print(" " + Math.round(emissionMatrix[i][m]*1000000.0)/1000000.0);
            }
        }
//        StringBuilder sb = new StringBuilder();
//        for (int step : mostLikelyPath) {
//            sb.append(step + " ");
//        }
//        System.out.print(sb.toString());
//        backwardsProb = new double[transisionMatrix.length][emissionSequence.length];

    }

    private static void calculateSequencePropForward(double[][] transisionMatrix, double[][] emissionMatrix, int[] emissionSequence, int currentState, int sequenceIndex, double ackProb) {
        for (int i = 0; i < transisionMatrix.length; i++) {
            double currentProb = ackProb * transisionMatrix[currentState][i] * emissionMatrix[i][emissionSequence[sequenceIndex]];
            if (currentProb > 0) {
                forwardProb[i][sequenceIndex] += currentProb;
                int newSequence = sequenceIndex + 1;
                if (newSequence < emissionSequence.length) {
                    calculateSequencePropForward(transisionMatrix, emissionMatrix, emissionSequence, i, sequenceIndex + 1, currentProb);
                }
            }
        }
    }

    //    private static double caculateSequenceProbBackward(double[][] transisionMatrix, double[][] emissionMatrix, int[] emissionSequence, int currentState, int sequenceIndex, double ackProb){
    //
    //    }

    private static double[] computeNextStateProb(double[] initialStatePDVector, double[][] transisionMatrix, double[][] emissionMatrix) {
        double[] nextStateVector = new double[initialStatePDVector.length];
        for (int transitionRow = 0; transitionRow < transisionMatrix.length; transitionRow++) {
            for (int transitionCol = 0; transitionCol < transisionMatrix.length; transitionCol++) {
                nextStateVector[transitionCol] += initialStatePDVector[transitionRow] * transisionMatrix[transitionRow][transitionCol];
            }
        }
        return nextStateVector;
    }

    private static double[] computeEmissionProb(double[][] emissionMatrix, double[] stateVector) {
        double[] emissionPVector = new double[emissionMatrix[0].length];

        for (int emissionIndex = 0; emissionIndex < emissionPVector.length; emissionIndex++) {
            for (int stateProbIndex = 0; stateProbIndex < stateVector.length; stateProbIndex++) {
                emissionPVector[emissionIndex] += stateVector[stateProbIndex] * emissionMatrix[stateProbIndex][emissionIndex];
            }
        }
        return emissionPVector;
    }

    public static String test() {
        StringBuilder sb = new StringBuilder();
        sb.append("Transission Matrix");
        for (int row = 0; row < transisionMatrix.length; row++) {
            sb.append('\n');
            for (int col = 0; col < transisionMatrix.length; col++) {
                sb.append(transisionMatrix[row][col]);
                sb.append('\t');
            }
        }
        sb.append('\n');
        sb.append('\n');
        sb.append("Emission Matrix");
        for (int row = 0; row < emissionMatrix.length; row++) {
            sb.append('\n');
            for (int col = 0; col < emissionMatrix[0].length; col++) {
                sb.append(emissionMatrix[row][col]);
                sb.append('\t');
            }
        }
        sb.append('\n');
        sb.append('\n');
        return sb.toString();
    }
}
