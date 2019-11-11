package FiducialAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.measure.Calibration;

import java.awt.*;
import java.util.ArrayList;

public class trackList {
    int[] ids;
    int[] lengths;
    int filterValue;
    int[] filteredIDs, filteredlengths;

    boolean Filtered = false;
    double[] meanX, meanY, filteredmeanXs, filteredmeanYs;
    int[][] trackLocIDs;

    public trackList(localisationList l) {
        ArrayList<Integer> trackIDs = l.getTrackIDs();
        ArrayList<Integer> trackLengths = l.getTrackLengths();
        ArrayList<Double> meanXs = l.getMeanX();
        ArrayList<Double> meanYs = l.getMeanY();
        ArrayList<Integer[]> trackLocIdList = l.getTrackLocIdsList();

        int size = trackIDs.size();

        ids = new int[size];
        lengths = new int[size];
        meanX = new double[size];
        meanY = new double[size];
        trackLocIDs = new int[size][];


        for (int i = 0; i < trackIDs.size(); i++) {
            ids[i] = trackIDs.get(i);
            lengths[i] = trackLengths.get(i);
            meanX[i] = meanXs.get(i);
            meanY[i] = meanYs.get(i);
            Integer[] tempIDs = trackLocIdList.get(i);
            trackLocIDs[i] = new int[tempIDs.length];

            for (int j = 0; j < tempIDs.length; j++) {
                trackLocIDs[i][j] = tempIDs[j];
            }

        }


    }


    public int[] getIds() {
        return this.ids;

    }

    public int[] getLengths() {
        return this.lengths;
    }

    public void FilterIds(int min) {
        int count = 0;

        for (int a:lengths) {
            if (a >= min) {
                count++;
            }
        }

        filteredIDs = new int[count];
        filteredlengths = new int[count];
        filteredmeanXs = new double[count];
        filteredmeanYs = new double[count];

        count = 0;

        for (int i = 0; i < lengths.length; i++) {
            if (lengths[i] >= min) {
                filteredIDs[count] = ids[i];
                filteredlengths[count] = lengths[i];
                filteredmeanXs[count] = meanX[i];
                filteredmeanYs[count] = meanY[i];

                count++;
            }
        }


        filterValue = min;
        Filtered = true;

    }

    public int getFilterValue() {
        return filterValue;
    }

    public int[] getFilteredIDs() {
        return filteredIDs;
    }

    public int[] getFilteredlengths() {
        return filteredlengths;
    }

    public double[] getFilteredlengths_double() {
        double[] d = new double[filteredlengths.length];

        for(int i=0;i<filteredlengths.length;i++){
            d[i] = filteredlengths[i];
        }

        return d;
    }

    public double[] getFilteredmeanXs() {
        return filteredmeanXs;
    }

    public double[] getFilteredmeanYs() {
        return filteredmeanYs;
    }

    public boolean isFiltered() {
        return Filtered;
    }

    public int[] getTrackLocs(int id) {

        boolean test = false;
        int[] out = null;
        for (int i = 0; i < ids.length; i++) {
            if (ids[i] == id) {
                out = trackLocIDs[i];
                test = true;
            }
        }

        if (!test) {
            out = null;
        }

        return out;

    }

    public ImagePlus addBeads(ImagePlus imp) {

        PointRoi pta = null;
        Calibration c = imp.getCalibration();

        for (int i = 0; i < filteredIDs.length; i++) {
            double x = filteredmeanXs[i] / c.pixelWidth;
            double y = filteredmeanYs[i] / c.pixelHeight;

            if (i == 0) {
                pta = new PointRoi(x, y);
            } else {
                pta.addPoint(x, y);
            }

        }

        imp.setRoi(pta);

        return imp;


    }
}
