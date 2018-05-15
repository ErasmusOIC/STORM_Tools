package FiducialAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

public class fiducial_Analysis implements PlugIn {

    public void run(String arg){
        ResultsTable rt = Analyzer.getResultsTable();
        int gap = 50;

        if(rt==null){
            IJ.showMessage("No Results table is open");
            return;
        }

        localisationList locList = new localisationList(rt);
        IJ.log("localisation loaded");

        locList.assignTracks(gap);
        IJ.log("tracks assigned");

        int frames = locList.getMaxF();

        locList.getTracks().FilterIds((int)(frames*0.5));

        trackList track = locList.getTracks();
        int[] ids = track.getFilteredIDs();
        int[] lengths = track.getFilteredlengths();
        double[] meanX = track.getFilteredmeanXs();
        double[] meanY = track.getFilteredmeanYs();


        for(int i=0;i<ids.length;i++){
            IJ.log(""+ids[i]+" "+lengths[i]+" "+(int) meanX[i]+" "+(int) meanY[i]);
        }

        ImagePlus imp = WindowManager.getCurrentImage();

        if(imp==null){
            double xSize = locList.getMaxX();
            double ySize = locList.getMaxY();

            Calibration c = new Calibration();
            c.pixelWidth = 5.0;
            c.pixelHeight = 5.0;
            c.setUnit("nm");

            ImageProcessor ip = new ByteProcessor((int) c.getRawX(xSize), (int)c.getRawY(ySize));
            imp = new ImagePlus("temp",ip);
            imp.setCalibration(c);

        }


        imp = track.addBeads(imp);

        int[] trackIDs = track.getTrackLocs(track.getFilteredIDs()[0]);
        localisation[] trackNumbers = locList.getLocList();

        for(int i=0;i<trackIDs.length;i++){
            IJ.log(trackIDs[i]+" "+trackNumbers[trackIDs[i]].getTrackID() +" "+ trackNumbers[trackIDs[i]].getX() );

        }






        imp.show();





        rt = locList.makeResultTable();
        rt.show("Results");

        IJ.log("finished");



    }
}
