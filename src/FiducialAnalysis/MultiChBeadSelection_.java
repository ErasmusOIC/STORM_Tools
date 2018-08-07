package FiducialAnalysis;

import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;

public class MultiChBeadSelection_ implements PlugIn {

    public void run(String arg){

        String primary,secondary;
        double[][] primary_xy = new double[3][];
        double[][] secondary_xy = new double[3][];

        //load two files

        OpenDialog od  =  new OpenDialog("Select Primary Data");
        primary = od.getPath();
        od = new OpenDialog("Select Secondary Data");
        secondary = od.getPath();


        Opener.openResultsTable(primary);

        ResultsTable primary_rt = Analyzer.getResultsTable();

        Opener.openResultsTable(secondary);

        ResultsTable secondary_rt = Analyzer.getResultsTable();


        //analyse two datasets

        localisationList primary_loc = new localisationList(primary_rt);
        localisationList secondary_loc = new localisationList(secondary_rt);

        primary_loc.assignTracks(50);
        secondary_loc.assignTracks(50);

        trackList primary_tracks = primary_loc.getTracks();
        trackList secondary_tracks = secondary_loc.getTracks();

        primary_tracks.FilterIds((int)(primary_loc.getMaxF()*0.3));
        secondary_tracks.FilterIds((int)(secondary_loc.getMaxF()*0.3));

        primary_xy[0] = primary_tracks.getFilteredmeanXs();
        primary_xy[1] = primary_tracks.getFilteredmeanYs();
        primary_xy[2] = primary_tracks.getFilteredlengths_double();

        secondary_xy[0] = secondary_tracks.getFilteredmeanXs();
        secondary_xy[1] = secondary_tracks.getFilteredmeanYs();
        secondary_xy[2] = secondary_tracks.getFilteredlengths_double();

        ResultsTable rt = new ResultsTable();
        Analyzer.setResultsTable(rt);

        int max;

       if(primary_xy[0].length>=secondary_xy[0].length){
           max = primary_xy[0].length;
        }else{
           max = secondary_xy[0].length;
       }

        for(int i=0;i<max;i++){
            rt.incrementCounter();
            if(i<primary_xy[0].length) {
                rt.addValue("X1", primary_xy[0][i]);
                rt.addValue("Y1", primary_xy[1][i]);
                rt.addValue("L1", primary_xy[2][i]);
            }
            if(i<secondary_xy[0].length) {
                rt.addValue("X2", secondary_xy[0][i]);
                rt.addValue("Y2", secondary_xy[1][i]);
                rt.addValue("L2", secondary_xy[2][i]);
            }
        }

        rt.show("Results");






    }
}
