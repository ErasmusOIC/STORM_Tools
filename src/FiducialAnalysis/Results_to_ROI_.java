package FiducialAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.PointRoi;
import ij.gui.Overlay;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import java.awt.*;

public class Results_to_ROI_ implements PlugIn {
    public void run(String arg) {

        ResultsTable rt = Analyzer.getResultsTable();
        if(rt==null){
            IJ.showMessage("No Results table is open");
            return;
        }
        double [] x1, y1, x2, y2, x3, y3;

            x1 = rt.getColumnAsDoubles(0);
            y1 = rt.getColumnAsDoubles(1);
            x2 = rt.getColumnAsDoubles(2);
            y2 = rt.getColumnAsDoubles(3);
            x3 = rt.getColumnAsDoubles(4);
            y3 = rt.getColumnAsDoubles(5);

        ImagePlus imp = WindowManager.getCurrentImage();
        Overlay ov = new Overlay();
        Calibration c = imp.getCalibration();
        String units;
        units = c.getUnit();
        //color 1
        PointRoi pta = null;
        String unit_img = "µm";
        if(units.equals(unit_img)){
            for (int i = 0; i < x1.length; i++) {
                double x = x1[i] / 1000.0 / c.pixelWidth;
                double y = y1[i] / 1000.0 / c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.red);
                ov.add(pta);

            }


            //color 2
            for (int i = 0; i < x2.length; i++) {
                double x = x2[i] / 1000.0/c.pixelWidth;
                double y = y2[i] / 1000.0/c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.blue);
                ov.add(pta);

            }
            //color 2
            for (int i = 0; i < x3.length; i++) {
                double x = x3[i] / 1000.0/ c.pixelWidth;
                double y = y3[i] / 1000.0/ c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.green);
                ov.add(pta);

            }
        } else {
            for (int i = 0; i < x1.length; i++) {
                double x = x1[i] / c.pixelWidth;
                double y = y1[i] / c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.red);
                ov.add(pta);

            }



            //color 2
            for (int i = 0; i < x2.length; i++) {
                double x = x2[i] / c.pixelWidth;
                double y = y2[i] / c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.blue);
                ov.add(pta);

            }
            //color 2
            for (int i = 0; i < x3.length; i++) {
                double x = x3[i] / c.pixelWidth;
                double y = y3[i] / c.pixelHeight;

                pta = new PointRoi(x, y);
                pta.setStrokeColor(Color.green);
                ov.add(pta);

            }
        }

        imp.setOverlay(ov);

        imp.show();
    }

}
