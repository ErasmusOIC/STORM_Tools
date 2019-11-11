package FiducialAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.*;
import java.util.concurrent.ExecutionException;

public class localisationList {

    localisation[] locList;
    trackList tracks;
    int maxX,maxY,maxF, maxTr;
    boolean tracksAssigned=false;

    ArrayList<Integer> trackIDs = new ArrayList<>();
    ArrayList<Integer> trackLengths = new ArrayList<>();
    ArrayList<Double> meanX = new ArrayList<>();
    ArrayList<Double> meanY = new ArrayList<>();
    ArrayList<Integer[]> trackLocIdsList = new ArrayList<>();


    public localisationList(localisation[] l){
        locList = l;
        this.calcMaxX();
        this.calcMaxY();
        this.calcMaxF();

    }


    public localisationList(localisationList l){

    }

    public localisationList(localisationList l, int[] ids){

    }



    public localisationList(ResultsTable rt){

        String[] headers_temp = rt.getHeadings();
        String[] headers = new String[headers_temp.length+1];
        for(int i=0;i<headers_temp.length;i++){
            headers[i] = headers_temp[i];
        }
        headers[headers_temp.length] = "none";

        GenericDialog gd = new GenericDialog("select column names");
        gd.addChoice("select_x", headers,headers[0]);
        gd.addChoice("select_y", headers,headers[0]);
        gd.addChoice("select_z *", headers,headers[0]);
        gd.addChoice("select_Frame Number", headers,headers[0]);
        gd.addChoice("select_Precision", headers,headers[0]);
        gd.addMessage("* select none for 2D localisations");
        gd.showDialog();

        if(gd.wasCanceled()){
            return;
        }

        String xlabel = gd.getNextChoice();
        String ylabel = gd.getNextChoice();
        String zlabel = gd.getNextChoice();
        String framelabel = gd.getNextChoice();
        String preclabel = gd.getNextChoice();

        locList = new localisation[rt.size()];

        for(int i=0;i<rt.size();i++){

             double x = rt.getValue(xlabel,i);
             double y = rt.getValue(ylabel,i);
             double z;

             if(!zlabel.equals("none")) {
                 z = rt.getValue(zlabel, i);
             }else{
                 z = 0.0;
             }
             int timeFrame = (int) rt.getValue(framelabel,i);
             double precision = rt.getValue(preclabel,i);
             locList[i] = new localisation(x,y,z, timeFrame, precision);



             IJ.showProgress(i ,rt.size());

        }

        this.calcMaxX();
        this.calcMaxY();
        this.calcMaxF();


    }

    public void assignTracks(int timeGap){


        int[] trackIDList = new int[locList.length];
        int[] trackCount = new int[locList.length];

        for(int i=0;i<trackIDList.length-1;i++){
            int counter = 0;
            Deque<Integer> stack = new ArrayDeque<>();

            if(trackIDList[i]==0){

                stack.addFirst(i);
                trackIDList[i] = i+1;


                int firstFrame = locList[i].getTimeFrame();
                double X_total = locList[i].getX();
                double Y_total = locList[i].getY();
                int firstID = i;
                int secondID = i+1;
                boolean loop;




                do{
                    loop = true;
                    int secondFrame = locList[secondID].getTimeFrame();

                    if(secondFrame-firstFrame > timeGap){
                        loop = false;
                    }
                    if(secondID>=locList.length-1){
                        loop=false;
                    }

                    if(overlap(locList[firstID],locList[secondID]) && loop && trackIDList[secondID]==0 && firstFrame!=secondFrame){
                        X_total = X_total+locList[secondID].getX();
                        Y_total = Y_total+locList[secondID].getY();
                        trackIDList[secondID] = i+1;
                        stack.addFirst(secondID);
                        firstID = secondID;
                        firstFrame = secondFrame;
                        counter++;
                        trackCount[secondID] = counter;

                    }

                    secondID++;


                }while(loop);

                if(counter==0){
                    stack.pop();
                }



                if(counter>0) {
                    Integer[] trackLocIds = new Integer[stack.size()];

                    //System.out.print(stack.size()+" ");

                    int stacksize = stack.size();

                    for(int j=0;j<stacksize;j++){
                        trackLocIds[j] = stack.pop();

                    }

                    //System.out.print(trackLocIds.length+" ");

                    //System.out.print(stack.size()+" ");


                    trackLocIdsList.add(trackLocIds);
                    trackIDs.add(i + 1);
                    trackLengths.add(counter);
                    meanX.add((X_total / (double) counter));
                    meanY.add((Y_total / (double) counter));

                    //System.out.println(trackLocIdsList.size()+" "+meanX.size());
                }

            }

            IJ.showProgress(i,trackIDList.length-1);
        }

        for(int i=0;i<locList.length;i++){
            locList[i].setTrackID(trackIDList[i]);
            locList[i].setTrackLength(trackCount[i]);
        }

        this.calcMaxTr();

        tracks = new trackList(this);
        tracksAssigned = true;

        trackIDs = null;
        trackLengths = null;
        meanX = null;
        meanY = null;
        trackLocIdsList = null;



    }

    private boolean overlap(localisation a, localisation b){

        double factor = 1.0;

        double x1 = a.getX();
        double y1 = a.getY();
        double p1 = a.getPrecision()*factor;

        double x2 = b.getX();
        double y2 = b.getY();
        double p2 = b.getPrecision()*factor;

        double left  = Math.max(x1-p1, x2-p2);
        double rigth = Math.min(x1+p1,x2+p2);
        double top = Math.max(y1-p1,y2-p2);
        double bottom = Math.min(y1+p1,y2+p2);

        boolean overlap;

        if(rigth > left && bottom > top){
            overlap=true;
        }else{
            overlap=false;
        }

        return overlap;


    }

    public void showLog(){
        for(int i=0;i<locList.length;i++){

            double x = this.locList[i].getX();
            double y = this.locList[i].getY();
            double z = this.locList[i].getZ();
            double p = this.locList[i].getPrecision();
            double f = this.locList[i].getTimeFrame();
            int trackID = this.locList[i].getTrackID();
            int trackCount = this.locList[i].getTrackLength();

            IJ.log(" "+(int)x+" "+(int)y+" "+(int)z+" "+(int)p+" "+(int)f+" ");
        }
    }

    public ResultsTable makeResultTable(){
        ResultsTable rt = new ResultsTable();

        for(int i=0;i<locList.length;i++){
            double x = this.locList[i].getX();
            double y = this.locList[i].getY();
            double z = this.locList[i].getZ();
            double p = this.locList[i].getPrecision();
            double f = this.locList[i].getTimeFrame();
            int trackID = this.locList[i].getTrackID();
            int trackCount = this.locList[i].getTrackLength();

            rt.incrementCounter();
            rt.addValue("X",x);
            rt.addValue("Y",y);
            rt.addValue("Z",z);
            rt.addValue("Precision",p);
            rt.addValue("Frame Number",f);
            rt.addValue("Track ID", trackID);
            rt.addValue("Track Length",trackCount);

        }

        return rt;
    }



    private void calcMaxX(){
        double max=0;
        for(int i=0;i<locList.length;i++){
            if(locList[i].getX()>max){
                max = locList[i].getX();
            }
        }

        maxX = (int)(max +1);


    }

    private void calcMaxY(){
        double max=0;
        for(int i=0;i<locList.length;i++){
            if(locList[i].getY()>max){
                max = locList[i].getY();
            }
        }

        maxY = (int)(max +1);

    }

    private void calcMaxF(){
        int max=0;
        for(int i=0;i<locList.length;i++){
            if(locList[i].getTimeFrame()>max){
                max = locList[i].getTimeFrame();
            }
        }

        maxF = max;

    }

    private void calcMaxTr(){
        int max=0;

        for(int i=0;i<locList.length;i++){
            if(locList[i].getTrackID()>max){
                max = locList[i].getTrackID();
            }
        }

        maxTr = max;

    }

    public int getMaxTr() {
        return this.maxTr;
    }

    public localisation[] getLocList(){
        return this.locList;
    }

    public trackList getTracks() {
        return this.tracks;
    }

    public ArrayList getTrackIDs(){

        return this.trackIDs;

    }

    public ArrayList getTrackLengths() {

        return this.trackLengths;
    }


    public ArrayList getMeanX() {
        return meanX;
    }

    public ArrayList getMeanY() {
        return meanY;
    }

    public int getMaxX() { return maxX; }

    public int getMaxY() { return maxY; }

    public boolean isTracksAssigned() {
        return tracksAssigned;
    }

    public int getMaxF() {
        return maxF;
    }

    public ArrayList<Integer[]> getTrackLocIdsList(){ return trackLocIdsList;}






}
