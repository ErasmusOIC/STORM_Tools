package FiducialAnalysis;

public class localisation {
    double x;
    double y;
    double z;
    int timeFrame;
    double precision;
    int trackID=0;
    int trackLength=0;
    int numFrames=0;

    public localisation(double x_,double y_, double z_, int timeFrame_, double precision_, int numFrames_){
        x = x_;
        y = y_;
        z = z_;
        timeFrame = timeFrame_;
        precision = precision_;
        numFrames = numFrames_;
    }

    public double getPrecision() {
        return precision;
    }

    public int getTimeFrame() {
        return timeFrame;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    public int getNumFrames(){ return numFrames;}

    public int getTrackID() {
        return trackID;
    }

    public int getTrackLength() {
        return trackLength;
    }

    public void setPrecision(double precision) {
        this.precision = precision;
    }

    public void setTimeFrame(int timeFrame) {
        this.timeFrame = timeFrame;
    }

    public void setX(double x) {
        this.x = x;
    }

    public void setY(double y) {
        this.y = y;
    }

    public void setZ(double z) {
        this.z = z;
    }

    public void setTrackID(int trackID) {
        this.trackID = trackID;
    }

    public void setTrackLength(int trackLength) {
        this.trackLength = trackLength;
    }
}
