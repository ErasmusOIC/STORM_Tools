///////////////////////////////////////////////////////////////
/// Script: "ROI_selection.ijm"                             ///
/// Authors: Lieke Koornneef & Johan Slotman                ///
/// Affiliation: Erasmus MC, Rotterdam, The Netherlands     ///
/// Contact: j.slotman@erasmusmc.nl  					    ///
/// License: LGPLv3                                         ///
/// Date: 13-01-2023                                        ///
///////////////////////////////////////////////////////////////

dir = getDirectory("Select directory");
files = getFileList(dir); 

tiffiles = newArray(1000);
count = 0;
for(i=0;i<files.length;i++){
	if(endsWith(files[i],".tif")){
    	tiffiles[count] = files[i];
        count++;
	}
}
tiffiles = Array.trim(tiffiles, count);

name = call("ij.Prefs.get","ST.Exp_name", "Exp1");
thr_r = call("ij.Prefs.get","ST.Thr_r", 500);
thr_d = call("ij.Prefs.get","ST.Thr_d", 300);

Dialog.create("Details of multicolor imaging");
Dialog.addString("Experiment name", name);
Dialog.addChoice("Multi-color image", tiffiles);
Dialog.addNumber("Threshold Rad51", thr_r); 
Dialog.addNumber("Threshold Dmc1", thr_d); 
Dialog.show();

name = Dialog.getString();
file = Dialog.getChoice();
thr_r = Dialog.getNumber();
thr_d = Dialog.getNumber();

call("ij.Prefs.set","ST.Exp_name", name);
call("ij.Prefs.set","ST.Thr_r", thr_r);
call("ij.Prefs.set","ST.Thr_d", thr_d);


//PREPARE IMAGE --------------------------------------------------------------------------------------------------------------

//Open image 
open(dir+file);
rename("storm image");
run("Set Scale...", "distance=1 known=5 pixel=1 unit=nm");
run("Duplicate...", "duplicate");
run("Split Channels");

//Select mask of chromosomal axes
selectWindow("C1-storm image-1");
run("Gaussian Blur...", "sigma=10");
setAutoThreshold("Huang dark");
run("Analyze Particles...", "size=2000-Infinity circularity=0.00-0.89 pixel show=Masks");
selectWindow("Mask of C1-storm image-1");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
roiManager("reset"); 
run("Invert LUT");

//Select RAD51 foci
selectWindow("C2-storm image-1");
setAutoThreshold("Huang dark");
run("Analyze Particles...", "size="+thr_r+"-Infinity pixel show=Masks");
selectWindow("Mask of C2-storm image-1");
run("Dilate");
run("Dilate");
run("Dilate");

//Select DMC1 foci 
selectWindow("C3-storm image-1");
setAutoThreshold("Huang dark");
run("Analyze Particles...", "size="+thr_d+"-Infinity pixel show=Masks");
selectWindow("Mask of C3-storm image-1");
run("Dilate");
run("Dilate");
run("Dilate");

//Combine masks 
selectWindow("Mask of C2-storm image-1");
run("Red");
run("8-bit");
selectWindow("Mask of C3-storm image-1");
run("Green");
run("8-bit");
run("Merge Channels...", "c1=[Mask of C2-storm image-1] c2=[Mask of C3-storm image-1] create keep");
rename("merge");
run("Duplicate...", "duplicate");

//Remove signal outside of chromosomal axes
imageCalculator("AND create", "Mask of C2-storm image-1","Mask of C3-storm image-1");
rename("overlapping_foci");
imageCalculator("add create", "overlapping_foci","Mask of C1-storm image-1");
setAutoThreshold("Huang dark");
run("Analyze Particles...", "size=210000-Infinity show=Masks exclude");
run("Create Selection");
roiManager("Add");

selectWindow("merge");
roiManager("Select", 0);
run("Clear Outside");
roiManager("reset"); 
close("Mask of C1-storm image-1");
close("Mask of Mask of C1-storm image-1");
close("Mask of C3-storm image-1");
close("Mask of C2-storm image-1");
close("C1-storm image-1");
close("C2-storm image-1");
close("C3-storm image-1");
close("overlapping_foci");
close("Result of overlapping_foci");
close("Mask of Result of overlapping_foci");
run("Select None");



//SELECTION ROIs--------------------------------------------------------------------------------------------------------------

selectWindow("merge")
run("Split Channels");
imageCalculator("Add create", "C1-merge","C2-merge");
rename("mask merge");

//Determine coordinates foci 
roiManager("reset");
run("Set Measurements...", "centroid redirect=None decimal=1");
close("Results"); 
run("Find Maxima...", "prominence=1 output=[Point Selection]");
roiManager("Add");
roiManager("Measure");
rn = nResults;
X = newArray(1000);
Y = newArray(1000);
newX = newArray(1000);
newY = newArray(1000);

for (i = 0; i < rn; i++) {
	X[i] = getResult("X", i);
	Y[i] = getResult("Y", i);	} 
X = Array.deleteValue(X, 0);
Y = Array.deleteValue(Y, 0);
close("Results");

//Combine adjacent foci
for (p = 0; p < rn ; p++) {
	threshold = 400; 		// in nm
	distance_new = 400;		// in nm
		for (i = 0; i < rn; i++) {
			if(p!=i){
			a = (Y[i])-(Y[p]) ;
			b = (X[i])-(X[p]) ;
			distance = sqrt((a*a) + (b*b));
			if(distance < distance_new) { 
				distance_new = distance ;
				wi = i;}}}			
	if(distance_new >= threshold) { 
		newX[p] = X[p];
		newY[p] = Y[p];}
		else {
			if(X[wi] < X[p]) { x_new = X[wi] + (abs(X[wi] - X[p])*0.5); } 
				else { x_new = X[p]+ (abs(X[wi]- X[p])*0.5); }
			if(Y[wi] < Y[p]) { y_new = Y[wi] + (abs(Y[wi] - Y[p])*0.5); } 
				else { y_new = Y[p] + (abs(Y[wi] - Y[p])*0.5); }
			newX[p] = x_new;
			newY[p] = y_new; }}

newX = Array.deleteValue(newX, 0);
newY = Array.deleteValue(newY, 0);

//Remove duplicates
for (p = 0; p < newX.length; p++) {
	p1 = (newX.length-1)-p;

	for (i = 0; i < newX.length; i++) {
		i1 = (newX.length-1)-i;
			
		if((p1!=i1) && (newX[i1]- newX[p1]) == 0){ 
			newX = Array.deleteIndex(newX, i1);
			newY = Array.deleteIndex(newY, i1);
			break; 
		}}}	
roiManager("reset");

//Create ROIs
d = 150; // in pixels 
for (i = 0; i < newX.length; i++) {
	xn = newX[i]-(0.5*d*5);
	yn = newY[i]-(0.5*d*5);
	makeOval(xn/5,yn/5,d,d);
	roiManager("add");}

//Change ROI numbers
for(i=0; i<roiManager("Count");i++) {
	roiManager("Select",i);
	roiManager("Rename",i+1);
}

//Manually adaptation
selectWindow("storm image");
roiManager("Show All");
roiManager("Select", Array.getSequence(roiManager("count")));
roiManager("Save", dir+name+"_allRois.zip");
waitForUser("Remove or add ROIs manually using -Manual ROI selection-");
numberofrois = roiManager("count");

//Change ROI numbers
roiManager("deselect");
for(i=0; i<roiManager("Count");i++) {
	roiManager("Select",i);
	roiManager("Rename",i+1);
}

//Save ROIs 
roiManager("Select", Array.getSequence(roiManager("count")));
roiManager("Save", dir+name+"_RoiSet.zip");

//Translate ROIs
open(dir+name+"_Roi_characteristics.csv");
X_roi = getResult("BX",0);
Y_roi = getResult("BY",0);
roiManager("Select", Array.getSequence(roiManager("count")));
RoiManager.translate(X_roi*200, Y_roi*200);

//Save ROIs
roiManager("Select", Array.getSequence(roiManager("count")));
roiManager("Save", dir+name+"_RoiSet_t.zip");

//Save ROI characteristics
close("Results"); 
newImage("test", "8-bit black", 10240, 10240, 1);
run("Set Scale...", "distance=1 known=5 unit=nm");
roiManager("Select", Array.getSequence(roiManager("count")));
run("Set Measurements...", "centroid bounding redirect=None decimal=1");
run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file save_column");
roiManager("Measure");
saveAs("Results", dir+name+"_RoiSet_t_characteristics.csv");
close("*");
close("Results");
selectWindow("ROI Manager"); 
run("Close"); 
close(name+"_RoiSet_t_characteristics.csv");
close(name+"_Roi_characteristics.csv");






