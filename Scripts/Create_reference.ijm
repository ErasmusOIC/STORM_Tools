///////////////////////////////////////////////////////////////
/// Script: "Create_reference.ijm"                          ///
/// Authors: Lieke Koornneef & Johan Slotman                ///
/// Affiliation: Erasmus MC, Rotterdam, The Netherlands     ///
/// License: LGPLv3                                         ///
/// Date: 13-01-2023                                        ///
///////////////////////////////////////////////////////////////

dir = getDirectory("Select directory");
files = getFileList(dir);
name = call("ij.Prefs.get","ST.Exp_name", "Exp1");
ch_nr = call("ij.Prefs.get","ST.Nr_Ch", 3);

Dialog.create("Details of multicolor imaging");
Dialog.addString("Experiment name", name);
Dialog.addChoice("Select Z-stack image", files, "none");
Dialog.addNumber("Number of channels", ch_nr);
Dialog.show();

name = Dialog.getString();
file = Dialog.getChoice();
ch_nr = Dialog.getNumber();

call("ij.Prefs.set","ST.Exp_name", name);
call("ij.Prefs.set","ST.Nr_Ch", ch_nr);

//Reference Image
open(dir+file);
rename("image");
roiManager("reset");
close("Results");

getDimensions(width, height, channels, slices, frames);
meanmax = 0;
id = 0;
for(i=1;i<=slices;i++){
	Stack.setPosition(1,i,1);
    getStatistics(area, mean, min, max, std, histogram);
    if(max>meanmax){
    	meanmax = max;
        id = i;
     }
}
channel1_selected = id;
Stack.setPosition(1,channel1_selected,1);
waitForUser("The slice with the highest maximum value is selected. Check if this is also the slice where the whole nucleus is in focus. If this is not the case, remember the number of the new selected slice.");
Dialog.create("Select most focused slice per channel");
Dialog.addNumber("Slice Ch1:", channel1_selected);
Dialog.show();
ch_1 = Dialog.getNumber();
run("Duplicate...", "duplicate channels=1 slices="+ch_1);

selectWindow("image");
if(ch_nr > 1){
	meanmax = 0;
	id = 0;
	for(i=1;i<=slices;i++){
		Stack.setPosition(2,i,1);
	    getStatistics(area, mean, min, max, std, histogram);
	    if(max>meanmax){
	    	meanmax = max;
	        id = i;
	     }
	}
	channel2_selected = id;
	Stack.setPosition(2,channel2_selected,1);
	
	waitForUser("The slice with the highest maximum value is selected. Check if this is also the slice where the whole nucleus is in focus. If this is not the case, remember the number of the new selected slice.");
	Dialog.create("Select most focused slice per channel");
	Dialog.addNumber("Slice Ch2:", channel2_selected);
	Dialog.show();
	ch_2 = Dialog.getNumber();
	run("Duplicate...", "duplicate channels=2 slices="+ch_2);
}

selectWindow("image");
if(ch_nr > 2){
	meanmax = 0;
	id = 0;
	for(i=1;i<=slices;i++){
		Stack.setPosition(3,i,1);
	    getStatistics(area, mean, min, max, std, histogram);
	    if(max>meanmax){
	    	meanmax = max;
	        id = i;
	     }
	}
	channel3_selected = id;
	Stack.setPosition(3,channel3_selected,1);
	
	waitForUser("The slice with the highest maximum value is selected. Check if this is also the slice where the whole nucleus is in focus. If this is not the case, remember the number of the new selected slice.");
	Dialog.create("Select most focused slice per channel");
	Dialog.addNumber("Slice Ch3:", channel3_selected);
	Dialog.show();
	ch_3 = Dialog.getNumber();
	run("Duplicate...", "duplicate channels=3 slices="+ch_3);
}

if(ch_nr == 1){
	selectWindow("image-1");
}

if(ch_nr == 2){
	run("Merge Channels...", "c1=image-1 c2=image-2 create");
}
if(ch_nr == 3){
	run("Merge Channels...", "c1=image-1 c2=image-2 c3=image-3 create");
}

saveAs("Tiff", dir+name+"_Composite.tif");

run("Flip Vertically", "stack");
run("Scale...", "x=20 y=20 z=1.0 width=10240 height=10240 depth=3 interpolation=Bilinear average create");
waitForUser("Select border of the nucleus using the Polygon selections tool");
roiManager("Add");
roiManager("Save", dir+name+"_Roi.roi");
run("Set Measurements...", "area mean min bounding redirect=None decimal=3");
roiManager("Measure");
saveAs("Results", dir+name+"_Roi_characteristics.csv");
run("Close All");
roiManager("reset");
close("Results");


