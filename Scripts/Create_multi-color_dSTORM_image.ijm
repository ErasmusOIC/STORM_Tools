///////////////////////////////////////////////////////////////
/// Script: "Create_multi-color_dSTORM_image.ijm"           ///
/// Authors: Lieke Koornneef & Johan Slotman                ///
/// Affiliation: Erasmus MC, Rotterdam, The Netherlands     ///
/// Contact: j.slotman@erasmusmc.nl  						///
/// License: LGPLv3                                         ///
/// Date: 23-01-2023                                        ///
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
none = newArray("none");
tiffiles = Array.concat(none,tiffiles);

name = call("ij.Prefs.get","ST.Exp_name", "Exp1");

Dialog.create("Details of multicolor imaging");
Dialog.addString("Experiment name", name);
Dialog.addChoice("Image of channel 1", tiffiles, "none");
Dialog.addChoice("Image of channel 2", tiffiles, "none");
Dialog.addChoice("Image of channel 3", tiffiles, "none");

Dialog.show();

name = Dialog.getString();
file1 = Dialog.getChoice();
file2 = Dialog.getChoice();
file3 = Dialog.getChoice();

call("ij.Prefs.set","ST.Exp_name", name);

if(!(file1 == "none")){
	run("Bio-Formats", "open=["+dir+file1+"] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("Ch1_original");
	run("Duplicate...", "duplicate channels=1");
	rename("Ch1_new");
	setMinAndMax(0,65535);
	run("16-bit");
	run("Magenta");
	run("Enhance Contrast", "saturated=0.35");
	close("Ch1_original");}

if(!(file2 == "none")){ 
	run("Bio-Formats", "open=["+dir+file2+"] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("Ch2_original");
	run("Duplicate...", "duplicate channels=1");
	rename("Ch2_new");
	setMinAndMax(0,65535);
	run("16-bit");
	run("Red");
	run("Enhance Contrast", "saturated=0.35");
	close("Ch2_original");}

if(!(file3 == "none")){
	run("Bio-Formats", "open=["+dir+file3+"] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("Ch3_original");
	run("Duplicate...", "duplicate channels=1");
	rename("Ch3_new");
	setMinAndMax(0,65535);
	run("16-bit");
	run("Green");
	run("Enhance Contrast", "saturated=0.35");
	close("Ch3_original");}
	
	
filenames = newArray(file1,file2,file3);
noneCount = 0;
chCount = 1;
mergeParam = "";
mergeName ="";
for(i=0;i<filenames.length;i++){
	if(filenames[i]=="none"){
		noneCount++;
	} else {
		mergeParam = mergeParam+"c"+chCount+"=Ch"+(i+1)+"_new ";
		mergeName = mergeName+"_Ch"+(i+1);
		chCount++;
	}
}
mergeParam = mergeParam+" create";
mergeName = mergeName+".tif";

if(noneCount<2){
	run("Merge Channels...", mergeParam);
 	saveAs("Tiff", dir+name+mergeName);
}else{
	showMessage("Only one image was selected, no merge needed");
}







