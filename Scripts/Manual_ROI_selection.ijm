///////////////////////////////////////////////////////////////
/// Script: "Manual_ROI_selection.ijm"                      ///
/// Authors: Maarten Paul & Lieke Koornneef                 ///
/// Affiliation: Erasmus MC, Rotterdam, The Netherlands     ///
/// Contact: j.slotman@erasmusmc.nl  				   	    ///
/// License: LGPLv3                                         ///
/// Date: 13-01-2023                                        ///
///////////////////////////////////////////////////////////////

//circle Clicker 
//Click to create cirle of 375 nm
//exit script by shift clicking when the point tool is on

//parameters
radius = 375; //nm
pixelsize = 5 //nm

Dialog.create("ROI clicker");
Dialog.addNumber("radius (nm)", 375);
Dialog.addNumber("pixel size (nm)", 5);

Dialog.show();
radius = Dialog.getNumber();
pixelsize = Dialog.getNumber();

run("Set Scale...", "distance=1 known=" + pixelsize +" pixel=1 unit=nm");
print("Clicker started");
setTool("point");
title=getTitle();
while(isOpen(title)){
while(IJ.getToolName()=="point" ){
	getCursorLoc(x,y,z,flags);
	if(flags==16){
	
	wait(50);
	makePoint(x, y);
	run("Enlarge...", "enlarge="+ radius);
	}
	if(flags==17){
		print("Clicker stopped");
		exit;
		}
}
}