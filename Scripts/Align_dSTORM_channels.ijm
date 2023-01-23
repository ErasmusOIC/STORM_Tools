/////////////////////////////////////////////////////////////////////////////////////////////
/// Script: "Align_dSTORM_channels.ijm"                                                   ///
/// Authors: Lieke Koornneef & Johan Slotman                                              ///
/// Affiliation: Erasmus MC, Rotterdam, The Netherlands                                   ///
/// License: LGPLv3                                                                       ///
/// Date: 13-01-2023                                                                      ///
///                                                                                       ///
/// This code is adapted from the DoM_Utrecht plugin for single-molecule analysis         ///
/// https://github.com/ekatrukha/DoM_Utrecht/                                             ///
/// Part of the class is reused: bead selection and filtering, b-spline alignment         ///
/// and application of alignment on a single-molecule localization data set.              ///
/////////////////////////////////////////////////////////////////////////////////////////////

dir = getDirectory("Select directory");
files = getFileList(dir);
none = newArray("none");
files = Array.concat(none,files);
maximum_gap_value = call("ij.Prefs.get","ST.maximum_gap_value",50);
track_length_value = call("ij.Prefs.get","ST.track_length_value",70);
name = call("ij.Prefs.get","ST.Exp_name", "Exp1");

Dialog.create("Details of multicolor imaging");
Dialog.addString("Experiment name", name);
Dialog.addMessage("Indicate your tables for channels, always starting at channel 1. If you imaged 2 channels, please indicate channel 3 as none.");
Dialog.addMessage("Channel 1 is used as reference for alignment");
Dialog.addChoice("Table for channel 1", files, "none");
Dialog.addChoice("Table for channel 2", files, "none");
Dialog.addChoice("Grouped table for channel 2", files, "none");
Dialog.addChoice("Table for channel 3", files, "none");
Dialog.addChoice("Grouped table for channel 3 ", files, "none");
Dialog.addNumber("Maximum Gap Value", maximum_gap_value);
Dialog.addNumber("Track Length Value", track_length_value);
Dialog.show();

name = Dialog.getString();
file1 = Dialog.getChoice();
file2 = Dialog.getChoice();
file2_g = Dialog.getChoice();
file3 = Dialog.getChoice();
file3_g = Dialog.getChoice();
maximum_gap_value = Dialog.getNumber();
track_length_value = Dialog.getNumber();

call("ij.Prefs.set","ST.maximum_gap_value",maximum_gap_value);
call("ij.Prefs.set","ST.track_length_value",track_length_value);
call("ij.Prefs.set","ST.Exp_name", name);

run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file save_column");

if(!(file2 == "none")){
	transformation_ch2 = name+"_Ch2_transformation.txt";
	file2_gt = name+"_Ch2_gt.txt";
	run("Autoselect fiducials", "maximum_gap="+maximum_gap_value+" track_length="+track_length_value+" select_primary="+dir+file1+" select_secondary="+dir+file2+" select_x=[Position X [nm]] select_y=[Position Y [nm]] select_z=none select_frame=[First Frame] select_precision=[Precision [nm]] select_x=[Position X [nm]] select_y=[Position Y [nm]] select_z=none select_frame=[First Frame] select_precision=[Precision [nm]]");
	run("Calculate Transformation", "select=1 select_0=2 save="+dir+transformation_ch2);
	run("Results... ", "open="+dir+file2_g);
	run("Apply Transformation", "open="+dir+transformation_ch2);
	run("Results... ", "open="+dir+file2_g);
	run("Apply Transformation", "open="+dir+transformation_ch2);
	saveAs("Results", dir+file2_gt);}

if(!(file3 == "none")){
	transformation_ch3 = name+"_Ch3_transformation.txt";
	file3_gt = name+"_Ch3_gt.txt";
	run("Autoselect fiducials", "maximum_gap="+maximum_gap_value+" track_length="+track_length_value+" select_primary="+dir+file1+" select_secondary="+dir+file3+" select_x=[Position X [nm]] select_y=[Position Y [nm]] select_z=none select_frame=[First Frame] select_precision=[Precision [nm]] select_x=[Position X [nm]] select_y=[Position Y [nm]] select_z=none select_frame=[First Frame] select_precision=[Precision [nm]]");
	run("Calculate Transformation", "select=1 select_0=2 save="+dir+transformation_ch3);
	run("Results... ", "open="+dir+file3_g);
	run("Apply Transformation", "open="+dir+transformation_ch3);
	run("Results... ", "open="+dir+file3_g);
	run("Apply Transformation", "open="+dir+transformation_ch3);
	saveAs("Results", dir+file3_gt);}





