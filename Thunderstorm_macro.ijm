//Run analysis
Channel1 = ;
Channel2 = ;

//Merge localization

//Apply molecule merge to detect fiducials, export
run("Show results table", "action=drift smoothingbandwidth=0.25 method=[Fiducial markers] ontimeratio=0.1 distancethr=40.0 save=false");
run("Export results", "floatprecision=5 filepath=[C:\\Thunderstorm test\\130628_data\\fiducials 488.csv] fileformat=[CSV (comma separated)] sigma=true intensity=true offset=true saveprotocol=true x=true y=true bkgstd=true id=true uncertainty_xy=true frame=true detections=true");


//Reset 
run("Show results table", "action=reset");


//Correct drift
run("Show results table", "action=drift smoothingbandwidth=0.25 method=[Fiducial markers] ontimeratio=0.2 distancethr=40.0 save=false");

//Correct drift save to file
run("Show results table", "action=drift smoothingbandwidth=0.25 path=[C:\\Thunderstorm test\\130628_data\\drift.json] method=[Fiducial markers] ontimeratio=0.1 distancethr=40.0 save=true");

