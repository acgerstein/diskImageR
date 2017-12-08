function alterImageSize(file) {
	selectWindow(getTitle);
	picWidth = getWidth();
	picHeight = getHeight();
	run("Size...", "width=1000 constrain interpolation=None");
	}

function makeLineE(centerX, centerY, length, angle) {
	angle = -angle * PI / 180;
	dX = cos(angle) * length;
	dY = sin(angle) * length;
	makeLine(centerX, centerY, centerX + dX, centerY + dY);
	}

function findDisk(file){
	run("Clear Results");
	selectWindow(getTitle);
	alterImageSize(getTitle);
	run("8-bit");
	setThreshold(181, 255);
	run("Convert to Mask");
	roiManager("reset");
	roiManager("Show All with labels");
	roiManager("Show All");
	print("Trying parameter set 1");
	run("Analyze Particles...", "size=700-2500 circularity=0.50-1.00 show=Outlines display exclude add");
	if (nResults !=16){
		print("Trying parameter set 2");
		run("Clear Results");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(150, 255);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=700-2500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults !=16){
		print("Trying parameter set 3");
		run("Clear Results");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(81, 255);
		run("Convert to Mask");
		run("Analyze Particles...", "size=700-2500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults !=16){
		print("Trying parameter set 4");
		run("Clear Results");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(200, 255);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=700-2500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Disks not identified, macro needs to be updated to account for photograph specifics.  Email Aleeza at acgerstein@gmail.com for assistance.");
	}
	if (nResults ==16){
		print("All 16 disks successfully identified");
	}

	if (nResults != 16){
		print("Unable to identify 16 disks, macro needs to be updated to account for photograph specifics. Email Aleeza at acgerstein@gmail.com for assistance.");
	}
}


//Actual work flow starts here:
print("Starting imageJ macro");
folders = getArgument;
delimiter = "*";
parts=split(folders, delimiter);
dir1 = parts[0];
dir2 = parts[1];
dir3 = parts[2];
knownDiam = parts[3];

print("Input directory: "+dir1);
print("Output directory: "+dir2);
print("Disk coordinates directory: "+dir3);
print("Disk diameter: "+knownDiam);
diam10 = 10/knownDiam;
fileList = getFileList(dir1);
print("Number of images: " + fileList.length);
setBatchMode(true);
for (i=0; i<fileList.length; i++){
//	print(i);
	showProgress(i+1, fileList.length);
	open(dir1 + fileList[i]);
	print("Current image: "+fileList[i]);
	outputFolder = dir2;
	//The filename is automatically set to be the title of the image (so title images accordingly)
	filename = substring(getTitle, 0, lengthOf(getTitle)-4);
	name = getTitle;

	run("Set Measurements...", "area mean centroid redirect=None decimal=0");
     run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

	findDisk(getTitle);
	selectWindow(getTitle);
	close();
	selectWindow(getTitle);
	run("Revert");

   //save the results table, will have to use this to figure out which disk is which based on X, Y coordinates
	saveAs("Results", dir3+filename+"_ResultsTable.txt");
	List.set("X0", getResult("X", 0));
	List.set("X1", getResult("X", 1));
	List.set("X2", getResult("X", 2));
	List.set("X3", getResult("X", 3));
	List.set("X4", getResult("X", 4));
	List.set("X5", getResult("X", 5));
	List.set("X6", getResult("X", 6));
	List.set("X7", getResult("X", 7));
	List.set("X8", getResult("X", 8));
	List.set("X9", getResult("X", 9));
	List.set("X10", getResult("X", 10));
	List.set("X11", getResult("X", 11));
	List.set("X12", getResult("X", 12));
	List.set("X13", getResult("X", 13));
	List.set("X14", getResult("X", 14));
	List.set("X15", getResult("X", 15));
	List.set("Y0", getResult("Y", 0));
	List.set("Y1", getResult("Y", 1));
	List.set("Y2", getResult("Y", 2));
	List.set("Y3", getResult("Y", 3));
	List.set("Y4", getResult("Y", 4));
	List.set("Y5", getResult("Y", 5));
	List.set("Y6", getResult("Y", 6));
	List.set("Y7", getResult("Y", 7));
	List.set("Y8", getResult("Y", 8));
	List.set("Y9", getResult("Y", 9));
	List.set("Y10", getResult("Y", 10));
	List.set("Y11", getResult("Y", 11));
	List.set("Y12", getResult("Y", 12));
	List.set("Y13", getResult("Y", 13));
	List.set("Y14", getResult("Y", 14));
	List.set("Y15", getResult("Y", 15));
	List.set("area", getResult("Area", 0));
     list = List.getList();


//
//	alterImageSize(getTitle);
//
	selectWindow(getTitle);
	picWidth = getWidth();
	picHeight = getHeight();
	run("Size...", "width=1000 constrain interpolation=None");
	run("8-bit");
	setMinAndMax(50, 250);

	area = getResult("Area", 0);
//	area = List.get("area");
	discDiam = 2*sqrt(area/3.1412);
	//the next line is just for debugging
     //	knownDiam = 6;
	//	print(discDiam);
	convert = discDiam/knownDiam;

	
//walk through each of the 16 disks

//
	for(m=0; m<16; m++) {
//		m =5;
//		print(m);
//		print(list);

		run("Clear Results");

		placeX = "X"+m;
		placeY = "Y"+m;
		centerX = List.get(placeX);
		centerY = List.get(placeY);
	
//		print(convert);
		makePoint(centerX, centerY);

//function makeLineE(centerX, centerY, length, angle) {
//	angle = -angle * PI / 180;
//	dX = cos(angle) * length;
//	dY = sin(angle) * length;
//	makeLine(centerX, centerY, centerX + dX, centerY + dY);
//	}


//		makeLineE(centerX, centerY, 25*convert, 5);


		Angle=0;
		while (Angle < 360){
			Angle = Angle + 2;
			makeLineE(centerX, centerY, 30*convert, Angle);

			// Get profile and display values in "Results" window
			profile = getProfile();
			for (j=0; j<profile.length; j++){
				k = nResults;
				setResult("X", k, j);
				setResult("Value", k, profile[j]);
			}
			updateResults();

			// Plot profile
			Plot.create("Profile", "X", "Value", profile);
		}
		n = m+1;
		//Save as spreadsheet compatible text file
		saveAs("Results", outputFolder+filename+"_"+n+".txt");
//		print("profile from disk "+ n+ " saved");
		}
	close();
}
