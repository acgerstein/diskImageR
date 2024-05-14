function alterImageSize(file) {
	selectWindow(file);
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
	run("Revert");
	alterImageSize(getTitle);
	run("8-bit");
	setThreshold(181, 255);
	run("Convert to Mask");
	roiManager("reset");
	roiManager("Show All with labels");
	roiManager("Show All");
	run("Analyze Particles...", "size=2500-4500 circularity=0.50-1.00 show=Outlines display exclude add");
	if (nResults ==0){
		print("Trying parameter set 2");
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
		run("Analyze Particles...", "size=2500-4500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying parameter set 3");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(81, 255);
		run("Convert to Mask");
		run("Analyze Particles...", "size=2500-4500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying parameter set 4");
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
		run("Analyze Particles...", "size=2500-4500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying parameter set 5");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(216, 255);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=2500-4500 circularity=0.50-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with less stringent circularity");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(181, 255);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=2000-5000 circularity=0.2-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with less stringent circularity, parameter set 2");
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
		run("Analyze Particles...", "size=2000-4000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with less stringent circularity, parameter set 3");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(81, 255);
		run("Convert to Mask");
		run("Analyze Particles...", "size=2000-4000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with less stringent circularity, parameter set 4");
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
		run("Analyze Particles...", "size=2000-4000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with different thresholding, parameter set 1");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(125, 162);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=2000-4000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with different thresholding, parameter set 2");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(113, 134);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=2000-50000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with different thresholding, parameter set 3");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(113, 173);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=2000-4000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		print("Trying with different thresholding, parameter set 4b");
		close();
		selectWindow(getTitle);
		run("Revert");
		alterImageSize(getTitle);
		run("8-bit");
		setThreshold(97, 129);
		run("Convert to Mask");
		roiManager("reset");
		roiManager("Show All with labels");
		roiManager("Show All");
		run("Analyze Particles...", "size=10000-50000 circularity=0.20-1.00 show=Outlines display exclude add");
	}
	if (nResults ==0){
		exit("Disk not identified, macro needs to be updated to account for photograph specifics.  Email Aleeza at gerst035@umn.edu for assistance");
	}
	if (nResults > 1){
		exit("More than one disk identified. Please ensure no other circles are present in the photograph (e.g., in labels) and rerun.");
	}
}

function findDiskLarge(file){
	run("Clear Results");
	print("Trying large disk parameter set 1");
	selectWindow(getTitle);
	run("Revert");
	alterImageSize(getTitle);
	run("8-bit");
	setThreshold(181, 255);
	run("Convert to Mask");
	roiManager("reset");
	roiManager("Show All with labels");
	roiManager("Show All");
	run("Analyze Particles...", "size=6000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  if (nResults ==0){
  	print("Trying parameter set 2");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(97, 129);
  	run("Convert to Mask");
  	roiManager("reset");
  	roiManager("Show All with labels");
  	roiManager("Show All");
  	run("Analyze Particles...", "size=8000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 3");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(97, 150);
  	run("Convert to Mask");
  	roiManager("reset");
  	roiManager("Show All with labels");
  	roiManager("Show All");
  	run("Analyze Particles...", "size=8000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 4");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(97, 200);
  	run("Convert to Mask");
  	run("Analyze Particles...", "size=10000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 5");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(97, 255);
  	run("Convert to Mask");
  	roiManager("reset");
  	roiManager("Show All with labels");
  	roiManager("Show All");
  	run("Analyze Particles...", "size=10000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 6");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(255, 255);
  	run("Convert to Mask");
  	roiManager("reset");
  	roiManager("Show All with labels");
  	roiManager("Show All");
  	run("Analyze Particles...", "size=8000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 7");
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
  	run("Analyze Particles...", "size=8000-20000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 8");
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
  	run("Analyze Particles...", "size=4000-25000 circularity=0.20-1.00 show=Outlines display exclude add");
  }
  if (nResults ==0){
  	print("Trying parameter set 9");
  	close();
  	selectWindow(getTitle);
  	run("Revert");
  	alterImageSize(getTitle);
  	run("8-bit");
  	setThreshold(243, 255);
  	run("Convert to Mask");
  	roiManager("reset");
  	roiManager("Show All with labels");
  	roiManager("Show All");
  	run("Analyze Particles...", "size=4000-25000 circularity=0.20-1.00 show=Outlines display exclude add");
  }


  if (nResults ==0){
  	exit("Disk not identified, macro needs to be updated to account for photograph specifics.  Email Aleeza at gerst035@umn.edu for assistance");
  }
  if (nResults > 1){
  	exit("More than one disk identified. Please ensure no other circles are present in the photograph (e.g., in labels) and rerun.");
  }
}

//Actual work flow starts here:
print("Starting imageJ macro");
parts=split(getArgument(), "*");

dir1 = parts[0];
dir2 = parts[1];
knownDiam = parts[2];

print("Input directory: "+dir1);
print("Output directory: "+dir2);
print("Disk diameter: "+knownDiam);
diam10 = 10/knownDiam;
list = getFileList(dir1);
print("Number of images: " + list.length);
setBatchMode(true);

//Append the appropriate file system separator, if required
if(!endsWith(dir1, File.separator)) {
	 dir1 = dir1 + File.separator;
}
if(!endsWith(dir2, File.separator)) {
	 dir2 = dir2 + File.separator;
}
outputFolder = dir2;
  
for (i=0; i<list.length; i++) {
  showProgress(i+1, list.length);

	open(dir1 + list[i]);
	print("Current image: "+list[i]);
	//setMinAndMax(50, 250);
	//The filename is automatically set to be the title of the image (so title images accordingly)
	filename = File.nameWithoutExtension;
	name = getTitle;

	run("Set Measurements...", "area mean centroid center perimeter redirect=None decimal=0");
  run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	//alterImageSize(getTitle);

	if(diam10>1.25){
		print("small disk");
		findDisk(getTitle);
	}
	if(diam10<=1.25){
		print("Large disk");
		findDiskLarge(getTitle);
	}
	close();
	run("Revert");
	alterImageSize(getTitle);
	centerX = getResult("X", 0);
	centerY = getResult("Y", 0);
	discDiam = 2*sqrt(getResult("Area")/3.1412);
	convert = discDiam/knownDiam;
	makePoint(centerX, centerY);
	run("Clear Results");

	setMinAndMax(50, 250);
	makeLineE(centerX, centerY, 40*convert, 5);

	Angle=0;
	while (Angle < 360){
		Angle = Angle + 5;
		makeLineE(centerX, centerY, 40*convert, Angle);

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

	//Save as spreadsheet compatible text file
	saveAs("Results", outputFolder+filename+".txt");
	close();
}
