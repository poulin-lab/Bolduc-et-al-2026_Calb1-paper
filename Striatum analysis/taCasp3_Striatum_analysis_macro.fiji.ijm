//***IMPORTANT NOTE BEFORE RUNNING THE MACRO***
//The file directory at row 11 must fit the directory where the chosen image is located, otherwise it won't work
//Furthermore, the following files must be generated within the working directory, otherwise it won't work: Noise, CSVfiles, ROIs

//Choose the image
path = File.openDialog("Select a File");

//Open the image
open(path);

//Extract the image name deleting the file name directory
image_name = getTitle();
output = path.replace(image_name, "")



// Split the channels
run("Split Channels");

//Ajust the contrast for each channels
selectImage("C1-"+image_name);
run("Enhance Contrast", "saturation=0.35");

selectImage("C2-"+image_name);
run("Enhance Contrast", "saturation=0.35");

// Tile the images
run("Tile");

waitForUser("Trace and add to the ROI manager the DV, ML and MV-DL arrows, respecting that order. Make three squares in the cortex for the noise. Then, Click OK");

//---------------------------------------------------------------------------------

//Add the DV, ML, and MV-DL axis to the images
selectImage("C1-"+image_name);
roiManager("Show All");

selectImage("C2-"+image_name);
roiManager("Show All");

waitForUser("Take a screen shot by pressing shift+command+3 and click OK");

//---------------------------------------------------------------------------------

//Save the ROIs
roiManager("save", output + "ROIs//" + image_name + "ROIset.zip");

//---------------------------------------------------------------------------------
//Extract the data

//For mCherry DV
//----------------
selectImage("C1-"+image_name);
roiManager("Select", 0);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output  + "CSVfiles" + File.separator + "mCherry_DV_" + image_name + ".csv");

//For mCherry ML
//----------------
selectImage("C1-"+image_name);
roiManager("Select", 1);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output  + "CSVfiles" + File.separator + "mCherry_ML_" + image_name + ".csv");

//For EYFP MV-DL
//----------------
selectImage("C1-"+image_name);
roiManager("Select", 2);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output  + "CSVfiles" + File.separator + "mCherry_MV-DL_" + image_name + ".csv");


//For TH DV
//----------------
selectImage("C2-"+image_name);
roiManager("Select", 0);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output  + "CSVfiles" + File.separator + "TH_DV_" + image_name + ".csv");


//For TH ML
//----------------
selectImage("C2-"+image_name);
roiManager("Select", 1);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output + "CSVfiles" + File.separator + "TH_ML_" + image_name + ".csv");

//For TH MV-DL
//----------------
selectImage("C2-"+image_name);
roiManager("Select", 2);

// Get profile and display values in "Results" window
run("Clear Results");
  profile = getProfile();
  for (i=0; i<profile.length; i++)
      setResult("Value", i, profile[i]);
  updateResults;
  
// Save as spreadsheet compatible text file
  saveAs("Results", output + "CSVfiles" + File.separator + "TH_MV-DL_" + image_name + ".csv");


//---------------------------------------------------------------------------------
//Make the ROI for the noise
//Close the arrows ROIs
roiManager("Select", 0);
roiManager("Delete");

roiManager("Select", 0);
roiManager("Delete");

roiManager("Select", 0);
roiManager("Delete");

//Then select the Noise ROIs and save them
selectImage("C1-"+image_name);
roiManager("Show All");
roiManager("Measure");

saveAs("Measurements", output + "CSVfiles" + File.separator + "Noise_mCherry_" + image_name + ".csv");
run("Clear Results");

selectImage("C2-"+image_name);
roiManager("Show All");
roiManager("Measure");

saveAs("Measurements", output + "CSVfiles" + File.separator + "Noise_TH_" + image_name + ".csv");
run("Clear Results");

//---------------------------------------------------------------------------------

//Close the session
//Close the images
selectImage("C1-"+image_name);
close();

selectImage("C2-"+image_name);
close();

//Close the arrows ROIs
roiManager("Select", 0);
roiManager("Delete");

roiManager("Select", 0);
roiManager("Delete");

roiManager("Select", 0);
roiManager("Delete");