original = getImageID();
run("Duplicate...", "duplicate"); //duplicating original
run("Z Project...", "projection=[Max Intensity]"); //make max intensity projection to get background
background = getImageID();
run("Subtract Background...", "rolling=70 light"); //improve background
run("Auto Threshold", "method=Huang2"); //threshold channels using the Huang method
//run("Fill Holes");
run("Maximum...", "radius=4"); // fill small holes in thresholded channels
run("Minimum...", "radius=4");
mask=getImageID(); //this is binary mask of channels
//Align mask
imageCalculator("Multiply create stack", "moving_bead_segmented.tif", mask);
//invert image and remove pixel scale
run("Invert", "stack");
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
//Run MTrack2 for tracking
run("MTrack2 ", "minimum=20 maximum=999999 maximum_=10 minimum_=1000 show_1");