# EyeDiagram
Better eye diagram generation code in MATLAB. Used a lot in my DPhil Thesis.

[![View Density Eye Diagram on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/116910-density-eye-diagram)

*How much better?*

A before and after is provided here:

![left is a typical MATLAB eyediagram() eye diagram, right is the new function](https://github.com/WillMatthews/EyeDiagram/blob/master/img/demo.png?raw=true)


Unlike MATLAB's built in `eyediagram()` function, the more data you provide `eyediagram_dense()` the prettier it looks, case in point (for a different communications system with a higher error rate):

![an even prettier eye diagram](https://github.com/WillMatthews/EyeDiagram/blob/master/img/demo2.png?raw=true)


*What does your code do?*
`eyediagram_dense` creates a raster image of an eye diagram, capturing signal density at each point on the image. This function creates a diagram, alongside returning all the important information captured within it as a MATLAB struct.

`eyediagram_dense()` also has the ability to detect zero crossings, and automatically label the x axis accordingly.

I mainly work with On-Off Keying, so for the y axis the system is able to detect two levels. In a future update I will be able to detect PAM.

With optional arguments, histograms can be placed along the x and y axis, for zero crossing distributions (a way to test if you have double-Dirac) and to sanity check your bit levels.

*What is coming next?*
My goals for this repo are as follows:
|Goal|Description  |
|--|--|
|Q factor estimation  |From the eye diagram we can  estimate Q factors, which are important quantities to know for any communications system. |
|BER estimation | From calculated Q factors we can estimate the BER. This is handy if we don't have all afternoon to resolve a 10^-12 BER.|
|Eye Opening Metrics|Use standard textbook metrics and return to the user.|
|Jitter|Quantify the jitter and return it to the user.|
|More...|I'm sure I've missed lots, but more eye diagram related features will be included|
