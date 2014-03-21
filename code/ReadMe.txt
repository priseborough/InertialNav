1) Build the file <main - closed loop test harness - single precision.cpp> and generate an executable
	a) By running 'make' (and 'make clean' after header edits)
	b) Or by running this line (e.g. on Windows):
	  g++ -o estimator_closed_loop_test estimator.cpp main_closed_loop_float.cpp

2) Run the executable in the same directory as the ATT, GPS, IMU, MAG, NTUN and timing.txt files:
	*nix: ./estimator_closed_loop_test
	Windows: estimator_closed_loop_test
3) It should generate the following output files (outputs have been provided for reference)

	CovDataOut.txt
	EulDataOut.txt
	MagFuse.txt
	RefVelPosDataOut.txt
	StateDataOut.txt
	TasFuse.txt
	VelPosFuse.txt
	
4) These can be plotted in Matlab using the PlotCcodeOutput.m script files