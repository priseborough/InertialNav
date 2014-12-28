1) Build the file <main - closed loop test harness - single precision.cpp> and generate an executable
	a) By running 'make' (and 'make clean' after header edits)
	b) Or by running this line (e.g. on Windows):
	  g++ -o estimator_closed_loop_test estimator.cpp main_closed_loop_float.cpp

2) Extract all files from the quadOptFlowLogData.zip archive

3) Run the executable from the same directory you placed the files in from step 2):
	*nix: ./estimator_closed_loop_test
	Windows: estimator_closed_loop_test

4) It should generate the following output files (Reference output data has been provided in ValidationOutputData.zip)

	CovDataOut.txt
	EulDataOut.txt
	FlowRawOut.txt
	GPSrawOut.txt
	MagFuse.txt
	OnboardVelPosDataOut
	OptFlowFuse.txt
	RefVelPosDataOut.txt
	StateDataOut.txt
	TasFuse.txt
	ValidationOut.txt
	VelPosFuse.txt
	
5) Generate plots by running the python scripts provided in the code directory from the command line, eg:

	python plot_position.py &
	python plot_optflow.py &
	python plot_rng.py

(Reference plots have been provided in ValidationOutputPlots.zip)
