clear all;
close all;
%LoadNavFilterTestDataStruct
LoadNavFilterTestData;
sim('NavFilterTestHarness24')
PlotNavFilterData24