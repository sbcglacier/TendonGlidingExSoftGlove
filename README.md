# TendonGlidingExSoftGlove
This repository includes the work done for my PhD 
The file Chung_ScaleFactors_F2.m calculates the scale factors for OpenSim for the index finger.
The file Chung_ScaleFactors_F3.m calculates the scale factors for OpenSim for the mid finger.
The file Chung_ScaleFactors_F4.m calculates the scale factors for OpenSim for the ring finger.
The file Chung_ScaleFactors_F5.m calculates the scale factors for OpenSim for the little finger.
The folder OpenSim_FingerMeshData contains the mesh files needed for the matlab files used to determine the OpenSim scale factors.

The file Apex_mesh_create_3Chambers.py is used to create the actuator in Apex.
The file Marc_setup_3Chambers.py is used to setup the chambers after the Apex mesh files have been created.

The SensitivityAnalysis.ipynb calculates the moments and bending angles from data gathered from Marc/Mentat.
Steps: Run the Marc_setup_3Chambers.py, obtain the ThreeActuator.proc file and run it. Sweep the nodes at the strain limiting layer; add the RBE elements; Correct pressure; select the 6 nodes for displacement then generate the Dx, Dz, and displacement report. 
Constrain the one side of the actuator for the reaction moments generate the Rx Rz plots and also the displacement moment report.
With the files created run the SensitivityAnalysis.ipynb file to calculate the respective moments and bending angle.

The BiPNA_cascadePlot.py uses affine transforms to cascade three actuators with varied parameters.
GenerateApex3ChamberMesh.py was used to create the mesh in Apex for the three chambers with altered width and length.
Marc_3ChamberSetupAutomation.py is used to setup the mesh in Marc after generating it in Apex.
