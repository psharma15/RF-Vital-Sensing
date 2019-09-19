# Lung Volume, BR, HR, HRV Analysis

*New and improved codes for doing various analysis of NCS vital signs data*

The codes updated in this section are the improved versions of all the previous codes before Fall 2019. 

## Measurement devices and data recording software:

### BIOPAC reference measurement (Biopac Acqknowledge)
* We have updated reference setup to BIOPAC MP36R system.
* Transducers include:
	* 3-lead Electrocardiogram (ECG)
	* Stress-based chest belts to be placed on thorax and abdomen.
	* Pre-calibrated pneumotach provides airflow information, that is converted to lung air volume by integration.
	* Finger Photoplethysmograph (PPG).
* Using Biopac Acqknowledge software to read the data and the start timestamp, accurate to seconds.

### USRP B210 or (B200 mini + B200) for NCS thorax (near heart) and NCS abdomen measurements (Labview)
* This has updated over the months - initially a B210 MIMO was implemented with same carrier frequency and different IF for thorax and abdomen NCS measurements. This MIMO is baseband synchronized.
* This setup also used a mic to record audio data (required by the study protocol, removed in following setup due to time constraint).
* The abdomen sensor was much affected by motion of the wire. So setup updated to work with more portable B200mini.
* The new Labview code was update of the previous version to software synchronize NCS thorax and abdomen.
* Now two-independent B200 mini with different carrier frequencies can be used for thorax and abdomen. 
* Labview code saves NCS data at 50kHz for each study routine separately, with time stamp accurate to milliseconds.




