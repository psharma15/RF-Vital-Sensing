# Lung Volume, BR, HR, HRV Analysis

*New and improved codes for doing various analysis of NCS vital signs data, Updated Fall 2019*

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
* Synchronized audio/ visual instructions for performing the routine along with viusal timer information.

## Data analysis codes:
The code hierarchy is: 
* filterLpHp.m (Filtering any waveforms)
* findMaxMin.m (Processing maxima and minima peaks in a waveform)
	* peakDet3.m (Finding the peaks using efficient, [Lu et al., "A semi-automatic method for peak and valley detection in free-breathing respiratory waveforms," Med Phys, 2006] algorithm)
* zeroCrossDet.m (Efficient zero-cross (ZC) detection along with slope information - for a noisy waveform)

--------------------------------------------------------------------------------------------------------------------------

* ncsBioProcess.m (Synchronizes and gets NCS-Biopac Lung Volume, can add BR and HR)
	* ncsBioSync.m (Synchronizing NCS and Biopac waveforms)
		* readNcsData (Reads and organizes NCS data)
		* readBioData (Reads and organizes Biopac data)
	* bioAirflowTV.m (Volume Estimation from airflow pneumotach: options to get TV/instantaneous vol.)
		* airflowToVol.m (Performs integration of airflow after detecting each cycle)
			* zeroCrossDet.m
	* bioBeltVolCalib.m (Calibrates Biopac chest belts using pneumotach volume data with different presets, using findMaxMin.m)
	* ncsVolCalib.m (Calibrates NCS respiration using pneumotach volume data with different presets, using findMaxMin.m)
	* bioBeltVol.m (Using calibrated coefficients to estimate volume from Biopac belts, using findMaxMin.m)
	* ncsVol (Using calibrated coefficients to estimate volume from NCS respiration, using findMaxMin.m)
	
--------------------------------------------------------------------------------------------------------------------------

* ncsAttnTestProcess.m (Synchronizes NCS-Biopac-Attention test observation and gets HRV information)
	* ncsBioSync.m (Synchronizing NCS and Biopac waveforms)
		* readNcsData (Reads and organizes NCS data)
		* readBioData (Reads and organizes Biopac data)
	* brEst.m (Estimates Breath Rate from input NCS or Biopac waveforms, using findMaxMin.m)
	* hrEst.m (Estimates Heart Rate from input NCS, using findMaxMin.m)
	* ecgHR.m (Estimates Heart Rate from Biopac ECG, using Matlab's function findpeaks.m)
	* hrvFeatureEst (Estimates Heart Rate Variation (HRV) features from NCS fundamental or harmonic heartbeat and ECG waveforms)
	
** Following is some more detail about different parts of code:**

### Volume Estimation

* So far we have estimated tidal volume from NCS and compared that with Hexosking pre-estimated tidal volume, which is defined as scaled version of difference between maximumk-inspiratory and minimum-expiratory respiration waveform for each respiration cycle.
* With Biopac airflow, we can get more accurate instantaneous volume information by integrating the airflow for each respiration cycle.
* Airflow waveform shows positive airflow during inspiration and negative airflow during expiration, with zero-crossing (ZC) at the inspire/expire transition. 
* This waveform is noisy and may have several ZC near the actual ZC, and needs accurate detection, this is implemented in zeroCrossDet.m
* Following the volume estimation, we fit this volume to calibrate Biopac belts and NCS respiration and we have several fitting equations to perform this. We have used Matlab Toolbox to do this fitting. 
* This calibration period ranges between 5-25 s. Once the calibration coefficients are estimated, they can be used to estimate the volume from belts and NCS, which are compared against each other.

### Peak Detection

* Several peak detection algorithms have been tried, with mimimum tuning parameters and efficient peak detection for non-stationary respiration waveforms, specially with the presence of different breathing (normal, slow-deep, fast breathing and breath hold durations).
* Earlier codes have used Automated Multiscaled Based Peak Detection (AMPD) algorith. 
* For this version, I've switched to a semi-automated method, with slightly more tuning parameters, but that stay fix for most cases (unless the frequency varation is huge, say breathing and heartbeat).
* Refer to above source for this algorithm, slight modifications have been made to implement this algorithm, which works well and provides both maxima and minima peaks with identifying information to distinguish maxima with minima.

### Heart Rate Variation

* We have estimated RR interval using both NCS and ECG and calculated HRV features including
	* mean RR (or NN) interval (mean of RR interval, ms)
	* SDNN (Standard deviation of NN interval, ms)
	* RMSSD (Square-root of the mean of the sum of the squares of differences between *adjacent* NN intervals)
	* SDSD (Standard-deviation of the differences between *adjacent* NN intervals)
	* pNN50 (Number of pairs of *adjacent* NN intervals differing by more than 50ms, divided by total number of all NN intervals)
	* LF power (Power in 0.04-0.15 Hz)
	* HF power (Power in 0.15-0.7 Hz)
* We have tried both fundamental and harmonic NCS, as peak from Harmonic is much more accurate, but with *low* SNR.







