## Measurement devices and data recording software
----------------------------------------------------------------------------------------------------------------------------------------

### [Hexoskin Smart Shirt](https://www.hexoskin.com/) reference 
- Provides the following sensors:
	- Fabric-electode Electrocardiogram (ECG)
	- Respiratory Inductance Plethysmography (RIP) chest belts
	- Accelerometer
- Also provides pre-calibrated Tidal Volume (TV) with limited accuracy.

### [BIOPAC](https://www.biopac.com/) reference   
- MP36R 4-channel system.
- Transducers include:
	-3-lead Electrocardiogram (ECG)
	- Stress-based chest belts to be placed on the thorax and abdomen.
	- Pre-calibrated pneumotach provides airflow information, that is converted to lung air volume by integration.
	- Finger Photoplethysmograph (PPG).
- Using Biopac Acqknowledge software to read the data and the start timestamp, accurate to seconds.

### [Ettus](https://www.ettus.com/) USRP for NCS thorax (near heart) and NCS abdomen measurements (Labview)
* This has updated over the months - initially a Ettus B210 MIMO was implemented with same carrier frequency and different IF for thorax and abdomen NCS measurements. This MIMO was baseband synchronized.
* This setup also used a mic to record audio data (required by the study protocol, removed in following setup due to time constraint).
* The abdomen sensor was much affected by motion of the wire. So setup updated to work with more portable B200mini.
* The new Labview code is updated version that includes software synchronization of NCS thorax and abdomen data coming from different SDR units (for eg., two B200 minis, or one B200 mini and one B200).
* These two-independent units are used with different carrier frequencies for thorax (1.82 GHz) and abdomen (1.9 GHz), with same IF (51 kHz). 
* Labview code saves NCS data at 50kHz for each study routine automatically when selected, with start time-stamp accurate to milliseconds.
* Synchronized audio and visual instructions for performing the routines along with viusal timer.


## Lung Volume, Breath Rate, Heart Rate, Heart Rate Variation Analysis
------------------------------------------------------------------------------------------------------------
### Volume Estimation

Lung volume is an important parameter to be measured for respiratory health monitoring. We can estimate instantaneous volume of air being exchanged by the lungs by calibrating the respiratory motion detected by the NCS.  
* Hexoskin provides pre-calibrated tidal volume, TV (lung volume only during normal breathing) estimate @1Hz sampling rate. 
* With Biopac airflow, we can get more accurate instantaneous volume information by integrating the airflow for each respiration cycle.
* Airflow waveform shows positive airflow during inspiration and negative airflow during expiration, with zero-crossing (ZC) at the inspire/expire transition. 
* This waveform is noisy and may have several ZC near the actual ZC, and needs accurate detection, this is implemented in zeroCrossDet.m. Following the volume estimation, we fit this volume to calibrate Biopac belts.
* NCS is calibrated (small calibration period of <10s) using the reference volume estimate (instantaneous volume with Biopac, and per breath volume with Hexoskin). Calibration involves fitting equations, with linear or quadrativ fitting (we have tried both, with linear being more accurate for most cases), using least-square fitting implemented by Matlab toolbox.

### BR and HR Estimation by Peak Detection

* Several peak detection algorithms have been tried, with mimimum tuning parameters and efficient peak detection for non-stationary respiration waveforms, specially with the presence of different breathing (normal, slow-deep, fast breathing and breath hold durations).
* Earlier codes have used Automated Multiscaled Based Peak Detection (AMPD) algorith. 
* For newer versions, I've switched to a semi-automated method, W. Lu, "A semi-automatic method for peak and valley detection in free-breathing respiratory waveforms.". It has slightly more tuning parameters, but they stay fix (unless the frequency varation is huge, say breathing and heartbeat).
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
	* 2D LF/HF graph
* We have tried both fundamental and harmonic NCS to estimate RR interval, as peak from Harmonic is much more accurate, but with *low* SNR.
