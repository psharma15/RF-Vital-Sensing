#Vital Sign Monitoring by Radio Frequency based Near-Field Coherent Sensing
==========================================================================================================================================
*Code for applications of radio frequency (RF) based near-field coherent sensing (NCS) in vital sign detection.*

## NCS Overview
NCS can record dielectric boundary movement of internal organs and body surfaces in the near-field region of the transmitter (Tx)
antenna. NCS can be implemented as either passive or active setup, with the Tx antenna on the chest, with optimal placement to get vital sign of interest (heart or breath or both). For the former setup, passive radiofrequency identification (RFID) tags can be put on the person's clothes to maximize the wearer comfort and minimize the tag cost, while receiver (Rx) can get the vital sign in the far-field. Mechanical movements that result in dynamic dielectric boundary changes are modulated onto the radio signals with unique digital identification (ID), which can be readily extended to monitor multiple tags and persons by a single RFID reader with good channel isolation. In the active tag approach, both Tx and Rx antennas are placed on the chest as a self-contained mobile unit without need of an external reference reader, which is then feasible for both indoor and outdoor applications. NCS is less sensitive to wearer movement and ambient motion which can be filtered out as the common-mode signal and is thus more feasible for continuous monitoring. 

As this setup provides comfortable non-invasive vital sign monitoring, it can be used for long-term monitoring. Among others, it can help improve diagnostics of respiratory diseases and sleep apnea, which can often be undetected and untreated due to lack of continous monitoring. With the ease of placing two independent sensors, we can easily differentiate thorax and abdomen breathing patterns to identify obstructive sleep apnea (OSA) by its thoracoabdominal asynchrony. 

We have estimated two key respiratory parameters breath rate (BR) and lung volume (LV) that can identify various respiratory dynamics, and also compared with reference signals derived from airflow pneumotach and calibrated chest belts. Additionally we can also get accurate heart rate variation (HRV) features from NCS that are compared with derived features from reference instruments.

## File organization
This folder orgainization is not in terms of code, but in terms of progress of the project. Initially individual components were focused (following are listed in the order oldest to current projects):
* **Motion detection in sleep [1]**
  * NCS with synchronized external ECG heartbeat waveform (reference instrument only added towards the end of this work, so the paper does not have refernce HR during motion corrected phase)
  * So far only Matlab codes are added, earlier Labview codes are in **EcgNcsCorrelation**.
* **Respiratory (Normal Breathing)**
  * A new reference: Hexoskin smart shirt is used for majority of this work. (https://www.hexoskin.com/)
  * Reference measurements from Hexoskin: ECG heart, thorax and abdomen chest belts respiration with calibrated lung volume estimate.
  * Focusing on normal breathing peak detection using modified Automated Multiscale Peak Detection (AMPD) algorithm, that is mostly automated, with no manually tuned threshold parameters.
  * Estimating lung volume by volume calibration, 
  * HR and BR estimation following the peak detection.
* **Respiratory (Breath Pattern, Coughing, Speaking) [2]**
  * Performing data collection and analysis with different breathing conditions (simulating various apnea and respiratory disorders), coughing and speaking.
  * Data analysis is similar coding as previous, but these abnormal breathing conditions require some more changes in the HR, BR estimation post-processing, as the peak-detection is prone to error.
  * Added Labview codes (improved versions are in recent folders).
* **Attention Test**
  * It includes a psychological test (Mackworth clock test)to detect attention and vigilance.
  * Aim is to compare HRV feature variation during relaxed and attention phases between NCS and reference instrument.
  * The test is written on PsyToolkit free-available online software (https://psytoolkit.org).
  * The script is based on https://www.psytoolkit.org/acknowledgements.html.
* **Mass Study 2019**
  * This project covers codes from all the previous projects with improved algorithm for peak detection etc.
  * The aim is to test the system for multiple people in multiple postures, hence many updates were needed:
    * Reference instrument is updated to Biopac system (https://www.biopac.com/).
    * Codes are automated as much as possible, with no tuning needed (once parameters are set) from person-to-person.
    * refer to the Readme of this for further details regarding hardware and code structure. 
 


## References
1. P. Sharma and E. C. Kan, “Sleep scoring with a UHF RFID tag by near field coherent sensing,” in IEEE MTT-S Int. Microw. Symp. Dig., 2018, pp. 1419–1422. (https://doi.org/10.1109/MWSYM.2018.8439216)
2. P. Sharma, X. Hui and E. C. Kan, "A wearable RF sensor for monitoring respiratory patterns," in IEEE Engineering in Medicine and Biology Society (EMBC), 2019 (Accepted).
3. X. Hui and E. C. Kan, “Monitoring vital signs over multiplexed radio by near-field coherent sensing,” Nat. Electron., vol. 1, pp. 74–78, 2018. (https://www.nature.com/articles/s41928-017-0001-0)
