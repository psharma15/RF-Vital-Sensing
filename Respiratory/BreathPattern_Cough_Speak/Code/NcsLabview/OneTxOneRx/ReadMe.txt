Refer to: D:\Research\SummerFall17Spring18\CnC\NCS\MotionDetectionSleep\MyExperiment\Labview

Edit 18 July: Major edit - Removed phase unwrapping - as it creates phase jumps - that after filtering, downsampling etc. are difficult to take care of automatically.


The above folder has codes for NCS data read.

ncs_breath == ncs_v4_timer in other files.

ncs_voice script is created, with different filter etc settings for high sample rate data collection. 

As this variation affects filtering data before downsampling, we have:
NCS_AmpPh - Original with fixed filter specifications. (for original ncs_breath)
NCS_AmpPh_FilterControl - With controllable filter. (for ncs_voice)