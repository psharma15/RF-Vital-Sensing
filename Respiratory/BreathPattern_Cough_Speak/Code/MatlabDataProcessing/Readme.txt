CodeStyle1 is originally copied from ... NCS/SleepStudySetup/CodeAndData/Code.

17 Jul 2018:
This is updated version of CodeStyle_Old, with less manual data conversions - leading to changes in NCS reading function etc. This should be the code in the following experiments and data collection.

Since data on 07/12 needs to be redone, both versions are kept. CodeStyle1 would go into archive after that.

First use saveTDMStoMAT.m to convert TDMS to mat. Then use ncsHxCompare to 
first sync NCS and Hexoskin and then get remaining results.
