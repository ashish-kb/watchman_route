#!bin/bash
#machines available 1-3 5-10 and 12-15
for i in {1..2}
do
     echo "running the script on system $i system 1= cvl05 2= cvl06 and so on"
	str=$( printf '%02d' $i ); #1 -> 01 2> 02
        cvlstr=$str
#$(($i+4))
     ssh -X ashishkb@cvl$str.ece.vt.edu "cd ~/project_code/GTSP$cvlstr/ && matlab -nodesktop -nosplash -r \"script_auto;exit;\"&"
     exit 


done



# nohup yourprocess 1>&2  | tee nohup.out &
# A=$( printf '%02d' $((A+1)) )
# cd ./gtsp_ashish/GTSP/ && matlab -nodesktop -nosplash -r "load('wrkspace49hr.mat');whos counter_struct;exit;"


