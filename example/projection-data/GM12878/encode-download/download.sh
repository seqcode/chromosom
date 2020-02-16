
xargs -L 1 curl -O -L < file-list1.txt
mv \?type\=Experiment\&files.output_type\=replicated+peaks\&files.file_type\=bed+narrowPeak histone-marks.data

xargs -L 1 curl -O -L < file-list2.txt
mv \?type\=Experiment\&files.output_type\=conservative+IDR+thresholded+peaks\&files.file_type\=bed+narrowPeak\&files.assembly\=hg19 tf.data

xargs -L 1 curl -O -L < file-list3.txt
mv \?type\=Experiment\&files.output_type\=peaks\&files.file_type\=bed+narrowPeak\&files.assembly\=hg19\&files.lab.title\=ENCODE+Processing+Pipeline\&files.status\=released dnase.data

gunzip *gz
