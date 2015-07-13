#!/usr/bin/env sh
TOOLS=/home/wuxiang/caffe/build/tools
mkdir log
export PATH=/usr/local/cuda-6.5/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-6.5/lib64:$LD_LIBRARY_PATH

#GLOG_logtostderr=0 GLOG_alsologtostderr=1 $TOOLS/caffe train --solver=./ocr_en_set001_solver.prototxt --gpu=1 
GLOG_logtostderr=0 GLOG_alsologtostderr=1 $TOOLS/caffe train --solver=./ocr_en_set001_solver_2.prototxt --gpu=1 --snapshot=./ocr_en_set001_iter_100000.solverstate