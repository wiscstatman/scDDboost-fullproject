#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import mpcontrol as mp
import os
import multiprocessing
import subprocess
import os


def RunD3E(input_dir, output_dir, label1, label2):
    args = ["python", "/ua/xiuyu/D3E-master/D3ECmd.py", input_dir, output_dir, label1, label2, "-m","1","-t",
        "0","-z","0","-n","0","-v"]
    subprocess.run(args)
    
def RunD3E_star(args):
    print(args[0])
    RunD3E(*args)
    
def main(argv):
    """Main program.

    Arguments:
    argv -- command line arguments
    """
    num_cores = 8


    mp_init = mp._multiprocessing_worker_init(os.getpid())
    pool = multiprocessing.Pool(num_cores, mp_init)
    ARGS = os.listdir("/ua/xiuyu/D3E_out/")
    
    label1 = "Gad2"
    label2 = "Cux2"
    proc_args = []
    for file in ARGS:
        input_dir = "/ua/xiuyu/D3E_out/" + file
        INDEX = file.split("_")[2].split(".")[0]
        output_dir = "/ua/xiuyu/D3E_res/res_" + str(INDEX) + ".txt"
        proc_args.append([input_dir,output_dir,label1,label2])

    pool.map(RunD3E_star, proc_args)
    
if __name__ == '__main__':
    main(sys.argv)

