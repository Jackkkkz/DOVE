import os
import numpy as np
from multiprocessing import Process
from multiprocessing import Pool
def getpdbname():
    irmsdresult=os.path.join(os.getcwd(),'iRMSD')
    listfile=os.listdir(irmsdresult)
    list1=[]
    for item in listfile:
        if(item[-4:]=='.txt'):
            list1.append(item[0:4])
    return list1
def check_goap(item1):
    goapdecoys=os.path.join(os.getcwd(),"goapdecoy")
    tempdir=os.path.join(goapdecoys,item1)
    if os.path.exists(tempdir):
        listfiles=os.listdir(tempdir)
        if len(listfiles)>20000:
            print('%s already generated'%item1)
            return False
    return True
from ops.os_operation import mkdir
def run_goap(item1):
    decoydataset=os.path.join(os.getcwd(),'decoys')
    goapset=os.path.join(os.getcwd(),'Goap')
    goapdecoys=os.path.join(os.getcwd(),"goapdecoy")
    pathgenerate=os.path.join(decoydataset,item1)
    listtemp=os.listdir(pathgenerate)
    listdecoy=[]
    for item in listtemp:
        if item[0:7]=='complex':
            listdecoy.append(item)
    os.chdir(goapset)
    #Copy files to the running directory
    os.system("cp fort.21_1.61_2 "+pathgenerate+"/fort.21_1.61_2")
    os.system("cp charge_inp.dat "+pathgenerate+"/charge_inp.dat")
    os.system("cp side_geometry.dat "+pathgenerate+"/side_geometry.dat")
    os.system("cp fort.31_g72_noshift5_new "+pathgenerate+"/fort.31_g72_noshift5_new")
    os.system("cp goap "+pathgenerate+"/goap")
    pathgoapdecoy=os.path.join(goapdecoys,item1)
    mkdir(pathgoapdecoy)
    os.chdir(pathgenerate)
    listrun=listdecoy
    if len(listrun)==0:
        print('no decoys avilable')
        return
    print('waiting dealing'+str(len(listrun)))
    file_object = open(str(item1)+'.inp','w')
    try:
        file_object.write(goapset+'\n')
        for item2 in listrun:
            file_object.write(str(item2)+'\n')
    finally:
        file_object.close()
    os.chdir(pathgenerate)
    os.system("./goap<"+str(item1)+".inp")
    allfiles=os.listdir(pathgenerate)
    list2=[]
    for item in allfiles:
        if item[-9:]=='_goap.pdb':
            list2.append(item)
    for item in list2:
        tmp_path=os.path.join(pathgoapdecoy,item)
        os.system("mv "+item+" "+tmp_path)


