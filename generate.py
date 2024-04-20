N_in = 207
N_hid_0 = 15
N_hid_1 = 80
Natoms=5
alpha=4
with open("test_read.f90","r")as f:
    s=f.read()
s=s.replace("N_in",str(N_in))
s=s.replace("N_hid_0",str(N_hid_0))
s=s.replace("N_hid_1",str(N_hid_1))
s=s.replace("Natoms*(Natoms-1)/2",str(Natoms*(Natoms-1)//2))
s=s.replace("\n     &        ","&\n     &        ")
s=s.replace("return","")
s=s.replace("real(wp)","double precision")
s=s.replace("bemsav(x,m,p)","bemsav(x,p)")
s=s.replace("Alpha",str(alpha))
with open("pipnn.f90","w")as f:
    f.write(s)

# with open("Makefile","a")as f:
#     f.write("""
# PES:
# 	python -m numpy.f2py -m NEW -h NEW.pyf pipnn.f90 --overwrite-signature
#     python -m numpy.f2py -c pipnn.f90 NEW.pyf 
# """)
import os
os.system("python -m numpy.f2py -m NEW -h NEW.pyf pipnn.f90 --overwrite-signature && python -m numpy.f2py -c pipnn.f90 NEW.pyf --compiler=intelem --fcompiler=intelem ")
