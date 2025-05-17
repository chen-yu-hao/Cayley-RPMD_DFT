import numpy as np
import matplotlib.pyplot as plt
# from PES import *
# initialize_potential()
import math
from scipy.constants import physical_constants
atomicMass = {
    "Mu":    0.113,
    "H":     1.00782503207,
    "D":     2.0141017778,
    "T":     3.0160492777,
    "He3":   3.0160293191,
    "He":    4.00260325415,
    "Li6":   6.015122795,
    "Li":    7.01600455,
    "Be":    9.0121822,
    "B5":    10.0129370,
    "B":     11.0093054,
    "C":     12.0,
    "C13":   13.0033548378,
    "N":     14.0030740048,
    "N15":   15.0001088982,
    "O":     15.99491461956,
    "O18":   17.9991610,
    "F":     18.99840322,
    "Ne":    19.9924401754,
    "Ne22":  21.991385114,
    "Na":    22.9897692809,
    "Mg":    23.985041700,
    "Mg25":  24.98583692,
    "Mg26":  25.982592929,
    "Al":    26.98153863,
    "Si":    27.9769265325,
    "Si29":  28.976494700,
    "Si30":  29.97377017,
    "P":     30.97376163,
    "S":     31.97207100,
    "S34":   33.96786690,
    "Cl":    34.96885268,
    "Cl37":  36.96590259,
    "Ar":    39.9623831225,
    "K":     38.96370668,
    "K41":   40.96182576,
    "Ca":    39.96259098,
    "Sc":    44.9559119,
    "Ti":    47.9479463,
    "V":     50.9439595,
    "Cr50":  49.9460442,
    "Cr":    51.9405075,
    "Cr53":  52.9406494,
    "Cr54":  53.9388804,
    "Mn":    54.9380451,
    "Fe54":  53.9396105,
    "Fe":    55.9349375,
    "Fe57":  56.9353940,
    "Co":    58.9331950,
    "Ni":    57.9353429,
    "Cu":    62.9295975,
    "Cu65":  64.9277895,
    "Zn64":  63.9291422,
    "Zn66":  65.9260334,
    "Zn67":  66.9271273,
    "Zn68":  67.9248442,
    "Ga":    68.9255736,
    "Ga71":  70.9247013,
    "Ge70":  69.9242474,       
    "Ge72":  71.9220758,
    "Ge73":  72.9234589,
    "Ge74":  73.9211778,
    "Ge76":  75.9214026,
    "As":    74.9215965,
    "Se74":  73.9224764,
    "Se76":  75.9192136,
    "Se77":  76.9199140,
    "Se78":  77.9173091,
    "Se80":  79.9165213,
    "Se82":  81.9166994,
    "Br":    78.9183371,
    "Br81":  80.9162906,
    "Kr":    83.911507,
    "Rb":    84.911789738,
    "Rb87":  86.909180527,
    "Sr":    87.9056121,
    "Y":     88.9058483,
    "Zr":    89.9047044,
    "Zr91":  90.9056458,
    "Zr92":  91.9050408,
    "Zr94":  93.9063152,
    "Zr96":  95.9082734,
}

global reactant1_atoms,reactant2_atoms,massfrac,Rinf0,number_of_transition_states,forming_bonds,breaking_bonds,forming_bond_lengths,breaking_bond_lengths,number_of_bonds,mass,get_potential

def reactant1_center_of_mass(position):
    # import np as np

    cm = np.zeros(3)
    for n in range(0, len(reactant1_atoms)):
        atom = reactant1_atoms[n]
        cm += massfrac[atom] * position[:, atom]

    return cm

def reactant2_center_of_mass(position):
    # import np as np

    cm = np.zeros(3)
    for n in range(0, len(reactant2_atoms)):
        atom = reactant2_atoms[n]
        cm += massfrac[atom] * position[:, atom]

    return cm
def value_s0(position):

    cm1 = reactant1_center_of_mass(position)
    cm2 = reactant2_center_of_mass(position)

    Rx = cm2[0] - cm1[0]
    Ry = cm2[1] - cm1[1]
    Rz = cm2[2] - cm1[2]
    R = np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

    s0 = Rinf0 - R

    return s0

def evaluate_all(position, Natoms):
    # import np as np

    values = np.zeros(breaking_bond_lengths.shape)

    m, n, atom1, atom2 = 0, 0, 0, 0
    Rx, Ry, Rz, R = 0, 0, 0, 0

    for n in range(0, number_of_transition_states):
        for m in range(0, number_of_bonds):
            # Forming bond
            atom1 = int(forming_bonds[n, m, 0])
            atom2 = int(forming_bonds[n, m, 1])
            Rx = position[0, atom1] - position[0, atom2]
            Ry = position[1, atom1] - position[1, atom2]
            Rz = position[2, atom1] - position[2, atom2]
            R = np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            # print(R)
            values[n, m] += (forming_bond_lengths[n, m] - R) #############################################################waiting for changing 
            
            # Breaking bond
            atom1 = int(breaking_bonds[n, m, 0])
            atom2 = int(breaking_bonds[n, m, 1])
            Rx = position[0, atom1] - position[0, atom2]
            Ry = position[1, atom1] - position[1, atom2]
            Rz = position[2, atom1] - position[2, atom2]
            R = np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            values[n, m] -= (breaking_bond_lengths[n, m] - R)
            # print(values)


    return values

def value_s1(position,Natoms):
    values=evaluate_all(position, Natoms)

    s1 = np.max(np.min(values, axis=1))
    return s1

def s1_gradient(position, Natoms):
    # import np as np

    ds1 = np.zeros(position.shape)

    values = evaluate_all(position, Natoms)
    n = np.argmin(np.max(values, axis=1))
    m = np.argmax(values[n, :])
    atom1=forming_bonds[n,m,0]
    atom2=forming_bonds[n,m,1]
    Rx = position[0,atom1] - position[0,atom2]
    Ry = position[1,atom1] - position[1,atom2]
    Rz = position[2,atom1] - position[2,atom2]
    Rinv = 1.0 / np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
    ds1[0,atom1] = ds1[0,atom1] - Rx * Rinv
    ds1[1,atom1] = ds1[1,atom1] - Ry * Rinv
    ds1[2,atom1] = ds1[2,atom1] - Rz * Rinv
    ds1[0,atom2] = ds1[0,atom2] + Rx * Rinv
    ds1[1,atom2] = ds1[1,atom2] + Ry * Rinv
    ds1[2,atom2] = ds1[2,atom2] + Rz * Rinv

    atom1=breaking_bonds[n,m,0]
    atom2=breaking_bonds[n,m,1]
    Rx = position[0,atom1] - position[0,atom2]
    Ry = position[1,atom1] - position[1,atom2]
    Rz = position[2,atom1] - position[2,atom2]
    Rinv = 1.0 / np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
    ds1[0,atom1] = ds1[0,atom1] + Rx * Rinv
    ds1[1,atom1] = ds1[1,atom1] + Ry * Rinv
    ds1[2,atom1] = ds1[2,atom1] + Rz * Rinv
    ds1[0,atom2] = ds1[0,atom2] - Rx * Rinv
    ds1[1,atom2] = ds1[1,atom2] - Ry * Rinv
    ds1[2,atom2] = ds1[2,atom2] - Rz * Rinv
    return ds1

def s0_gradient(position):
    cm1 = reactant1_center_of_mass(position)
    cm2 = reactant2_center_of_mass(position)

    Rx = cm2[0] - cm1[0]
    Ry = cm2[1] - cm1[1]
    Rz = cm2[2] - cm1[2]
    Rinv = 1.0/np.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

    ds0=np.zeros(position.shape)

    for n in range(0, len(reactant1_atoms)):
        atom = reactant1_atoms[n]
        ds0[0,atom] = Rx * Rinv * massfrac[(atom)]
        ds0[1,atom] = Ry * Rinv * massfrac[(atom)]
        ds0[2,atom] = Rz * Rinv * massfrac[(atom)]
    for n in range(0, len(reactant2_atoms)):
        atom = reactant2_atoms[n]
        ds0[0,atom] = -Rx * Rinv * massfrac[(atom)]
        ds0[1,atom] = -Ry * Rinv * massfrac[(atom)]
        ds0[2,atom] = -Rz * Rinv * massfrac[(atom)]
    return ds0

def xi(position):
    Natoms=position.shape[1]
    s0=value_s0(position)
    s1=value_s1(position,Natoms)
    return s0 / (s0 - s1)
def dxi(position):
    Natoms=position.shape[1]
    s0=value_s0(position)
    s1=value_s1(position,Natoms)
    ds0=s0_gradient(position)
    ds1=s1_gradient(position,Natoms)
    return (s0 * ds1 - s1 * ds0) / ((s0 - s1) * (s0 - s1))

def SHAKE_xi(position,xi_current):
    Natoms=position.shape[1]
    s0=value_s0(position)
    s1=value_s1(position,Natoms)
    return xi_current * s1 + (1 - xi_current) * s0

def SHAKE_dxi(position,xi_current):
    Natoms=position.shape[1]
    s0=value_s0(position)
    s1=value_s1(position,Natoms)
    ds0=s0_gradient(position)
    ds1=s1_gradient(position,Natoms)
    return xi_current * ds1 + (1 - xi_current) * ds0


#calculate Hessian from function dvdq
def cal_HESSIAN(position):
    Natoms=position.shape[0]
    dx=1e-4
    HESSIAN=np.zeros([Natoms,3,Natoms,3])
    for i in range(Natoms):
        for j in range(3):
            q_1=np.copy(position)
            q_2=np.copy(position)
            q_1[i,j]=position[i,j]+dx
            q_2[i,j]=position[i,j]-dx
            HESSIAN[i,j,:,:]=((get_potential(q_1.T)[1][:,:,0].T-
                               get_potential(q_2.T)[1][:,:,0].T)/(2*dx))
    
    return HESSIAN

#calculate Hessian from function dvdq
def cal_HESSIAN_T(position):
    Natoms=position.shape[1]
    dx=1e-4
    HESSIAN=np.zeros([3,Natoms,3,Natoms])
    for i in range(3):
        for j in range(Natoms):
            q_1=np.copy(position)
            q_2=np.copy(position)
            q_1[i,j]=position[i,j]+dx
            q_2[i,j]=position[i,j]-dx
            HESSIAN[i,j,:,:]=((get_potential(q_1)[1][:,:,0]-
                               get_potential(q_2)[1][:,:,0])/(2*dx))
    
    return HESSIAN

def constrain_to_dividing_surface(position,momentum,xi_current,dt,maxiter=100000):
    mult=0
    qctemp=np.zeros(position.shape)
    for iter in range(maxiter):
        coeff = mult*dt**2
        for i in range(3):
            for j in range(qctemp.shape[1]):
                qctemp[i,j]=position[i,j]+coeff * SHAKE_dxi(position,xi_current)[i,j]
        sigma=SHAKE_xi(qctemp,xi_current)
        dsigma = 0.0
        dxi_new=SHAKE_dxi(qctemp,xi_current)
        for i in range(3):
            for j in range(qctemp.shape[1]):
                dsigma = dsigma + dt**2*dxi_new[i,j] * SHAKE_dxi(position,xi_current)[i,j]
        dx=sigma/dsigma
        mult=mult - dx
        if (abs(dx) < 1e-6) or abs(sigma)< 1e-8:
            break
        if iter==maxiter-1:
            print("Warning: maximum number of iterations reached")

        # print(SHAKE_xi(qctemp,xi_current),dx,dsigma,mult,coeff)
    momentum=momentum+mult*dt*SHAKE_dxi(position,xi_current)
    return qctemp,momentum

def constrain_momentum_to_dividing_surface(momentum,value_dxi):
    coeff1=np.sum(value_dxi*momentum)
    coeff2=np.sum(value_dxi*value_dxi)
    value_lambda=-coeff1/coeff2
    momentum=momentum+value_lambda*value_dxi
    return momentum

def reactants(
  atoms = ['O', 'O', 'C', 'H'],
  reactant1Atoms = [1,4],
  reactant2Atoms = [2,3],
  Rinf = (15.0,"angstrom"),
):
    global reactant1_atoms,reactant2_atoms,Rinf0,massfrac,mass,atoms_label
    atoms_label=atoms
    reactant1_atoms=np.array(reactant1Atoms)-1
    reactant2_atoms=np.array(reactant2Atoms)-1
    mass=[]
    for i in atoms:
        mass.append(atomicMass[i])
    Rinf0=Rinf[0]
    if Rinf[1]=="angstorm":
        Rinf0=Rinf0*	1.8897259886
    massfrac=np.zeros(len(atoms))
    masss1=0
    for i in reactant1_atoms:
        masss1+=mass[i]
    for i in reactant1_atoms:
        massfrac[i]=mass[i]/masss1
    masss2=0
    for i in reactant2_atoms:
        masss2+=mass[i]
    for i in reactant2_atoms:
        massfrac[i]=mass[i]/masss2

        mass=np.array(mass)
    

def cal_MEP(position,xi_list,maxin=20000,dt=10,maxiter=10000,dxi=0.0001,dt_xi=1):
    E=[]
    position_list=[]
    # number_xi=(xi_list[0]-xi_list[1])//dxi
    # print(number_xi)
    dt0=dt
    for xi_current in xi_list:
        print(xi_current,end="\r")
        # for i in range(abs(int(number_xi))):
        #     # print(xi_current+(number_xi-(xi_list[0]-xi_list[1])/abs((xi_list[0]-xi_list[1]))*i)*dxi)
        #     position,_= constrain_to_dividing_surface(position,0,xi_current+(number_xi-(xi_list[0]-xi_list[1])/abs((xi_list[0]-xi_list[1]))*i)*dxi,dt_xi,maxiter)
        count=1
        
        # dvdx=1
        position0=position
        dt=dt0
        E1=get_potential(position)[0][0]
        E0=-10000.0
        dtnum=0

        while ((E1-E0>=0) or (E1-E0<-1e-7)) :
                
            position0=position
            E0,dvdx,info=get_potential(position)
            dvdx=dvdx[:,:,0]

            # position_HESS=cal_HESSIAN_T(position).reshape([3*len(mass),3*len(mass)])
            # try:
            #     HESS_I=np.linalg.inv(position_HESS)
            # except:
            #     HESS_I=np.linalg.pinv(position_HESS)
            # dvdx=np.dot(HESS_I,dvdx.reshape(-1)).reshape([3,-1])
            # dvdx=dvdx*np.random.random([3,len(mass)])

            dvdx=constrain_momentum_to_dividing_surface(dvdx,SHAKE_dxi(position,xi_current))
            if (np.abs(dvdx*dt).max())<5e-6:
                if dt<1e-5:
                    position=position-(5e-6*dvdx/np.abs(dvdx).max()+1e-5*(0.5-np.random.random(position.shape)))
                else:
                    position=position-5e-6*dvdx/np.abs(dvdx).max()
            elif (np.abs(dvdx*dt).max())>5e-4:
                position=position-5e-4*dvdx/np.abs(dvdx).max()
            else:
                position=position-dvdx*dt
            # print(np.abs(dvdx*dt).max())
            
            position,_= constrain_to_dividing_surface(position,0,xi_current,dt_xi,maxiter)
            count+=1

            if count>maxin:
                print("maxin reached",xi_current,E0-E1,dt)
                break

            E1,dvdx,info=get_potential(position)
            E1=E1[0]
            E0=E0[0]
            if E0-E1<0:
                # print(E0,E1-E0,dt,xi_current,"true")    
                # print((E1-E0>0) or (E1-E0<-1e-6))
                dtnum+=1
                # if dtnum==10:
                dt=dt/2
                #     dtnum=0
                E1=E0
                position=position0
                
                 
            else:
                # print(E0,E1-E0,dt,xi_current)
                dt=dt0
                # print((E1-E0>0) or (E1-E0<-1e-6))
        # if type(E0)==float:
        if (type(E0)==np.ndarray):
            E.append(E0[0])
        else:

        # else:
            E.append(E0)
        

        position_list.append(position)
    return E,position_list

def TS_optimiser(TS,grid_max=1e-6,v=0.1,randomc="True"):
    if randomc:
        TS_0=TS
        g=1
        g0=2
        count=1
        g_counter=0
        import random
        while g>grid_max:
            TS_list=[]
            g=[]

            g.append(get_potential(TS_0)[1].max())
            TS_list.append(TS_0)
            TS_HESS=np.matrix(cal_HESSIAN_T(TS_0).reshape([TS.shape[1]*3,TS.shape[1]*3]))
            try:
                inverse=np.linalg.inv(TS_HESS)
            except:
                inverse=np.linalg.pinv(TS_HESS)
            grid=get_potential(TS)[1][:,:,0]
            grid_HESS=np.dot(TS_HESS.I,get_potential(TS)[1].reshape([TS.shape[1]*3,1])).reshape(3,TS.shape[1])
            grid0=get_potential(TS)[1].max()
            for i in range(1000):
                TS_0=TS_0-np.array(grid_HESS)*v*(1*np.random.random([3,TS.shape[1]])-0.5)*random.random()
                # TS_0=TS_0-v*(1*np.random.random([3,5])-0.5)
                # TS_0=TS_0-np.dot(inverse,grid0*v*(1*np.random.random(15)-0.5)).reshape([3,-1])
                # TS_0=TS_0-np.array(grid_HESS)*v*random.random()

                TS_list.append(TS_0)
                g.append(get_potential(TS_0)[1].max())
            g=np.array(g)
            TS_0=TS_list[np.argmin(g)]
            
            g=get_potential(TS_0)[1].max()
            if g==g0:
                g_counter+=1
                if g_counter==20:
                    v=v/2
                    g_counter=0
            else:
                g_counter=0
            print("NO. %i g: %e            %e"%(count,g,v),end="\r")
            count+=1
            g0=g
        TS_0=np.array(TS_0)
        # print("transition state frequency:")
        # print(" ".join([str(i) for i in cal_freq_cm(TS_0)]))
        return TS_0
    else:
        TS_0=TS
        g=1
        g0=2
        count=1
        g_counter=0
        import random
        while g>grid_max:
            TS_list=[]
            g=[]

            g.append(get_potential(TS_0)[1].max())
            TS_list.append(TS_0)
            TS_HESS=np.matrix(cal_HESSIAN_T(TS_0).reshape([TS.shape[1]*3,TS.shape[1]*3]))
            try:
                inverse=np.linalg.inv(TS_HESS)
            except:
                inverse=np.linalg.pinv(TS_HESS)
            grid=get_potential(TS)[1][:,:,0]
            grid_HESS=np.dot(TS_HESS.I,get_potential(TS)[1].reshape([TS.shape[1]*3,1])).reshape(3,TS.shape[1])
            grid0=get_potential(TS)[1].max()
            if np.abs(np.array(grid_HESS)*v).max()>1e-3:
                grid_HESS=np.array(grid_HESS)/np.abs(np.array(grid_HESS)*v).max()*1e-3
            TS_0=TS_0-np.array(grid_HESS)*v
            TS_list.append(TS_0)
            g.append(get_potential(TS_0)[1].max())
            g=np.array(g)
            TS_0=TS_list[np.argmin(g)]
            g=get_potential(TS_0)[1].max()
            if g==g0:
                g_counter+=1
                if g_counter==20:
                    v=v/2
                    g_counter=0
            else:
                g_counter=0
            print("NO. %i g: %e            %e"%(count,g,v),end="\r")
            count+=1
            g0=g
        TS_0=np.array(TS_0)
        # print("transition state frequency:")
        # print(" ".join([str(i) for i in cal_freq_cm(TS_0)]))
        return TS_0

def transitionState(
  geometry = (
   [[ 0.00353956,  0.0, 1.34353059],
    [ 2.26544198,  0.0, 1.90551190],
    [ 1.09957058,  0.0, 1.85551665],
    [ 0.00324387,  0.0, 0.00412347]],
   "angstrom",
  ),
  formingBonds = [(1,3)],
  breakingBonds = [(4,1)],
):
    global TS,forming_bond_lengths,breaking_bond_lengths,forming_bonds,breaking_bonds
    TS=np.array(geometry[0])
    if geometry[1]=="angstorm":
        TS=TS*	1.8897259886
    forming_bond_lengths=[[]]
    breaking_bond_lengths=[[]]
    forming_bonds=[]
    breaking_bonds=[]
    forming_bonds.append (np.array(formingBonds )-1)
    breaking_bonds.append(np.array(breakingBonds)-1)
    for i in formingBonds:
        dis=np.sqrt(np.sum((TS[i[0]-1]-TS[i[1]-1])**2))
        forming_bond_lengths[0].append(dis)
    for i in breakingBonds:
        dis=np.sqrt(np.sum((TS[i[0]-1]-TS[i[1]-1])**2))
        breaking_bond_lengths[0].append(dis)
def xyz_print(xyz):
    global atoms_label
    xyz=xyz*0.529177
    for i in range(xyz.shape[1]):
        print(atoms_label[i],end=" ")
        print(xyz[0,i],end=" ")
        print(xyz[1,i],end=" ")
        print(xyz[2,i],end=" \n")
def equivalentTransitionState( 
    formingBonds,
    breakingBonds ):
    
    global TS,forming_bond_lengths,breaking_bond_lengths,forming_bonds,breaking_bonds
    # forming_bond_lengths.append([])
    # breaking_bond_lengths.append([])
    forming_bonds. append(np.array(formingBonds )-1)
    breaking_bonds.append(np.array(breakingBonds)-1)
    # for i in formingBonds:
    #     dis=np.sqrt(np.sum((TS[i[0]-1]-TS[i[1]-1])**2))
    #     forming_bond_lengths[-1].append(dis)
    # for i in breakingBonds:
    #     dis=np.sqrt(np.sum((TS[i[0]-1]-TS[i[1]-1])**2))
    #     breaking_bond_lengths[-1].append(dis)
    Natoms=TS.T.shape[1]
    mapping = {}
    for bond1, bond2 in zip(forming_bonds[0]+1, formingBonds):
        atom11, atom12 = bond1
        atom21, atom22 = bond2
        if atom11 in mapping and mapping[atom11] != atom21:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom11, mapping[atom11], atom21))
        elif atom21 in mapping and mapping[atom21] != atom11:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom21, mapping[atom21], atom11))
        else:
            mapping[atom11] = atom21
            mapping[atom21] = atom11
        if atom12 in mapping and mapping[atom12] != atom22:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom12, mapping[atom12], atom22))
        elif atom22 in mapping and mapping[atom22] != atom12:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom22, mapping[atom22], atom12))
        else:
            mapping[atom12] = atom22
            mapping[atom22] = atom12
    for bond1, bond2 in zip(breaking_bonds[0]+1, breakingBonds):
        atom11, atom12 = bond1
        atom21, atom22 = bond2
        if atom11 in mapping and mapping[atom11] != atom21:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom11, mapping[atom11], atom21))
        elif atom21 in mapping and mapping[atom21] != atom11:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom21, mapping[atom21], atom11))
        else:
            mapping[atom11] = atom21
            mapping[atom21] = atom11
        if atom12 in mapping and mapping[atom12] != atom22:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom12, mapping[atom12], atom22))
        elif atom22 in mapping and mapping[atom22] != atom12:
            raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom22, mapping[atom22], atom12))
        else:
            mapping[atom12] = atom22
            mapping[atom22] = atom12
    for atom in range(1, Natoms+1):
        if atom not in mapping:
            mapping[atom] = atom
    
    geometry = np.zeros_like(TS.T)
    # print(geometry.shape,TS.shape)
    for atom in range(Natoms):
        geometry[:,atom] = TS.T[:,mapping[atom+1]-1]
    # geometry = geometry.T
    
    formingBonds = np.array(formingBonds, int)
    breakingBonds = np.array(breakingBonds, int)
    
    
    Nforming_bonds = formingBonds.shape[0]
    Nbreaking_bonds =breakingBonds.shape[0]
    
    # print(formingBonds)
    # print(breakingBonds)
    formingBondLengths = np.empty(Nforming_bonds)
    breakingBondLengths = np.empty(Nbreaking_bonds)
    
    for m in range(Nforming_bonds):
        atom1 = formingBonds[m,0] - 1
        atom2 = formingBonds[m,1] - 1
        Rx = geometry[0,atom1] - geometry[0,atom2]
        Ry = geometry[1,atom1] - geometry[1,atom2]
        Rz = geometry[2,atom1] - geometry[2,atom2]
        R = math.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        formingBondLengths[m] = R
        # print(R,atom1,atom2,geometry.T)
    forming_bond_lengths.append(formingBondLengths)
    
    for m in range(Nbreaking_bonds):
        atom1 = breakingBonds[m,0] - 1
        atom2 = breakingBonds[m,1] - 1
        Rx = geometry[0,atom1] - geometry[0,atom2]
        Ry = geometry[1,atom1] - geometry[1,atom2]
        Rz = geometry[2,atom1] - geometry[2,atom2]
        R = math.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        breakingBondLengths[m] = R
    breaking_bond_lengths.append(breakingBondLengths)





        

def init(get_potential0):
    global TS,forming_bond_lengths,breaking_bond_lengths,forming_bonds,breaking_bonds,number_of_transition_states,number_of_bonds
    forming_bond_lengths =np.array(forming_bond_lengths )
    breaking_bond_lengths=np.array(breaking_bond_lengths)
    forming_bonds=        np.array(forming_bonds        )
    breaking_bonds=       np.array(breaking_bonds       )
    number_of_transition_states=forming_bonds.shape[0]
    number_of_bonds=forming_bonds.shape[1]
    global get_potential
    get_potential=get_potential0

def cal_freq_cm(position):
    

    E_h = physical_constants["Hartree energy"][0]
    a_0 = physical_constants["Bohr radius"][0]
    N_A = physical_constants["Avogadro constant"][0]
    c_0 = physical_constants["speed of light in vacuum"][0]
    e_c = physical_constants["elementary charge"][0]
    e_0 = physical_constants["electric constant"][0]
    mu_0 = physical_constants["mag. constant"][0]
    Natoms=position.shape[1]
    position=position.T
    HESSIAN=cal_HESSIAN(position)
#     print(HESSIAN.shape)
    for i in range(Natoms):
        for j in range(Natoms):
            HESSIAN[i,:,j,:]=HESSIAN[i,:,j,:]/np.sqrt(mass[i]*mass[j])

    HESSIAN = HESSIAN.reshape(3 * Natoms, 3 * Natoms)

    e, q = np.linalg.eigh(HESSIAN)
    center_coord = (position * mass[:,None]).sum(axis=0) / mass.sum()
    centered_coord = position - center_coord


    rot_tmp = np.zeros((Natoms, 3, 3))
    rot_tmp[:, 0, 0] =   centered_coord[:, 1]**2 + centered_coord[:, 2]**2
    rot_tmp[:, 1, 1] =   centered_coord[:, 2]**2 + centered_coord[:, 0]**2
    rot_tmp[:, 2, 2] =   centered_coord[:, 0]**2 + centered_coord[:, 1]**2
    rot_tmp[:, 0, 1] = - centered_coord[:, 0] * centered_coord[:, 1]
    rot_tmp[:, 1, 2] = - centered_coord[:, 1] * centered_coord[:, 2]
    rot_tmp[:, 2, 0] = - centered_coord[:, 2] * centered_coord[:, 0]
    rot_tmp[:, 1, 0] = - centered_coord[:, 0] * centered_coord[:, 1]
    rot_tmp[:, 2, 1] = - centered_coord[:, 1] * centered_coord[:, 2]
    rot_tmp[:, 0, 2] = - centered_coord[:, 2] * centered_coord[:, 0]
    rot_tmp = (rot_tmp * mass[:, None, None]).sum(axis=0)
    _, rot_eig = np.linalg.eigh(rot_tmp)

    rot_coord = np.einsum("At, ts, rw -> Asrw", centered_coord, rot_eig, rot_eig)
    proj_scr = np.zeros((Natoms, 3, 6))
    proj_scr[:, (0, 1, 2), (0, 1, 2)] = 1
    proj_scr[:, :, 3] = (rot_coord[:, 1, :, 2] - rot_coord[:, 2, :, 1])
    proj_scr[:, :, 4] = (rot_coord[:, 2, :, 0] - rot_coord[:, 0, :, 2])
    proj_scr[:, :, 5] = (rot_coord[:, 0, :, 1] - rot_coord[:, 1, :, 0])
    proj_scr *= np.sqrt(mass)[:, None, None]
    proj_scr.shape = (-1, 6)
    proj_scr /= np.linalg.norm(proj_scr, axis=0)
    e_tr, _ = np.linalg.eigh(np.dot(proj_scr.T ,np.dot(HESSIAN, proj_scr)))
    proj_inv = np.zeros((Natoms * 3, Natoms * 3))
    proj_inv[:, :6] = proj_scr
    cur = 6
    for i in range(0, Natoms * 3):
        vec_i = np.einsum("Ai, i -> A", proj_inv[:, :cur], proj_inv[i, :cur])
        vec_i[i] -= 1
        if np.linalg.norm(vec_i) > 1e-11:
            proj_inv[:, cur] = vec_i / np.linalg.norm(vec_i)
            cur += 1
        if cur >= Natoms * 3:
            break
    proj_inv = proj_inv[:, 6:]
    e, q = np.linalg.eigh(proj_inv.T @ HESSIAN @ proj_inv)
    freq_cm_1 = np.sqrt(np.abs(e * E_h * 1000 * N_A / a_0**2)) / (2 * np.pi * c_0 * 100) * ((e > 0) * 2 - 1)
    return freq_cm_1

def cal_Hessian_q(position):
    

    E_h = physical_constants["Hartree energy"][0]
    a_0 = physical_constants["Bohr radius"][0]
    N_A = physical_constants["Avogadro constant"][0]
    c_0 = physical_constants["speed of light in vacuum"][0]
    e_c = physical_constants["elementary charge"][0]
    e_0 = physical_constants["electric constant"][0]
    mu_0 = physical_constants["mag. constant"][0]
    Natoms=position.shape[1]
    position=position.T
    HESSIAN=cal_HESSIAN(position)
#     print(HESSIAN.shape)
    for i in range(Natoms):
        for j in range(Natoms):
            HESSIAN[i,:,j,:]=HESSIAN[i,:,j,:]/np.sqrt(mass[i]*mass[j])

    HESSIAN = HESSIAN.reshape(3 * Natoms, 3 * Natoms)

    e, q = np.linalg.eigh(HESSIAN)
    center_coord = (position * mass[:,None]).sum(axis=0) / mass.sum()
    centered_coord = position - center_coord


    rot_tmp = np.zeros((Natoms, 3, 3))
    rot_tmp[:, 0, 0] =   centered_coord[:, 1]**2 + centered_coord[:, 2]**2
    rot_tmp[:, 1, 1] =   centered_coord[:, 2]**2 + centered_coord[:, 0]**2
    rot_tmp[:, 2, 2] =   centered_coord[:, 0]**2 + centered_coord[:, 1]**2
    rot_tmp[:, 0, 1] = - centered_coord[:, 0] * centered_coord[:, 1]
    rot_tmp[:, 1, 2] = - centered_coord[:, 1] * centered_coord[:, 2]
    rot_tmp[:, 2, 0] = - centered_coord[:, 2] * centered_coord[:, 0]
    rot_tmp[:, 1, 0] = - centered_coord[:, 0] * centered_coord[:, 1]
    rot_tmp[:, 2, 1] = - centered_coord[:, 1] * centered_coord[:, 2]
    rot_tmp[:, 0, 2] = - centered_coord[:, 2] * centered_coord[:, 0]
    rot_tmp = (rot_tmp * mass[:, None, None]).sum(axis=0)
    _, rot_eig = np.linalg.eigh(rot_tmp)

    rot_coord = np.einsum("At, ts, rw -> Asrw", centered_coord, rot_eig, rot_eig)
    proj_scr = np.zeros((Natoms, 3, 6))
    proj_scr[:, (0, 1, 2), (0, 1, 2)] = 1
    proj_scr[:, :, 3] = (rot_coord[:, 1, :, 2] - rot_coord[:, 2, :, 1])
    proj_scr[:, :, 4] = (rot_coord[:, 2, :, 0] - rot_coord[:, 0, :, 2])
    proj_scr[:, :, 5] = (rot_coord[:, 0, :, 1] - rot_coord[:, 1, :, 0])
    proj_scr *= np.sqrt(mass)[:, None, None]
    proj_scr.shape = (-1, 6)
    proj_scr /= np.linalg.norm(proj_scr, axis=0)
    e_tr, _ = np.linalg.eigh(np.dot(proj_scr.T ,np.dot(HESSIAN, proj_scr)))
    proj_inv = np.zeros((Natoms * 3, Natoms * 3))
    proj_inv[:, :6] = proj_scr
    cur = 6
    for i in range(0, Natoms * 3):
        vec_i = np.einsum("Ai, i -> A", proj_inv[:, :cur], proj_inv[i, :cur])
        vec_i[i] -= 1
        if np.linalg.norm(vec_i) > 1e-11:
            proj_inv[:, cur] = vec_i / np.linalg.norm(vec_i)
            cur += 1
        if cur >= Natoms * 3:
            break
    proj_inv = proj_inv[:, 6:]
    e, q = np.linalg.eigh(proj_inv.T @ HESSIAN @ proj_inv)
    freq_cm_1 = np.sqrt(np.abs(e * E_h * 1000 * N_A / a_0**2)) / (2 * np.pi * c_0 * 100) * ((e > 0) * 2 - 1)
    q_unnormed = np.einsum("AtQ, A -> AtQ", (proj_inv @ q).reshape(Natoms, 3, (proj_inv @ q).shape[-1]), 1 / np.sqrt(mass))
    q_unnormed = q_unnormed.reshape(-1, q_unnormed.shape[-1])
    q_normed = q_unnormed / np.linalg.norm(q_unnormed, axis=0)
    return q_normed,HESSIAN
