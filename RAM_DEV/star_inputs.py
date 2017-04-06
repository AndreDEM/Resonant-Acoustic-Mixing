import numpy as np
import scipy.optimize as opt

def Pi(Par, Target):
    """Par needs to include at least 1 or more parameters one would like to keep constant.
    One of these variables should atleast be able to describe the amount of particle types to add, eg. d,N,...,C,mu"""

# Check what variables is in the parameters
    s = ['H','D','A','f','Ts','d','N','p','E','nu','Rn','Rt','Wn','Wt','mu']
    sl = len(s)
    var = {}
    check = np.ones(10)
    for I in range(5,sl):
        if s[I] in Par:
            L = len(Par[s[I]])
            if I >= 10:
                L = L-1
            break
    
    for I in range(0,sl):
        if s[I] in Par:
            globals()[s[I]] = Par[s[I]]
#####################################################
    def Pi_1(var):
        S = ['f','Ts']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var[k]
                k = k+1
        return f*Ts - Target[0]
    #Determine missing variables in Pi_1 domain
    if 'Ts' in Par:
        Par['f'] = float(opt.fsolve(Pi_1,5))
    else:
        Par['Ts'] = float(opt.fsolve(Pi_1,0.5))

    check[0] = Pi_1(1)
#####################################################
    def Pi_2(var):
        S = ['H','D']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var[k]
                k = k+1
        return H/D - Target[1]
    if 'H' in Par:
        Par['D'] = float(opt.fsolve(Pi_2,0.02))
    else:
        Par['H'] = float(opt.fsolve(Pi_2,0.05))
        
    check[1] = Pi_2(1)
#####################################################
    def Pi_3(var):
        S = ['A','f']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var[k]
                k = k+1
        g = 9.81
        return (A*(np.pi*f*2)**2)/g - Target[2]
    Par['A'] = float(opt.fsolve(Pi_3,0.01))
    
    check[2] = Pi_3(1)
#####################################################
    param = []
    def Pi_0n(var):
        S = ['H','D','d','N']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var
        #print(((d**3)*(N))/((D**2)*H) - Target[3])
        #print('d',d)
        return ((d**3)*(N))/((D**2)*H) - Target[3]
    if 'N' not in Par:
        change = 'N'
        param = opt.fsolve(Pi_0n,10*np.ones(L),xtol = 1E-04)
        param = param.astype(int)
        globals()['N'] = param
    else:
        change = 'd'
        param = opt.fsolve(Pi_0n,0.001*np.ones(L),xtol = 1E-04)
        globals()['d'] = param
    Par[change] = param
    
    check[3] = np.linalg.norm(Pi_0n(np.ones(L)))
#####################################################
    param = []
    def Pi_1n(var):
        S = ['f','d','p','E']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var
        return ((E/1E09)/((d)*(d)*(f**2)*p)) - Target[4]
    if 'E' not in Par:
        change = 'E'
        param = opt.fsolve(Pi_1n,80000*np.ones(L))
        globals()['E'] = param
    else:
        change = 'p'
        param = opt.fsolve(Pi_1n,400*(1/Target[4][0])*np.ones(L))
        globals()['p'] = param
    Par[change] = param
    
    check[4] = np.linalg.norm(Pi_1n(np.ones(L)))
#####################################################
    param = []
    def Pi_2n(var):
        S = ['nu']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var
        return nu - Target[5]    
    change = 'nu'
    param = opt.fsolve(Pi_2n,0.001*np.ones(L))
    globals()['nu'] = param
    Par[change] = param
    
    check[5] = np.linalg.norm(Pi_2n(Par))
#####################################################
    param = []
    def Pi_3n(var):
        S = ['Rn','Rt']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = np.reshape(var,[L,L])
        return np.reshape(Rn/Rt - Target[6],[L**2])
    
    if 'Rn' not in Par:
        change = 'Rn'
        param = np.reshape(opt.fsolve(Pi_3n,0.1*np.ones([L,L])),[L,L])
        globals()['Rn'] = param
    else:
        change = 'Rt'
        param = np.reshape(opt.fsolve(Pi_3n,0.1*np.ones([L,L])),[L,L])
        globals()['Rt'] = param
    Par[change] = param
    
    check[6] = np.linalg.norm(Pi_3n(np.ones([L,L])))
#####################################################
    param = []
    def Pi_4n(var):
        S = ['Wn','Wt']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var
                
        return Wn/Wt - Target[7]
    
    if 'Wn' not in Par:
        change = 'Wn'
        param = opt.fsolve(Pi_4n,0.1*np.ones([L]))
        globals()['Wn'] = param
    else:
        change = 'Wt'
        param = opt.fsolve(Pi_4n,0.1*np.ones([L]))
        globals()['Wt'] = param
    Par[change] = param
    
    check[7] = np.linalg.norm(Pi_4n(np.ones([L])))
#####################################################
    param = []
    def Pi_5n(var):
        S = ['mu']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = np.reshape(var,[L,L])
                
        return np.reshape(mu - Target[8],[L**2])
    change = 'mu'
    param = np.reshape(opt.fsolve(Pi_5n,0.1*np.ones([L,L])),[L,L])
    globals()['mu'] = param
    Par[change] = param
    
    check[8] = np.linalg.norm(Pi_5n(np.ones([L,L])))
#####################################################
    param = []
    def Pi_6n(var):
        S = ['muW']
        k = 0
        for i in range(0,len(S)):
            if S[i] not in Par:
                globals()[S[i]] = var
                
        return muW - Target[9]
    change = 'muW'
    param = opt.fsolve(Pi_6n,0.1*np.ones([L]))
    globals()['muW'] = param
    Par[change] = param
    
    check[9] = np.linalg.norm(Pi_6n(np.ones([L])))
    
    return Par, check

    
def writer(Par,DIR,name):
    import write_macro as wm
    wm.write_header(name,DIR)

    H = Par['H']
    R = Par['D']/2
    mD = Par['D']/7.5

    wm.volume_set(name,H,R,mD,DIR)
    
    L = len(Par['d'])
    for n in range(0,L):
        wm.add_particles(name,'P'+str(n+1),n,Par['p'][n],Par['E'][n],200E9,DIR)
    
    k = 0
    for n in range(0,L):
        for m in range(n,L):
            mu = Par['mu'][n,m]
            Rn = Par['Rn'][n,m]
            Rt = Par['Rt'][n,m]
            #ASSUME COHESION AND ADHESION EFFECTS ARE 0
            W = 0
            
            wm.particle_interaction(name,'P'+str(n+1),'P'+str(m+1),mu,Rn,Rt,W,k,DIR)
            k = k+1
            
    for n in range(0,L):
        mu = Par['muW'][n]
        Rn = Par['Wn'][n]
        Rt = Par['Wt'][n]
        #ASSUME COHESION AND ADHESION EFFECTS ARE 0
        W = 0
            
        wm.particle_interaction(name,'P'+str(n+1),'Wall',mu,Rn,Rt,W,k,DIR)
        k = k+1
    
    k = 0
    for n in range(0,L):
        F = Par['N'][n]*4
        wm.flowrate(name,F,k,DIR)
        k = k+1
    
    k = 0
    for n in range(0,L):
        d = Par['d'][n]
        wm.injector(name,H,D,d,'P'+str(n+1),k,DIR)
        k = k+1

    A = Par['A']
    f = Par['f']
    start = 0.5 + 0.25*L
    TS = Par['Ts']
    
    wm.speed(name,A,f,start,TS,DIR)
    wm.add_table(name,5,start,L-1,DIR)
    wm.save_and_run(name,DIR)


