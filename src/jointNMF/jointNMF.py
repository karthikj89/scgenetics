import numpy as np
from time import time
from sklearn.utils import check_random_state
from sklearn.decomposition import NMF

# constants
_largenumber = 1E100
_smallnumber = 1E-06
_smallnumber2 = 1E-04

class JointNMF:
    def __init__(self, Xh, Xd, Wh=None, Hh=None, Wshh=None, Hshh=None, Wshd=None, Hshd=None, Wd=None, Hd=None, 
                 nh_components=10, nsh_components=10, nd_components=10, gamma=1, mu=0.1):
        self.Xh = np.copy(Xh)
        self.Xd = np.copy(Xd)

        # initializations
        self.nh_components = nh_components
        self.nsh_components = nsh_components
        self.nd_components = nd_components
        
        self.maxiters = 1000
        self.tol = _smallnumber
        
        self.mu = mu
        self.gamma = gamma
        
        
        # initialize the matrix using result from standard NMF
        nmfh = NMF(n_components = nsh_components + nh_components)
        nmfd = NMF(n_components = nsh_components + nd_components)
        Wh_init = nmfh.fit_transform(Xh)
        Hh_init = nmfh.components_
        Wd_init = nmfd.fit_transform(Xd)
        Hd_init = nmfd.components_
        
        # healthy specific programs
        if Wh is None:
            self.Wh = Wh_init[:,self.nsh_components:]
        else:
            if (Wh.shape != (self.Xh.shape[0], self.nh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wh = np.copy(Wh)
            
        if Hh is None:
            self.Hh = Hh_init[self.nsh_components:,:]
        else:
            if (Hh.shape != (self.nh_components, self.Xh.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hh = np.copy(Hh)
            
        # healthy derived shared programs
        if Wshh is None:
            self.Wshh = Wh_init[:,:self.nsh_components]
        else:
            if (Wshh.shape != (self.Xh.shape[0], self.nsh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wshh = np.copy(Wshh)
            
        if Hshh is None:
            self.Hshh = Hh_init[:self.nsh_components,:]
        else:
            if (Hshh.shape != (self.nsh_components, self.Xh.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hshh = np.copy(Hshh)
            
        # disease derived shared programs
        if Wshd is None:
            self.Wshd = Wd_init[:,:self.nsh_components]
        else:
            if (Wshd.shape != (self.Xd.shape[0], self.nsh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wshd = np.copy(Wshd)
            
        if Hshd is None:
            self.Hshd = Hd_init[:self.nsh_components,:]
        else:
            if (Hshd.shape != (self.nsh_components, self.Xd.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hshd = np.copy(Hshd)  
        
        # disease specific programs
        if Wd is None:
            self.Wd = Wd_init[:,nsh_components:]
        else:
            if (Wd.shape != (self.Xd.shape[0], self.nd_components)):
                raise ValueError("Initial Wd has wrong shape.")
            self.Wd = np.copy(Wd)
        
        if Hd is None:
            self.Hd = Hd_init[nsh_components:,:]
        else:
            if (Hd.shape != (self.nd_components, self.Xd.shape[1])):
                raise ValueError("Initial Wd has wrong shape.")
            self.Hd = np.copy(Hd)
            
    
    @property
    def cost(self):
        diff1 = 0.5*np.linalg.norm(self.Xh - np.dot(self.Wh, self.Hh) - np.dot(self.Wshh, self.Hshh), ord='fro')**2
        diff2 = 0.5*np.linalg.norm(self.Xd - np.dot(self.Wd, self.Hd) - np.dot(self.Wshd, self.Hshd), ord='fro')**2
        diff3 = (self.mu/2)*np.linalg.norm(self.Wshh, ord='fro')**2 + (self.mu/2)*np.linalg.norm(self.Wh, ord='fro')**2
        diff4 = (self.mu/2)*np.linalg.norm(self.Wshd, ord='fro')**2 + (self.mu/2)*np.linalg.norm(self.Wd, ord='fro')**2        
        diff5 = -1*np.sum(self.gamma*np.dot(self.Wshh.T,self.Wshd))
        #diff3 = self.mu/2*np.linalg.norm(self.Wh, ord='fro')**2
        #diff4 = self.mu/2*np.linalg.norm(self.Wd, ord='fro')**2
        #diff5 = self.gamma/2*np.linalg.norm(self.Wshh-self.Wshd, ord='fro')**2
        chi2 = diff1 + diff2 + diff3 + diff4 + diff5
        return chi2
    
    def solve(self, W_only=False, H_only=False, sparsemode=False, maxiters=None, tol=None):
        
        
        # initialize the state
        t0 = time()        
        niter = 0
        chi2 = self.cost
        oldchi2 = _largenumber
        
        while (niter < self.maxiters) and np.abs((oldchi2-chi2)/chi2) > self.tol: #((oldchi2-chi2)/chi2 > self.tol):
            # update Wh
            Wh_up = np.dot(self.Xh, self.Hh.T)
            Wh_down = np.dot(self.Wh, np.dot(self.Hh, self.Hh.T)) + self.mu*self.Wh + _smallnumber2
            self.Wh = self.Wh*(Wh_up/Wh_down)
            
            # update Hh
            Hh_up = np.dot(self.Wh.T, self.Xh)
            Hh_down = np.dot(np.dot(self.Wh.T, self.Wh), self.Hh) + _smallnumber2
            self.Hh = self.Hh*(Hh_up/Hh_down)
            
            # udpate Wshh
            Wshh_up = np.dot(self.Xh, self.Hshh.T) + self.gamma*self.Wshd
            Wshh_down = np.dot(self.Wshh, np.dot(self.Hshh, self.Hshh.T)) + self.mu*self.Wshh + _smallnumber2
            self.Wshh = self.Wshh*(Wshh_up/Wshh_down)
            
            # update Hshh
            Hshh_up = np.dot(self.Wshh.T, self.Xh)
            Hshh_down = np.dot(np.dot(self.Wshh.T, self.Wshh), self.Hshh) + _smallnumber2
            self.Hshh = self.Hshh*(Hshh_up/Hshh_down)
            
            # update Wd
            Wd_up = np.dot(self.Xd, self.Hd.T)
            Wd_down = np.dot(self.Wd, np.dot(self.Hd, self.Hd.T)) + self.mu*self.Wd + _smallnumber2
            self.Wd = self.Wd*(Wd_up/Wd_down)
            
            # update Hd
            Hd_up = np.dot(self.Wd.T, self.Xd)
            Hd_down = np.dot(np.dot(self.Wd.T, self.Wd), self.Hd) + _smallnumber2
            self.Hd = self.Hd*(Hd_up/Hd_down)
            
            # update Wshd
            Wshd_up = np.dot(self.Xd, self.Hshd.T) + self.gamma*self.Wshh
            Wshd_down = np.dot(self.Wshd, np.dot(self.Hshd, self.Hshd.T)) + self.mu*self.Wshd + _smallnumber2
            self.Wshd = self.Wshd*(Wshd_up/Wshd_down)
            
            # update Hshd
            Hshd_up = np.dot(self.Wshd.T, self.Xd)
            Hshd_down = np.dot(np.dot(self.Wshd.T, self.Wshd), self.Hshd) + _smallnumber2
            self.Hshd = self.Hshd*(Hshd_up/Hshd_down)
            
            # chi2
            oldchi2 = chi2
            chi2 = self.cost
            
            #print(niter, oldchi2, chi2)

            # Some quick check. May need its error class ...
            if (not np.isfinite(chi2)):
               raise ValueError("NMF construction failed, likely due to missing data")

            if (np.mod(niter, 20)==0):
                print("Current Chi2={0:.4f}, Previous Chi2={1:.4f}, \
                      Change={2:.4f}% @ niters={3}".format(chi2,oldchi2,((oldchi2-chi2)/oldchi2)*100.,niter), flush=True)

            niter += 1
            if (niter == self.maxiters):
                print("Iteration in re-initialization reaches maximum number = {0}".format(niter), flush=True)
                
        time_used = (time()-t0)/60.
        print("Took {0:.3f} minutes to reach current solution.".format(time_used), flush=True)

        return (chi2, time_used)