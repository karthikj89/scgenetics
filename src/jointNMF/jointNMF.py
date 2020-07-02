import numpy as np
import scipy

from scipy import sparse
from time import time
from sklearn.utils import check_random_state
from sklearn.decomposition import NMF
from sklearn.utils.extmath import safe_sparse_dot

# constants
_largenumber = 1E100
_smallnumber = 1E-06
_smallnumber2 = 1E-04


class JointNMF:
    def __init__(self, Xh, Xd, Wh=None, Hh=None, Wshh=None, Hshh=None, Wshd=None, Hshd=None, Wd=None, Hd=None, 
                 nh_components=10, nsh_components=10, nd_components=10, gamma=1, mu=None):
        self.Xh = scipy.sparse.csr_matrix(Xh).copy()
        self.Xd = scipy.sparse.csr_matrix(Xd).copy()

        # initializations
        self.nh_components = nh_components
        self.nsh_components = nsh_components
        self.nd_components = nd_components
        
        self.maxiters = 1000
        self.tol = _smallnumber
        
        self.gamma = gamma
        
        
        # initialize the matrix using result from standard NMF
        nmfh = NMF(n_components = nsh_components + nh_components)
        nmfd = NMF(n_components = nsh_components + nd_components)
        
        # healthy programs
        if Wh is None:
            self.Wh = scipy.sparse.csr_matrix(nmfh.fit_transform(Xh))
        else:
            if (Wh.shape != (self.Xh.shape[0], self.nh_components + self.nsh_components)):
                raise ValueError("Initial Wh has wrong shape.")
            self.Wh = np.copy(Wh)
            
        if Hh is None:
            self.Hh = scipy.sparse.csr_matrix(nmfh.components_)
        else:
            if (Hh.shape != (self.nh_components + self.nsh_components, self.Xh.shape[1])):
                raise ValueError("Initial Wh has wrong shape.")
            self.Hh = np.copy(Hh) 
        
        # disease programs
        if Wd is None:
            self.Wd = scipy.sparse.csr_matrix(nmfd.fit_transform(Xd))
        else:
            if (Wd.shape != (self.Xd.shape[0], self.nd_components + self.nsh_components)):
                raise ValueError("Initial Wd has wrong shape.")
            self.Wd = np.copy(Wd)
        
        if Hd is None:
            self.Hd = scipy.sparse.csr_matrix(nmfd.components_)
        else:
            if (Hd.shape != (self.nd_components + self.nsh_components, self.Xd.shape[1])):
                raise ValueError("Initial Wd has wrong shape.")
            self.Hd = np.copy(Hd)

        # option for user input mu or estimated mu 
        if mu:
            self.mu = mu
        else:
            healthy_diff = 0.5*sparse.linalg.norm(self.Xh - safe_sparse_dot(self.Wh, self.Hh), ord='fro')**2
            disease_diff = 0.5*sparse.linalg.norm(self.Xd - safe_sparse_dot(self.Wd, self.Hd), ord='fro')**2
            denominator = sparse.linalg.norm(self.Wh, ord='fro')**2 + sparse.linalg.norm(self.Wd, ord='fro')**2
            self.mu = (healthy_diff + disease_diff)/denominator
    
    @property
    def cost(self):
        self.Wh = self.Wh.tocsr()
        self.Wd = self.Wd.tocsr()
        self.Wshh = self.Wh[:,:self.nsh_components]
        self.Wshd = self.Wd[:,:self.nsh_components]
        
        diff1 = 0.5*sparse.linalg.norm(self.Xh - safe_sparse_dot(self.Wh, self.Hh), ord='fro')**2
        diff2 = 0.5*sparse.linalg.norm(self.Xd - safe_sparse_dot(self.Wd, self.Hd), ord='fro')**2
        diff3 = (self.mu/2)*sparse.linalg.norm(self.Wh, ord='fro')**2 + (self.mu/2)*sparse.linalg.norm(self.Wd, ord='fro')
        diff4 = (self.gamma/2)*sparse.linalg.norm(self.Wshh-self.Wshd, ord='fro')**2
        chi2 = diff1 + diff2 + diff3 + diff4
        return chi2
    
    def solve(self, W_only=False, H_only=False, sparsemode=False, maxiters=None, tol=None):
        
        
        # initialize the state
        t0 = time()        
        niter = 0
        chi2 = self.cost
        oldchi2 = _largenumber
        maxiters = maxiters if maxiters else self.maxiters
        while (niter < maxiters) and np.abs((oldchi2-chi2)/chi2) > self.tol: #((oldchi2-chi2)/chi2 > self.tol):
            # update Wh
            scale2=np.append((self.gamma+self.mu)*np.ones(self.nsh_components), (self.mu)*np.ones(self.nh_components))
            self.Wshd = self.Wd[:,:self.nsh_components]
            Wh_up1 = safe_sparse_dot(self.Xh, self.Hh.T)
            Wshd_transform = self.Wshd.multiply(scipy.sparse.csr_matrix(self.gamma*np.ones((self.Xh.shape[0], self.nsh_components))))
            zeros = scipy.sparse.csr_matrix(np.zeros((self.Xh.shape[0], self.nh_components)))
            Wh_up2 = scipy.sparse.hstack((Wshd_transform, zeros)) #+ _smallnumber2
            Wh_down = safe_sparse_dot(self.Wh, safe_sparse_dot(self.Hh, self.Hh.T)) + safe_sparse_dot(self.Wh, np.diag(scale2)) #+ _smallnumber2
            self.Wh = self.Wh.multiply((Wh_up1 + Wh_up2)/Wh_down).tocsr()
            
            # update Hh
            Hh_up = safe_sparse_dot(self.Wh.T, self.Xh) #+ _smallnumber2
            Hh_down = safe_sparse_dot(safe_sparse_dot(self.Wh.T, self.Wh), self.Hh) #+ _smallnumber2
            self.Hh = self.Hh.multiply(Hh_up/Hh_down).tocsr()
            
            # update Wd
            scale2=np.append((self.gamma+self.mu)*np.ones(self.nsh_components), self.mu*np.ones(self.nd_components))
            self.Wshh = self.Wh[:,:self.nsh_components]
            Wd_up1 = safe_sparse_dot(self.Xd, self.Hd.T)
            Wshh_transform = self.Wshh.multiply(scipy.sparse.csr_matrix(self.gamma*np.ones((self.Xd.shape[0], self.nsh_components))))
            zeros = np.zeros((self.Xd.shape[0], self.nd_components))
            Wd_up2 = scipy.sparse.hstack((Wshh_transform, zeros)) #+ _smallnumber2
            Wd_down = safe_sparse_dot(self.Wd, safe_sparse_dot(self.Hd, self.Hd.T)) + safe_sparse_dot(self.Wd, np.diag(scale2)) #+ _smallnumber2 
            self.Wd = self.Wd.multiply((Wd_up1 + Wd_up2)/Wd_down).tocsr()
                        
            # update Hd
            Hd_up = safe_sparse_dot(self.Wd.T, self.Xd) #+ _smallnumber2
            Hd_down = safe_sparse_dot(safe_sparse_dot(self.Wd.T, self.Wd), self.Hd) #+ _smallnumber2
            self.Hd = self.Hd.multiply(Hd_up/Hd_down).tocsr()
                        
            # chi2
            oldchi2 = chi2
            chi2 = self.cost
            
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