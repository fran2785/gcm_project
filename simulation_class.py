class System: 
    
    def __init__(self,filename,seed):
        self.filename = filename  # contains the informations about the system
        self.seed = seed
        self.N = 0 
        self.box_lenght = 0 
        self.rho = 0 
        self.box = [] 
        self.position = [] 
        self.eo = 0
        self.wo = 0 
    
    def initialize(self):
        self.N, self.box_lenght, pos = read_cnf_atoms(f'./initial/{self.filename}')
        self.position = pos.transpose()
        self.box = np.array([self.box_lenght, self.box_lenght, self.box_lenght])
        self.rho = self.N/(self.box_lenght)**3
        self.eo, self.wo = mc.montecarlo.interaction(self.position,self.box)
        
        
class Simulation: 
    
    def __init__(self,particle_system,nsteps,dr):
        self.system = particle_system   # connection to the system class
        self.nsteps = nsteps            # total number of steps
        self.dr = dr
        self.count = np.array([0])      # accepted configuration count
        self.position = self.system.position
        self.position_list = []
        
    def simulation(self,T): 
        mc.montecarlo.put_seed(self.system.seed)
        for i in tqdm(range(self.nsteps)):
            mc.montecarlo.mc_step(self.position,self.count,self.dr,self.system.box,T)
            if i%10 == 0: 
                temporary_position = np.copy(self.position)
                self.position_list.append(temporary_position)
        
class Analysis: 
    
    def __init__(self,particle_system,simulation,nbin,dr):
        self.system = particle_system 
        self.sim = simulation 
        self.nbin = nbin
        self.max_distance = dr
        self.e = [] 
        self.w = [] 
        self.bins = []
        self.ghisto = []
        
    def get_energy_virial(self): 
        for i in range(len(self.sim.position_list)):
            pos = self.sim.position_list[i]
            ep, wp = mc.montecarlo.interaction(pos,self.system.box)
            self.e.append(ep) ; self.w.append(wp)
        return np.array(self.e), np.array(self.w)

    def get_gr(self,nequil):
        histo_list = []
        for i in range(len(self.sim.position_list)):
            pos = self.sim.position_list[i]
            bins = np.zeros(self.nbin) ; histo = np.zeros(self.nbin,dtype=int)
            gr.distribution_function.gr_self(pos,self.system.box,self.max_distance,histo,bins)
            omega =2.0*np.pi*self.max_distance*(bins**2)*self.system.rho*self.system.N
            histo_list.append(histo/omega)
        self.bins = bins
        self.ghisto = np.mean(histo_list[nequil:],axis=0)
        return self.bins, self.ghisto
