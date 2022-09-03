import numpy as np
DELTA_T = 1

class Protein():

  def __init__(self,beta=0,dillution=0,start=0,degradation= lambda t : 0):
    """
    beta - rate of production (constant or a function to compute production rate)
    start - initial level
    dillution - constant
    degradation - function

    TODO: give degradation rate + activator class together?
    """
    self.time = 0
    self.level = start
    self.beta = beta
    self.dillution = dillution # how often does it happen? how does it update?
    self.degradation = degradation
    self.alpha = self.degradation(self.time) + self.dillution
    self.level_record = self.levels_vs_time()

  def levels_vs_activator(self,X,K,n):
      # n is the steepness of activation function
      # X is the current level of activator
      # K is activation coefficient (hardwired, depends on affinity)
      return self.beta*(X**n/(K**n+X**n))

  def levels_vs_time(self,time=10):
      lvls = {0:0}
      for j,t in enumerate(np.arange(0,time+DELTA_T,DELTA_T)):
          if t == 0:
              continue
          else:
            lvls[t] = lvls[t-DELTA_T] +self.dx_dt() * DELTA_T
      return lvls

  def level_at_time(self,time):
      return self.level_record[time]

  def update_state(self):
    self.time += DELTA_T
    # try this formula for level too: Y = Y_st(1-e^(-alpha*t)) - only for simple response genes
    self.level = self.level + self.dx_dt()*DELTA_T
    return self.time,max(self.level,0)

  def dx_dt(self):
    production = 0
    alpha = self.dillution + self.degradation(self.time)
    degradation = np.sum(alpha)*self.level
    production = self.beta
    # commented out because activators / repressors are removed from this version

    # if self.activators or self.repressors:
    #   for a in self.activators:
    #     production += self.positive_reg_dxdt(a)
    #   for r in self.repressors:
    #     production -= self.negative_reg_dxdt(r)
    # else:
    #     production = self.beta
    return production-degradation

  def simple_reg_lvl(self,t):
      Y = Y_st(1 - np.exp(-self.alpha * t))

  def positive_reg_dxdt(self, Y,K=0.1):
      """
      Y - activator (Protein class)
      K - dissociation constant between X-Y TODO: constant or a function
      kmax - maximum transcription rate of X  - (should be beta?) TODO: check this statement
      """
      kmax=self.beta
      return kmax*(Y.level / (K+Y.level))

  def negative_reg_dxdt(self, Y,K=0.1):
      """
      Y - repressor (Protein class)
      K - dissociation constant between X-Y TODO: constant or a function
      kmax - maximum transcription rate of X - (should be beta?)
      """
      kmax=self.beta
      return kmax*(K / (K+Y.level))