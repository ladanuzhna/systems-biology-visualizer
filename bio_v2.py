

class ProteinRepressilator():

    def __init__(self, start, max_production_rate,degradation_rate, leaky_production_rate,
                 repressor=None, K=None, n=None):
        self.start = start
        self.concentration = start
        self.degradation_rate = degradation_rate
        self.production_rates = {'max': max_production_rate, 'leaky': leaky_production_rate}
        self.r_kinetics = {'r': repressor, 'K': K, 'n': n}

    def add_repressor(self, repressor, K, n):
        self.r_kinetics = {'r': repressor, 'K': K, 'n': n}

    def reset(self):
        self.concentration = self.start

    def update(self):
        self.concentration = self.concentration + self.rate_of_change()
        return self.concentration

    def rate_of_change(self):
        regulated_production = self.production_rates['max'] * (self.r_kinetics['K'] ** self.r_kinetics['n'] / (
                    self.r_kinetics['K'] ** self.r_kinetics['n'] + self.r_kinetics['r'].concentration **
                    self.r_kinetics['n']))
        production = self.production_rates['leaky'] + regulated_production
        degradation = self.concentration * self.degradation_rate
        return production - degradation

class ProteinNegativeAuto():

    def __init__(self, start, max_production_rate,degradation_rate, leaky_production_rate,
                 repressor=None, K=None, n=None):
        self.start = start
        self.concentration = start
        self.alpha = degradation_rate
        self.production_rates = {'max': max_production_rate, 'leaky': leaky_production_rate}
        self.r_kinetics = {'r': repressor, 'K': K, 'n': n}

    def add_repressor(self, repressor, K, n):
        self.r_kinetics = {'r': repressor, 'K': K, 'n': n}

    def reset(self):
        self.concentration = self.start

    def update(self,concentration=None):
        if not concentration:
            concentration = self.concentration
        self.concentration = max(concentration + self.rate_of_change(concentration),0)
        return self.concentration

    def rate_of_change(self, concentration=None):
        if not concentration:
            concentration = self.concentration
        return self.production_rate(concentration) - self.degradation_rate(concentration)

    def production_rate(self, concentration=None):
        if not concentration:
            concentration = self.concentration
        regulated_production = self.production_rates['max'] * (self.r_kinetics['K'] ** self.r_kinetics['n'] / (
                    self.r_kinetics['K'] ** self.r_kinetics['n'] + concentration **
                    self.r_kinetics['n']))
        return self.production_rates['leaky'] + regulated_production

    def degradation_rate(self, concentration=None):
        if not concentration:
            concentration = self.concentration
        return concentration * self.alpha

class ProteinSimpleReg():

    def __init__(self, start, max_production_rate,degradation_rate):
        self.start = start
        self.concentration = start
        self.alpha = degradation_rate
        self.production_rates = {'max': max_production_rate}

    def reset(self):
        self.concentration = self.start

    def update(self,concentration=None):
        if not concentration:
            concentration = self.concentration
        self.concentration = max(concentration + self.rate_of_change(),0)
        return self.concentration

    def rate_of_change(self,concentration=None):
        if not concentration:
            concentration = self.concentration
        return self.production_rate(concentration) - self.degradation_rate(concentration)

    def production_rate(self,concentration=None):
        if not concentration:
            concentration = self.concentration
        return self.production_rates['max']

    def degradation_rate(self,concentration=None):
        if not concentration:
            concentration = self.concentration
        return concentration * self.alpha