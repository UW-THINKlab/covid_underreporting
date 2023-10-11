#!/usr/bin/env python3

####################################################
# filename: InterventionSIR.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: specifies the SIR models
# with and without or peak-minimizing
# intervention according to our
# model of intervention functions b(t)
####################################################

import numpy as np 
import pandas as pd
from scipy.integrate import odeint

import optimize_interventions as oi

class InterventionSIR():
    """
    class for handling our SIR
    model with interventions
    """
    count = 0

    def __init__(self,
                 b_func = None,
                 R0 = None,
                 gamma = None,
                 inits = None
    ):
        if b_func is None:
            # use the class for handeling no interventions
            b_func = Intervention()
        self.b_func = b_func     # assign a intervention type
        self.R0 = R0
        self.gamma = gamma
        self.inits = inits       # save initial states
        self.reset()
        
    def reset(self):
        self.state = self.inits                                      # assign SIR states
        self.time = 0
        self.time_ts = np.array([])                                  # empty the timeframe
        self.state_ts = np.array([[], [], []]).reshape((-1, 6))      # empty the state frame
        self.state_obs = self.inits
        self.report_rate = 0.15
        self.counter = 0

    

    def deriv(self, state, time):
        # pass in state and time
        S, I, R = state
        self.counter += 1
        beta = self.R0 * self.gamma                                  # caclculate transmission rate from R0 and gamma
        b = self.b_func(time, beta, self.gamma, S, I)                # what should be the next step transmission rate given a duration(time),
                                                                    # transmission rate (beta), recovery rate (self.gamma), States (S and I)
                                                                    # possibly add underreporting here
        
        dS = -b * beta * S * I                                        # change in S I R
        dI = b * beta * S * I - self.gamma * I                 
        dR = self.gamma * I

        # c = self.c_func()
        # dI_obs = (I - I_obs) * report_rate
        # dR_obs = gamma * (I - I_obs) * report_rate * c
        # dS_obs = 0 - (I - I_obs) * report_rate - gamma * (I - I_obs) * report_rate * c

        return np.array([dS, dI, dR])

    def deriv_new(self, state, time):
        # pass in state and time
        S, I, R, S_obs, I_obs, R_obs = state
        beta = self.R0 * self.gamma 
        self.b_func.report_rate = self.r_func(time)                                  # caclculate transmission rate from R0 and gamma
        b = self.b_func(time, beta, self.gamma, S,I)                # what should be the next step transmission rate given a duration(time),
                                                                    # transmission rate (beta), recovery rate (self.gamma), States (S and I)
        report_rate = self.r_func(time)                                                            # possibly add underreporting here
        dS = -b * beta * S * I                                        # change in S I R
        dI = b * beta * S * I - self.gamma * I                 
        dR = self.gamma * I

        dI_obs = dI * (report_rate) 
        dR_obs = self.gamma * I_obs

        dS_obs = 0 - dI_obs - dR_obs

        return np.array([dS, dI, dR, dS_obs, dI_obs, dR_obs])

    def deriv_null(self, state, time):  # without any intervention so transmission rate is the unchanged
        S, I, R = state
        beta = self.R0 * self.gamma
        dS = -beta * S * I
        dI = beta * S * I - self.gamma * I
        dR = self.gamma * I

        return np.array([dS, dI, dR])

    def r_func(self, time):
        if self.case == 0: # increasing
            report_rate = 0.8 / (self.b_func.t_i + self.b_func.tau * self.b_func.f) * time + 0.15
        elif self.case == 1: # decreasing
            report_rate = 0.1 - 0.1/ (self.b_func.t_i + self.b_func.tau * self.b_func.f) * time + 0.05
        elif self.case == 2: # suppose there are two periods within tau*f
            report_rate = (np.sin(4 * np.pi / (self.b_func.t_i + self.b_func.tau * self.b_func.f) * time) + 1) * 0.9 /2 + 0.05
        elif self.case == 3: # fixed rate
            report_rate = 0.15
            #report_rate should use 0.15
            
        else: # no underreporting
            report_rate = 1   

        return report_rate
    
    def integrate(self, final_time, fineness = 10000):

        # Below is the original code
        times = np.linspace(self.time,
                            final_time,
                            fineness)

        results = odeint(self.deriv_new, self.state, times)
        self.state = results[-1]
        self.time = final_time
        self.time_ts = np.concatenate([self.time_ts,
                                       times])
        self.state_ts = np.concatenate([self.state_ts,
                                        results])
        # Print the I_max value
        I_max_achieved = np.max(self.state_ts[:,1])
        print('The I_max achieved is: {}'.format(I_max_achieved))
        # print the cumulatively infected individuals
        total_infected = self.state_ts[:,2][-1]
        print('The cumulatively infected individuals are: {}'.format(total_infected))
        df = pd.DataFrame(self.state_ts, columns = ['S', 'I','R', 'S_obs','I_obs','R_obs'])
        
        for i in range(0, len(df)):
            I = df.iloc[i]['I']
            S = df.iloc[i]['I']
            
            I_used = np.minimum(S+I, I * self.b_func.report_rate / self.b_func.correct_rate)
            S_used = S + (I - I_used) 
            df.iloc[i]["I_obs"] = I_used
            df.iloc[i]["S_obs"] = S_used
            df.iloc[i]["R_obs"] = 1 - I_used - S_used
        
        df.to_csv(self.outpath + 'states for case {} tau {} cc {}.csv'.format(self.case, self.b_func.tau, self.b_func.correct_rate), index=False)

        return (times, results)

    def integrate_null(self, final_time, fineness = 10000):
        times = np.linspace(self.time,
                            final_time,
                            fineness)
        results = odeint(self.deriv_null, self.state, times)
        return (times, results)
    
    def I_max_SI(self, S_x, I_x):
        """
        get the maximum value of I(t) 
        in the window from t s.t. S = S_x,
        I = I_x to t = infinity
        """
        return (S_x + I_x - 
                (1/self.R0) * np.log(S_x) - 
                (1/self.R0) + 
                (1/self.R0) * np.log(1/self.R0))

    def I_of_S(self, S):
    # wholly susceptible population without intervention, express I in terms of S
    # Eq 14 15
        S0, I0, Rec0 = self.inits
        return (I0 + S0 - (1/self.R0) * np.log(S0) -
                S + (1/self.R0) * np.log(S))

    def t_of_S(self, S_target):
        S0, I0, Rec0 = self.inits
        if np.isnan(S_target):
            raise ValueError("Cannot find time "
                             "for non-numeric/nan S\n\n"
                             "check that S is being "
                             "calculated correctly")
        def deriv(t, S_val):
            I = self.I_of_S(S_val)
            return -1 / (self.R0 * self.gamma * S_val * I)
        return odeint(deriv, 0, np.linspace(S0, S_target, 2))[-1]

    def get_I_max(self,
                  allow_boundary_max = True):
        last_timestep_error = (
            "Max at last timestep. "
            "You likely need to "
            "increase integration "
            "max time. If this was expected, "
            "set allow_boundary_max = True")
        first_timestep_error = (
            "Max at first timestep. "
            "Your model may be misspecified. "
            "If this was expected, "
            "set allow_boundary_max = True")
        ## check that we didn't get a boundary soln
        wheremax = np.argmax(self.state_ts[:, 1])
        if not allow_boundary_max:
            if wheremax == self.state_ts[:, 1].size: 
                raise ValueError(last_timestep_error)
            elif wheremax == 0:
                raise ValueError(first_timestep_error)
        return self.state_ts[wheremax, 1]
    
    def get_t_peak(self):
        return self.t_of_S(1 / self.R0)

    def __repr__(self):
        return ("InterventionSIR with R0 = {}, "
                "gamma = {}, and an intervention "
                "function {}".format(
                    self.R0,
                    self.gamma,
                    self.b_func))



class Intervention():
    """
    class for defining intervention
    functions b(t)
    """

    def __init__(self,
                 tau = None,
                 t_i = None,
                 sigma = None,
                 f = None,
                 S_i_expected = None,
                 I_i_expected = None,
                 strategy = None):
        self.tau = tau
        self.t_i = t_i
        self.sigma = sigma
        self.f = f
        self.S_i_expected = S_i_expected
        self.I_i_expected = I_i_expected
        self.strategy = strategy
        # self.report_rate = 0.15
        self.I_thres = 1
        self.lasttime = -1
        self.correct_rate = 1

        self.repertoire = {
            "fixed": self.fixed_b,
            "maintain-suppress-time": self.maintain_suppress_time,
            "maintain-suppress-state": self.maintain_suppress_state,
            "full-suppression": self.fixed_b}

    def __call__(self,
                 time,
                 beta,
                 gamma,
                 S,
                 I):
        return self.repertoire[self.strategy](
            time,
            beta,
            gamma,
            S,
            I)

    def fixed_b(self,
                time,
                beta,
                gamma,
                S,
                I):
        """
        Fixed intervention of strictness
        sigma
        """
        if time >= self.t_i and time < self.t_i + self.tau:
            result = self.sigma
        else:
            result = 1

        return result

    def maintain_suppress_time(self,
                              time,
                              beta,
                              gamma,
                              S,
                              I):
        """
        Variable maintain/suppress
        intervention tuned by 
        current time
        """
        
        
        if time >= self.t_i and time < self.t_i + self.tau * self.f:
            # print(time)
            # print(self.report_rate)

            I_used = np.minimum(S+I, I * self.report_rate / self.correct_rate)
            S_used = S + (I - I_used)               


            S_expected = (S_used -
                         gamma *  (360/9999)*
                         I_used)       
            # if self.counter % 100 == 0:
            #     print(self.counter)
            #     print(time)
                # print(self.I_i_observed)
                # print(self.S_i_observed) 
            # orginally 

            # S_expected = (S_used -
            #               gamma * (360/9999) *
            #               I_used)
            # result = np.maximum(np.minimum(gamma / (beta * S_expected), 1),0)
            result = np.maximum(np.minimum(gamma / (beta * S_expected), 1),0)
            

        elif (time >= self.t_i + self.tau * self.f and
              time < self.t_i + self.tau):
            result = 0
        else:
            result = 1
            
        self.lasttime = time


        """ R_i_observed = 1 - self.I_i_observed - self.S_i_observed

        # if the report_rate does not vary with time,
        # self.I_i_observed = np.minimum((I - self.I_i_observed) * self.report_rate + self.I_i_observed, 1.0)
        # if the report_rate varies with time,
        # self.report_rate = np.maximum(0.15 + (0.05 - 0.15) * time / (self.t_i + self.tau * self.f),0.05)
        # print(self.report_rate)
        self.I_i_observed = np.minimum((I - self.I_i_observed) * (self.report_rate) + self.I_i_observed, 1.0)

        # correct_rate = (I - self.I_i_observed) * self.report_rate / self.I_i_observed
        # print(correct_rate)
        correct_rate = 0.75
        # self.report_rate = 0.15 + (0.05 - 0.15) * time / (self.t_i + self.tau * self.f)
        # print(self.report_rate)

        if self.I_i_observed == 1.0:
            print('Warning, I exceeds 1')
        R_i_observed = R_i_observed + gamma * (I - self.I_i_observed) * (self.report_rate ) * (360/9999)
        self.S_i_observed = 1 - self.I_i_observed - R_i_observed

        R_i_corrected = 1 - self.I_i_corrected -self.S_i_corrected

        self.I_i_corrected = np.minimum(self.I_i_observed / correct_rate,0.5)
        R_i_corrected = np.minimum(R_i_observed / correct_rate,1)
        # R_used = R_used + self.I_used * gamma
        self.S_i_corrected = 1 - self.I_i_corrected - R_i_corrected """


        return result

    def maintain_suppress_state(self,
                               time,
                               beta,
                               gamma,
                               S,
                               I):
        """
        Variable maintain/suppress
        intervention tuned by 
        current state of the system
        (S(t), I(t))
        """
        # need to recalculate the implementation time



        if time >= self.t_i and time < self.t_i + self.tau * self.f:
            result = gamma / (beta * S)
        elif (time >= self.t_i + self.tau * self.f and
              time < self.t_i + self.tau):
            result = 0
        else:
            result = 1
        return result
        



## helper functions for the above:
def get_SI_expected(
        t_i_expected,
        R0,
        gamma,
        I0 = 1e-6,
        Rec0 = 0,
        integration_fineness = 1000000):
    """
    What S_i and I_i do we
    expect in a time tuned
    intervention that we
    plan to start at a time
    t_i_expected?
    """
    S0 = 1 - I0 - Rec0
    def deriv(state, time):
        beta = R0 * gamma
        S, I = state
        dS = -beta * S * I
        dI = beta * S * I - gamma * I
        return np.array([dS, dI])

    state = odeint(deriv, [S0, I0],
                   np.linspace(0, t_i_expected, integration_fineness))
    expected_S_i, expected_I_i = state[-1]
    return (expected_S_i, expected_I_i)

def make_state_tuned_variable_b_func(tau, t_i, f):
    """
    create a function to execute
    the variable-b intervention
    with parameters t_i, tau and f
    where we operate based on the current
    state
    """
    return Intervention(
        tau = tau,
        t_i = t_i,
        f = f,
        strategy = "maintain-suppress-state")

def make_time_tuned_variable_b_func(tau,
                                    t_i,
                                    f,
                                    S_i_expected,
                                    I_i_expected):
    """
    create a function to execute
    the variable-b intervention
    with parameters t_i, tau and f
    where we operate based on the current
    time
    """                                          
    return Intervention(
        tau = tau,
        f = f,
        t_i = t_i,
        S_i_expected = S_i_expected,
        I_i_expected = I_i_expected,
        strategy = "maintain-suppress-time")



def make_fixed_b_func(tau, t_i, sigma):
    """
    create a function to execute
    the fixed intervention
    with parameters t_i, tau and sigma 
    """
    return Intervention(
        tau = tau,
        t_i = t_i,
        sigma = sigma,
        strategy = "fixed")


def check_gamma(gamma,
                model,
                offset,
                strategy = None,
                R0 = None,
                tau = None,
                verbose = False):
    """
    Check that a gamma 
    value produces a viable
    t_i given the model
    and its intervention
    function
    """
    if strategy is None:
        strategy = model.b_func.strategy
    if R0 is None:
        R0 = model.R0
    if tau is None:
        tau = model.b_func.tau

    model.gamma = gamma
        
    if strategy in ["maintain-suppress-time", "maintain-suppress-state"]:
        opt_S_i, opt_f = oi.calc_Sf_opt(
            R0,
            gamma * tau)
        t_i = model.t_of_S(opt_S_i)[0]
    elif strategy in ["fixed"]:
        fixed_S_i, fixed_sigma = oi.calc_Sb_opt(
            R0,
            gamma,
            tau)
        t_i = model.t_of_S(fixed_S_i)[0]
    elif strategy in ["full-suppression"]:
        full_supp_S_i = oi.calc_S_var_opt(
            R0,
            gamma * tau,
            0)
        t_i = model.t_of_S(full_supp_S_i)[0]
    if verbose:
        print("strategy:", strategy)
        print("gamma:", gamma)
        print("t_i:", t_i)
        print("offset:", offset)
    return min(t_i + offset, t_i) > 0

