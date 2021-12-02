#BLOCKING PROBABILITY SIMULATION
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import scipy.stats as st
import sympy as sym
import scipy
import operator
import queue
import queueNoDonation as bq
import math
from scipy.stats import norm

#DEFINE GLOBAL VARIABLES
b = 200
#uniform dist: 100
#beta dist: 200 #expected bail amount
b_sq = 400*400/12 + 200*200
# uniform dist 200*200/3 + 100*100
# beta dist 400*400/12 + 200*200 
#Var[Beta] + E[Beta]^2
runs = 1000#0
err = 0.05 #should be 
M01 =  8281 #beta(1,1)
#my guess 8775 
#for lambda = 10: 66780
#8629.82-- beta by Pender
#1619.82 #1600 - 100*(np.random.binomial(1, 0.46))
#initial = 1700 #np.random.binomial(1, 0.2) #is 1 with 0.7 probability
#initial = 3606807.74 #for weibull dist
#initial = 1835.13649 #for exponential dist
#initial = 1736.64688 #1624.99246 #for gamma dist
#
#M02 = solve_MO(err, 6000, b)  # ADD EQUATION
inf_funding = -1000000000
time_period = 200
inc = 0.1
lam_b = 1
G_bar = 30



colors1 = ["r-", "y-", "g-", "c-", "b-", "m-"]
colors2 = ["pink", "violet", "lime", "dodgerblue", "darkblue", "indigo"]
color_pair = [["tomato", "red", "lightcoral"], 
              ["yellow", "gold", "lightyellow"], 
              ["palegreen", "lawngreen", "honeydew"], 
              ["aqua", "deepskyblue", "lightcyan"], 
              ["violet", "blueviolet", "thistle"]]
bail_amt = ["200", "500", "750", "1000", "1500", "2000"] #CHANGED FROM 2500



#BAIL GENERATORS
def gen_deterministic_bail():
    return 100

def gen_uni_bail():
    return np.random.uniform(0,200)

def gen_weibull_bail():
    return 50+np.random.weibull(0.223637)

def gen_exp_bail():
    return np.random.exponential(100)

def gen_gamma_bail():
    return np.random.gamma(2,50)

def gen_beta_bail():
    return 400*np.random.beta(1,1)
#2000*np.random.beta(1,9)

def gen_semiuni_bail():
    np.random.choice([500,750,1000,1500,2000,2500],p = [0.181818, 0.181818,0.181818,0.181818,0.181818,0.0909090909])

def initial_deterministic():
    return np.random.binomial(1, 0.2)*100 + 1600

def mo1_eq(x, err, lam_b, b, G_bar, b_sq):
    return x + lam_b*b*G_bar + norm.ppf(1-err)*math.sqrt(lam_b*b_sq*G_bar)
    
#def mo2_eq(x, err, avgMO, b):
#    return ((1 - err)*(0.5 - math.erf((avgMO - x)/(200*math.sqrt(30)))) - (0.5 -  math.erf((avgMO + b - x)/(200*math.sqrt(30)))))


#def solve_MO2(err, avgMO, b):
#    stu = -1000
#    stu_eq = 100000
#    for i in range(9000):
#        if stu_eq >= np.abs(mo2_eq(i, err, avgMO, b)):
#            stu = i
#            stu_eq = np.abs(mo2_eq(stu, err, avgMO, b))
#    return stu
def mo2_eq(M02, e, b):
    return ((1 - e)*(0.5 - math.erf((-1*M02 + lam_b*b*G_bar)/(math.sqrt(lam_b*b_sq*G_bar)))) - (0.5 - math.erf((b - M02 + lam_b*b*G_bar)/(math.sqrt(lam_b*b_sq*G_bar)))))
    

def solve_MO2(e, b):
    stu = -1000
    stu_eq = 100000
    for i in range(9000):
        #print(mo2_eq(i, err, b))
        if abs(stu_eq) >= abs((mo2_eq(i, e, b))):
            stu = i
            stu_eq = np.abs(mo2_eq(stu, e, b))
            #print("happened" + str(stu))
    return stu

def solve_mu_sd(time_period):
    form_E_M = np.zeros(time_period)
    form_Var_M = np.zeros(time_period)
    for time in range(time_period):
        G = 0
        if (time <= 10):
            G = 0
        elif (time <= 20):
            G = 0.2*(time - 10)
        elif (time <= 30):
            G = 2 + 0.4*(time - 20)
        elif (time <= 40):
            G = 6 + 0.6*(time - 30)
        elif (time <= 50):
            G = 12 + 0.8*(time - 40)
        else:
            G = time - 30
        form_E_M[time] = M01 - lam_b*(time)*b + lam_b*b*(G)
        form_Var_M[time] = time*(lam_b*b_sq) - lam_b*b_sq*G
    return form_E_M, form_Var_M

#def mo2_eq(x, err, avgMO, b):
#    return ((1 - err)*(0.5 - math.erf((avgMO - x)/(math.sqrt(lam_b*b_sq*30)))) - (0.5 -  math.erf((avgMO + b - x)/(math.sqrt(lam_b*b_sq*30)))))
#((1 - err)*(0.5 - math.erf((avgMO - x)/(316*math.sqrt(2)))) - (0.5 -  math.erf((avgMO + b - x)/(316*math.sqrt(2)))))
    

#def solve_MO2(err, avgMO, b):
#    stu = -1000
#    stu_eq = 100000
#    for i in range(9000):
#        if stu_eq >= np.abs(mo2_eq(i, err, avgMO, b)):
#            stu = i
#            stu_eq = np.abs(mo2_eq(stu, err, avgMO, b))
#    return stu
    

def bail_simulation(x, runs, err, initial, min_amt, time_period, inc, lam_b, bail_generator):
    avg_M = 0
    var_M = []
    avg_O = 0
    avg_Q = 0
    avg_block = 0
    queues = []
    tb = 0
    report = runs/10
    bl = 0
    for i in range(runs):
        if ((i % report) == 0):
            print(str(i*100 / runs) + "%")
        t, M, O, Q, queue, block, total_blocked = bq.sim(initial, min_amt, 0, time_period, inc, lam_b, bail_generator)
        #initial, min_funding, expense, time_period, inc, lam_b, ranking, gen_bail
        avg_M += np.array(M)/runs
        var_M.append(M)
        avg_O += np.array(O)/runs
        avg_Q += np.array(Q)/runs
        tb += total_blocked/runs
        b_idx = int(b/100 - 1)
        bl += np.array(block[b_idx])/runs
        blocking = (np.array(M))
        #print(blocking)
        blocking[blocking < x] = 1
        blocking[blocking >= x] = 0
        #print(blocking)
        #if (abs(sum(blocking - block[0])) > 0):
            #print(sum(blocking - block[0]))
        avg_block = blocking/runs  
        queues.append(np.array(queue))
    avg_queue = np.mean(queues, axis = 0).flatten()
    var_queue = np.var(queues, axis = 0).flatten()
    var_M_t =  np.var(var_M, axis = 0).flatten()
    #print(sum(avg_block))
    return t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M_t, tb


bail_gen = gen_beta_bail #gen_uni_bail#


#GRAPH GENERATORS:
def gen_block_prob_for_e():
    fig, axs = plt.subplots(1, 1)
    e = np.array([0.01, 0.05, 0.1, 0.25, 0.5])
    for i in range(len(e)):
        MO1_form = mo1_eq(b, e[i], lam_b, b, G_bar, b_sq)
        MO2_form = solve_MO2(e[i], 6000, b)
        t, avg_block, var_block, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(runs, err, MO1_form, inf_funding, time_period, inc, lam_b, bail_gen)
        axs.plot(t, avg_block[0], color_pair[i][0], label = "Average block probability by simulation for M1 solved by formula for epsilon = " + str(e[i]))
        t, avg_block, var_block, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(runs, err, MO2_form, inf_funding, time_period, inc, lam_b, bail_gen)
        axs.plot(t, avg_block[0], color_pair[i][1], label = "Average block probability by simulation for M2 solved by formula for epsilon = " + str(e[i]))
        axs.plot(np.arange(10, time_period), np.ones(time_period-10)*e[i], color_pair[i][2])#, label = "Block probability epsilon")
    axs.legend()
    plt.show()
    
def gen_block_prob_for_e1():
    fig, axs = plt.subplots(1, 1)
    axs.set_ylim([0, 0.6])
    e = np.array([0.01, 0.05, 0.1, 0.25, 0.5])
    for i in range(len(e)):
        MO1_form = mo1_eq(b, e[4-i], lam_b, b, G_bar, b_sq)
        t, avg_block, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(b, runs, e[4-i], MO1_form, inf_funding, time_period, inc, lam_b, bail_gen)
        axs.plot(t, avg_block, color_pair[4-i][1], label = "Average block probability by simulation for M1 solved by formula for epsilon = " + str(e[4-i]))
        axs.plot(np.arange(10, time_period), np.ones(time_period-10)*e[4-i], color_pair[4-i][1], linestyle='dashed')
    axs.legend(loc='upper left')
    plt.show()
    
def gen_block_prob_for_e2():
    fig, axs = plt.subplots(1, 1)
    axs.set_ylim([0, 0.55])#[0, 0.15])
    e = np.array([0.01, 0.05, 0.1, 0.25, 0.5])#[0.005, 0.01, 0.025, 0.05, 0.1])#
    for i in range(len(e)):
        MO2_form = solve_MO2(e[4-i], b)
        #print(MO2_form)
        t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(b, runs, e[4-i], MO2_form, 0, time_period, inc, lam_b, bail_gen)
        axs.plot(t, bl, color_pair[4-i][1], label = "Average block probability by simulation for M2 solved by formula for epsilon = " + str(e[4-i]))
        axs.plot(np.arange(10, time_period), np.ones(time_period-10)*e[4-i], color_pair[4-i][1], linestyle='dashed')
    axs.legend(loc='upper left')
    plt.show()
    
    
def gen_err_v_m():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("blocking probability, " + r'$\epsilon$')
    axs.set_ylabel("initial M0 by  formula")
    e = np.arange(0.01, 0.99, 0.01)
    MO1_form = np.zeros((len(e)))
    for i in range(len(e)):
        MO1_form[i] = mo1_eq(b, e[i], lam_b, b, G_bar, b_sq)
    axs.plot(e, MO1_form, "indigo", label = "M1 calculated by formula for different error values")
    MO2_form = np.zeros((len(e)))
    for i in range(len(e)):
        MO2_form[i] = solve_MO2(e[i], b)
    axs.plot(e, MO2_form, "forestgreen", label = "M2 calculated by formula for different error values")
    axs.legend()
    plt.show()
    
    
def gen_err_v_m1():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("blocking probability, " + r'$\epsilon$')
    axs.set_ylabel("initial M0 by  formula")
    e = np.arange(0.01, 0.99, 0.01)
    MO1_form = np.zeros((len(e)))
    for i in range(len(e)):
        MO1_form[i] = mo1_eq(b, e[i], lam_b, b, G_bar, b_sq)
    axs.plot(e, MO1_form, "indigo", label = "M1 calculated by formula for different error values")
    axs.legend()
    plt.show()
    
def gen_err_v_m2():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("blocking probability, " + r'$\epsilon$')
    axs.set_ylabel("initial M0 by  formula")
    e = np.arange(0.01, 0.99, 0.01)
    MO2_form = np.zeros((len(e)))
    for i in range(len(e)):
        MO2_form[i] = solve_MO2(e[i], b)
    axs.plot(e, MO2_form, "indigo", label = "M2 calculated by formula for different error values")
    axs.legend()
    plt.show()
    
def gen_avg_var_sim():
    fig, axs = plt.subplots(2, 1)
    form_E_M, form_Var_M = solve_mu_sd(time_period)
    t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(b, runs, err, M01, inf_funding, time_period, inc, lam_b, bail_gen)
    axs[0].set_xlabel("time, t")
    axs[0].set_ylabel("Expected Value of M(t)")
    axs[1].set_xlabel("time, t")
    axs[1].set_ylabel("Variance of M(t)")
    axs[0].plot(t, avg_M, "darkblue", label = "Average M[t] by simulation")
    axs[1].plot(t, var_M, "g-", label = "Variance of M[t] by simulation")
    axs[0].plot(np.arange(time_period), form_E_M, "dodgerblue", label = "Average M[t] by formula")
    axs[1].plot(np.arange(time_period), form_Var_M, "lime", label = "Variance of M[t] by formula")
    axs[0].legend()
    axs[1].legend()
    plt.show()

def search_m2():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("blocking probability, " + r'$\epsilon$')
    axs.set_ylabel("initial M0")
    m2 = np.arange(0, 6000, 5)#0.01, 0.99, 0.01)
    e = np.zeros(len(m2))
    for i in range(len(m2)):
        avg_bl = 0
        run = 100
        for j in range(run):
            t, M, O, Q, queue, block, total_blocked = bq.sim(m2[i], 0, 0, time_period, inc, lam_b, bail_gen)
            avg_bl += np.array(block[0])/run
        e[i] = np.mean(avg_bl[50:])
    axs.plot(e, m2, "forestgreen", label = "Initial capital for blocking model searched by simulation to get "+ r'$\epsilon$' +" blocking probability")
    axs.legend()
    plt.show()
    
def form_and_sim_m2():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("blocking probability, " + r'$\epsilon$')
    axs.set_ylabel("initial M0")
    m2 = np.arange(0, 6000, 5)#0.01, 0.99, 0.01)
    e = np.zeros(len(m2))
    for i in range(len(m2)):
        #print(i)
        avg_bl = 0
        run = 25
        for j in range(run):
            t, M, O, Q, queue, block, total_blocked = bq.sim(m2[i], 0, 0, time_period, inc, lam_b, bail_gen)
            avg_bl += np.array(block[0])/run
        e[i] = np.mean(avg_bl[50:])
    axs.plot(e, m2, "forestgreen", label = "Initial capital for blocking model searched by simulation to get "+ r'$\epsilon$' +" blocking probability")
    e = np.arange(0.01, 0.99, 0.01)
    MO2_form = np.zeros((len(e)))
    for i in range(len(e)):
        MO2_form[i] = solve_MO2(e[i], b)
    axs.plot(e, MO2_form, "indigo", label = "Initial capital for blocking model calculated by formula for different " + r'$\epsilon$')
    axs.legend()
    plt.show()
    
        
#form_and_sim_m2()
#form_and_sim_m2()
#gen_block_prob_for_e2()
#gen_avg_var_sim()

#bail_simulation input: runs, err, initial, min_amt, time_period, inc, lam_b
#t, avg_block, var_block, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(runs, err, M01, inf_funding, time_period, inc, lam_b, bail_gen)

#t, avg_block, var_block, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(runs, err, initial_deterministic, inf_funding, time_period, inc, lam_b, bail_gen)

#M1_IC_t, M1_IC_avg_block, M1_IC_var_block, M1_IC_avg_queue, M1_IC_var_queue, M1_IC_avg_M, M1_IC_avg_Q,  M1_IC_var_M, M1_IC_total_blocked = bail_simulation(runs, err, M01, inf_funding, time_period, inc, lam_b, bail_gen)

#M1_B_t, M1_B_avg_block, M1_B_var_block, M1_B_avg_queue, M1_B_var_queue, M1_B_avg_M, M1_B_avg_Q,  M1_B_var_M, M1_B_total_blocked = bail_simulation(runs, err, M01, 0, time_period, inc, lam_b, bail_gen)

#M2_IC_t, M2_IC_avg_block, M2_IC_var_block, M2_IC_avg_queue, M2_IC_var_queue, M2_IC_avg_M, M2_IC_avg_Q,  M2_IC_var_M, M2_IC_total_blocked = bail_simulation(runs, err, M02, inf_funding, time_period, inc, lam_b, bail_gen)

#M2_B_t, M2_B_avg_block, M1_B_var_block, M2_B_avg_queue, M2_B_var_queue, M2_B_avg_M, M2_B_avg_Q,  M2_B_var_M, M2_B_total_blocked = bail_simulation(runs, err, M02, 0, time_period, inc, lam_b, bail_gen)



                    
                    
                    
#UNCOMMENT TO SHOW AVG BLOCK PROBABILITY
#axs.plot(t, avg_block[0], "m-", label = "Average block probability by simulation")
#axs.plot(np.arange(10, time_period), np.ones(time_period-10)*err, "violet", label = "Block probability epsilon")

#axs.legend()
#UNCOMMENT FOR AVG AND VAR OF SIM vs FORMULA 





#UNCOMMENT FOR M1 VS M2 graphs
#axs.plot(M1_IC_t, M1_IC_avg_block[j], colors1[j], label = "M1 with infinite credit")
#print(np.mean(M1_IC_avg_block[j][50:]))
#axs.plot(M1_B_t, M1_B_avg_block[j], colors2[j], label = "M1 with blocking")
#axs.legend()
#plt.show()

#axs.plot(M2_IC_t, M2_IC_avg_block[j], colors1[j], label = "M2 with infinite credit for bail "+str(bail_amt[j]) )
#axs.plot(M2_B_t, M2_B_avg_block[j], colors2[j], label = "M2 with blocking for bail "+str(bail_amt[j]))
#axs.legend()

#if (time <= 10):
#        G = time
#    elif (time <= 20):
#        G = 10 + 0.8*(time - 10)
#    elif (time <= 30):
#        G = 18 + 0.6*(time - 20)
#    elif (time <= 40):
#        G = 24 + 0.4*(time - 30)
#    elif (time <= 50):
#        G = 28 + 0.2*(time - 40)
#    else:
#        G = 30


#MO1_norm_dist = (b - M1_avg_M)/np.sqrt(MO1_var_M)
#MO2_norm_dist = (b - MO2_avg_M)/np.sqrt(MO2_var_M)

#axs[2].plot(MO1_t, st.norm.cdf(MO1_norm_dist), "deepskyblue", label = "M01 Normal CDF of blocking prob for b = 100")
#axs[2].plot(MO1_t, MO1_avg_block[5], "chartreuse", label = "M01 Average Block for b = 100")
#axs[2].legend()

#axs[3].plot(MO2_t, st.norm.cdf(MO2_norm_dist), "royalblue", label = "M02 Normal CDF of blocking prob for b = 100")
#axs[3].plot(MO2_t, MO2_avg_block[5], "springgreen", label = "M02 Average Block for b = 100")
#axs[3].legend()

gen_avg_var_sim()



