#BLOCKING PROBABILITY SIMULATION
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import scipy.stats as st
import sympy as sym
import scipy
import operator
import queue
import queueNoDonation_thrm_3_4 as bq_3_4
import queueNoDonation_thrm_3_5 as bq_3_5
import math
from scipy.stats import norm

#DEFINE GLOBAL VARIABLES
b = 200
#uniform dist 1: 100
#beta, deterministic, uniform 2 dist: 200 #expected bail amount
b_sq = 400*400/12 + 200*200 #400*400/12 + 200*200 #400*400/12 +
#deterministic: 200*200
# uniform 1 dist 200*200/12 + 100*100
# uniform 2 dist 400*400/12 + 200*200
# beta dist 400*400/12 + 200*200 
#Var[Beta] + E[Beta]^2
d = 0#50
d_sq = d*d
runs = 500#0#10000
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
time_period = 1000#0#0
inc = 1
lam_b = 2
lam_d = 1
#G_bar = 30
p = 0#.03
c = 50


#TRY:
#smaller incriments 


#lam_b*b*(1-p) - lam_d*d
#2*200*0.97-50*1

colors1 = ["r-", "y-", "g-", "c-", "b-", "m-"]
colors2 = ["pink", "violet", "lime", "dodgerblue", "darkblue", "indigo"]
color_pair = [["tomato", "red", "lightcoral"], 
              ["yellow", "gold", "lightyellow"], 
              ["palegreen", "lawngreen", "honeydew"], 
              ["aqua", "deepskyblue", "blue"], 
              ["violet", "blueviolet", "thistle"]]
bail_amt = ["200", "500", "750", "1000", "1500", "2000"] #CHANGED FROM 2500



#BAIL GENERATORS
def gen_deterministic_bail():
    return 200

def gen_uni_bail():
    return np.random.uniform(0,200)

def gen_uni2_bail():
    return np.random.uniform(0,400)

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

    
def solve_G(time):
    if (time < 0):
        print("ERROR")
        return -1000000
    elif (time <= 10):
        return 0
    elif (time <= 20):
        return 0.2*(time - 10)
    elif (time <= 30):
        return 2 + 0.4*(time - 20)
    elif (time <= 40):
        return 6 + 0.6*(time - 30)
    elif (time <= 50):
        return 12 + 0.8*(time - 40)
    else:
        return time - 30
    
def solve_G_bar(time):
    G = 0
    if (time <= 10):
        G = time
    elif (time <= 20):
        G = 10 + 0.8*(time - 10)
    elif (time <= 30):
        G = 18 + 0.6*(time - 20)
    elif (time <= 40):
        G = 24 + 0.4*(time - 30)
    elif (time <= 50):
        G = 28 + 0.2*(time - 40)
    else:
        G = 30
    return G
 
    
def solve_mu_sig(time_period):
    mu_t = np.zeros((time_period))
    sig_t = np.zeros((time_period))
    temp_sig_p = np.zeros((time_period))
    temp_sig_p2 = np.zeros((time_period))
    for time in range(1, time_period):
        G = solve_G(time)
        G_bar = solve_G_bar(time)
        mu_t[time] = lam_d*d*time - lam_b*b*time + (1-p)*lam_b*b*G
        #lam_d*d*time-lam_b*b*time + (1 - p)*lam_b*b*G
        sig_t[time] = math.sqrt(lam_d*d_sq*time + lam_b*b_sq*G_bar + lam_b*b_sq*(p**2)*G) #- 2*p
        #temp_sig_p2[time] = math.sqrt((lam_d*d_sq+lam_b*b_sq)*time + lam_b*b_sq*(p**2 - 1)*G)#- 2*p
    return mu_t, sig_t#, temp_sig_p#, temp_sig_p2
    
    
def solve_c_t(time_period, initial):
    c_t = np.zeros((time_period))
    mu_t, sig_t = solve_mu_sig(time_period)
    for time in range(1, time_period):
        #SHOULD CT BE ABLE TO BE NEGATIVE?
        c_t_temp = (b - initial - mu_t[time] - sig_t[time]*norm.ppf(1-err))/time
        print(c_t_temp)
        c_t[time] = max(0, c_t_temp) 
        # - lam_d*d*time + lam_b*b*time - lam_b*b*p*G
    return(c_t)

def solve_alt_c_t(time_period, initial):
    c_t = np.zeros((time_period))
    mu_t, sig_t = paper_mu_sig(time_period)
    for time in range(1, time_period):
        #SHOULD CT BE ABLE TO BE NEGATIVE?
        c_t_temp = (mu_t[time-1] - mu_t[time]) + (sig_t[time-1]- sig_t[time])*norm.ppf(1-err)#)/time
        print(c_t_temp)
        c_t[time] = max(0, c_t_temp) 
        # - lam_d*d*time + lam_b*b*time - lam_b*b*p*G
    return(c_t)

def solve_lam_d_t(time_period, initial):
    lam_d_t = np.zeros((time_period))
    report = time_period/10
    for time in range(1, time_period):
        if ((time % report) == 0):
            print(str(time*100 / time_period) + "%")
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
        temp_lam_d = np.arange(0.1, 10.1, 0.1)
        temp_diff = np.ones(len(temp_lam_d))*1000000
        #diff_not_change = 0
        for i in range(len(temp_lam_d)):
            sig_t = math.sqrt((temp_lam_d[i]*d_sq+lam_b*b_sq)*time + lam_b*b_sq*(p*p)*G)
            l_d = (b - initial - (sig_t*norm.ppf(1-err)) - c*time + lam_b*b*time - lam_b*b*p*G)/(d*time)
            #print(norm.ppf(err))
            temp_diff[i] =  abs(l_d - (temp_lam_d[i]))
            #if (abs(diff_ld) < temp_diff):
                #    temp_diff = abs(diff_ld)`1
            #    temp_lam_d = i
            #    diff_not_change = 0
            #if(time ==  2): 
            #    print("time = " + str(time) + " lam d = " + str(temp_lam_d) + " and diff " + str(temp_diff))
        #if (temp_diff == 10000000):
        #    lam_d_t[time] = 0
        #else:
        #fig, axs = plt.subplots(1, 1)
        #axs.set_xlabel("lam d")
        #axs.set_ylabel("the difference between rhs and lhs")
        #axs.plot(temp_lam_d, temp_diff, color_pair[0][0])
        #plt.show()
        min_l = np.argmin(temp_diff)
        #print(temp_diff[max(min_l - 3, 0):min_l + 3]) #np.hstack((temp_lam_d, 
        #print(temp_diff[min_l])
        lam_d_t[time] = temp_lam_d[min_l]
        #print(lam_d_t[time])
        #print(temp_diff)
    return(lam_d_t)
    
 #c approaches lam_b*b*(1-p) - lam_d*d
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
    

def bail_simulation(x, initial, min_amt, bail_generator, input_c_lam_d, theorem):
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
        if (theorem == 4):
            lam_d = 2
            t, M, O, Q, queue, block, total_blocked = bq_3_4.sim(initial, min_amt, 0, time_period, inc, lam_b, lam_d, bail_generator, input_c_lam_d, d, p)
        else: #theorem = 5            
            t, M, O, Q, queue, block, total_blocked = bq_3_5.sim(initial, min_amt, 0, time_period, inc, lam_b, input_c_lam_d, bail_generator, c, p)
        #initial, min_funding, expense, time_period, inc, lam_b, lam_d, gen_bail, c_t, p=0, a=0
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
    #print(avg_M)
    return t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M_t, tb


bail_gen = gen_beta_bail#gen_deterministic_bail#gen_beta_bail #gen_uni_bail#


#GRAPH GENERATORS:
initial = 1200#

def gen_avg_var_sim():
    fig, axs = plt.subplots(2, 1)
    form_E_M, form_Var_M = solve_mu_sig(time_period)
    form_E_M = form_E_M + np.ones(len(form_E_M))*initial        
    form_Var_M = form_Var_M**2
    #form_Var_M_p = form_Var_M_p**2
    print(form_Var_M)
    ct = np.zeros(time_period+1)
    t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M, tb = bail_simulation(b, initial, inf_funding, bail_gen, ct, 4)
    axs[0].set_xlabel("time, t")
    axs[0].set_ylabel("Expected Value of M(t)")
    axs[1].set_xlabel("time, t")
    axs[1].set_ylabel("Variance of M(t)")
    axs[0].plot(t, avg_M, "darkblue", label = "Average M[t] by simulation")
    axs[1].plot(t, var_M, "g-", label = "Variance of M[t] by simulation")
    axs[1].plot(np.arange(time_period), form_Var_M, "lime", label = "Variance of M[t] by formula")
    axs[0].plot(np.arange(time_period), form_E_M, "indigo", label = "Average M[t] by formula")
    #axs[1].plot(np.arange(time_period), form_Var_M_p, "lawngreen", label = "Variance of M[t] by my formula")
    axs[0].legend()
    axs[1].legend()
    plt.show()

def graph_c_t_test():
    fig, axs = plt.subplots(2, 2)
    t = np.arange(0, time_period+1)
    mu_t, sig_t = solve_mu_sig(time_period+1)
    ct = solve_c_t(time_period+1, initial)
    ct_alt = solve_alt_c_t(time_period+1, initial)
    cumulative_ct = np.cumsum(ct)
    cumulative_alt_ct = np.cumsum(ct_alt)
    #print(cumulative_ct)
    to_graph = np.zeros(time_period+1)
    to_consider = np.zeros(time_period+1)
    for time in range(1, len(to_graph)):
        print((b - initial - cumulative_ct[time] - mu_t[time])/sig_t[time])
        to_graph[time] = (b - initial - cumulative_ct[time] - mu_t[time])/sig_t[time]
        to_consider[time] = (b - initial - cumulative_alt_ct[time] - mu_t[time])/sig_t[time]
    #print(cumulative_ct)
    axs[0][0].set_xlabel("mu")
    axs[0][0].plot(t, mu_t, color_pair[0][0])
    #limit_c = np.ones(time_period+1)*norm.ppf(1-err)
    axs[1][0].plot(t, sig_t, color_pair[0][0])
    axs[1][0].set_xlabel("sig")
    #axs[1].plot(t, limit_c, color_pair[0][1], linestyle='dashed')
    axs[0][1].set_xlabel("ct")
    axs[0][1].plot(t, ct, color_pair[0][0])
    axs[1][1].plot(t, to_graph, color_pair[2][0])
    #axs[1][1].plot(t, to_consider, color_pair[1][1])
    print(to_consider)
    #print(to_graph)
    plt.show()
    
def graph_c_t():
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("time, t")
    axs.set_ylabel("cash drift, c(t)")
    t = arange(0, time_period+1)
    ct = solve_c_t(time_period+1, initial)
    limit_c = np.ones(time_period+1)*(lam_b*b*(1-p)-lam_d*d)
    axs.plot(t, ct, color_pair[0][0])
    axs.plot(t, limit_c, color_pair[0][1], linestyle='dashed')
    plt.show()

def graph_c_t_sim():
    ct = solve_c_t(time_period+1, initial)
    print(ct)
    print(sum(ct))
    t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M_t, tb = bail_simulation(b, initial, inf_funding, bail_gen, ct, 4)
    limit_c = np.ones(time_period+1)*(lam_b*b*(1-p)-lam_d*d)
    err_t = np.ones(time_period+1)*err
    d_err_t = np.ones(time_period+1)*err*2
    fig, axs = plt.subplots(3, 1)
    axs[0].set_xlabel("time, t")
    axs[0].set_ylabel("cash drift, c(t)")
    axs[1].set_xlabel("time, t")
    axs[1].set_ylabel("blocking probability")
    axs[2].set_xlabel("time, t")
    axs[2].set_ylabel("Money in fund, M(t)")
    axs[0].plot(t[:2000], ct[:2000], color_pair[3][0])
    axs[0].plot(t[:2000], limit_c[:2000], color_pair[4][1], linestyle='dashed')
    axs[1].plot(t[:2000], bl[:2000], color_pair[3][1])
    axs[1].plot(t[:2000], err_t[:2000], color_pair[4][1], linestyle='dashed')
    axs[2].plot(t[:2000], avg_M[:2000], color_pair[3][2])
    plt.show()
    fig, axs = plt.subplots(3, 1)
    #plt.title(str(runs) + " runs over " + str(time_period) + "units of time with " + str(bail_gen) + " bail distribution and initial capital "+str(initial), pad=10)
    axs[0].set_xlabel("time, t")
    axs[0].set_ylabel("cash drift, c(t)")
    axs[1].set_xlabel("time, t")
    axs[1].set_ylabel("blocking probability")
    axs[2].set_xlabel("time, t")
    axs[2].set_ylabel("Money in fund, M(t)")
    axs[0].plot(t, ct, color_pair[3][0])
    axs[0].plot(t, limit_c, color_pair[4][1], linestyle='dashed')
    axs[1].plot(t, bl, color_pair[3][1])
    axs[1].plot(t, err_t, color_pair[4][1], linestyle='dashed')
    axs[1].plot(t, d_err_t, color_pair[4][1], linestyle='dashed')
    axs[2].plot(t, avg_M, color_pair[3][2])
    
    plt.show()

def graph_lam_d_t():
    print((lam_b*b*(1-p)-c)/d)
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("time, t")
    axs.set_ylabel("arrival rate of donations, lam_d(t)")
    t = arange(0, time_period+1)
    lam_d_t = solve_lam_d_t(time_period+1, initial)
    print(lam_d_t)
    limit_lam_d = np.ones(time_period+1)*(lam_b*b*(1-p)-c)/d
    axs.plot(t, lam_d_t, color_pair[0][0])
    axs.plot(t, limit_lam_d, color_pair[0][1], linestyle='dashed')
    plt.show()

    
def graph_lam_d_sim():
    lam_d_t = solve_lam_d_t(time_period+1, initial)
    print(lam_d_t)
    t, avg_block, bl, avg_queue, var_queue, avg_M, avg_Q, var_M_t, tb = bail_simulation(b, initial, inf_funding, bail_gen, lam_d_t, 5)
    err_t = np.ones(time_period+1)*err
    fig, axs = plt.subplots(2, 1)
    axs[0].set_xlabel("time, t")
    axs[0].set_ylabel("blocking probability")
    axs[1].set_xlabel("time, t")
    axs[1].set_ylabel("Money in fund, M(t)")
    axs[0].plot(t, bl, color_pair[3][1])
    axs[0].plot(t, err_t, color_pair[4][1], linestyle='dashed')
    axs[1].plot(t, avg_M, color_pair[3][2])
    plt.show()
    
def cool_feature():
    ct_1 = np.zeros((time_period+1))
    ct_2 = np.zeros((time_period+1))
    for time in range(1, time_period+1):
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
        sig_t = math.sqrt((lam_d*d_sq+lam_b*b_sq)*time + lam_b*b_sq*(p*p-1)*G)
        ct_1[time] = (b - initial - sig_t*norm.ppf(1-err) - lam_d*d*time + lam_b*b*time - lam_b*b*p*G)/time
        ct_2[time] = (b - initial - sig_t*norm.ppf(1-err) - lam_d*(d - 20)*time + lam_b*b*time - lam_b*b*p*G)/time
        #if we have 20 less per time period, we should need difference*lam_d = 20* 2 = 40 more from constant funding:
    diff = ct_2-ct_1
    fig, axs = plt.subplots(1, 1)
    axs.set_xlabel("time, t")
    axs.set_ylabel("money")
    t = arange(0, time_period+1)
    axs.plot(t, ct_1, color_pair[3][0], label = "cash drift, c(t), with donation = 50")
    axs.plot(t, ct_2, color_pair[3][1], label = "cash drift, c(t), with donation = 30")
    axs.plot(t, diff, color_pair[4][0], label = "difference between cash drifts")
    axs.plot(t, np.ones(time_period+1)*40, color_pair[4][1], linestyle='dashed')
    axs.legend()
    print(diff)
    plt.show()
#THIS MEANS WE ARE ALSO ESSENTIALLY SOLVING FOR THE DONATION AMOUNT: it is just a function of lam_d and the difference so it is linearly related; duh by equation

#bail_simulation(x, initial, min_amt, bail_generator, input_c_lam_d, 5)
gen_avg_var_sim()
#graph_c_t_test()
#graph_c_t()
#graph_c_t_sim()
#cool_feature()

#graph_lam_d_t()
#graph_lam_d_sim()
#graph_lam_d_t()
    