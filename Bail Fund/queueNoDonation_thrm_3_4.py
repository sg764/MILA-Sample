import matplotlib.pyplot as plt
import numpy as np
import operator
import queue
import sys

#CLASSES
class donation:
    kind = "d"
    def __init__(self, time, donation_amt):
        self.t = time
        self.d = donation_amt

class bail_request:
    kind = "br"
    def __init__(self, time, t_to_trial, wait_time, bail_amt):
        self.t = time
        self.s = time + t_to_trial #saves overall court time 
        self.w = 0
        self.b = bail_amt
        self.accepted = 'f'
        #f = False, never reached, w = Waiting in queue, l = loaned money t, t = True, loaned and returned, r = waited but was never helped 

#PRINT CLASSES:
def print_d(d):
    print(" Donation at time " + str(d.t) + " of amount " + str(d.d))

def print_b(b):
    if(b.accepted == 'f'): 
        print(" Bail Request at time " + str(b.t) + " with trial at time " + str(b.s) + " of amount " + str(b.b) + " was not accepted" )
    elif (b.accepted == 'w'): 
        print(" Bail Request at time " + str(b.t) + " with trial at time " + str(b.s) + " of amount " + str(b.b) + " is in the queue and will wait until " + str(b.w))
    elif (b.accepted == 'l'): 
        print(" Bail Request at time " + str(b.t) + " with trial at time " + str(b.s) + " of amount " + str(b.b) + " is currently loaned the bail money")
    elif (b.accepted == 't'): 
        print(" Bail Request at time " + str(b.t) + " with trial at time " + str(b.s) + " of amount " + str(b.b) + " was loaned the bail money and returned it")    
    else:
        print(" Bail Request at time " + str(b.t) + " with trial at time " + str(b.s) + " with wait time " + str(b.w) + " of amount " + str(b.b) + " waited but didn't get helped in time")

def print_events(e):
    for i in range(len(e)):
        if (e[i].kind == "d"):
            print_d(e[i])
        else:
            print_b(e[i])

#GENERATORS:
def get_event(events, e):
    i = 0
    if (e.kind == "d"):
        while (i < len(events)):
            if (events[i].kind == e.kind and events[i].t == e.t and events[i].d == e.d):
                return i
            i+=1
    else:
        while (i < len(events)):
            if (events[i].kind == e.kind and events[i].t == e.t and events[i].b == e.b and events[i].s == e.s):
                return i
            i+=1
    return -1

def gen_time_to_trial(prob = [0.2, 0.2, 0.2, 0.2, 0.2]):
    #t = np.random.uniform(30,180)
    #if (t < 60):
    #    t += 50 #makes the P(80 < t < 110) is twice P(110 < t < 180) 
    return np.random.choice([10, 20, 30, 40, 50], p=prob)

def gen_wait_time(t, s):
    #w = (580 + t - s)/500
    #return min(t*np.random.normal(w, 0.2, 1), t)
    return s

    
def gen_poisson_process(lam, max):
    p = []
    p.append(np.random.exponential(float(1/float(lam))))
    t = []
    t.append(p[0])
    temp_p = np.random.exponential(float(1/float(lam)))
    while (t[len(t)-1] + temp_p < max):
        p.append(temp_p)
        t.append(t[len(t)-1] + p[len(p)-1])
        temp_p = np.random.exponential(float(1/float(lam)))
    return p, t

#SIMULATION:
def sim(initial, min_funding, expense, time_period, inc, lam_b, lam_d, gen_bail, c_t, don_amt, p=0, a=0):
    #initial = initial()
#GENERATE BAIL REQUESTS
    #The poisson distribution of the time between bail requests
    Na_b, Nt_b = gen_poisson_process(lam_b, time_period)
    Na_d, Nt_d = gen_poisson_process(lam_d, time_period)
    br = []
    for i in range(len(Na_b)):
        s = gen_time_to_trial()
        br.append(bail_request(Nt_b[i],  s, gen_wait_time(Nt_b[i], s), gen_bail()))
    events = br #sorted(br, key=operator.attrgetter('t'))
    donations = []
    for j in range(len(Na_d)):
        donations.append(donation(Nt_d[j], don_amt))
        #print(donations[j].d)
#SIMULATE TIME
    MT = [initial] + [None]*int(time_period/inc)
    QT = [0] + [None]*int(time_period/inc)
    OT = [0] + [None]*int(time_period/inc)
    Queue = [0] + [None]*int(time_period/inc)
    blocked = np.zeros((6, 1+int(time_period/inc))) #200, 500, 750, 1000, 1500, 2000
    times = np.arange(0, time_period+inc, inc) 
    cur_index = 1
    cur_bail = []
    bail_queue = [] #queue.Queue()
    total_funds = initial
    time = inc 
    i = 0
    d = 0
    c_total = 0
    donations_total = 0
    total_blocked = 0
    def queue_sort(e):
        #a = 0, acts like STT, a = 1, acts like LBRA
        return (e.b**(2*a))*(e.s-e.t)**(2-2*a)
    while (time <= time_period):
        #UPDATE FUNDS:
        total_funds += c_t[time]
        c_total += c_t[time]
        #UPDATE CURRENT BAIL
        while(0 < len(cur_bail) and cur_bail[0].s < time):
            ev = cur_bail.pop(0)
            #print(ev.b)
            total_funds += (ev).b*(1-p)
            events[get_event(events, ev)].accepted = 't'
            #temp += 1
        #UPDATE BY NORMAL EVENTS
        #ADD DONATION:
        while(d < len(donations) and donations[d].t < time):
            #print(donations[d].d)
            total_funds += donations[d].d
            
            donations_total += donations[d].d
            d += 1
        while (i < len(events) and events[i].t < time):
            if (total_funds - events[i].b >= min_funding) and (len(bail_queue) == 0):
                events[i].accepted = 'l'
                total_funds -= events[i].b
                cur_bail = sorted(cur_bail + [events[i]], key= operator.attrgetter('s'))
            else:
                #print("REJECTED BAIL " + str(events[i].b))
                events[i].accepted = 'w'
                total_blocked+= 1
                bail_queue = bail_queue + [events[i]]
                if (a >= 0):
                    bail_queue = sorted(bail_queue, key = queue_sort)
            i+=1
        
        #UPDATE BAIL QUEUE
        temp = 0
        #CHECKING IF ANY PEOPLE ARE DONE WAITING
        while(temp < len(bail_queue)):
            if (bail_queue[temp].w < time):
                ev = bail_queue.pop(temp)
                events[get_event(events, ev)].accepted = 'r'
            else:
                temp += 1
        #CHECK IF PEOPLE CAN BE LOANED
        while (len(bail_queue) != 0) and (total_funds - bail_queue[0].b >= min_funding):
            events[get_event(events, bail_queue[0])].accepted = 'l'
            total_funds -= bail_queue[0].b
            cur_bail.append(bail_queue[0])
            ev = bail_queue.pop(0)
        #print("at time " + str(events[i].t) + " fund has $" + str(total_funds))
        #SAVE FUNCTION VALUES
        MT[round(time/inc)] = total_funds
        QT[round(time/inc)] = len(cur_bail)
        Queue[round(time/inc)] = len(bail_queue)
        OT[round(time/inc)] = 0
        for j in range(len(cur_bail)):
            OT[round(time/inc)] += cur_bail[j].b
        blocked[0, round(time/inc)] = 1 if (total_funds < 100) else 0 
        blocked[1, round(time/inc)] = 1 if (total_funds <  200) else 0
        blocked[2, round(time/inc)] = 1 if (total_funds <  750) else 0
        blocked[3, round(time/inc)] = 1 if (total_funds < 1000) else 0
        blocked[4, round(time/inc)] = 1 if (total_funds < 1500) else 0
        blocked[5, round(time/inc)] = 1 if (total_funds < 2000) else 0
        
        time+=inc   
    #print(donations_total)
    #print(c_total)
    return times, MT, OT, QT, Queue, blocked, total_blocked/len(Na_b)

#initial, min_funding, expense, time_period, inc, d_c, d_t, lam_b, a=0.5
#times, MT, OT, QT, queue, b = sim(2000, 0, 0, 20, 0.1, 500, 1, 1, 0.5)
#np.set_printoptions(threshold=sys.maxsize)
#\print(b[3])
#for i in range(5):
#    print(b[i+1]-b[i])
#print(MT)
#print(queue)
#print(times)
#for i in range(len(queue)):
#    if (queue[i] != 0):
#        print(times[i])
#plt.step(times, QT)
#plt.show()
#AVERAGE: 
#avg_O = 0
#avg_Q = 0
#avg_queue = 0
#for i in range(1000):
#    t, M, O, Q, queue = sim(2000, 0, 0, 100, 0.1, 3000, 1, 0.5, 1)
    #initial, min_funding, expense, time_period, inc, d_c, d_t,, lam_b
#    avg_O += 0.0001*np.array(O)
#    avg_Q += 0.0001*np.array(Q)
#    avg_queue += 0.0001*np.array(queue)
#plt.step(t, avg_Q)
#plt.show()
#print(avg_O)    
#print(avg_queue)



