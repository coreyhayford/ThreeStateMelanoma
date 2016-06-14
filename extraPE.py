# ic_result = []
# for i,p in enumerate(model.parameters[:3]):
#     ic_result.append([])
#     print p,
#     for n in range(len(ICs())):
#         print result.x[i+n*3+7],
#         ic_result[-1].append(result.x[i+n*3+7])
#     print
# xlist = [1,2,3,4,5]
# x_array = np.array(xlist)
# ylist = [6,7,8,9,10,11,12]
# y_array = np.array(ylist)
# 
# x_y_list = [x_array, y_array]
# print x_y_list
# print type(x_y_list)
# x_y_array = np.array(x_y_list)
# print x_y_array
# print type(x_y_array)

 
#print xopt
#print type(xopt)



### testing ideas ###


# print x_array + y_array
# print (x_array + y_array) / 2
# print 3/2
# print [(0,None)]*10


#def mean(all_lists):
#    return sum(all_lists) / len(all_lists)
#all_lists = [nl2_rep1, nl2_rep2, nl2_rep3]
#average_lists = map(mean, zip(*all_lists))
#plt.plot(time, average_lists, ms = 12, mfc = "0.50")

#for i in range(len(time)):
#    print time[i], solver.yobs["Obs_All"][i]



#sum_lists = [sum(i) for i in zip(*all_lists)]
#print all_lists
#print sum_lists
#array = np.array([[nl2_rep1], [nl2_rep2], [nl2_rep3]])
#rep_average = np.mean(array, axis = 0)
#print rep_average
#plt.plot(time, rep_average, ms = 12, mfc = "0.01")

#def sublist(list, n):
#    for i in range(1, len(list), n):
#        yield list[i: i+n]

#print sublist(time, 10)

#exp_data = np.array([[time[i], nl2[i]] for i in range(1, len(time)) if time[i] < 10])
#print exp_data

#for i in range(1, len(time), 10):
#    time_new = time[i, i+10],
#    print time_new

#def sublist(exp_data,10):
#    sub = []; result = []
#    for i in exp_data:
#        sub+=[i]
#        if len(sub)==n: result+=[sub]; sub =[]
#    if sub:result += [sub]
#    return result    
        
#cell_count = data_parental_treated["Cell.Nucleus"]
#print cell_count

# def eq(par,initial_cond,start_t,end_t,incr):
#     # Create time space
#     t  = np.linspace(start_t, end_t,incr)
#     # Create ODE system
#     def funct(y,t):
#        A=y[0]
#        B=y[1]
#        C=y[2]
#        par = k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC
#        # the model equations from PySB model
#        f0 = -A*k_AB - A*k_death_A + A*k_div_A + B*k_BA
#        f1 = A*k_AB - B*k_BA - B*k_BC - B*k_death_B + B*k_div_B + C*k_CB
#        f2 = B*k_BC - C*k_CB -C*k_death_C + C*k_div_C
#        return [f0, f1, f2]
     # Integrate over the equations
#     ds = integrate.odeint(funct,initial_cond,t)
#     print ds
#     return (ds[:,0],ds[:,1],ds[:,2],t)

#rates = (k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC)

#### Define Initial Conditions ###
#A_0 = 1000
#B_0 = 1000
#C_0 = 1000
#y_0 = [A_0, B_0, C_0]

#### Redefine model steps - encompassed in PySB model ###
#start_time=0.0
#end_time=120.0
#intervals=10
#mt=np.linspace(start_time,end_time,intervals)

#### Model index to compare to data ###
#findindex=lambda x:np.where(mt>=x)[0][0]
#mindex=map(findindex,time)
#print mindex

# time_sim = [time[0]]
# 
# for i in range(1, len(time)):
#     if time[i] < time[i-1]:
#         break
#     time_sim.append(time[i])

# solver.run # saved in array - species y, obs in yobs, 