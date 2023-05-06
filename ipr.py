# coding=utf-8
from __future__ import division 
import math
import numpy as np
from numpy import power
# import psapy.FluidProps as FluidProps
import matplotlib.pyplot as plt
import pandas as pd

def Darcy_IPR(k,h,visc, re,rw, s, P, OilFVF, nPoints):
    """Function to calculate IPR using Darcy's Equation.  It returns a list with a pair of Pressure and rates"""
    PwfList=[]
    QList=[]
    QList.append(0)
    PwfList.append(P)

    mStep=P/nPoints
    i=1

    while (i<=nPoints):
        
        Pwf=PwfList[i-1]-mStep
        Q= (k*h/visc)*(P-Pwf)/(141.2*OilFVF*visc*(math.log(re/rw)-0.75+s))
        
        QList.append(Q)
        PwfList.append(Pwf)

        i=i+1

    DarcyList=[QList,PwfList]

    return DarcyList

def VogelIPR(P, Pb, Pwf, Qo, nPoints):
    """Function to calculate IPR using Vogel's Equation.  It returns a list with a pair of Pressure and rates"""
    
    PwfList=[]
    QList=[]
    QList.append(0)
    PwfList.append(P)
    VogelList=[]
    mStep=P/nPoints
    i=1

    if Pwf>=Pb:
        J=Qo/(P-Pwf)
    else:
        J=Qo/((P-Pb)+((Pb/1.8)*(1-0.2*(Pwf/Pb)-0.8*(Pwf/Pb)**2)))

    while (i<=nPoints):   
        Pwfs=PwfList[i-1]-mStep
        
        if Pwfs>=Pb:
            Q=J*(P-Pwfs)
        else:
            
            Qb=J*(P-Pb)
            Q=Qb+(J*Pb/1.8)*(1-0.2*(Pwfs/Pb)-0.8*(Pwfs/Pb)**2)

        QList.append(Q)
        PwfList.append(Pwfs)

        i=i+1

    VogelList=[QList,PwfList]

    return VogelList

def Vogel_DarcyIPR(P, k,h,visc, re,rw, s, OilFVF, Pb, nPoints):
    """Function to calculate IPR using Vogel's Equation.  It returns a list with a pair of Pressure and rates"""
    
    PwfList=[]
    QList=[]
    QList.append(0)
    PwfList.append(P)
    VogelList=[]
    mStep=P/nPoints
    i=1
    # get j
    J= (k*h/visc)/(141.2*OilFVF*visc*(math.log(re/rw)-0.75+s))       
    # loop over the pressure points
    while (i<=nPoints):             
        Pwfs=PwfList[i-1]-mStep
        # abore the pupble point pressrue
        if Pwfs>=Pb:
            Q=J*(P-Pwfs)
        # blow ' two phase flow
        else:
            Qb=J*(P-Pb)
            Q=Qb+(J*Pb/1.8)*(1-0.2*(Pwfs/Pb)-0.8*(Pwfs/Pb)**2)
            
        QList.append(Q)
        PwfList.append(Pwfs)
        i=i+1
    VogelList=[QList,PwfList]
    return VogelList

# ----------------------- vogel with flow efficiensy ---------------------------
def vogel_fe_saturated(P, Pwf,Qo,fe,step):
    c = (1- (Pwf/P))
    b = (1.8 * fe * c ) - (0.8 * np.power(fe,2) * np.power(c,2) )
    Qmax = Qo /b
    PwfList=[]
    QList=[]
    QList.append(0)
    PwfList.append(P)
    last_pwf = P * (1 - (1/fe))
    last_pwf = 0 if last_pwf <0 else last_pwf
    for pwf in range(P,int(last_pwf),-step):
        cc = (1- (pwf/P))
        Q = Qmax*(1.8 * fe * (cc)  - 0.8 * np.power(fe,2) * np.power(cc,2) )
        QList.append(Q)
        PwfList.append(pwf)
    # the last point
    ccc = (1- (last_pwf/P))
    Q = Qmax*(1.8 * fe * (ccc)  - 0.8 * np.power(fe,2) * np.power(ccc,2) )
    PwfList.append(last_pwf)
    QList.append(Q)

    return QList , PwfList ,Qmax

def vogel_fe_saturated_fes(P,Pwf,Qo,fe,fes,step):
    # get the qmax first for the present fe
    c = (1- (Pwf/P))
    b = (1.8 * fe * c ) - (0.8 * np.power(fe,2) * np.power(c,2) )
    Qmax = Qo /b
    # constract the equation
    def Q(pwf,fe):
        ccc = (1- (pwf/P))
        Q = Qmax*(1.8 * fe * (ccc)  - 0.8 * np.power(fe,2) * np.power(ccc,2) )
        return Q 
    # loop over all fe having the pwf and get the corresponding Q
    Ps=[]
    Qs=[]
    for fe  in [fe,*fes]:
        PwfList=[]
        QList=[]
        QList.append(0)
        PwfList.append(P)
        last_pwf = P * (1 - (1/fe))
        last_pwf = 0 if last_pwf <0 else last_pwf
        for pwf in range(int(P),int(last_pwf),-step):
            cc = (1- (pwf/P))
            Q = Qmax*(1.8 * fe * (cc)  - 0.8 * np.power(fe,2) * np.power(cc,2) )
            QList.append(Q)
            PwfList.append(pwf)
        # the last point
        ccc = (1- (last_pwf/P))
        Q = Qmax*(1.8 * fe * (ccc)  - 0.8 * np.power(fe,2) * np.power(ccc,2) )
        PwfList.append(last_pwf)
        QList.append(Q)
        # append the final P and Q 
        Ps.append(PwfList)
        Qs.append(QList)

    return [Qs,Ps]

# ----------------------------------- voggel the general form --------------------------------------
def vogel_fe_undersaturated_fes(P,Pb,Pwf,Qo,fe,fes,step):
    c = (1-(Pwf/P))
    # get j for the present pwf and q
    if Pwf>=Pb:
        j = Qo/(P-Pwf)
    else:
        j= Qo/((P-Pb)+((Pb/1.8)*( 1.8 * c - .8* fe * np.power(c,2))))
        
    FEs = [fe,*fes]
    Js = [j]
    
    # get all the js for the FEs_test
    for i in range(1,len(FEs)):
        Js.append( (Js[i-1] * FEs[i]) / FEs[i-1] )
    
    # now we have all the Js and FEs
    # loop over all fe having the pwf and get the corresponding Q
    Ps=[]
    Qs=[]
    for j ,fe  in zip(Js,FEs):
        
        PwfList=[]
        QList=[]
        QList.append(0)
        PwfList.append(P)
        last_pwf = P * (1 - (1/fe))
        last_pwf = 0 if last_pwf <0 else last_pwf
        
        # iterate over the pressures
        for pwf in range(int(P),int(last_pwf),-step):
            # if we are in the undersaturated regrion
            if pwf>=Pb:
                Q=j*(P-pwf)
            else:
                Qb = j*(P-Pb)
                a = (j*Pb / 1.8)
                d = 1- (pwf / Pb )
                Q = Qb +a * (1.8*d -.8 * fe * d**2)
            QList.append(Q)
            PwfList.append(pwf)
            
        # the last point
        Qb = j*(P-Pb)
        a = (j*Pb / 1.8)
        d = 1- (last_pwf / Pb )
        Q = Qb +a * (1.8*d -.8 * fe * d**2)
        PwfList.append(last_pwf)
        QList.append(Q)
        
        
        # append the final P and Q 
        Ps.append(PwfList)
        Qs.append(QList)
        
    return [Qs,Ps]

#---------------------- futrue ipr ----------------------

def vogel_future_saturated(Pp,Pf,vis_pres,beta_pres,per_pres,vis_f,beta_f,per_f, Pwf, Qo, step):
    """Function to calculate IPR using Vogel's Equation.  It returns a list with a pair of Pressure and rates"""
    
    PwfList_pres=[]
    PwfList_f = []
    QList_presnt=[]
    Qlist_future = []
    QList_presnt.append(0)
    Qlist_future.append(0)
    PwfList_pres.append(Pp)
    PwfList_f.append(Pf)
    present_pressure_function = per_pres / (vis_pres * beta_pres)
    future_pressure_function = per_f / (vis_f * beta_f)
    
    Qmax_present =Qo/(1-0.2*(Pwf/Pp)-0.8*(Pwf/Pp)**2)
    Qmax_future = Qmax_present * ((Pf * future_pressure_function)/(Pp * present_pressure_function))

    for pwfp in range(int(Pp-step),-1,-step):
        PwfList_pres.append(pwfp)
        Q= Qmax_present * (1-0.2*(pwfp/Pp)-0.8*(pwfp/Pp)**2)
        QList_presnt.append(Q)
    for pwff in range(Pf-step,-1,-step):
        PwfList_f.append(pwff)
        Q= Qmax_future * (1-0.2*(pwff/Pf)-0.8*(pwff/Pf)**2)
        Qlist_future.append(Q)
        
    return [QList_presnt,PwfList_pres,"Present_IPR"],[Qlist_future,PwfList_f,"Future_IPR"]

# for undersaturated
def vogel_future_undersaturated(Pp,Pb,Pf,vis_pres,beta_pres,per_pres,vis_f,beta_f,per_f, Pwf, Qo, step):
    """Function to calculate IPR using Vogel's Equation.  It returns a list with a pair of Pressure and rates"""
    
    PwfList_pres=[]
    PwfList_f = []
    QList_presnt=[]
    Qlist_future = []
    QList_presnt.append(0)
    Qlist_future.append(0)
    PwfList_pres.append(Pp)
    PwfList_f.append(Pf)
    present_pressure_function = per_pres / (vis_pres * beta_pres)
    future_pressure_function = per_f / (vis_f * beta_f)
    if Pwf>=Pb:
        j_present = Qo/(Pp-Pwf)
    else:
        j_present = Qo/((Pp-Pb)+(1-0.2*(Pwf/Pp)-0.8*(Pwf/Pp)**2))
        
    j_future = j_present * (future_pressure_function / present_pressure_function)
    # get the present IPR
    for pwfp in range(Pp-step,-1,-step):
        if pwfp >=Pb :
            Q = j_present*(Pp-pwfp)
        else:
            Qb = j_present*(Pp-Pb)
            Q = Qb + (j_present*Pb/1.8) * (1-0.2*(pwfp/Pb)-0.8*(pwfp/Pb)**2)  
        PwfList_pres.append(pwfp)
        QList_presnt.append(Q)
        
    # get the future IPR  
    # check if the new pr in saturated
    if Pf <= Pb :
        # so for saturated
        for pwff in range(Pf-step,-1,-step):
            Q= (j_future * Pf / 1.8) * (1-0.2*(pwff/Pf)-0.8*(pwff/Pf)**2)
            # add the values
            Qlist_future.append(Q)
            PwfList_f.append(pwff)
    else: 
        # so it is under saturated
        for pwff in range(Pf-step,-1,-step):
            if pwff >= Pb:
                Q = j_future * (Pp - pwff)
            else:
                Qb = j_future * (Pf - Pb)
                Q = Qb + (j_future * Pb / 1.8) * (1 - 0.2 * (pwff / Pb) - 0.8 * (pwff / Pb) ** 2)
            # add the values
            PwfList_f.append(pwff)
            Qlist_future.append(Q)
            
    return [QList_presnt,PwfList_pres,"Present_IPR"],[Qlist_future,PwfList_f,"Futrue_IPR"]

# -------------------------- future for well with FE ---------------------------------------------
def vogel_future_fe_saturated(Pp,Pf,fe,vis_pres,beta_pres,per_pres,vis_f,beta_f,per_f, Pwf, Qo, step):
    c = (1- (Pwf/Pp))
    b = (1.8 * fe * c ) - (0.8 * np.power(fe,2) * np.power(c,2) )
    # calculate the Qmax_present , Qmax_future
    Qmax_present = Qo /b  # present Qmax (at fe=1)
    present_pressure_function = per_pres / (vis_pres * beta_pres)
    future_pressure_function = per_f / (vis_f * beta_f)
    Qmax_future = Qmax_present * ((Pf * future_pressure_function)/(Pp * present_pressure_function))
    # get current j
    j_present = Qo / ((Pp - Pwf) + ((Pwf / 1.8) * (1.8 * c - .8 * fe * np.power(c, 2))))
    # get the future j
    j_future = j_present * (future_pressure_function / present_pressure_function)
    # get fe future
    fe_f = fe * ( j_future / j_present)

    
    # the new Qmax for the future 
    PwfList_present=[Pp]
    QList_present=[0]
    PwfList_future = [Pf]
    QList_future = [0]
    last_pwf_presnt = Pp * (1 - (1/fe))
    last_pwf_future = Pf * (1 - (1/fe))
    last_pwf_presnt = 0 if last_pwf_presnt <0 else last_pwf_presnt 
    last_pwf_future = 0 if last_pwf_future <0 else last_pwf_future 
    
    # get values for present IPR
    for pwf in range(Pp,int(last_pwf_presnt),-step):
        cc = (1- (pwf/Pp))
        Q = Qmax_present*(1.8 * fe * (cc)  - 0.8 * np.power(fe,2) * np.power(cc,2) )
        QList_present.append(Q)
        PwfList_present.append(pwf)
    # the last point
    ccc = (1- (last_pwf_presnt/Pp))
    Q = Qmax_present*(1.8 * fe * (ccc)  - 0.8 * np.power(fe,2) * np.power(ccc,2) )
    PwfList_present.append(last_pwf_presnt)
    QList_present.append(Q)
    
    # get values for the future IPR
    for pwf in range(Pf,int(last_pwf_future),-step):
        cc = (1- (pwf/Pf))
        Q = Qmax_future*(1.8 * fe_f * (cc)  - 0.8 * np.power(fe_f,2) * np.power(cc,2) )
        QList_future.append(Q)
        PwfList_future.append(pwf)
    # the last point
    ccc = (1- (last_pwf_future/Pf))
    Q = Qmax_future*(1.8 * fe_f * (ccc)  - 0.8 * np.power(fe_f,2) * np.power(ccc,2) )
    PwfList_future.append(last_pwf_future)
    QList_future.append(Q)

    return [QList_present, PwfList_present,"Present_IPR"] ,[QList_future,PwfList_future,"Future_IPR"]

#---------------------------- future IPR FE!=1 and undersaturated
def vogel_future_fe_undersaturated(Pp,Pf,Pb,fe,vis_pres,beta_pres,per_pres,vis_f,beta_f,per_f, Pwf, Qo, step):
    PwfList_pres = [Pp]
    PwfList_f = [Pf]
    QList_presnt = [0]
    Qlist_future = [0]
    present_pressure_function = per_pres / (vis_pres * beta_pres)
    future_pressure_function = per_f / (vis_f * beta_f)
    c = (1 - (Pwf / Pp))
    # get present j
    if Pwf >= Pb:
        j_present = Qo / (Pp - Pwf)
    else:
        j_present = Qo / ((Pp - Pb) + ((Pb / 1.8) * (1.8 * c - .8 * fe * np.power(c, 2))))
    # get the future j
    j_future = j_present * (future_pressure_function / present_pressure_function)
    # get fe future
    fe_f = fe * (j_future / j_present)

    # get the last point of the pwf
    last_pwf_presnt = Pp * (1 - (1 / fe))
    last_pwf_future = Pf * (1 - (1 / fe))
    last_pwf_presnt = 0 if last_pwf_presnt < 0 else last_pwf_presnt
    last_pwf_future = 0 if last_pwf_future < 0 else last_pwf_future

    # get values for the present IPR
    for pwf in range(Pp, int(last_pwf_presnt), -step):
        # if we are in the undersaturated regrion
        if pwf >= Pb:
            Q = j_present * (Pp - pwf)
        else:
            Qb = j_present * (Pp - Pb)
            a = (j_present * Pb / 1.8)
            d = 1 - (pwf / Pb)
            Q = Qb + a * (1.8 * d - .8 * fe * d ** 2)
        QList_presnt.append(Q)
        PwfList_pres.append(pwf)
    # the last point
    Qb = j_present * (Pp - Pb)
    a = (j_present * Pb / 1.8)
    d = 1 - (last_pwf_presnt / Pb)
    Q = Qb + a * (1.8 * d - .8 * fe * d ** 2)
    PwfList_pres.append(last_pwf_presnt)
    QList_presnt.append(Q)

    # check first if the new Pr above or below the bubble point
    if Pf <= Pb :
        Qmax = (j_future * Pf ) / 1.8
        # get values for the future IPR
        for pwf in range(Pf,int(last_pwf_future),-step):
            cc = (1 - (pwf/Pf))
            Q = Qmax * (1.8 * fe_f * (cc)  - 0.8 * np.power(fe_f,2) * np.power(cc,2) )
            Qlist_future.append(Q)
            PwfList_f.append(pwf)
        # the last point
        ccc = (1- (last_pwf_future/Pf))
        Q = Qmax*(1.8 * fe_f * (ccc)  - 0.8 * np.power(fe_f,2) * np.power(ccc,2) )
        PwfList_f.append(last_pwf_future)
        Qlist_future.append(Q)
            
    else:
        # so it is under saturated
        for pwf in range(Pf, int(last_pwf_future), -step):
            # if we are in the undersaturated regrion
            if pwf >= Pb:
                Q = j_future * (Pf - pwf)
            else:
                Qb = j_future * (Pf - Pb)
                a = (j_future * Pf / 1.8)
                d = 1 - (pwf / Pb)
                Q = Qb + a * (1.8 * d - .8 * fe_f * d ** 2)
            Qlist_future.append(Q)
            PwfList_f.append(pwf)
        # the last point
        Qb = j_future * (Pf - Pb)
        a = (j_future * Pf / 1.8)
        d = 1 - (last_pwf_future/ Pb)
        Q = Qb + a * (1.8 * d - .8 * fe_f * d ** 2)
        PwfList_f.append(last_pwf_future)
        Qlist_future.append(Q)
    return [ QList_presnt , PwfList_pres ,"Present_IPR"] , [Qlist_future ,PwfList_f ,"Future_IPR"]

# processing functions 
def process_vlp(data):
    q ,p = data[0],data[1]
    p_min_index =0
    for i ,v in enumerate(p):
        if v == min(p):
            p_min_index = i 
    q=q[p_min_index:] 
    p=p[p_min_index:] 
    return q,p

# ---------------------------------------------------------- fetkovich ----------------------------------------
# class Fetkovic:
#     def __init__(self, Pr,Pb, test_points: list):
#         """
#         test_points : list of tuples of each is (q,pwf)
#         """
#         self.pb = Pb
#         self.pr = Pr
#         self.test_points = test_points
#         self.n = None
#         self.c = None
#         self.result = None

#     # get the n value from regression
#     def get_n(self):
#         # check number of test points
#         if len(self.test_points) == 1:
#             self.n = 1
#         else:
#             # separate the data
#             q = np.array([i[0] for i in self.test_points])
#             pwf = np.array([i[1] for i in self.test_points])
#             # square the difference
#             pr_pwf_2 = np.power(self.pr, 2) - np.power(pwf, 2)
#             # log values
#             log_q = np.log10(q)
#             self.log_q = log_q
#             log_p = np.log10(pr_pwf_2)
#             result = lr(x=log_q, y=log_p)
#             self.n = round((1 / result.slope), 4)
#             self.result = result

#     def get_c(self):
#         if not self.n:
#             raise Exception("Get the (n) value first ")
#         if len(self.test_points) == 1:
#             p = self.test_points[0][1]
#             q = self.test_points[0][0]
#             pp = (np.power(self.pr, 2) - np.power(p, 2))
#             self.c = q / pp
#         else:
#             a = self.result.slope
#             b = self.result.intercept
#             q = self.log_q.mean()
#             p = b + q * a
#             Pwf = (10 ** p)
#             q = (10 ** q)
#             c = q / (Pwf ** (self.n))
#             self.c = c

#     # draw the ipr
#     def draw(self):
#         Pwf = [self.pr]
#         Q = [0]
#         step = 20
#         for p in range(int(self.pr - step), 0, -20):
#             Pwf.append(p - 15)
#             pp = (np.power(self.pr, 2) - np.power(p, 2))
#             q = self.c * np.power(pp, self.n)
#             Q.append(q)
#         return [[Q, Pwf,"Present"]]

#     # main funcion
#     def main(self):
#         self.get_n()
#         self.get_c()
#         data = self.draw()
#         return data

# # --------------------------------------------- future for Fetkovich -----------------------------------
# class Future_Fetkovic(FetkovicSaturated):
#     def __init__(self, Pr, test_points: list,Prf):
#          # initalize the orginal model
#         super().__init__(Pr, test_points)
#         self.prf = Prf
#         self.new_c = None
#     def get_new_c(self):
#         if not self.c :
#             raise Exception("You must get c_present first")
#         self.new_c = self.c * (self.prf / self.pr)
#     def draw_future(self):
#         if not self.new_c:
#             raise Exception("Get the c_present value first")
#         Pwf = [self.prf]
#         Q = [0]
#         step = 20
#         for p in range(int(self.prf - step), 0, -20):
#             Pwf.append(p - 15)
#             pp = (np.power(self.prf, 2) - np.power(p, 2))
#             q = self.new_c * np.power(pp, self.n)
#             Q.append(q)
#         return [Q, Pwf]
#     def main(self):
#         self.get_n()
#         self.get_c()
#         self.get_new_c()
#         data1 = self.draw()
#         data2= self.draw_future()
#         return [[*data1[0],"Present"],[*data2,"Future"]]
    
# # ------------------- jones etal -----------------------------------------------
# class Jones:
#     def __init__(self, Pr, test_points: list):
#         """
#         test_points : list of tuples of each is (q,pwf)
#         """
#         self.pr = Pr
#         self.test_points = test_points
#         self.a = None
#         self.b = None
#         self.result = None

#     # get the n value from regression
#     def get_a_b(self):
#         # separate the data
#         q = np.array([i[0] for i in self.test_points])
#         pwf = np.array([i[1] for i in self.test_points])
#         # get ( pr - pwf) / q
#         y = (self.pr - pwf) / q
#         # fit the data
#         self.result = lr(x=q, y=y)
#         self.a = self.result.intercept
#         self.b = self.result.slope

#     # draw the ipr
#     def draw(self):
#         Pwf = [self.pr]
#         Q = [0]
#         step = 20
#         for p in range(int(self.pr - step), 0, -20):
#             Pwf.append(p - 15)
#             aa = np.power(self.a, 2) + 4 * self.b * (self.pr - p)
#             aa = np.power(aa, .5)
#             q = (-self.a + aa) / (2 * self.b)
#             Q.append(q)
#         return [Q, Pwf]

#     # main funcion
#     def main(self):
#         self.get_a_b()
#         data = self.draw()
#         return data
