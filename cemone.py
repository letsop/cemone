#!/usr/bin/env python
""" computes the outputs of CEMONE model by the characteristic method"""
from math import *
from numpy import *
__author__ = "Marie Postel, Sorbonne Universite"
__copyright__ = "Copyright 2018, LJLL,SU"
__credits__ = ["Marie Postel", "Frederique Clement", "Sylvie Schneider-Maunoury",
               "Alice Karam"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marie Postel"
__email__ = "marie.postel@sorbonne_universite.fr"
__status__ = "Production"
def sigmoid(t,a,b,c,d=None):
    """ return the value of sigmoid for t """
    # t, time
    # a,b,c,d sigmoid parameters
    if d==None:
        d=0.
    return a+(d-a)/(1.+exp(-b*(t-c)))

def sigabcd(t,a,b,c,d):
    return array([sigmoid(ti,a,b,c,d) for ti in t])

def fap(t,kap,tp,sp,tm,sm):
    return kap*(1.-sigmoid(t,1.,sp,tp,0.))*sigmoid(t,1.,sm,tm,0.)

def fapt(t,kap,tp,sp,tm,sm):
    return array([fap(ti,kap,tp,sp,tm,sm) for ti in t])

def compute_IP_sums(IP1a,IP2a,ca_2,da_1,da_2,Na_1,Na_2):
    """ integrate age-time solution over cell cycle to obtain time cell counts"""
    IP1=da_1*(0.5*IP1a[0]+0.5*IP1a[Na_1]+sum(IP1a[1:Na_1]))
    IP2=da_2*(0.5*IP2a[0]+0.5*IP2a[Na_2]+sum(IP2a[1:Na_2])+0.5*(1.-ca_2)*IP2a[Na_2]-0.5*ca_2*IP2a[Na_2-1])
    
    return IP1,IP2


def prepare_arrays(T1_cycle,T1M,T2_cycle,Nm,Tstart,Tend):
    
    Na_1=int(floor(Nm*T1_cycle/T1M))
    da_1=T1_cycle/Na_1
    T = Tend-Tstart # simulation length (day)
    dt=da_1
    Nt=int(floor(T*24./dt))+1
    t_ = Tstart*ones(Nt+1)+dt*arange(Nt+1)/24.

    da_2=dt
    Na_2=int(floor(T2_cycle/da_2))+1
    ca_2=(Na_2*da_2-T2_cycle)/da_2
    
    a1_=da_1*arange(Na_1+1)
    a2_=da_2*arange(Na_2+1)
    a2_[Na_2]=T2_cycle
    
    IP1  = zeros(Nt+1)
    IP1a = zeros(Na_1+1)
    IP1w = zeros(Na_1+1)
    
    
    IP2  = zeros(Nt+1)
    IP2a = zeros(Na_2+1)
    IP2w = zeros(Na_2+1)
    
    N  = zeros(Nt+1)
    N_AP  = zeros(Nt+1)
    alpha_  = zeros(Nt+1)
    beta_   = zeros(Nt+1)
    gamma_  = zeros(Nt+1)
    delta_  = zeros(Nt+1)
    
    return dt,da_1,da_2,Nt,Na_1,Na_2,ca_2,t_,alpha_,beta_,gamma_,delta_,a1_,a2_,IP1,IP1a,IP1w,IP2,IP2a,IP2w,N,N_AP

def prepare_param_and_arrays(val_arg):
    [a_alpha,b_alpha,c_alpha,d_alpha,a_beta,b_beta,c_beta,d_beta,a_delta,b_delta,c_delta,d_delta,a_gamma,b_gamma,c_gamma,d_gamma,d_dAPb,Tstart,Tend,T1G1,T1S,T1G2,T1M,T2G1,T2S,T2G2,T2M,Nm]=val_arg
    Nm = int(Nm)   # number of mesh points in the shorter phase (usually M)
    T1_cycle =T1G1+T1S+T1G2+T1M
    T2_cycle =T2G1+T2S+T2G2+T2M
    
    list_arrays=prepare_arrays(T1_cycle,T1M,T2_cycle,Nm,Tstart,Tend)
    return list_arrays


def compute_time_solution(val_arg,list_arrays):
    """ compute solution with characteristic method (See Algorithm in [PKPSMC])"""
    [alpha_0,s_alpha,t_alpha,alpha_1,beta_0,s_beta,t_beta,beta_1,delta_0,s_delta,t_delta,delta_1,gamma_0,s_gamma,t_gamma,gamma_1,d_dAPb,Tstart,Tend,T1G1,T1S,T1G2,T1M,T2G1,T2S,T2G2,T2M,Nm]=val_arg
    dt,da_1,da_2,Nt,Na_1,Na_2,ca_2,t_,alpha_,beta_,gamma_,delta_,a1_,a2_,IP1,IP1a,IP1w,IP2,IP2a,IP2w,N,N_AP=list_arrays
    alpha_=sigabcd(t_,alpha_0,s_alpha,t_alpha,alpha_1)
    beta_=sigabcd(t_,beta_0,s_beta,t_beta,beta_1)
    gamma_=sigabcd(t_,gamma_0,s_gamma,t_gamma,gamma_1)
    delta_=sigabcd(t_,delta_0,s_delta,t_delta,delta_1)

    flow_AP = d_dAPb*delta_[0]*(1.-alpha_[0])*beta_[0]
    f_AP1 = gamma_[0]*flow_AP
    f_AP2 = (1.-gamma_[0])*flow_AP
    i=0
    IP1a[0]=f_AP1
    IP2a[0]=f_AP2+2*IP1a[Na_1]
    N[0]=(1.-alpha_[0])*(1.-beta_[0])*delta_[0]*d_dAPb
    N_AP[0]=(1.-alpha_[0])*(1.-beta_[0])*delta_[0]*d_dAPb
    
    for i in range(Nt):
        [IP1[i],IP2[i]]=compute_IP_sums(IP1a,IP2a,ca_2,da_1,da_2,Na_1,Na_2)
        t=t_[i+1]
        gamma = gamma_[i+1]
        beta  = beta_[i+1]
        alpha = alpha_[i+1]
        delta = delta_[i+1]
        delta_IP=1.
        
        flow_AP = d_dAPb*delta*(1.-alpha)*beta
        f_AP1 = gamma*flow_AP
        f_AP2 = (1.-gamma)*flow_AP
        
        IP1w[1:Na_1+1]=IP1a[0:Na_1]
        IP1w[0]=f_AP1
        
        IP2w[0]=f_AP2+2*IP1w[Na_1]
        IP2w[1:Na_2+1]=IP2a[0:Na_2]
        IP2w[Na_2]=ca_2*IP2a[Na_2-1]+(1.-ca_2)*IP2a[Na_2-2]
        
        IP1a=IP1w
        IP2a=IP2w
        N[i+1]=N[i]+dt*((1.-alpha)*(1.-beta)*delta*d_dAPb+2*IP2a[Na_2]*delta_IP)
        N_AP[i+1]=N_AP[i]+dt*(1.-alpha)*(1.-beta)*delta*d_dAPb

    i=Nt
    [IP1[i],IP2[i]]=compute_IP_sums(IP1a,IP2a,ca_2,da_1,da_2,Na_1,Na_2)
    
    return [IP1,IP2,N]
def compute_sol(input_arg):
    T1G1=20.9
    T1S=6.4
    T1G2=1.6
    T1M=0.5
    
    T2G1=21.3
    T2S=2.8
    T2G2=1.6
    T2M=0.5

    alpha_0=1.
    alpha_1=0.
    delta_0=1.
    delta_1=0.
    
    [K_AP,t_alpha,s_alpha,t_delta,s_delta,gamma_0,s_gamma,t_gamma,gamma_1,beta_0,s_beta,t_beta,beta_1]=input_arg
    
    Nm=5
    Tstart=9.
    Tend=20.
    val_arg=[alpha_0,s_alpha,t_alpha,alpha_1,beta_0,s_beta,t_beta,beta_1,delta_0,s_delta,t_delta,delta_1,gamma_0,s_gamma,t_gamma,gamma_1,K_AP,Tstart,Tend,T1G1,T1S,T1G2,T1M,T2G1,T2S,T2G2,T2M,Nm]
    T1C=T1G1+T1S+T1G2+T1M
    T2C=T2G1+T2S+T2G2+T2M
    list_arrays=prepare_param_and_arrays(val_arg)
    t_=list_arrays[7]
    dt=list_arrays[0]
    alpha_=list_arrays[8]
    delta_=list_arrays[11]
    
    [IP1,IP2,N]=compute_time_solution(val_arg,list_arrays)
    
    ntim=len(alpha_)
    t_=linspace(Tstart,Tend,ntim)
    alpha_=sigabcd(t_,alpha_0,s_alpha,t_alpha,alpha_1)
    delta_=sigabcd(t_,delta_0,s_delta,t_delta,delta_1)
    F_AP_=K_AP*(1.-alpha_)*delta_
    F_AP_int=[0]
    for i in range(len(t_)-1):
        F_AP_int.append(F_AP_int[i]+dt*F_AP_[i])
    return t_,IP1,IP2,N,F_AP_int







