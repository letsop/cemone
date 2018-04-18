#!/usr/bin/env python
""" Creates a 6 panel interactive figure to display outputs of CEMONE model"""
from cemone import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import interactive
import math
__author__ = "Marie Postel, Sorbonne Universite"
__copyright__ = "Copyright 2018, LJLL,SU"
__credits__ = ["Marie Postel", "Frederique Clement", "Sylvie Schneider-Maunoury",
               "Alice Karam"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marie Postel"
__email__ = "marie.postel@sorbonne_universite.fr"
__status__ = "Production"
# fix figure dimensions
fig, ax = plt.subplots(figsize=(10,6))
left_gamma=0.05
bottom_gamma=0.5
width_gamma=0.25
height_gamma=0.35

left_fap=0.35
bottom_fap=0.5
width_fap=0.25
height_fap=0.35

left_ip=0.05
bottom_ip=0.1
width_ip=0.25
height_ip=0.35

left_ipp=0.35
bottom_ipp=0.1
width_ipp=0.25
height_ipp=0.35

left_n=0.66
bottom_n=0.1
width_n=0.25
height_n=0.35

left_int=0.66
bottom_int=0.5
width_int=0.25
height_int=0.35

plt.subplots_adjust(left=left_gamma, bottom=bottom_gamma+0.4*height_gamma,right=left_gamma+width_gamma, top=bottom_gamma+height_gamma)

#fix default parameter values
Tmax=20
Tmin=10
t = np.arange(Tmin, Tmax, 0.1)
gamma_0 = 1
gamma_1 = 1
t_gamma= 14.65
s_gamma=2.53
beta_0 = 0.9
beta_1 = 0.9
t_beta= 14.
s_beta=1.
kapmax=40.
kap0=6.29
tp0=11.55
sp0=1.09
tm0=16.65
sm0=6.3

#Display gamma,beta and F_AP curves with default parameter values
(a0,b0,c0,d0)=(gamma_0,s_gamma,t_gamma,gamma_1)
sg = sigabcd(t,gamma_0,s_gamma,t_gamma,gamma_1)
plt.plot(t, sg, lw=2, color='red',label=r'$\gamma$ ref')
lgamma, = plt.plot(t, sg, lw=2, color='blue',label=r'$\gamma$')
(a0,b0,c0,d0)=(beta_0,s_beta,t_beta,beta_1)
sb = sigabcd(t,beta_0,s_beta,t_beta,beta_1)
plt.plot(t, sb, '--',lw=2, color='red',label=r'$\beta$ ref')
lbeta, = plt.plot(t, sb, '--',lw=2, color='blue',label=r'$\beta$')
plt.text(9.,1.05,'A', fontsize=10, fontweight='bold')
plt.text(21.,1.05,'B', fontsize=10, fontweight='bold')
plt.text(33.,1.05,'C', fontsize=10, fontweight='bold')
plt.text(9.,-0.9,'D', fontsize=10, fontweight='bold')
plt.text(21.,-0.9,'E', fontsize=10, fontweight='bold')
plt.text(33.,-0.9,'F', fontsize=10, fontweight='bold')
plt.legend(loc=6,prop={'size': 6})
plt.axis([Tmin, Tmax, 0., 1.])

s = fapt(t,kap0,tp0,sp0,tm0,sm0)
FAP = plt.axes([left_fap, bottom_fap+0.35*height_fap, width_fap, 0.65*height_fap], facecolor='white')
FAP.plot(t, s, lw=2, color='red', label=r'$F_{AP}$ ref')
lfap, = FAP.plot(t, s, lw=2, color='blue', label=r'$F_{AP}$')
FAP.legend(loc=2,prop={'size': 6})
FAP.axis([Tmin, Tmax, 0., kapmax])




# compute model output for default parameter values
var_arg=[kap0,tp0,sp0,tm0,sm0,gamma_0,s_gamma,t_gamma,gamma_1,beta_0,s_beta,t_beta,beta_1]
t_,IP1,IP2,N,F_AP_int=compute_sol(var_arg)

# display model output IP,IPP,N and F_AP_int for default parameter values
plt_IP=plt.axes([left_ip, bottom_ip, width_ip, height_ip], facecolor='white')
plt_IP.plot(t_, IP1+IP2, lw=2, color='red',label='IP ref')
l_IP, = plt_IP.plot(t_, IP1+IP2, lw=2, color='blue',label='IP')
plt.legend(loc=2,prop={'size': 6})
plt_IP.set_xlim([Tmin, Tmax])
maxIP=max(IP1+IP2)

plt_IPP=plt.axes([left_ipp, bottom_ipp, width_ipp, height_ipp], facecolor='white')
plt_IPP.plot(t_, IP1, lw=2, color='red',label='IPP ref')
l_IPP, = plt_IPP.plot(t_, IP1, lw=2, color='blue',label='IPP')
plt.legend(loc=2,prop={'size': 6})
plt_IPP.set_xlim([Tmin, Tmax])
maxIPP=max(IP1)

plt_N=plt.axes([left_n, bottom_n, width_n, height_n], facecolor='white')
plt_N.plot(t_, N, lw=2, color='red',label='N ref')
l_N, = plt_N.plot(t_, N, lw=2, color='blue',label='N')
plt.legend(loc=2,prop={'size': 6})
plt_N.set_xlim([Tmin, Tmax])
maxN=max(N)

plt_FAP_int=plt.axes([left_int, bottom_int, width_int, height_int], facecolor='white')
plt_FAP_int.plot(t_, F_AP_int, lw=2, color='red',label=r'$\int F_{AP}$ ref')
l_int, = plt_FAP_int.plot(t_, F_AP_int, lw=2, color='blue',label=r'$\int F_{AP}$')
plt.legend(loc=2,prop={'size': 6})
plt_FAP_int.set_xlim([Tmin, Tmax])
maxfapint=max(F_AP_int)

# display 4 sliders to modify the 4 parameters of gamma
xgamma=left_gamma+width_gamma/2
axcolor = 'lightgoldenrodyellow'
agxa0 = plt.axes([xgamma, bottom_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
agxb0 = plt.axes([xgamma, bottom_gamma+0.08*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
agxc0 = plt.axes([xgamma, bottom_gamma+0.16*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
agxd0 = plt.axes([xgamma, bottom_gamma+0.24*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)

sga0 = Slider(agxa0, r'$\gamma_0$', 0., 1., valinit=gamma_0)
sgb0 = Slider(agxb0, r'$s_\gamma$', 0., 10.0, valinit=s_gamma)
sgc0 = Slider(agxc0, r'$t_\gamma$', 10., 20.0, valinit=t_gamma)
sgd0 = Slider(agxd0, r'$\gamma_1$', 0., 1., valinit=gamma_1)

# display 4 sliders to modify the 4 parameters of beta
axcolor = 'lightgoldenrodyellow'
xbeta=0.02
abxa0 = plt.axes([xbeta, bottom_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
abxb0 = plt.axes([xbeta, bottom_gamma+0.08*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
abxc0 = plt.axes([xbeta, bottom_gamma+0.16*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)
abxd0 = plt.axes([xbeta, bottom_gamma+0.24*height_gamma, width_gamma/3, 0.03*height_gamma], facecolor=axcolor)

sba0 = Slider(abxa0, r'$\beta_0$', 0., 1., valinit=beta_0)
sbb0 = Slider(abxb0, r'$s_\beta$', 0., 10.0, valinit=s_beta)
sbc0 = Slider(abxc0, r'$t_\beta$', 10., 20.0, valinit=t_beta)
sbd0 = Slider(abxd0, r'$\beta_1$', 0., 1., valinit=beta_1)

# display 5 sliders to modify the 5 parameters of F_AP
axcolor = 'lightgoldenrodyellow'
xfap=left_fap
fxkap0 = plt.axes([xfap, bottom_fap, width_fap/4,0.03*height_fap], facecolor=axcolor)
fxtp0 = plt.axes([xfap, bottom_fap+0.1*height_fap, width_fap/4,0.03*height_fap], facecolor=axcolor)
fxsp0 = plt.axes([xfap, bottom_fap+0.2*height_fap, width_fap/4,0.03*height_fap], facecolor=axcolor)
fxtm0 = plt.axes([xfap+width_fap/2, bottom_fap+0.05*height_fap, width_fap/4,0.03*height_fap], facecolor=axcolor)
fxsm0 = plt.axes([xfap+width_fap/2, bottom_fap+0.15*height_fap, width_fap/4,0.03*height_fap], facecolor=axcolor)

fkap0 = Slider(fxkap0, r'$K_{AP}$', 1., kapmax, valinit=kap0)
ftp0 = Slider(fxtp0, r'$t_+$', 10., 14.0, valinit=tp0)
fsp0 = Slider(fxsp0, r'$s_+$', 0.1, 10.0, valinit=sp0)
ftm0 = Slider(fxtm0, r'$t_-$', 14., 20., valinit=tm0)
fsm0 = Slider(fxsm0, r'$s_-$', 0.1, 10., valinit=sm0)

# callback function (called if one parameter value is changed)
def update_param(val):
    gamma_0 = sga0.val
    s_gamma = sgb0.val
    t_gamma = sgc0.val
    gamma_1 = sgd0.val
    lgamma.set_ydata(sigabcd(t,gamma_0,s_gamma,t_gamma,gamma_1))
    beta_0 = sba0.val
    s_beta = sbb0.val
    t_beta = sbc0.val
    beta_1 = sbd0.val
    lbeta.set_ydata(sigabcd(t,beta_0,s_beta,t_beta,beta_1))
    kap= fkap0.val
    tp = ftp0.val
    sp = fsp0.val
    tm = ftm0.val
    sm = fsm0.val
    lfap.set_ydata(fapt(t,kap,tp,sp,tm,sm))
    var_arg=[kap,tp,sp,tm,sm,gamma_0,s_gamma,t_gamma,gamma_1,beta_0,s_beta,t_beta,beta_1]
    t_,IP1,IP2,N,F_AP_int=compute_sol(var_arg)
    l_IP.set_ydata(IP1+IP2)
    l_IPP.set_ydata(IP1)
    l_N.set_ydata(N)
    l_int.set_ydata(F_AP_int)
    plt_FAP_int.set_ylim([0, max(maxfapint,max(F_AP_int))])
    plt_N.set_ylim([0, max(max(N),maxN)])
    plt_IP.set_ylim([0, max(max(IP1+IP2),maxIP)])
    plt_IPP.set_ylim([0, max(max(IP1),maxIPP)])
    fig.canvas.draw_idle()

sga0.on_changed(update_param)
sgb0.on_changed(update_param)
sgd0.on_changed(update_param)
sgc0.on_changed(update_param)
sba0.on_changed(update_param)
sbb0.on_changed(update_param)
sbd0.on_changed(update_param)
sbc0.on_changed(update_param)
fkap0.on_changed(update_param)
ftp0.on_changed(update_param)
fsp0.on_changed(update_param)
ftm0.on_changed(update_param)
fsm0.on_changed(update_param)


# if the reset button is pressed, all parameter values are set back to their default values
resetax = plt.axes([0.8, 0.015, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sba0.reset()
    sbb0.reset()
    sbc0.reset()
    sbd0.reset()
    sga0.reset()
    sgb0.reset()
    sgc0.reset()
    sgd0.reset()
    fkap0.reset()
    ftp0.reset()
    fsp0.reset()
    ftm0.reset()
    fsm0.reset()
button.on_clicked(reset)

plt.show()
