from numpy import pi, sin, cos, sqrt
from math import atan2, atan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


def opt_tau_anal(e,w,W):    
    r  = sqrt(w*W*(1.0+e**2))
    th = atan(-sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    return r*cos(th) + 1j*(r*sin(th))

def J(e, w, W, tau=False):
    if not tau:
        tau = opt_tau_anal(e,w,W)
    r     = 0.5*np.sqrt(1.0 + (tau.real/tau.imag)**2)
    c1_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*w)/((w-tau.real)**2+(e*w+tau.imag)**2)
    cN_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*W)/((W-tau.real)**2+(e*W+tau.imag)**2)
    R     = np.sqrt(tau.real**2+tau.imag**2)*np.sqrt((e**2+1.0))/(2.0*abs(tau.real*e+tau.imag))
    C_im  = e*(tau.real**2+tau.imag**2)/(2.0*tau.imag*(tau.real*e+tau.imag))
    return np.sqrt(r**2/(R**2-C_im**2+2.0*C_im*c1_im))
  
def calc_circles(freq, tau, e):
    om  = 2.0*pi*freq*(1.0-1j*e)
    eta = om/(om-tau)
    #e   = -om[0].imag/om[0].real
    C   = 0.0 + 1j*( (e*abs(tau)**2)/(2.0*tau.imag*(tau.imag+e*tau.real)) )
    R   = sqrt( abs(tau)**2*(e**2+1.0)/(4.0*(tau.imag+e*tau.real)**2) )
    ck  = np.zeros((len(om),), dtype=complex)
    for k in range(0,len(om)):
        ck[k] = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
    r  = abs(tau/(tau-np.conj(tau)))
    return C, R, ck, r

def draw_circles(tau_re, tau_im, eps, freq):
    om  = 2.0*pi*freq*(1.0-1j*eps)
    NOP = 1000
    th  = np.linspace(0.0, 2.0*pi, NOP)
    
    C, R, c, r = calc_circles(freq, tau_re+1j*tau_im, eps)
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag

    ax.plot(X, Y, 'k')
    ax.plot(C.real, C.imag, 'kx', markersize=10)
    
    for k in range(0,len(om)):
        x = r*np.cos(th)+c[k].real
        y = r*np.sin(th)+c[k].imag
        ax.plot(x, y, col[k]+'--')
        ax.plot(c[k].real, c[k].imag, color=col[k], marker='x', markersize=10)
        
    ax.axhline(linewidth=0.5, color='k')
    ax.axvline(linewidth=0.5, color='k')
    

axis_color = 'lightgoldenrodyellow'
fig = plt.figure()

# Draw the plot
Nom  = 5
fmin = 1.0
fmax = 9.0
eps  = 0.5
freq = np.linspace(fmin,fmax,Nom)
om   = 2.0*np.pi*freq*(1.0-1j*eps)
tau = opt_tau_anal(eps,om[0].real,om[-1].real) 

cc  = list('grcmy')
col = list('b')
j = -1
for k in range(1,Nom-1):
    j=j+1
    if (j>4):
        j=0
    col.append(cc[j])
col.append('b')
    
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.25, bottom=0.25)
draw_circles(tau.real, tau.imag, eps, freq)
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.axis('equal')


# Add two sliders for tweaking the parameters
eps_slider_ax     = fig.add_axes([0.25, 0.05, 0.65, 0.03], axisbg=axis_color)
eps_slider        = Slider(eps_slider_ax, r'$\epsilon$', 0.0, 1.0, valinit=eps)
tau_re_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
tau_re_slider     = Slider(tau_re_slider_ax, r'Re($\tau$)', 0.0, om[-1].real, valinit=tau.real)
tau_im_slider_ax  = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axis_color)
tau_im_slider     = Slider(tau_im_slider_ax, r'Im($\tau$)', -om[-1].real, 0.0, valinit=tau.imag)
def sliders_on_changed(val):
    draw_circles(tau_re_slider.val, tau_im_slider.val, eps_slider.val, freq)
    fig.canvas.draw_idle()
eps_slider.on_changed(sliders_on_changed)
tau_re_slider.on_changed(sliders_on_changed)
tau_im_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.25, 0.1, 0.04])
reset_button = Button(reset_button_ax, r"Go to $\mathbf{\tau^\ast}$", color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    eps_slider.reset()
    tau_re_slider.reset()
    tau_im_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.7, 0.15, 0.2], axisbg=axis_color)
color_radios = RadioButtons(color_radios_ax, ('fmin = ', 'fmax = ', r'$N_f = $'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

#J_ax = fig.add_axes([0.25, 0.7, 0.15, 0.2], axisbg=axis_color)

plt.show()