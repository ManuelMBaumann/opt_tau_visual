from numpy import pi, sin, cos, sqrt
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
  
def calc_circles(om, tau):
    eta = om/(om-tau)
    e   = -om[0].imag/om[0].real
    C   = 0.0 + 1j*( (e*abs(tau)**2)/(2.0*tau.imag*(tau.imag+e*tau.real)) )
    R   = sqrt( abs(tau)**2*(e**2+1.0)/(4.0*(tau.imag+e*tau.real)**2) )
    ck  = np.zeros((len(om),))
    for k in range(0,len(om)):
        ck[k] = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
    r  = abs(tau/(tau-np.conj(tau)))
    return C, R, ck, r

def signal(amp, freq):
    return amp * sin(2 * pi * freq * t)

axis_color = 'lightgoldenrodyellow'
fig = plt.figure()

# Draw the plot
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(-5, 5, 0.01)
amp_0 = 5
freq_0 = 3
[line] = ax.plot(t, signal(amp_0, freq_0), linewidth=2, color='red')
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.axis('equal')

# Add two sliders for tweaking the parameters
amp_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
amp_slider = Slider(amp_slider_ax, r'Re($\tau$)', 0.1, 10.0, valinit=amp_0)
freq_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axis_color)
freq_slider = Slider(freq_slider_ax, r'Im($\tau$)', 0.1, 30.0, valinit=freq_0)
def sliders_on_changed(val):
    line.set_ydata(signal(amp_slider.val, freq_slider.val))
    fig.canvas.draw_idle()
amp_slider.on_changed(sliders_on_changed)
freq_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, r"Go to $\mathbf{\tau^\ast}$", color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    freq_slider.reset()
    amp_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.7, 0.15, 0.2], axisbg=axis_color)
color_radios = RadioButtons(color_radios_ax, (r'$\epsilon = $', 'fmin = ', 'fmax = ', r'$N_f = $'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

#J_ax = fig.add_axes([0.25, 0.7, 0.15, 0.2], axisbg=axis_color)

plt.show()