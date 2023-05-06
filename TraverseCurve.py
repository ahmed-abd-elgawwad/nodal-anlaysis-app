import BeegsAndBrell as BB
import numpy as np
import matplotlib.pyplot as plt

def pressure_traverse(liquid_rate,depth,tht=150,twf=150,glr=0,wc=0,gas_grav =0.65,oil_grav = 35,wtr_grav = 1.07,diameter=1.922,angle=90.0,
                      margin=30,thp=0,sample_size =100):
    depths = np.linspace(0, depth, sample_size)
    t_grad = temp_gradient(tht, twf, depth)
    temps = tht + t_grad * depths
    thp+=margin
    water_rate = liquid_rate * wc
    oil_rate = liquid_rate * (1 - wc)
    gor = glr / (1 - wc)
    p = []
    dpdz = []
    for i in range(len(depths) - 1):
        if i == 0:
            p.append(thp)

        if len(dpdz) == 0:
            pavg = p[0] + .2 * (depths[1] - depths[0])
        else:
            pavg = p[i] + dpdz[i - 1] * .5 * (depths[i + 1] - depths[i])

        while (True):
            dpdz_step = BB.Pgrad(pavg, temps[i], oil_rate, water_rate, gor, gas_grav, oil_grav, wtr_grav, diameter,
                                  angle)
            if (p[i] + .5 * dpdz_step * (depths[i + 1] - depths[i]) - pavg < .1):
                break
            else:
                pavg = p[i] + .5 * dpdz_step * (depths[i + 1] - depths[i])

        dpdz.append(dpdz_step)
        dz = (depths[i + 1] - depths[i])
        pressure = p[i] + dz * dpdz[i - 1]
        p.append(pressure)
    p[:] = [number - margin for number in p]
    return p,depths

def temp_gradient(t0,t1, depth):
    if depth==0:
        return 0  
    else:
        return abs(t0-t1)/depth

def plot_traverse(liquid_rate,depth,ax=None,figure=None,tht=150,twf=150,glr=0,wc=0,gas_grav =0.65,oil_grav = 35,wtr_grav = 1.07,diameter=1.922,angle=90.0,
                      margin=30,thp= 0,sample_size =100,rate=False,G=False,Wc=False):
    p, depths = pressure_traverse(liquid_rate=liquid_rate, depth=depth, tht=tht, twf=twf,
                                        glr=glr, wc=wc, gas_grav=gas_grav, oil_grav=oil_grav,
                                        wtr_grav=wtr_grav, diameter=diameter, angle=angle, margin=margin,
                                        thp=thp, sample_size=sample_size)

    plt.ion()
    if rate:
        ax.plot(p, depths, label=f'Liquid Rate = {liquid_rate} STB\day \n')
    elif G:
        ax.plot(p, depths, label=f'GLR =  {glr} scf\stb \n')
    elif Wc:
        ax.plot(p, depths, label=f'WC =  {glr} ')

    else:
        ax.plot(p, depths, label=f'D = {diameter} inch \n'
                                 f'Liquid Rate = {liquid_rate} STB\day \n'
                                 f'GLR =  {glr} scf\stb \n'
                                 f'wc = {wc}\nGas gravity = {gas_grav} \n'
                                 f'Oil gravity=  {oil_grav} API\n'
                                 f'Water Gravity = {wtr_grav}\n\n')
    ax.set_xlabel('Pressure psi')
    ax.set_ylabel('Depth ft')
    plt.grid(which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.ylim(0, depths[-1])
    plt.legend()
    plt.gca().invert_yaxis()
    figure.suptitle("Traverse Curve", size=20)
    return figure

#------------------------------------- vlp --------------------------
def vlp(min_q, max_q, points, depth, tht, twf, GLR, wc, gas_grav, oil_grav, wtr_grav, diameter, angle, thp,
        sample_size=100, margin=30):
    # showing the graph
    rates = np.linspace(min_q, max_q, points)
    bhps = []
    for q in rates:
        p, depths = pressure_traverse(liquid_rate=q, depth=depth, tht=tht, twf=twf,
                                            glr=GLR, wc=wc, gas_grav=gas_grav, oil_grav=oil_grav,
                                            wtr_grav=wtr_grav, diameter=diameter, angle=angle, margin=margin,
                                            thp=thp, sample_size=sample_size)
        bhp = p[-1]
        bhps.append(bhp)
    # figure
    fig2, ax2 = plt.subplots(1, 1)
    plt.grid(which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.ylim(0, depths[-1])
    plt.plot(rates, bhps, '-')
    ax2.set_xlabel('Rate STB/day')
    ax2.set_ylabel('Pressure psi')
    fig2.suptitle("VLP Curve", size=18)
    return fig2