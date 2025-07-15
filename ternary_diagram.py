import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

theta = np.pi/3
allowed_err = 0.1

##Data on Curie Temperatures (in deg C) for given TMX(0.4,0.7,1.0) and oxidation, Z. Digitised from Readman and OReilly 1972 
#Readman, P.W. & O’Reilly, W., 1972. Magnetic Properties of Oxidized (Cation-Deficient) Titanomagnetites (Fe, Ti, )
#3O4. J. geomagn. geoelec, 24, 69–90. doi:10.5636/jgg.24.69

data_RO1972 = {
    0.4: {'Z': [0.00, 0.22, 0.50, 0.62, 0.73, 0.81, 0.92, 0.98], 'Tc': [334, 352, 390, 402, 439, 472, 495, 519]},
    0.7: {'Z': [0.00, 0.39, 0.60, 0.72, 0.85, 0.97], 'Tc': [73, 125, 209, 253, 308, 381]},
    1.0: {'Z': [0.00, 0.33, 0.39, 0.46, 0.63, 0.75, 0.87], 'Tc': [-154, -120, -100, -52, 32, 81, 162]}
}

# Replot figure from Readman and O'Reilly
fig1, ax1 = plt.subplots(figsize=(8, 6))
for TMX, data in data_RO1972.items():
    # Plot original data
    ax1.scatter(data['Z'], data['Tc'], label=f'x={TMX}')

    # Fit a polynomial of degree 3 (cubic polynomial) to the data
    coeffs = np.polyfit(data['Z'], data['Tc'], 3)
    poly = np.poly1d(coeffs)
    
    # Generate a smoother curve using the polynomial
    Z_smooth = np.linspace(0, max(data['Z']), 100)
    Tc_smooth = poly(Z_smooth)
    
    plt.plot(Z_smooth, Tc_smooth, label=f'x={TMX}')

    # Print polynomial coefficients
    print(f"x={TMX}: Polynomial Coefficients: {coeffs}")

ax1.set_xlabel('Z')
ax1.set_ylabel('Tc')
ax1.set_title('Curie Temperature vs. Oxidation Value')
ax1.legend()
ax1.grid(True)
fig1.tight_layout()
plt.show()




def ternary_to_binary(a,b,c):
    if a+b+c > 1-allowed_err and a+b+c < 1+allowed_err :
        print("Components sum to 1. All good.")
    else:
        print("Components do not sum to 1. Check components.")
        return 0,0
    
    y = a*np.sin(theta)
    x = c + y/np.tan(theta)
    return x, y

def binary_to_ternary(x,y):
    a=y/np.sin(theta)
    c= -y/np.tan(theta) + x
    b = 1 -a -c
    return a,b,c


def ternary_plot(fig,ax):
    # Ensure the sum of the three components is 1
    #assert np.allclose(np.sum(data, axis=1), 1.0), "Sum of components must be 1"

    # Plot the triangle outline
    ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3)/2, 0], 'k-')

    # Plot grid lines parallel to the sides of the triangle
    theta = np.pi/3
    for i in np.linspace(0.1, 1, 10):
        print(i)
        ax.plot([(1-i), (1-i)*np.cos(theta)], [0, (1-i) *np.sin(theta)], 'g--', alpha=0.2)
        ax.plot([i, 1-(1-i)*np.cos(theta)], [0, (1-i) *np.sin(theta)], 'r--', alpha=0.2)
        ax.plot([i*np.cos(theta), 1-i*np.cos(theta)], [i*np.sin(theta), i*np.sin(theta)], 'b--', alpha=0.2)

        # Annotate the ends of the lines with i
        ax.annotate(f'{i:.1f}', xy=((1-i)*np.cos(theta), (1-i) *np.sin(theta)), xytext=(-10, 5), textcoords='offset points', ha='center', va='center', color='g')
        ax.annotate(f'{i:.1f}', xy=(i, 0), xytext=(-10, -10), textcoords='offset points', ha='center', va='center', color='r')
        ax.annotate(f'{i:.1f}', xy=(1-i*np.cos(theta), i*np.sin(theta)), xytext=(10, 0), textcoords='offset points', ha='center', va='center', color='b')
    
    ax.annotate('A', xy=(0.5, np.sin(theta)), xytext=(0, 20), textcoords='offset points', ha='center', va='center', color='b', fontsize=16)
    ax.annotate('$Ti^{4+}$', xy=(0.5, np.sin(theta)), xytext=(0, 40), textcoords='offset points', ha='center', va='center', color='b', fontsize=16)

    ax.annotate('B', xy=(0, 0), xytext=(-20, -20), textcoords='offset points', ha='center', va='center', color='g', fontsize=16)
    ax.annotate('$Fe^{2+}$', xy=(0, 0), xytext=(-40, -40), textcoords='offset points', ha='center', va='center', color='g', fontsize=16)

    ax.annotate('C', xy=(1, 0), xytext=(20, -20), textcoords='offset points', ha='center', va='center', color='r', fontsize=16)
    ax.annotate('$Fe^{3+}$', xy=(1, 0), xytext=(40, -40), textcoords='offset points', ha='center', va='center', color='r', fontsize=16)


    for TMX in [20,40,60,80]:
        a,b,c = titanomagnetite(TMX)
        x,y = ternary_to_binary(a,b,c)
        print(f'TM{TMX} a:{a:.2f}, b:{b:.2f}, c:{c:.2f}')
        ax.plot(x,y,'k|', label=f'TM{TMX} a:{a:.2f}, b:{b:.2f}, c:{c:.2f}')
        ax.annotate(f'TM{TMX}', xy=(x, y), xytext=(-25, 0), textcoords='offset points', ha='center', va='center', color='k')

    #Magnetite
    TMX=0
    a,b,c = titanomagnetite(TMX)
    xTM0,yTM0 = ternary_to_binary(a,b,c)
    ax.plot(xTM0,yTM0,'k*')
    ax.annotate(f'Magnetite', xy=(xTM0, yTM0), xytext=(0, -20), textcoords='offset points', ha='center', va='center', color='k')
    ax.annotate('$1/3Fe_3O_4$', xy=(xTM0, yTM0), xytext=(0, -35), textcoords='offset points', ha='center', va='center', color='k')

    #Ulvospinel
    TMX=100
    a,b,c = titanomagnetite(TMX)
    xTM100,yTM100 = ternary_to_binary(a,b,c)
    ax.plot(xTM100,yTM100,'k*')
    ax.annotate(f'Ulvospinel', xy=(xTM100, yTM100), xytext=(-30, 10), textcoords='offset points', ha='center', va='center', color='k')
    ax.annotate(r'$\frac{1}{3} Fe_2TiO_4$', xy=(xTM100, yTM100), xytext=(-35, 25), textcoords='offset points', ha='center', va='center', color='k')
    
    ax.plot([xTM100, xTM0], [yTM100, yTM0], 'k-', alpha=0.8)

    ##Titanium oxides
    a,b,c = 1,0,0
    x,y = ternary_to_binary(a,b,c)
    ax.annotate(f'Rutile, Anatase, Akaogiite, Brookite', xy=(x,y), xytext=(30, 10), textcoords='offset points', ha='left', va='center', color='k')
    ax.annotate('$TiO_2$', xy=(x,y), xytext=(35, 25), textcoords='offset points', ha='center', va='center', color='k')

    ##Hematite/Maghemite
    a,b,c = 0,0,1
    xTH0,yTH0 = ternary_to_binary(a,b,c)
    ax.annotate(f'Hematite, Maghemite', xy=(xTH0,yTH0), xytext=(30, 10), textcoords='offset points', ha='left', va='center', color='k')
    ax.annotate(r'$\frac{1}{2} \alpha Fe_2O_3$, $\frac{1}{2} \gamma Fe_2O_3$', xy=(xTH0,yTH0), xytext=(35, 25), textcoords='offset points', ha='left', va='center', color='k')

    ##Ilmenite
    a,b,c = 0.5,0.5,0
    xTH100,yTH100 = ternary_to_binary(a,b,c)
    ax.annotate(f'Ilmenite', xy=(xTH100,yTH100), xytext=(-30, 10), textcoords='offset points', ha='right', va='center', color='k')
    ax.annotate(r'$\frac{1}{2} FeTiO_3$', xy=(xTH100,yTH100), xytext=(-35, 25), textcoords='offset points', ha='right', va='center', color='k')

    ##Wustite
    a,b,c = 0,1,0
    x,y = ternary_to_binary(a,b,c)
    ax.annotate(f'Wüstite', xy=(x,y), xytext=(-30, 10), textcoords='offset points', ha='right', va='center', color='k')
    ax.annotate(r'$FeO$', xy=(x,y), xytext=(-35, 25), textcoords='offset points', ha='right', va='center', color='k')

    ##Titanohematite / Hemoilmenite
    for THY in [20,40,60,80]:
        a,b,c = titanohematite(THY)
        x,y = ternary_to_binary(a,b,c)
        print(f'TH{THY} a:{a:.2f}, b:{b:.2f}, c:{c:.2f}')
        ax.plot(x,y,'k|', label=f'TH{THY} a:{a:.2f}, b:{b:.2f}, c:{c:.2f}')
        #ax.annotate(f'TH{THY}', xy=(x, y), xytext=(20, 0), textcoords='offset points', ha='center', va='center', color='k')
    ax.plot([xTH100, xTH0], [yTH100, yTH0], 'k-', alpha=0.8)

    for z in [0.2,0.4,0.6,0.8,1]:
        TMX = 60 ## Use TM60 as illustrative
        a,b,c = oxidise(TMX, z)
        x,y = ternary_to_binary(a,b,c)
        ax.plot(x,y,'k.')
        ax.annotate(f'z={z}', xy=(x, y), xytext=(0, -10), textcoords='offset points', ha='center', va='center', color='k')
    a,b,c = oxidise(60, 0)
    xTM60, yTM60 = ternary_to_binary(a,b,c)
    a,b,c = oxidise(60,1)
    xTM60ox, yTM60ox = ternary_to_binary(a,b,c)
    ax.plot([xTM60, xTM60ox], [yTM60, yTM60ox], 'k',linestyle ='dotted', alpha=0.8)


    # Define a colormap and normalization
    colormap = plt.get_cmap('hot')
    normalize = plt.Normalize(vmin=-271, vmax=1000)

    # Loop through data_RO1972
    xList=[]
    yList=[]
    TcList=[]
    for TMX, data in data_RO1972.items():
        TMX = TMX * 100  # Convert to percentage
        Z_values = data['Z']
        Tc_values = data['Tc']
        for Z, Tc in zip(Z_values, Tc_values):
            a, b, c = oxidise(TMX, Z)
            x, y = ternary_to_binary(a, b, c)
            xList.append(x)
            yList.append(y)
            TcList.append(Tc)
            color_value = Tc  # Using Tc as the color value
            color = colormap(normalize(color_value))  # Map value to color using colormap and normalization
            ax.plot(x, y, marker='o', linestyle='', markersize=5, alpha=1, 
                    color=color, markeredgecolor='black', label=f'TM{TMX} Z:{Z:.2f}, Tc:{Tc:.2f}')
            
    ##include the point for pure magnetite (Tc =580c from Tauxe)
    Z=0 
    TMX=0
    Tc=580
    a, b, c = oxidise(TMX, Z)
    x, y = ternary_to_binary(a, b, c)
    xList.append(x)
    yList.append(y)
    TcList.append(Tc)
    color_value = Tc  # Using Tc as the color value
    color = colormap(normalize(color_value))  # Map value to color using colormap and normalization
    ax.plot(x, y, marker='o', linestyle='', markersize=5, alpha=1, 
            color=color, markeredgecolor='black', label=f'TM{TMX} Z:{Z:.2f}, Tc:{Tc:.2f}')
    

    ##include the point for pure hematite (Tc =680c from Petersen and Bleil, 1982)
    Z=1 
    TMX=0
    Tc=680
    a, b, c = oxidise(TMX, Z)
    x, y = ternary_to_binary(a, b, c)
    xList.append(x)
    yList.append(y)
    TcList.append(Tc)
    color_value = Tc  # Using Tc as the color value
    color = colormap(normalize(color_value))  # Map value to color using colormap and normalization
    ax.plot(x, y, marker='o', linestyle='', markersize=5, alpha=1, 
            color=color, markeredgecolor='black', label=f'TM{TMX} Z:{Z:.2f}, Tc:{Tc:.2f}')
    
    
    # Convert lists to 1-D arrays
    xList = np.array(xList)
    yList = np.array(yList)
    TcList = np.array(TcList)

    # Add color bar indicating the value of Curie Temperature
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=normalize)
    sm.set_array([])  # Dummy array for the color bar
    cbar = fig.colorbar(sm, ax=ax, label='Curie Temperature (°C)')

    # Create a grid of x and y values
    X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    # Interpolate z values on the grid
    Z = griddata((xList, yList), TcList, (X, Y), method='cubic')
    # Contour plot
    contour = ax.contour(X, Y, Z, levels =[-100,0,100,200,300,400,500,600], cmap=colormap, norm=normalize)
    # Add labels to the contour lines
    def fmt_func(value):
        return f"{value:.0f}°C"
    ax.clabel(contour, inline=True, fontsize=12, fmt=fmt_func)


    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)

    # Remove axes and bounding box
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #ax.legend() ##comment this in to add all labels


##Titanomagnetites
TMX=60 #The titanium composition in the titanomagnetite solid solution as a %. 100 for ulvospinel, 0 for magnetite
def titanomagnetite(TMX):
    TMX = TMX/100 #convert from percentage to fraction
    c=2/3*(1-TMX)
    a=TMX/3
    b=1-a-c
    return a,b,c 

def titanohematite(THY):
    THY = THY/100 #convert from percentage to fraction
    c =1-THY
    a = THY/2
    b=1-a-c
    return a,b,c

def oxidise (TMX, z):
    """Calculates the coordinates in a, b, c space for a given titanomagnetite and degree of oxidation"""
    ##a (Ti content) is kept fixed 
    ##x_oxy is the x coordinate of the oxidised mineral
    a_int, b_int, c_int = titanomagnetite(TMX)
    x_int, y_int = ternary_to_binary(a_int,b_int,c_int)
    x_max, x_min = limits(y_int)
    x_oxy = z*(x_max-x_int) + x_int 
    y_oxy = y_int
    a_oxy, b_oxy, c_oxy = binary_to_ternary(x_oxy, y_oxy)
    return(a_oxy, b_oxy, c_oxy)

def limits(y):
    """Returns the maximum and minimum allowed value of x for a given y"""
    h = 3**(1/2)/2#height of the triangle
    ##when y=h, x_max=0.5
    ##when y=0, x_max=1
    x_max = 1 - 0.5*y/h
    x_min = 0.5*(1-x_max)
    return x_max, x_min





# Plotting
fig2, ax2 = plt.subplots(figsize=(8, 8))
ternary_plot(fig2, ax2)
ax2.set_title('Ternary Phase Diagram')
ax2.grid(False)  # Turn off the default grid
fig2.savefig('Ternary_Phase_Diagram.svg', format='svg')
plt.show()
