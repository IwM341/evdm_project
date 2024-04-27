import numpy as np
import matplotlib.pyplot as plt

def GFunc(x0,x1,F,N = 100):
    grid = np.linspace(x0,x1,N)
    values = F(grid)
    return (grid,values)
def GFunc1(x0,x1,F1,F2,N = 100):
    grid = np.linspace(x0,x1,N)
    values1 = F1(grid)
    values2 = F2(grid)
    return (grid,values1,values2)

Emin = -2
E1_pad = -1.1
E1 = -1
E2_pad = -0.7
E2 = -0.5

m_line_grid = np.linspace(-1,-0.5,100)
def m_line_func(x):
    return np.sqrt(2*(1+x))/(1.5+x)

m_line_grid1 = np.linspace(E1_pad,E2_pad,100)
def m_line_func1(x):
    return np.select([x<=E2_pad], [np.sqrt((-E1_pad+x)/(E2_pad-E1_pad))/(1 +1.2*(-E2_pad+x))],1)

plt.plot(m_line_grid,m_line_func(m_line_grid),color = 'black')
plt.plot(m_line_grid1,m_line_func1(m_line_grid1),linestyle='dashed')



plt.xlim((Emin,0))
plt.ylim((0,1))

E_intersec_blue = -0.99
plt.plot([Emin,0],[0.8,0.8],color = '#e219e5',linestyle='dashed')
plt.fill_between([Emin,E1_pad],[0.8,0.8],color = '#eedf9c') #Zone 1
plt.fill_between([Emin,E_intersec_blue],[0.8,0.8],[1 ,1],color = '#9ceeee') #Zone 2

E_b_grd = np.linspace(E1_pad,E_intersec_blue)
E_b_grd1 = np.linspace(E_intersec_blue,E2_pad)

plt.fill_between(E_b_grd1,m_line_func1(E_b_grd1),1+E_b_grd1*0,color = '#9ceeee') #Zone 2
plt.fill_between(E_b_grd,m_line_func1(E_b_grd),E_b_grd*0+0.8,color = '#eedf9c') #Zone 1

E_intersec_black = -0.88

plt.fill_between(*GFunc(E1_pad,E_intersec_blue,m_line_func1),color = '#a7e6c4')

m_grid1 = np.linspace(E_intersec_blue,E_intersec_black,100)
m_grid2 = np.linspace(E_intersec_black,-0.5,100)

L08_1 = m_grid1*0+0.8
black_1 = m_line_func(m_grid1)
blue_1 = m_line_func1(m_grid1)

black_2 = m_line_func(m_grid2)
blue_2 = m_line_func1(m_grid2)

plt.fill_between(m_grid1,black_1,L08_1,color = '#a7e6c4') #Zone 3

plt.fill_between(m_grid1,L08_1,blue_1,color = '#b28be0') #Zone 4
plt.fill_between(m_grid2,black_2 ,blue_2,color = '#b28be0') #Zone 4

plt.fill_between(*GFunc(-1,E_intersec_black,m_line_func),color = '#7c8976') #Zone 5
plt.fill_between([E_intersec_black,0],[0.8,0.8],color = '#7c8976') #Zone 5

plt.fill_between(*GFunc1(E_intersec_black,
                         -0.5,lambda x: x*0+0.8,
                         m_line_func),color = '#dce1bf') #Zone 6
plt.fill_between([-0.5,0],[0.8,0.8],[1,1],color = '#dce1bf') #Zone 6

plt.xlabel('E')
plt.ylabel('l')

#plt.plot([-1.1,-1.1],[0,1],linestyle='dashed')

plt.text(-1.75,0.4,r"$I_1(u_0, u_1)$") #Zone1
plt.text(-1.75,0.9,r"$I_2(u_0 \pm u_1)$") #Zone2

plt.text(-1.07,0.4,r"$I_3(\delta_1)$") #Zone3

plt.text(-0.98,0.9,r"$I_4(\delta_1,\delta_0)$") #Zone4

plt.text(-0.98,0.9,r"$I_4(\delta_1,\delta_0)$") #Zone4

plt.text(-0.6,0.4,r"$I_5(u_1 - exact,u_0)$") #Zone5

plt.text(-0.6,0.9,r"$I_6(u_1 - exact,\delta_0)$") #Zone6

plt.savefig('diagramm1.png')