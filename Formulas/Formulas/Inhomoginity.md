# Учет неоднородностей

## 1

### Время полета от $r_{min}$ до $r_{max}$

$$
T = \int_{r_{min}}^{r_{max}} {\cfrac{dr}{\sqrt{\phi(r) - e - l^2/r^2}}} = |u = r^2| =
$$
$$
\frac{1}{2}\int_{r_{min}^2}^{r_{max}^2} {\cfrac{du}{\sqrt{u(\phi(u) - e) - l^2}}} =
| u = \cfrac{u_{max}+u_{min}}{2} - \cfrac{u_{max}-u_{min}}{2} \cos{\psi}|
$$
$$
\frac{1}{2}\int_{0}^{\pi} {\cfrac{d\psi}{\sqrt{S(u(\psi))}}} =
$$
где
$$
S(u(\psi)) = \cfrac{u(\phi(u) - e) - l^2}{(u-u_{min})(u_{max}-u)}
$$

### Вращение траектории

* Изменение угла равно

$$
\varphi_{T} = \int {\frac{l}{r^2} dt} =  \frac{1}{2}\int_{0}^{\pi} {
    \cfrac{1}{u(\psi)} \cfrac{d\psi}{\sqrt{S(u(\psi))}}} 
$$

* Координаты рассматриваются так, что:

$$
    x = r \sin{\theta} \cos{\varphi}
$$
$$
    y = r \sin{\varphi}
$$
$$
    z = r \cos{\theta} \cos{\varphi}
$$

* Тогда изменение угла на траектории равно:

$$
\delta \varphi  =  \frac{1}{2}\int_{0}^{\pi} {
    \cfrac{d\psi}{u(\psi)} \left(\cfrac{1}{\sqrt{S(u(\psi))}} -\cfrac{1}{\sqrt{\rho_0}}   \right)}
$$

считаем, что $\rho_0 \cdot u_{max}\cdot u_{min} = l^2$

Если $S(u(\psi)) = const$, то $\varphi_{T} = \pi/2$, а если $\phi(r) ~ 1/r$ то  $\varphi_{T} = \pi$. В остальных случаях траектория будет апериодической, поэтому **влиянием неоднородностей на $\delta \varphi$ можно пренебречь**

### Изменение углового момента

Пусть $V$ – малая поправка к $\phi(r)$. Тогда

$$
    \frac{d}{dt} \vec{L} = [\vec{r},\dot{\vec{v}}] = -[\vec{r},\nabla V]
$$

В случае, если $V$ осесимметричный

$$
    \nabla V = \vec{n}_r \cfrac{\partial V}{\partial r} + \vec{n}_{\theta} \cfrac{1}{r_K} \cfrac{\partial V}{\partial {\theta}_K}
$$
$$
    \frac{d}{dt} \vec{L} = -\cfrac{\partial V}{\partial \theta_K} \vec{n}_{\varphi}
$$

Далее, единичный вектор в напревлении $\vec{L}$ равен

$$
    \vec{n}_{L} = (-\sin{\theta},0,-\cos{\theta})
$$

а единичный вектор $\vec{n}_{\varphi}$ равен

$$
    \vec{n}_{\varphi} = \cfrac{(x,y,0)}{\sqrt{x^2+y^2}} = \cfrac{(\sin{\theta} \cos{\varphi},\sin{\varphi},0)}{\sqrt{\sin^2{\theta} \cos^2{\varphi}+\sin^2{\varphi}}}
$$

* Поперечные изменения $\vec{L}$ неважны, так как приводят лишь к вращению вокруг $z$
плоскости траектории (что несущественно). Поэтому важно лишь изменение вдоль $x$

$$
    \delta L^{x}_{T} = \frac{1}{2}\int_{0}^{\pi} {\cfrac{d\psi}{\sqrt{S(u(\psi))}}}
    \cfrac{\partial V}{\partial \theta_K}(r[u(\psi)], \theta)
    \cfrac{\sin{\theta} \cos{\varphi}}{\sqrt{\sin^2{\theta} \cos^2{\varphi}+\sin^2{\varphi}}}
$$

В итоге, у нас имеется следующая параметризация:

$$
    (L_z,L_x,\varphi_T)
$$

$$
    \cfrac{d}{dt}L_x = X(L_z,L_x,\varphi_T)
$$

### Далее

Скорее всего $\varphi_T$ будет менятся достаточно быстро, поэтому мы можем
усреднить по этому параметру производную $\dot{L}_x$. Как итог получим некое движение (медленное) с уравнением:

$$
    \cfrac{d L_x}{dt} = \langle X(L_z,L_x,\varphi_T)\rangle _{\varphi_T} = X(L_z,L_x)
$$

**И тут тревога: уравнение первого порядка имеет только один аттрактор – точку, а значит – что-то не так. А не так что-то при $\theta = \pi/2$.**

В любом случае, можно считать, что одна траектория на самом деле рапределена по $L_x$ с
плотостью обратно пропорциональной скорости

$$
    \rho(L_x) \sim  \cfrac{1}{|dL_x/dt|} =  \cfrac{1}{|X(L_z,L_x)|}
$$
