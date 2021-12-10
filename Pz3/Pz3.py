'''
Моделирование распространения ЭМ волны, падающей на границу
вакуум - идеальный диэлектрик.
Используются граничные условия ABC первой степени.
'''
import math
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt
import numpy

import tools


class Sampler:
    def __init__(self, discrete: float):
        self.discrete = discrete

    def sample(self, x: float) -> int:
        return math.floor(x / self.discrete + 0.5)


class GaussianPlaneWave:
    ''' Класс с уравнением плоской волны для гауссова сигнала в дискретном виде
    d - определяет задержку сигнала.
    w - определяет ширину сигнала.
    Sc - число Куранта.
    eps - относительная диэлектрическая проницаемость среды, в которой расположен источник.
    mu - относительная магнитная проницаемость среды, в которой расположен источник.
    '''

    def __init__(self, d, w, Sc=1.0, eps=1.0, mu=1.0):
        self.d = d
        self.w = w
        self.Sc = Sc
        self.eps = eps
        self.mu = mu

    def getE(self, m, q):
        '''
        Расчет поля E в дискретной точке пространства m
        в дискретный момент времени q
        '''
        return numpy.exp(-(((q - m * numpy.sqrt(self.eps * self.mu) / self.Sc) - self.d) / self.w) ** 2)


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Скорость света в вакууме
    c = 299792458.0

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime_s = 4e-9

    # Размер области моделирования в метрах
    maxSize_m = 0.5


    
    # Дискрет по пространству в м
    dx = 1e-4

    v=c/numpy.sqrt(2)

    # Переход к дискретным отсчетам
    # Дискрет по времени
    dt = dx * Sc / v

    sampler_x = Sampler(dx)
    sampler_t = Sampler(dt)

    #Время расчета в отсчетах
    maxTime = sampler_t.sample(maxTime_s)

    #Размер области моделирования в отсчетах
    maxSize = sampler_x.sample(maxSize_m)


    # Положение источника в метрах
    sourcePos_m = 0.25
    sourcePos = math.floor(sourcePos_m / dx + 0.5) # Положение источника в отсчетах

  
    probesPos_m = 0.4
    # Датчики для регистрации поля
    probesPos = [math.floor( probesPos_m / dx + 0.5)]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]


    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[:] = 2.0

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)
    source = GaussianPlaneWave(300, 75, Sc, eps[sourcePos], mu[sourcePos])


    # Ez[-2] в предыдущий момент времени
    oldEzRight = Ez[-2]

    # Расчет коэффициентов для граничных условий


    tempRight = Sc / numpy.sqrt(mu[-1] * eps[-1])
    koeffABCRight = (tempRight - 1) / (tempRight + 1)

    # Параметры отображения поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -1.1
    display_ymax = 1.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    #display = tools.AnimateFieldDisplay(maxSize,
                                        #display_ymin, display_ymax,
                                       #display_ylabel)
    display = tools.AnimateFieldDisplay(dx, dt,
                                        maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel)

    display.activate()
    display.drawProbes(probesPos)
    display.drawSources([sourcePos])

    for q in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= Sc / (W0 * mu[sourcePos - 1]) * source.getE(0, q)

        # Расчет компоненты поля E
        Hy_shift = Hy[: -1]
        Ez[1:-1] = Ez[1: -1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1: -1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += (Sc / (numpy.sqrt(eps[sourcePos] * mu[sourcePos])) *
                          source.getE(-0.5, q + 0.5))

        # Граничные условия слева
        Ez[0] = 0
        # Граничные условия ABC первой степени
        Ez[-1] = oldEzRight + koeffABCRight * (Ez[-2] - Ez[-1])
        oldEzRight = Ez[-2]

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 100 == 0:
            display.updateData(display_field, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    #tools.showProbeSignals(probes, -1.1, 1.1)
EzSpec = fftshift(numpy.abs(fft(probe.E)))
df = 1.0 / (maxTime * dt)
freq = numpy.arange(-maxTime / 2 , maxTime / 2 , 1)*df
tlist = numpy.arange(0, maxTime * dt, dt)

# Вывод сигнала и спектра
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.set_xlim(8e-10, 1.6e-9)
ax1.set_ylim(-0.2, 1.1)
ax1.set_xlabel('t, с')
ax1.set_ylabel('Ez, В/м')
ax1.plot(tlist, probe.E/numpy.max(probe.E))
ax1.minorticks_on()
ax1.grid()
ax2.set_xlim(-1e11, 1e11)
ax2.set_ylim(0, 1.1)
ax2.set_xlabel('f, Гц')
ax2.set_ylabel('|S| / |Smax|, б/р')
ax2.plot(freq, EzSpec / numpy.max(EzSpec))
ax2.minorticks_on()
ax2.grid()
plt.show()
