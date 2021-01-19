import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.axes_grid1
import matplotlib.widgets

class Player(FuncAnimation):
    def __init__(self, ims, frames=None, init_func=None, fargs=None,
                 save_count=None, mini=0, maxi=100, pos=(0.125, 0.92),
                 interval=25, **kwargs):


        self.fig, self.ax = plt.subplots()

        self.point = self.ax.imshow(ims[:, :, 0])

        self.i = 0
        self.min = 0
        self.max = ims.shape[-1]
        self.runs = True
        self.forwards = True

        self.func = self.update
        self.setup(pos)
        FuncAnimation.__init__(self, self.fig, self.update,
                               frames=self.play(), init_func=init_func,
                               fargs=fargs, save_count=save_count,
                               interval=interval,
                               **kwargs )

    def update(self, i):
        self.time_slider.eventson = False
        self.time_slider.set_val(i)
        self.time_slider.eventson = True

        self.point.set_data(ims[:, :, i])
    def update_playback_speed(self, i):
        FuncAnimation.__init__(self, self.fig, self.update,
                               frames=self.play(), init_func=init_func,
                               fargs=fargs, save_count=save_count,
                               interval=interval,
                               **kwargs )

    def play(self):
        while self.runs:
            self.i = self.i+self.forwards-(not self.forwards)
            if self.i > self.min and self.i < self.max:
                yield self.i
            else:
                self.stop()
                yield self.i

    def start(self):
        self.runs=True
        self.event_source.start()

    def stop(self, event=None):
        self.runs = False
        self.event_source.stop()

    def forward(self, event=None):
        self.forwards = True
        self.start()
    def backward(self, event=None):
        self.forwards = False
        self.start()
    def oneforward(self, event=None):
        self.forwards = True
        self.onestep()
    def onebackward(self, event=None):
        self.forwards = False
        self.onestep()

    def onestep(self):
        if self.i > self.min and self.i < self.max:
            self.i = self.i+self.forwards-(not self.forwards)
        elif self.i == self.min and self.forwards:
            self.i+=1
        elif self.i == self.max and not self.forwards:
            self.i-=1
        self.func(self.i)
        self.fig.canvas.draw_idle()

    def set_time(self, i):
        self.i = int(i)
        self.stop()
        self.func(self.i)
        self.fig.canvas.draw_idle()

    def setup(self, pos):
        playerax = self.fig.add_axes([pos[0],pos[1], 0.22, 0.04])
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(playerax)
        bax = divider.append_axes("right", size="80%", pad=0.05)
        sax = divider.append_axes("right", size="80%", pad=0.05)
        fax = divider.append_axes("right", size="80%", pad=0.05)
        ofax = divider.append_axes("right", size="100%", pad=0.05)
        self.button_oneback = matplotlib.widgets.Button(playerax, label=u'$\u29CF$')
        self.button_back = matplotlib.widgets.Button(bax, label=u'$\u25C0$')
        self.button_stop = matplotlib.widgets.Button(sax, label=u'$\u25A0$')
        self.button_forward = matplotlib.widgets.Button(fax, label=u'$\u25B6$')
        self.button_oneforward = matplotlib.widgets.Button(ofax, label=u'$\u29D0$')
        self.button_oneback.on_clicked(self.onebackward)
        self.button_back.on_clicked(self.backward)
        self.button_stop.on_clicked(self.stop)
        self.button_forward.on_clicked(self.forward)
        self.button_oneforward.on_clicked(self.oneforward)
        slider_ax = self.fig.add_axes((0.1, 0.0, 0.85, 0.04))
        self.time_slider = matplotlib.widgets.Slider(slider_ax,
                                                     label='Frame',
                                                     valmin=0,
                                                     valmax=self.max,
                                                     valinit=0.0)
        self.time_slider.on_changed(self.set_time)



### using this class is as easy as using FuncAnimation:



ims = np.random.random((300, 300, 30))
ani = Player(ims)

plt.show()
