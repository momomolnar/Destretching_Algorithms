"""
Flicker routine for visual inspection of the differences between
two images. (following flicker in IDL)

Author: Momchil Molnar (momo@nso.edu)
Date: 31/12/2020
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def replot():
    """

    :return:
    """


def flicker(im1, im2, wait_time=.25, repetitions=100):
    """
    The main flicker function.

    : im1: ndarray [nx, ny]
        Image 1 to be plotted
    : im2: ndarray [nx, ny]
         Image 2 to be plotted
    : wait_time: int
        delay between replots [seconds]

    :return: success code (0==success!)
    """

    fig = plt.figure()
    sz_im = im1.shape[0]

    ims = []
    im = plt.imshow(im1, animated=True)
    ims.append([im])
    im = plt.imshow(im2, animated=True)
    ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=wait_time*1e3, blit=True,
                                    repeat_delay=wait_time*1e3)

    plt.show()

    return 0


def test_flicker():
    im1 = np.random.random((1000, 1000))
    im2 = np.random.random((1000, 1000))*.0
    flicker(im1, im2)
    # assert (flicker(im1, im2) == 0)


if __name__ == '__main__':
    test_flicker()
