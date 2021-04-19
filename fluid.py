"""
Based on the Jos Stam paper https://www.researchgate.net/publication/2560062_Real-Time_Fluid_Dynamics_for_Games
and the mike ash vulgarization https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

https://github.com/Guilouf/python_realtime_fluidsim
"""
import numpy as np
import math


class Fluid:

    def __init__(self):
        self.rotx = 1
        self.roty = 1
        self.cntx = 1
        self.cnty = -1

        self.size = 60  # map size
        self.dt = 0.2  # time interval
        self.iter = 2  # linear equation solving iteration number

        self.diff = 0.0000  # Diffusion
        self.visc = 0.0000  # viscosity

        self.s = np.full((self.size, self.size), 0,
                         dtype=float)        # Previous density
        self.density = np.full((self.size, self.size), 0,
                               dtype=float)  # Current density

        # array of 2d vectors, [x, y]
        self.velo = np.full((self.size, self.size, 2), 0, dtype=float)
        self.velo0 = np.full((self.size, self.size, 2), 0, dtype=float)

    def step(self):
        self.diffuse(self.velo0, self.velo, self.visc)

        # x0, y0, x, y
        self.project(self.velo0[:, :, 0], self.velo0[:, :, 1],
                     self.velo[:, :, 0], self.velo[:, :, 1])

        self.advect(self.velo[:, :, 0], self.velo0[:, :, 0], self.velo0)
        self.advect(self.velo[:, :, 1], self.velo0[:, :, 1], self.velo0)

        self.project(self.velo[:, :, 0], self.velo[:, :, 1],
                     self.velo0[:, :, 0], self.velo0[:, :, 1])

        self.diffuse(self.s, self.density, self.diff)

        self.advect(self.density, self.s, self.velo)

    def lin_solve(self, x, x0, a, c):
        """Implementation of the Gauss-Seidel relaxation"""
        c_recip = 1 / c

        for iteration in range(0, self.iter):
            # Calculates the interactions with the 4 closest neighbors
            x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (x[2:, 1:-1] +
                                                   x[:-2, 1:-1] + x[1:-1, 2:] + x[1:-1, :-2])) * c_recip

            self.set_boundaries(x)

    def set_boundaries(self, table):
        """
        Boundaries handling
        :return:
        """

        if len(table.shape) > 2:  # 3d velocity vector array
            # Simulating the bouncing effect of the velocity array
            # vertical, invert if y vector
            table[:, 0, 1] = - table[:, 0, 1]
            table[:, self.size - 1, 1] = - table[:, self.size - 1, 1]

            # horizontal, invert if x vector
            table[0, :, 0] = - table[0, :, 0]
            table[self.size - 1, :, 0] = - table[self.size - 1, :, 0]

        table[0, 0] = 0.5 * (table[1, 0] + table[0, 1])
        table[0, self.size - 1] = 0.5 * \
            (table[1, self.size - 1] + table[0, self.size - 2])
        table[self.size - 1, 0] = 0.5 * \
            (table[self.size - 2, 0] + table[self.size - 1, 1])
        table[self.size - 1, self.size - 1] = 0.5 * table[self.size - 2, self.size - 1] + \
            table[self.size - 1, self.size - 2]

    def diffuse(self, x, x0, diff):
        if diff != 0:
            a = self.dt * diff * (self.size - 2) * (self.size - 2)
            self.lin_solve(x, x0, a, 1 + 6 * a)
        else:  # equivalent to lin_solve with a = 0
            x[:, :] = x0[:, :]

    def project(self, velo_x, velo_y, p, div):
        # numpy equivalent to this in a for loop:
        # div[i, j] = -0.5 * (velo_x[i + 1, j] - velo_x[i - 1, j] + velo_y[i, j + 1] - velo_y[i, j - 1]) / self.size
        div[1:-1, 1:-1] = -0.5 * (
            velo_x[2:, 1:-1] - velo_x[:-2, 1:-1] +
            velo_y[1:-1, 2:] - velo_y[1:-1, :-2]) / self.size
        p[:, :] = 0

        self.set_boundaries(div)
        self.set_boundaries(p)
        self.lin_solve(p, div, 1, 6)

        velo_x[1:-1, 1:-1] -= 0.5 * (p[2:, 1:-1] - p[:-2, 1:-1]) * self.size
        velo_y[1:-1, 1:-1] -= 0.5 * (p[1:-1, 2:] - p[1:-1, :-2]) * self.size

        self.set_boundaries(self.velo)

    def advect(self, d, d0, velocity):
        dtx = self.dt * (self.size - 2)
        dty = self.dt * (self.size - 2)

        for j in range(1, self.size - 1):
            for i in range(1, self.size - 1):
                tmp1 = dtx * velocity[i, j, 0]
                tmp2 = dty * velocity[i, j, 1]
                x = i - tmp1
                y = j - tmp2

                if x < 0.5:
                    x = 0.5
                if x > (self.size - 1) - 0.5:
                    x = (self.size - 1) - 0.5
                i0 = math.floor(x)
                i1 = i0 + 1.0

                if y < 0.5:
                    y = 0.5
                if y > (self.size - 1) - 0.5:
                    y = (self.size - 1) - 0.5
                j0 = math.floor(y)
                j1 = j0 + 1.0

                s1 = x - i0
                s0 = 1.0 - s1
                t1 = y - j0
                t0 = 1.0 - t1

                i0i = int(i0)
                i1i = int(i1)
                j0i = int(j0)
                j1i = int(j1)

                try:
                    d[i, j] = s0 * (t0 * d0[i0i, j0i] + t1 * d0[i0i, j1i]) + \
                        s1 * (t0 * d0[i1i, j0i] + t1 * d0[i1i, j1i])
                except IndexError:
                    # tmp = str("inline: i0: %d, j0: %d, i1: %d, j1: %d" % (i0, j0, i1, j1))
                    # print("tmp: %s\ntmp1: %s" %(tmp, tmp1))
                    raise IndexError
        self.set_boundaries(d)

    def turn(self):
        self.cntx += 1
        self.cnty += 1
        if self.cntx == 3:
            self.cntx = -1
            self.rotx = 0
        elif self.cntx == 0:
            self.rotx = self.roty * -1
        if self.cnty == 3:
            self.cnty = -1
            self.roty = 0
        elif self.cnty == 0:
            self.roty = self.rotx
        return self.rotx, self.roty


if __name__ == "__main__":
    try:
        import matplotlib.pyplot as plt
        from matplotlib import animation
        from matplotlib.colors import LinearSegmentedColormap
        import re

        # We get velocity data from line and return it to update_im function
        def getVelocityDataFromLine(line):
            temp = line[8:].split("|")
            temp[0] = re.sub(r'[()]', '', temp[0])
            x, y = temp[0].split(",")
            velocity = temp[1].split("=")
            velocity[1] = re.sub(r'[()]', '', velocity[1])
            a, b = velocity[1].split(",")
            return int(x), int(y), int(a), int(b)

        # We get density data from line and return it to update_im function
        def getDensityDataFromLine(line):
            temp = line[15:].split("|")
            temp[0] = re.sub(r'[()]', '', temp[0])
            x, y = temp[0].split(",")
            x1, x2 = x.split(":")
            y1, y2 = y.split(":")
            density = temp[1].split("=")
            return int(x1), int(x2), int(y1), int(y2), int(density[1])

        # We get figure data from line and return it to update_im function
        def getFigureData(line):
            temp = re.sub(r'[()]', '', line[15:])
            x, y = temp.split(",")
            x1, x2 = x.split(":")
            y1, y2 = y.split(":")
            return int(x1), int(x2), int(y1), int(y2)

        inst = Fluid()

        # We open input.txt file
        inputFile = open('input.txt', 'r')
        # We save all lines of input.txt file
        lines = inputFile.readlines()

        def update_im(i):
            # We read input.txt line by line
            for line in lines:
                # If first velocity animation (1V) is detected in line
                if(line[0:2] == "1V"):
                    x, y, a, b = getVelocityDataFromLine(line)
                    # First velocity animation equation
                    inst.velo[x, y] = [a, b]

                # If second velocity animation (2V) is detected in line
                elif(line[0:2] == "2V"):
                    x, y, a, b = getVelocityDataFromLine(line)
                    # Second velocity animation equation
                    inst.velo[x, y] = [4*math.cos(a), 3*math.sin(3 * b)]

                # If third velocity animation (3V) is detected in line
                elif(line[0:2] == "3V"):
                    x, y, a, b = getVelocityDataFromLine(line)
                    # Third velocity animation equation
                    inst.velo[x, y] = [2*math.cos(a), 2*math.sin(b)]

                # If density (D) is detected in line
                elif(line[0] == "D"):
                    x1, x2, y1, y2, density = getDensityDataFromLine(line)
                    # Density equation
                    inst.density[x1:x2, y1:y2] += density

                # If figure (F) is detected in line
                elif(line[0] == "F"):
                    x1, x2, y1, y2 = getFigureData(line)
                    for x in range(x1, x2):
                        for y in range(y1, y2):
                            inst.density[x, y] = 0
                            inst.velo[x, y] = 0

            inst.step()
            im.set_array(inst.density)
            q.set_UVC(inst.velo[:, :, 1], inst.velo[:, :, 0])
            print(f"Density sum: {inst.density.sum()}")
            im.autoscale()

        fig = plt.figure()

        # We create new colormap
        cmapRGB = {'red':   [(0.0,  0.0, 0.0),
                             (0.5,  1.0, 1.0),
                             (1.0,  1.0, 1.0)],

                   'green': [(0.0,  0.0, 0.0),
                             (0.25, 0.0, 0.0),
                             (0.75, 1.0, 1.0),
                             (1.0,  1.0, 1.0)],

                   'blue':  [(0.0,  0.0, 0.0),
                             (0.5,  0.0, 0.0),
                             (1.0,  1.0, 1.0)]}

        # We register new colormap
        plt.register_cmap(cmap=LinearSegmentedColormap(
            'cmap', cmapRGB))

        # plot density
        im = plt.imshow(inst.density, cmap="cmap",
                        vmax=100, interpolation='bilinear')

        # plot vector field
        q = plt.quiver(inst.velo[:, :, 1],
                       inst.velo[:, :, 0], scale=10, angles='xy')
        anim = animation.FuncAnimation(
            fig, update_im, interval=1, save_count=1000)
        #anim.save("movie.mp4", fps=30, extra_args=['-vcodec', 'libx264'])
        plt.show()

    except ImportError:
        import imageio

        frames = 30

        flu = Fluid()

        video = np.full((frames, flu.size, flu.size), 0, dtype=float)

        for step in range(0, frames):
            flu.density[4:7, 4:7] += 100  # add density into a 3*3 square
            flu.velo[5, 5] += [1, 2]

            flu.step()
            video[step] = flu.density

        imageio.mimsave('./video.gif', video.astype('uint8'))
