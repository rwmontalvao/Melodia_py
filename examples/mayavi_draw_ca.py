import glob

import seaborn as sns

from mayavi import mlab


def get_color(idt):
    pal = sns.color_palette('Set1', 7)
    colours = [(p[0], p[1], p[2]) for p in pal]
    if idt > 0:
        color = colours[idt - 1]
    else:
        color = (0.5, 0.5, 0.5)

    return color


mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0), size=(800 * 2, 600 * 2))

files = glob.glob("*a.csv")

for in_file in files:
    x = []
    y = []
    z = []
    g = []
    with open(in_file) as f:
        for line in f:
            tokens = list(map(float, line.split(",")))
            x.append(tokens[1])
            y.append(tokens[2])
            z.append(tokens[3])
            g.append(int(tokens[4]))

    ini = 0
    prev = g[ini]
    blocks = []
    for end in range(ini + 1, len(g)):
        curr = g[end]
        if curr != prev or end == len(g) - 1:
            blocks.append((ini, end, prev))
            prev = curr
            ini = end

    for i, block in enumerate(blocks):
        ini, end, c = block
        px = [v for v in x[ini:end]]
        py = [v for v in y[ini:end]]
        pz = [v for v in z[ini:end]]
        mlab.plot3d(px, py, pz, tube_radius=0.25, color=get_color(c))
        mlab.points3d(px, py, pz, scale_factor=1.0, color=get_color(c))

        if i + 1 < len(blocks):
            # end current block
            x1 = x[end - 1]
            y1 = y[end - 1]
            z1 = z[end - 1]
            c1 = c

            # init next block        
            ini, end, c = blocks[i + 1]
            x2 = x[ini]
            y2 = y[ini]
            z2 = z[ini]
            c2 = c

            xc = (x1 + x2) / 2.0
            yc = (y1 + y2) / 2.0
            zc = (z1 + z2) / 2.0

            mlab.plot3d([x1, xc], [y1, yc], [z1, zc], tube_radius=0.25, color=get_color(c1))
            mlab.plot3d([xc, x2], [yc, y2], [zc, z2], tube_radius=0.25, color=get_color(c2))

mlab.show()
