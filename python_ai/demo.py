import random
import colorsys
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from voronoi import generateL1Voronoi


def random_normal(sharpness):
    return sum(random.random() for _ in range(sharpness)) / sharpness


def get_cell_color(index, total):
    hue = (index * 0.618033988749895) % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.5, 0.95)
    return (r, g, b)


def get_color(size):
    return {
        4: '#4286f4',
        8: '#44f453',
        16: '#931d78',
        32: '#ff3c35',
        64: '#f4ad42',
        128: '#009182',
        256: '#993300',
        512: '#669999',
        1024: '#800000',
        2048: '#333300',
    }.get(size, '#000000')


def main():
    width = 400
    height = 400
    num_points = 89

    random.seed(1)
    raw = [[int(random_normal(2) * width), int(random_normal(2) * height)] for _ in range(num_points)]
    sites = [list(p) for p in raw]

    vector_points = generateL1Voronoi(sites, width, height)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    drawn_bisectors = set()
    for i, site in enumerate(vector_points):
        poly_pts = [(x, height - y) for x, y in site['polygonPoints']]
        poly = Polygon(poly_pts, closed=True, edgecolor='black',
                       facecolor=get_cell_color(i, len(vector_points)), linewidth=0.5, alpha=0.8)
        ax1.add_patch(poly)
        ax1.plot(site['site'][0], height - site['site'][1], 'ko', markersize=3)

        for bisector in site['bisectors']:
            bid = id(bisector)
            if bid in drawn_bisectors:
                continue
            drawn_bisectors.add(bid)
            pts = [(x, height - y) for x, y in bisector['points']]
            color = get_color(bisector.get('mergeLine', 0))
            ax2.plot([p[0] for p in pts], [p[1] for p in pts],
                     color=color, linewidth=0.8)
        ax2.plot(site['site'][0], height - site['site'][1], 'ko', markersize=2)

    ax1.set_xlim(0, width)
    ax1.set_ylim(0, height)
    ax1.set_aspect('equal')
    ax1.set_title('Diagram')

    ax2.set_xlim(0, width)
    ax2.set_ylim(0, height)
    ax2.set_aspect('equal')
    ax2.set_title('Merge Process')

    plt.tight_layout()
    plt.savefig('python_ai/voronoi_demo.png', dpi=150)
    plt.close()
    print("Diagram saved to python_ai/voronoi_demo.png")
    print("Sites:", len(vector_points))


if __name__ == "__main__":
    main()
