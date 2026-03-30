import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def display(pngFile):
    image = mpimg.imread(pngFile)
    plt.imshow(image)
    plt.show()


pngs = [
        'defrost.png',
       ]

for png in pngs:
    display(png)

