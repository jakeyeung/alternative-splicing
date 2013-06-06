'''
Created on 2013-06-05

@author: jyeung
Copied from 
http://kitchingroup.cheme.cmu.edu/software/python/matplotlib/interactive-annotations-in-matplotlib
'''


import math

import pylab
import matplotlib


class AnnoteFinder:
    """
    callback for matplotlib to display an annotation when points are clicked on.    The
    point which is closest to the click and within xtol and ytol is identified.
        
    Register this function like this:
        
    scatter(xdata, ydata)
    af = AnnoteFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
        self.data = zip(xdata, ydata, annotes)
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if axis is None:
            self.axis = pylab.gca()
        else:
            self.axis= axis
        self.drawnAnnotations = {}
        self.links = []

    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return(math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 ))

    def __call__(self, event):
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            if (self.axis is None) or (self.axis==event.inaxes):
                annotes = []
                for x,y,a in self.data:
                    if    (clickX-self.xtol < x < clickX+self.xtol) and    (clickY-self.ytol < y < clickY+self.ytol) :
                        annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
                if annotes:
                    annotes.sort()
                    _, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)

    def drawAnnote(self, axis, x, y, annote):
        """
        Draw the annotation on the plot
        """
        if self.drawnAnnotations.has_key((x,y)):
            markers = self.drawnAnnotations[(x,y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.axis.figure.canvas.draw()
        else:
            # t = axis.text(x,y, " - %s"%(annote), )
            t = axis.text(x,y, "%s"%(annote), horizontalalignment='left', verticalalignment='bottom')
            m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x,y)] =(t,m)
            self.axis.figure.canvas.draw()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
        for x,y,a in annotesToDraw:
            self.drawAnnote(self.axis, x, y, a)


if __name__ == '__main__':
    '''
    Example code...
    '''
    import matplotlib.pylab as plt

    x = range(10)
    y = range(10)
    annotes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

    plt.scatter(x,y)
    af =  AnnoteFinder(x,y, annotes)
    plt.connect('button_press_event', af)

    plt.show()