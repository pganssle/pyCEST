"""
Handles the GUI creation functions
"""
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon

from matplotlib.colors import colorConverter

import numpy as np

from collections import namedtuple

Point = namedtuple('Point', ('x', 'y'), verbose=False)

LEFTMOUSE = 1
MIDDLEMOUSE = 2
RIGHTMOUSE = 3


class PolygonSelector(object):
    CLEAN = 1
    DRAWING = 2
    DRAWN = 4

    def __init__(self, ax, onselect=None, polygonprops=None):
        self.state = self.CLEAN

        ### Line-drawing stuff
        self.polygon = None

        self.poly_vertices = []
        self.polygonprops = polygonprops or {}

        ### Figure-related stuff
        self.ax = ax
        self.fig = self.ax.get_figure()
        self.canvas = self.fig.canvas

        ### Callback-related stuff
        self._cids = {}
        self.onselect = onselect

        self.connect_default_events()

    def press(self, event):
        if event.inaxes and self.state & (self.CLEAN | self.DRAWING):
            point = self._get_event_position(event)

            completing_polygon = ((event.button == LEFTMOUSE and event.dblclick == True) or
                                  (event.button == RIGHTMOUSE and self.state & self.DRAWING))
            adding_point = event.button == LEFTMOUSE or completing_polygon

            if adding_point:
                self._add_point_to_poly(point)

            if completing_polygon:
                self._complete_polygon()

            self.canvas.draw()

    def move(self, event):
        if event.inaxes:
            ax = event.inaxes
            
            if (self.state & self.DRAWING and
                (event.button == None or event.button == 1) and
                self.polygon != None):
                point = self._get_event_position(event)

                self._update_poly(point)

                self.canvas.draw()

    def _add_point_to_poly(self, point):
        # Lock our last line in place
        self._update_poly(point)
        self.poly_vertices.append(point)

        # Create a new line
        if self.polygon is None:
            self._new_polygon(point)

    def _get_event_position(self, event):
        return Point(event.xdata, event.ydata)

    def _new_polygon(self, point):
        self.state = self.DRAWING
        self.polygon = plt.Polygon([point], closed=False, **self.polygonprops)

        self.ax.add_patch(self.polygon)

    def _update_poly(self, new_point):
        if self.polygon is not None:
            self.polygon.set_xy(self.poly_vertices + [new_point])

    def _complete_polygon(self):
        self.polygon.set_closed(True)
        self.fig.canvas.draw()
        self.state = self.DRAWN

        if self.onselect:
            self.onselect(self.poly_vertices)

    def connect_default_events(self):
        self.connect_event('motion_notify_event', self.move, 'move')
        self.connect_event('button_press_event', self.press, 'press')

    def connect_event(self, event, callback, alias=None):
        """ Connect callback with an event       """
        alias = alias or event
        cid = self.canvas.mpl_connect(event, callback)
        self._cids[alias] = cid

        return cid

    def disconnect_events(self):
        for k, cid in self._cids.items():
            self.canvas.mpl_disconnect(cid)

        self.cids = {}


class ROIHandler:
    def __init__(self, axis, **poly_kwargs):
        """
        Supply an axis on which to plot the boxes
        """
        self.ax = axis
        self.rois = []
        self.selector = None

        dflt_kwargs = dict(fc='#0000b3',
                           fa=0.25,
                           ec='k',
                           ea=1.0,
                           lw=1)

        dflt_kwargs.update(poly_kwargs)
        poly_kwargs = dflt_kwargs

        # Allow 'face_alpha' / 'fa' and 'edge_alpha' / 'ea'
        face_alpha = poly_kwargs.pop('fa', None)
        face_alpha = poly_kwargs.pop('face_alpha', face_alpha)

        face_color = poly_kwargs.pop('fc', None)
        face_color = poly_kwargs.pop('face_color', face_color)
        poly_kwargs['fc'] = colorConverter.to_rgba(face_color, alpha=face_alpha)

        edge_alpha = poly_kwargs.pop('ea', None)
        edge_alpha = poly_kwargs.pop('edge_alpha', edge_alpha)

        edge_color = poly_kwargs.pop('ec', None)
        edge_color = poly_kwargs.pop('edge_color', edge_color)
        poly_kwargs['ec'] = colorConverter.to_rgba(edge_color, alpha=edge_alpha)

        self.poly_kwargs = poly_kwargs

    def connect(self):
        self.selector = PolygonSelector(self.ax, self.onselect,
                                        polygonprops=self.poly_kwargs)
        self.selector.connect_default_events()


    def disconnect(self):
        self.selector.disconnect_events()
        self.selector = None

    def onselect(self, polygon):
        self.rois.append(polygon)
        self.disconnect()

