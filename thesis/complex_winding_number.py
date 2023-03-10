"""
This file implements a routine to perform the "point in polygon" inclusion test
"""

###### Complex winding number ############

#!/usr/bin/env python
#
# routine for performing the "point in polygon" inclusion test

# Copyright 2001, softSurfer (www.softsurfer.com)
# This code may be freely used and modified for any purpose
# providing that this copyright notice is included with it.
# SoftSurfer makes no warranty for this code, and cannot be held
# liable for any real or imagined damage resulting from its use.
# Users of this code must verify correctness for their application.

# translated to Python by Maciej Kalisiak <mac@dgp.toronto.edu>

#   a Point is represented as a tuple: (x,y)

# ===================================================================

# is_left(): tests if a point is Left|On|Right of an infinite line.

#   Input: three points P0, P1, and P2
#   Return: >0 for P2 left of the line through P0 and P1
#           =0 for P2 on the line
#           <0 for P2 right of the line
#   See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"


def is_left(P0, P1, P2):
    return (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])


# ===================================================================

# wn_PnPoly(): winding number test for a point in a polygon
#     Input:  P = a point,
#             V[] = vertex points of a polygon
#     Return: wn = the winding number (=0 only if P is outside V[])


def wn_PnPoly(P, V):
    wn = 0  # the winding number counter

    # repeat the first vertex at end
    V = tuple(V[:]) + (V[0],)

    # loop through all edges of the polygon
    for i in range(len(V) - 1):  # edge from V[i] to V[i+1]
        if V[i][1] <= P[1]:  # start y <= P[1]
            if V[i + 1][1] > P[1]:  # an upward crossing
                if is_left(V[i], V[i + 1], P) > 0:  # P left of edge
                    wn += 1  # have a valid up intersect
        else:  # start y > P[1] (no test needed)
            if V[i + 1][1] <= P[1]:  # a downward crossing
                if is_left(V[i], V[i + 1], P) < 0:  # P right of edge
                    wn -= 1  # have a valid down intersect
    return wn


def translate(L):
    """
    take a list of complex values
    return the list of the coordinates of each complex in the R^2 plan
    """
    V = []
    for z in L:
        V.append((z.real, z.imag))
    return V


def Indice(L):
    """
    take a list of (x_i,x_j) representing a polygon
    return the winding number of the origin with respect to the curve
    """
    return wn_PnPoly((0, 0), translate(L))
