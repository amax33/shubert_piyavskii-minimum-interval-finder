import numpy as np
import matplotlib.pyplot as plt
import math
from collections import namedtuple

def example_function(x):
    return math.pow(x-2,2) * math.sin(3*x)

Pt = namedtuple('Pt', ['x', 'y'])

def _get_sp_intersection(A, B, l):
    t = ((A.y - B.y) - l*(A.x - B.x)) / (2*l)
    return Pt(A.x + t, A.y - t*l)

def shubert_piyavskii(f, a, b, l, ϵ, δ=0.01):
    m = (a + b) / 2
    A, M, B = Pt(a, f(a)), Pt(m, f(m)), Pt(b, f(b))
    pts = [A, _get_sp_intersection(A, M, l), M, _get_sp_intersection(M, B, l), B]
    Δ = float('inf')
    while Δ > ϵ:
        i = min(range(len(pts)), key=lambda i: pts[i].y)
        P = Pt(pts[i].x, f(pts[i].x))
        Δ = P.y - pts[i].y
        P_prev = _get_sp_intersection(pts[i-1], P, l)
        P_next = _get_sp_intersection(P, pts[i+1], l)
        del pts[i]
        pts[i:i] = [P_next, P, P_prev]
    intervals = []
    min_index = 2 * (min(range(0, len(pts), 2), key=lambda i: pts[i].y)) - 1
    i = min_index
    for j in range(1, len(pts), 2):
        if pts[j].y < pts[i].y:
            dy = pts[i].y - pts[j].y
            x_lo = max(a, pts[j].x - dy / l)
            x_hi = min(b, pts[j].x + dy / l)
            if intervals and intervals[-1][1] + δ >= x_lo:
                intervals[-1] = (intervals[-1][0], x_hi)
            else:
                intervals.append((x_lo, x_hi))
    return intervals

intervals = shubert_piyavskii(example_function, 0, 4, 6, 0.01)
print("Intervals where the function might contain its minimum:")
for interval in intervals:
    print(interval)
