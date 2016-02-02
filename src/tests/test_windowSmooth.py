from unittest import TestCase

import numpy as np

from aquarium.geom.smooth import WindowSmooth


class TestWindowSmooth(TestCase):
    simple_data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
    simple_output = np.array([[5., 6., 7., ], [5.125, 6.125, 7.125], [5.125, 6.125, 7.125], [5.125, 6.125, 7.125]])  # w3 r2 mean

    def test_shape(self):
        for w in range(1, 20):
            for r in range(1, 3):
                ws = WindowSmooth(window=w, recursive=r)
                output = ws(self.simple_data)
                self.assertEquals(output.shape, self.simple_data.shape)



    def test_value(self):
        ws = WindowSmooth(window=3, recursive=2)
        for r,o in zip(ws(self.simple_data),self.simple_output):
            for e,l in zip(r,o):
                self.assertEquals(e,l)