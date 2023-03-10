"""
This file aims at testing utils functions
"""

import numpy as np
from time import time

from thesis.utils import sort


if __name__ == "__main__":
    random_list = np.random.uniform(0, 10, 1000)
    random_list2 = random_list.copy()

    start = time()
    a = sort(random_list)
    end = time()
    print(f"Pierre sort time: {end - start:.3f} s")

    start = time()
    b = np.sort(random_list2)
    end = time()
    print(f"Numpy  sort time: {end - start:.3f} s")

    assert a.all() == b.all()
