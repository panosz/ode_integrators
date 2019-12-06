if __name__ == "__main__":
    from dynamic_analysis import harmonics as hr

    a = hr.HarmonicSeries([1, 2, 3, 4, 5])

    print(a.nth_amplitude(3))
